#!/usr/bin/env python

import os, re
import argparse
from Bio import pairwise2
from Bio import SeqIO
from Bio.SubsMat.MatrixInfo import blosum62
from functools import partial
import json
import math
from multiprocessing import Pool
import shutil
import subprocess
from tqdm import tqdm

# Import from JASPAR tools
from jaspartools import JTglobals

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line using argparse.
    """

    parser = argparse.ArgumentParser(description="compares the DBD sequence(s) of the given TF(s) to those of homologous TFs stored in JASPAR, and infers the TF(s) binding profile(s) from the best compared JASPAR homologous TF(s) as potentially recognized by the given TF(s).")

    parser.add_argument("fasta", help="FASTA file")
    parser.add_argument("files", help="files directory (from make_files.py)")

    # Optional args
    parser.add_argument("--dummy", metavar="DUMMYDIR", default="/tmp/",
        help="dummy directory (default = /tmp/)")
    parser.add_argument("-l", "--latest", default=False, action="store_true",
        help="use the lastest JASPAR version (default = False)")
    parser.add_argument("-n", default=5, type=int, 
        help="n parameter for the Rost's curve (e.g. n=5 corresponds to ~99%% correctly assigned homologs; default = 5)")
    parser.add_argument("-o", metavar="OUTFILE",
        help="output file (default = STDOUT)")
    parser.add_argument("--taxon", default=JTglobals.taxons,
        choices=JTglobals.taxons, nargs="*", metavar="TAXON",
        help="taxonomic group (default = \"%s\")" % " ".join(JTglobals.taxons))
    parser.add_argument("--threads", default=1, type=int,
        help="number of threads to use (default = 1)")

    return parser.parse_args()

def main():

    # Parse arguments
    args = parse_args()

    # Infer profiles
    infer_profiles(args.fasta, args.files, args.dummy, args.latest, args.n, args.o,
        args.taxon, args.threads)

def infer_profiles(fasta_file, files_dir, dummy_dir="/tmp/", latest=False, n=5,
    output_file=None, taxons=JTglobals.taxons, threads=1):

    # Initialize
    base_name = os.path.basename(__file__)
    dummy_dir = os.path.join(dummy_dir, "%s.%s" % (base_name, os.getpid()))
    dummy_file = os.path.join(dummy_dir, "inferred_profiles.tsv")

    # Create dummy dir
    if not os.path.exists(dummy_dir):
        os.makedirs(dummy_dir)

    # Get sequences as SeqRecords
    # Note: https://biopython.org/wiki/SeqRecord
    seq_records = []
    for seq_record in SeqIO.parse(JTglobals._get_file_handle(fasta_file), "fasta"):
        seq_records.append(seq_record)

    # Load JSON files
    global domains, jaspar
    domains, jaspar = _load_json_files(files_dir)

    # Write
    JTglobals.write(dummy_file,
        "Query\tTF Name\tTF Matrix\tE-value\tQuery Start-End\tTF Start-End\tDBD %ID")
    # Infer SeqRecord profiles
    pool = Pool(threads)
    parallelization = partial(infer_SeqRecord_profiles, files_dir=files_dir,
        dummy_dir=dummy_dir, latest=latest, n=n, taxons=taxons)
    for inference_results in tqdm(pool.imap(parallelization, iter(seq_records)),
        desc="Profile inference", total=len(seq_records)):
        # Sort by E-value, TF Name and Matrix
        if latest:
            inference_results.sort(key=lambda x: (x[3], x[1], -float(x[2][2:])))
        else:
            inference_results.sort(key=lambda x: (x[3], x[1], float(x[2][2:])))
        # For each inference...
        for i in range(len(inference_results)):
            # Use the lastest version of JASPAR
            if latest and i > 0:
                if inference_results[i][2][:6] == inference_results[i - 1][2][:6]:
                    continue
            # Write
            JTglobals.write(dummy_file, "\t".join(map(str, inference_results[i])))
    pool.close()
    pool.join()

    # Write
    if output_file:
        shutil.copy(dummy_file, output_file)
    else:
        with open(dummy_file) as f:
            # For each line...
            for line in f:
                JTglobals.write(None, line.strip("\n"))

    # Remove dummy dir
    shutil.rmtree(dummy_dir)

def _load_json_files(files_dir):

    with open(os.path.join(files_dir, "domains.json")) as f:
        domains = json.load(f)
    with open(os.path.join(files_dir, "jaspar.json")) as f:
        jaspar = json.load(f)

    return domains, jaspar

def infer_SeqRecord_profiles(seq_record, files_dir, dummy_dir="/tmp/", latest=False,
    n=5, taxons=JTglobals.taxons):

    # Initialize
    inference_results = []

    # Homology search
    blast_results = _SeqRecord_blast_search(seq_record, files_dir, dummy_dir, taxons)
    blast_homologs = _filter_results_under_the_Rost_seq_id_curve(blast_results, n)

    # For each result...
    for result in blast_homologs:
        # Initialize
        (query, target, query_start_end, target_start_end, e_value, score) = result
        # For each result...
        for result in _SeqRecord_profile_inference(seq_record, target, files_dir):
            # Initialize
            (gene_name, matrix, identities) = result
            # Add result
            inference_results.append([seq_record.id, gene_name, matrix, e_value,
                query_start_end, target_start_end, identities])

    return inference_results

def _SeqRecord_blast_search(seq_record, files_dir, taxons=JTglobals.taxons):

    # Initialize
    blast_results = set()
    outfmt = "6 sseqid pident length qstart qend sstart send evalue bitscore ppos qlen slen"

    # For each taxon...
    for taxon in taxons:

        # Taxon db
        taxon_db = os.path.join(files_dir, "%s.fa" % taxon)

        # Homology search
        try:

            # Run BLASTP
            
            process = subprocess.Popen([
                    "blastp",
                    "-db", taxon_db,
                    "-outfmt", outfmt
                ],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE)
            fasta_sequence = ">%s\n%s" % (seq_record.id, seq_record.seq)
            process.stdin.write(fasta_sequence.encode())
            (blast_records, blast_errors) = process.communicate()

            # For each BLAST+ record...
            for blast_record in blast_records.decode("utf-8").split("\n"):

                # Custom BLAST+ record:
                # (1) identifier of target sequence;
                # (2) percentage of identical matches;
                # (3) alignment length;
                # (4-5, 6-7) start and end-position in query and in target;
                # (8) E-value;
                # (9) bit score;
                # (10) percentage of positive-scoring matches; and
                # (4-7, 11, 12) joint coverage (i.e. square root of the coverage
                # on the query and the target).
                blast_record = blast_record.split("\t")

                # Skip if not a BLAST+ record
                if len(blast_record) != 9: continue

                # Get BLAST+ record
                target = blast_record[0]
                percent_identities = float(blast_record[1])
                alignment_length = int(blast_record[2])
                query_start_end = "%s-%s" % (blast_record[3], blast_record[4])
                target_start_end = "%s-%s" % (blast_record[5], blast_record[6])
                e_value = float(blast_record[7])
                score = float(blast_record[8])
                percent_similarity = float(blast_record[9])
                query_coverage = (int(blast_record[4]) - int(blast_record[3])) + 1) * 100 / float(blast_record[10])
                target_coverage = (int(blast_record[6]) - int(blast_record[5])) + 1) * 100 / float(blast_record[11])
                joint_coverage = math.sqrt(query_coverage * target_coverage)

                # Add BLAST+ record to search results
                blast_results.add((seq_record.id, target, query_start_end, target_start_end, e_value, score, percent_identities, alignment_length, percent_similarity, joint_coverage))

        except:
            raise ValueError("Could not exec BLAST+!")

    # Return results sorted by score
    return list(sorted(blast_results, key=lambda x: x[-1], reverse=True))

def _filter_results_under_the_Rost_seq_id_curve(blast_results, n=5):

    # Initialize
    blast_homologs = []

    # For each result...
    for result in blast_results:

        # If homologs...
        if _is_alignment_over_Rost_seq_id_curve(percent_identities, alignment_length, n):

            # Add homolog
            blast_homologs.append(result)

def _is_alignment_over_Rost_seq_id_curve(percent_identities, L, n=5):
    """
    This function returns whether an alignment is over the Rost's pairwise sequence
    identity curve or not.
    """
    return percent_identities >= _get_Rost_cutoff_percent_identities(L, n)

def _get_Rost_cutoff_percent_identities(L, n=5):
    """
    This function returns the Rost's cut-off percentage of identical residues for an
    alignment of length "L".
    """
    return n + (480 * pow(L, -0.32 * (1 + pow(math.e, float(-L) / 1000))))

def _SeqRecord_profile_inference(seq_record, uniacc, files_dir):

    # Initialize
    inference_results = {}

    # Load JSON files
    global domains, jaspar
    try:
        domains, jaspar
    except NameError:
        domains, jaspar = _load_json_files(files_dir)

    # If domains...
    if uniacc in domains:
        # For each domain...
        for domain in domains[uniacc][0]:
            # For each pairwise alignment...
            for alignment in _pairwise_alignment(
                seq_record.seq, domain):
                # If alignment does not satisfy the threshold...
                identities = _get_alignment_identities(
                    alignment[0], alignment[1]) / float(len(domain))
                if identities >= float(domains[uniacc][1]):
                    # For each JASPAR matrix... #
                    for matrix, gene_name in jaspar[uniacc]:
                        # Infer matrix
                        inference_results.setdefault((gene_name, matrix), identities)
                        if identities > inference_results[(gene_name, matrix)]:
                            inference_results[(gene_name, matrix)] = identities

    return [[i[0], i[1], inference_results[i]] for i in inference_results]

def _pairwise_alignment(A, B):
    """
    This function returns the alignments between two sequences "A" and "B" using
    dynamic programming.
    """
    try:
        # Parameters from EMBOSS needle
        return pairwise2.align.globalds(A, B, blosum62, -10.0, -0.5)
    except:
        return []

def _get_alignment_identities(A, B):
    """
    This function returns the number of identities between two aligned sequences
    "A" and "B". If "A" and "B" have different lengths, returns None.
    """
    if len(A) == len(B):
        return len([i for i in range(len(A)) if A[i] == B[i]])

    return None

#-------------#
# Main        #
#-------------#

if __name__ == "__main__": main()