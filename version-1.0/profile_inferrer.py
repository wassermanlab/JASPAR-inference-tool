#!/usr/bin/env python

"""
tool:    jaspartools infer (i.e. profile_inferrer.py)
version: 0.0.1
summary: compares the DBD sequence of the given TF to those of homologous
         TFs in JASPAR, and infers the TF binding profiles from the best
         compared JASPAR homologous TFs as potentially recognized by the
         given TF.

usage:   jaspartools infer [-h] [--dummy] [-n] [-o] [--threads]
                           [--fungi] [--insects] [--nematodes]
                           [--plants] [--vertebrates] [-l] fasta
"""

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

# Import my functions
import functions

#-------------#
# Functions   #
#-------------#

def parse_args(extra=""):
    """
    This function parses arguments provided via the command line using argparse.
    """

    parser = argparse.ArgumentParser(
        prog="jaspartools infer",
        usage=argparse.SUPPRESS,
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("file", help="input file in FASTA format")
    parser.add_argument("files", help="files directory (from make_files.py)")

    # Optional args
    parser.add_argument("--dummy", default="/tmp/",
        help="dummy directory (default = /tmp/)", metavar="")
    parser.add_argument("-n", "--n-param", default=5, type=int, metavar="",
        help="\"n\" parameter for the Rost's curve (default = 5)")
    parser.add_argument("-o", "--out-file", metavar="",
        help="output file (default = stdout)")
    parser.add_argument("--threads", default=1, type=int, metavar="",
        help="number of threads to use (default = 1)")

    # JASPAR args
    jaspar_group = parser.add_argument_group("jaspar options")
    jaspar_group.add_argument("--fungi", action="store_true",
        help="use profiles from the JASPAR CORE fungi collection")
    jaspar_group.add_argument("--insects", action="store_true",
        help="use profiles from the JASPAR CORE insects collection")
    jaspar_group.add_argument("--nematodes", action="store_true",
        help="use profiles from the JASPAR CORE nematodes collection")
    jaspar_group.add_argument("--plants", action="store_true",
        help="use profiles from the JASPAR CORE plants collection")
    jaspar_group.add_argument("--vertebrates", action="store_true",
        help="use profiles from the JASPAR CORE vertebrates collection")
    jaspar_group.add_argument("-l", "--latest", default=False, action="store_true",
        help="use profiles from the lastest JASPAR release")

    return parser.parse_args()

def main():

    # Parse arguments
    args = parse_args()

    # Taxons
    taxons = []
    if args.fungi: taxons.append("fungi")
    if args.insects: taxons.append("insects")
    if args.nematodes: taxons.append("nematodes")
    if args.plants: taxons.append("plants")
    if args.vertebrates: taxons.append("vertebrates")
    if not taxons: taxons = ["fungi", "insects", "nematodes", "plants", "vertebrates"]

    # Infer profiles
    infer_profiles(args.file, args.files, args.dummy, args.latest,
        args.n_param, args.out_file, taxons, args.threads)

def infer_profiles(fasta_file, files_dir, dummy_dir="/tmp/", latest=False, n=5,
    output_file=None, taxons=["fungi", "insects", "nematodes", "plants", "vertebrates"],
    threads=1):

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
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        seq_records.append(seq_record)

    # Load JSON files
    global domains, jaspar
    domains, jaspar = _load_json_files(files_dir)

    # Write
    functions.write(dummy_file,
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
            functions.write(dummy_file, "\t".join(map(str, inference_results[i])))
    pool.close()
    pool.join()

    # Write
    if output_file:
        shutil.copy(dummy_file, output_file)
    else:
        with open(dummy_file) as f:
            # For each line...
            for line in f:
                functions.write(None, line.strip("\n"))

    # Remove dummy dir
    shutil.rmtree(dummy_dir)

def _load_json_files(files_dir):

    with open(os.path.join(files_dir, "domains.json")) as f:
        domains = json.load(f)
    with open(os.path.join(files_dir, "jaspar.json")) as f:
        jaspar = json.load(f)

    return domains, jaspar

def infer_SeqRecord_profiles(seq_record, files_dir, dummy_dir="/tmp/", latest=False,
    n=5, taxons=["fungi", "insects", "nematodes", "plants", "vertebrates"]):

    # Initialize
    inference_results = []

    # Homology search
    homology_search_results = _SeqRecord_homology_search(seq_record, files_dir,
        dummy_dir, n, taxons)

    # For each result...
    for result in homology_search_results:
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

def _SeqRecord_homology_search(seq_record, files_dir, dummy_dir="/tmp/", n=5,
    taxons=["fungi", "insects", "nematodes", "plants", "vertebrates"]):

    # Initialize
    search_results = set()

    # For each taxon...
    for taxon in taxons:
        # Taxon db
        taxon_db = os.path.join(files_dir, "%s.fa" % taxon)
        # Homology search
        try:
            process = subprocess.Popen([
                "blastp",
                "-db", taxon_db,
                "-outfmt", "6"],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE)
            fasta_sequence = ">%s\n%s" % (seq_record.id, seq_record.seq)
            process.stdin.write(fasta_sequence.encode())
            (blast_records, blast_errors) = process.communicate()
            # For each BLAST+ record...
            for blast_record in blast_records.decode("utf-8").split("\n"):
                # A BLAST+ record is formatted as a tab-separated list w/ 12 columns:
                # (1,2) identifiers for query and target sequences;
                # (3) percentage sequence identity
                # (4) alignment length;
                # (5) number of mismatches;
                # (6) number of gap openings;
                # (7-8, 9-10) start and end-position in query and in target;
                # (11) E-value; and
                # (12) bit score.
                blast_record = blast_record.split("\t")
                # Skip if not a BLAST+ record
                if len(blast_record) != 12: continue
                # Initialize
                target = blast_record[1]
                perc_identity = float(blast_record[2])
                alignment_length = int(blast_record[3])
                query_start_end = "%s-%s" % (blast_record[6], blast_record[7])
                target_start_end = "%s-%s" % (blast_record[8], blast_record[9])
                e_value = float(blast_record[10])
                score = float(blast_record[11])
                # If homologs...
                if _is_alignment_over_Rost_sequence_identity_curve(
                    round(alignment_length * perc_identity/100), alignment_length, n):
                    # Add homolog to search results
                    search_results.add((seq_record.id, target, query_start_end,
                        target_start_end, e_value, score))
        except:
            raise ValueError("Could not exec BLAST+!")

    # Return results sorted by score
    return list(sorted(search_results, key=lambda x: x[-1], reverse=True))

def _is_alignment_over_Rost_sequence_identity_curve(identities, alignment_length, n=5):
    """
    This function evaluates whether an alignment is over the Rost's sequence identity
    curve or not.
    """
    return identities >= _get_Rost_ID_threshold(alignment_length, n)

def _get_Rost_ID_threshold(L, n=5):
    """
    This function returns the Rost sequence identity threshold for a given alignment of
    length "L".
    """    
    return n + (480 * pow( L, float('-0.32') *\
        (1 + pow(float(repr(math.e)), float(repr(float(-L) / 1000))))))

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
    This function returns the number of identities between two aligned sequences "A"
    and "B". If "A" and "B" have different lengths, returns None.
    """
    if len(A) == len(B):
        return len([i for i in range(len(A)) if A[i] == B[i]])

    return None

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()