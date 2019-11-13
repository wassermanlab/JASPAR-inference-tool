#!/usr/bin/env python

import os, re
import argparse
from Bio import SearchIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SubsMat.MatrixInfo import blosum62
from functools import partial
import json
import math
from multiprocessing import Pool
import numpy as np
import shutil
import string
import subprocess
import sys
from tqdm import tqdm

# Append JASPAR-profile-inference to path
root_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(root_dir)

# Import globals
from __init__ import Jglobals

#-------------#
# Functions   #
#-------------#

usage_msg = """
usage: %s --fasta-file FILE --files-dir DIR
                        --models-dir DIR
""" % os.path.basename(__file__)

help_msg = """%s
  --fasta-file FILE   one or more sequences in FASTA format
  --files-dir DIR     output directory from get_files.py
  --models-dir DIR    output directory from regression.py

optional arguments:
  -h, --help          show this help message and exit
  --dummy-dir DIR     dummy directory (default = /tmp/)
  --output-file FILE  output file (default = STDOUT)
  --threads INT       threads to use (default = 1)

inference arguments:
  -l, --latest        use the latest version of each profile
  --rost INT          n parameter for the Rost's sequence
                      identity curve (default = 5; i.e. ~99%%
                      of correctly assigned homologs)
  --taxon [STR ...]   taxon(s) to use (default = all)
""" % usage_msg

def parse_args():
    """
    This function parses arguments provided via the command line and returns an {argparse} object.
    """

    # Initialize
    parser = argparse.ArgumentParser(add_help=False)

    # Mandatory args
    parser.add_argument("--fasta-file")
    parser.add_argument("--files-dir")
    parser.add_argument("--models-dir")

    # Optional args
    optional_group = parser.add_argument_group("optional arguments")
    optional_group.add_argument("-h", "--help", action="store_true")
    optional_group.add_argument("--dummy-dir", default="/tmp/")
    optional_group.add_argument("--output-file")
    optional_group.add_argument("-l", "--latest", action="store_true")
    optional_group.add_argument("--rost", default=5)
    optional_group.add_argument("--taxon", nargs="*", default=Jglobals.taxons)
    optional_group.add_argument("--threads", default=1)

    args = parser.parse_args()

    check_args(args)

    return(args)

def check_args(args):
    """
    This function checks an {argparse} object.
    """

    # Print help
    if args.help:
        print(help_msg)
        exit(0)

    # Check mandatory arguments
    if not args.fasta_file or not args.files_dir or not args.models_dir:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "arguments \"--fasta-file\" \"--files-dir\" \"--models-dir\" are required\n"]
        print(": ".join(error))
        exit(0)

    # Check "--threads" argument
    try:
        args.threads = int(args.threads)
    except:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"--threads\"", "invalid value", "\"%s\"\n" % args.threads]
        print(": ".join(error))
        exit(0)

    # Check "--rost" argument
    try:
        args.rost = int(args.rost)
    except:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"--rost\"", "invalid value", "\"%s\"\n" % args.rost]
        print(": ".join(error))
        exit(0)

    # Check "--taxon" argument
    for taxon in args.taxon:
        if taxon not in Jglobals.taxons:
            error = [os.path.basename(__file__), "error", "argument \"--taxon\"", "invalid value", "\"%s\"" % taxon]
            print(": ".join(error))
            exit(0)

    return(args) 

def main():

    # Parse arguments
    args = parse_args()

    # Infer profiles
    infer_profiles(args.fasta_file, args.files_dir, args.models_dir,
        args.dummy_dir, args.output_file, args.threads, args.latest,
        args.rost, args.taxon)

def infer_profiles(fasta_file, files_dir, models_dir, dummy_dir="/tmp/",
    output_file=None, threads=1, latest=False, n=5, taxons=Jglobals.taxons):

    # Initialize
    base_name = os.path.basename(__file__)
    pid = os.getpid()
    
    # Create dummy dir
    dummy_dir = os.path.join(dummy_dir, "%s.%s" % (base_name, pid))
    dummy_file = os.path.join(dummy_dir, "inferred_profiles.tsv")
    if not os.path.exists(dummy_dir):
        os.makedirs(dummy_dir)

    # Get sequences as SeqRecords
    # Note: https://biopython.org/wiki/SeqRecord
    seq_records = []
    for seq_record in Jglobals.parse_fasta_file(fasta_file):
        seq_records.append(seq_record)

    # Load JSON files
    global domains, jaspar, models, profiles
    domains, jaspar, models = _load_json_files(files_dir, models_dir)

    # Write
    columns = ["Query", "TF Name", "TF Matrix", "E-value", "Query Start-End",
        "TF Start-End", "DBD %ID", "Similarity Regression"]
    Jglobals.write(dummy_file, "\t".join(columns))

    # Infer SeqRecord profiles
    pool = Pool(threads)
    parallelized = partial(infer_SeqRecord_profiles, files_dir=files_dir,
        dummy_dir=dummy_dir, latest=latest, n=n, taxons=taxons)
    for inference_results in tqdm(pool.imap(parallelized, seq_records),
        total=len(seq_records)):
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
            Jglobals.write(dummy_file, "\t".join(map(str, inference_results[i])))
    pool.close()
    pool.join()

    # Write
    if output_file:
        shutil.copy(dummy_file, output_file)
    else:
        with open(dummy_file) as f:
            # For each line...
            for line in f:
                Jglobals.write(None, line.strip("\n"))

    # Remove dummy dir
    shutil.rmtree(dummy_dir)

def _load_json_files(files_dir, models_dir, taxons=Jglobals.taxons):

    # Initialize
    jaspar = {}
    pfams = {}
    profiles = {}
    uniprots = {}

    with open(os.path.join(files_dir, "pfam-DBDs.json")) as f:
        domains = json.load(f)

    for taxon in taxons:
        with open(os.path.join(files_dir, "%s.pfam.json" % taxon)) as f:
            for key, values in json.load(f).items():
                pfams.setdefault(key, values)
        with open(os.path.join(files_dir, "%s.profiles.json" % taxon)) as f:
            for key, values in json.load(f).items():
                profiles.setdefault(key, values)
        with open(os.path.join(files_dir, "%s.uniprot.json" % taxon)) as f:
            for key, values in json.load(f).items():
                uniprots.setdefault(key, values)

    for uniprot in pfams:

        # Add Pfam domains
        jaspar.setdefault(uniprot, {})
        jaspar[uniprot].setdefault("pfam", pfams[uniprot])

        # Add profiles
        jaspar[uniprot].setdefault("profiles", [])
        for profile in uniprots[uniprot][0]:
            jaspar[uniprot]["profiles"].append([profile, profiles[profile]])

    handle = Jglobals._get_file_handle(os.path.join(models_dir, "models.json.gz"))
    models = json.load(handle)
    handle.close()

    return(domains, jaspar, models)

def infer_SeqRecord_profiles(seq_record, files_dir, dummy_dir="/tmp/",
    latest=False, n=5, taxons=Jglobals.taxons):

    # Initialize
    cutoffs = {}
    inference_results = []

    # BLAST+ search
    blast_results = BLAST(seq_record, files_dir, taxons)

    # Get Pfam DBDs
    SeqRecord_Pfam_alignments = _get_SeqRecord_Pfam_alignments(seq_record,
        files_dir, dummy_dir)
    SeqRecord_DBDs = [v[0] for v in SeqRecord_Pfam_alignments]
    SeqRecord_alignments = [v[1] for v in SeqRecord_Pfam_alignments] 
    pfam_results = _get_results_Pfam_alignments(blast_results)

    # Skip if no Pfam DBDs
    if len(SeqRecord_Pfam_alignments) == 0:
        return(inference_results)

    # Filter results
    filtered_results = _filter_results(blast_results, pfam_results,
        SeqRecord_DBDs, n)

    # Get cut-offs on the percentage of sequence identity
    for pfam_ac in domains:
        if domains[pfam_ac][0] in SeqRecord_DBDs:
            cutoffs.setdefault(domains[pfam_ac][0], domains[pfam_ac][1])

    # Get similarity, coefficients and Y cut-off
    DBDs = "+".join(SeqRecord_DBDs)
    if DBDs in models:
        similarity, coeffs, Y = models[DBDs][:3]
        coeffs = np.array(coeffs)
    else:
        similarity = None

    # For each result...
    for result in filtered_results:

        # Inference: percentage of sequence identity
        pids = []
        pid_cutoffs = []
        for a in range(len(SeqRecord_alignments)):
            s1 = _removeLowercase(SeqRecord_alignments[a])
            s2 = _removeLowercase(pfam_results[result[1]][1][a])
            pids.append(sum(_fetchXs(s1, s2))/float(len(s1)))
            pid_cutoffs.append(pids[-1] >= cutoffs[SeqRecord_DBDs[a]])
        if True in pid_cutoffs:
            identities = np.mean(np.array(pids))
        else:
            identities = None

        # Inference: similarity regression
        Xs = []
        similarity_regression = None
        if similarity is not None:
            sr = []
            sr_cutoff = False
            for a in range(len(SeqRecord_alignments)):
                s1 = _removeLowercase(SeqRecord_alignments[a])
                s2 = _removeLowercase(pfam_results[result[1]][1][a])
                Xs.extend(_fetchXs(s1, s2, similarity=similarity))
            similarity_regression = sum(np.array(Xs) * coeffs)
            if similarity_regression < Y:
                similarity_regression = None

        # Add result
        if identities is not None or similarity_regression is not None:
            for matrix, gene_name in jaspar[result[1]]["profiles"]:
                inference_results.append([result[0], gene_name, matrix,
                    result[4], result[2], result[3], identities,
                    similarity_regression])

    return(inference_results)

def BLAST(seq_record, files_dir, taxons=Jglobals.taxons):

    # Initialize
    blast_results = set()
    outfmt = "sseqid pident length qstart qend sstart send evalue bitscore ppos qlen slen"

    # For each taxon...
    for taxon in taxons:

        # Taxon db
        taxon_db = os.path.join(files_dir, "%s.fa" % taxon)

        # Run BLAST+
        cmd = "blastp -db %s -outfmt \"6 %s\"" % (taxon_db, outfmt)
        process = subprocess.Popen([cmd], shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
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
            if len(blast_record) != 12: continue

            # Get BLAST+ record
            target_id = blast_record[0]
            percent_identities = float(blast_record[1])
            alignment_length = int(blast_record[2])
            query_start_end = "%s-%s" % (blast_record[3], blast_record[4])
            target_start_end = "%s-%s" % (blast_record[5], blast_record[6])
            e_value = float(blast_record[7])
            score = float(blast_record[8])
            percent_similarity = float(blast_record[9])
            query_aligned_residues = int(blast_record[4]) - int(blast_record[3]) + 1
            query_length = float(blast_record[10])
            target_aligned_residues = int(blast_record[6]) - int(blast_record[5]) + 1
            target_length = float(blast_record[11])
            query_coverage = query_aligned_residues * 100 / query_length
            target_coverage = target_aligned_residues * 100 / target_length
            joint_coverage = math.sqrt(query_coverage * target_coverage)

            # Add BLAST+ record to search results
            blast_results.add((seq_record.id, target_id, query_start_end, target_start_end, e_value, score, percent_identities, alignment_length, percent_similarity, joint_coverage))

    # Return results sorted by score
    return(list(sorted(blast_results, key=lambda x: x[-1], reverse=True)))

def _get_SeqRecord_Pfam_alignments(seq_record, files_dir, dummy_dir="/tmp/"):

    # Initialize
    alignments = []
    hmm_db = os.path.join(files_dir, "pfam-DBDs", "all_DBDs.hmm")

    # Make seq file
    seq_file = os.path.join(dummy_dir, ".seq.fasta")
    _makeSeqFile(seq_record, seq_file)

    # For each DBD...
    for pfam_id_std, start, end, evalue in hmmScan(seq_file, hmm_db, dummy_dir,
        non_overlapping_domains=True):

        # Initialize
        hmm_file = os.path.join(files_dir, "pfam-DBDs", "%s.hmm" % pfam_id_std)

        # Make seq file
        sub_seq_record = SeqRecord(seq_record.seq[start:end], id=seq_record.id,
            name=seq_record.name, description=seq_record.description)
        _makeSeqFile(sub_seq_record, seq_file)

        # Add DBDs
        alignment = hmmAlign(seq_file, hmm_file)
        alignments.append((pfam_id_std, alignment, start+1, end, evalue))

    return(alignments)

def _makeSeqFile(seq_record, file_name=".seq.fa"):

    # Remove seq file if exists...
    if os.path.exists(file_name):
        os.remove(file_name)

    # Write
    Jglobals.write(file_name, seq_record.format("fasta"))

def hmmScan(seq_file, hmm_file, dummy_dir="/tmp/", non_overlapping_domains=False):

    # Initialize
    out_file = os.path.join(dummy_dir, ".out.txt")

    # Scan
    cmd = "hmmscan --domtblout %s %s %s" % (out_file, hmm_file, seq_file)
    process = subprocess.run([cmd], shell=True, stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL)

    # Read domains
    domains = _readDomainsTab(out_file)

    # Remove output file
    if os.path.exists(out_file):
        os.remove(out_file)

    # Filter overlapping domains
    if non_overlapping_domains:
        domains = _getNonOverlappingDomains(domains)

    # Yield domains one by one
    for pfam_ac, start, end, evalue in sorted(domains, key=lambda x: x[1]):

        yield(pfam_ac, start, end, evalue)

def _readDomainsTab(file_name):
    """
    From PMID:22942020;
    A hit has equal probability of being in the same clan as a different clan
    when the E-value is 0.01 (log 10 = -2). When the E-value is 10-5, the pro-
    bability that a sequence belongs to the same clan is >95%.

    From CIS-BP paper;
    We scanned all protein sequences for putative DNA-binding domains (DBDs)
    using the 81 Pfam (Finn et al., 2010) models listed in (Weirauch and
    Hughes, 2011) and the HMMER tool (Eddy, 2009), with the recommended de-
    tection thresholds of Per-sequence Eval < 0.01 and Per-domain conditional
    Eval < 0.01.
    """

    # Initialize
    domains = []
    cutoff_mod = 1e-5
    cutoff_dom = 0.01

    # For each result...
    for res in SearchIO.parse(file_name, "hmmscan3-domtab"):

        # For each model...
        for mod in res.iterhits():

            # Skip poor models
            if mod.evalue > cutoff_mod:
                continue

            # For each domain...
            for dom in mod.hsps:

                # Skip poor domains
                if dom.evalue_cond > cutoff_dom:
                    continue

                # Append domain
                domains.append((mod.id, dom.query_start, dom.query_end,
                    dom.evalue_cond))

    return(domains)

def _getNonOverlappingDomains(domains):
    """
    Do domains 1 & 2 overlap?
    ---------1111111---------
    -------22222-------------  True
    ----------22222----------  True
    -------------22222-------  True
    -----22222---------------  False
    ---------------22222-----  False
    """

    # Initialize
    nov_domains = []

    # Sort domains by e-value
    for domain in sorted(domains, key=lambda x: x[-1]):

        # Initialize
        domains_overlap = False

        # For each non-overlapping domain...
        for nov_domain in nov_domains:

            if domain[1] < nov_domain[2] and domain[2] > nov_domain[1]:
                domains_overlap = True
                break

        # Add non-overlapping domain
        if not domains_overlap:
            nov_domains.append(domain)

    return(nov_domains)

def hmmAlign(seq_file, hmm_file):

    # Align
    cmd = "hmmalign --outformat PSIBLAST %s %s" % (hmm_file, seq_file)
    process = subprocess.check_output([cmd], shell=True, universal_newlines=True)

    return(_readPSIBLASToutformat(process))

def _readPSIBLASToutformat(psiblast_alignment):

    # Initialize
    alignment = ""

    # For each chunk...
    for chunk in psiblast_alignment.split("\n"):

        # If alignment substring...
        m = re.search("\s+(\S+)$", chunk)
        if m:
            alignment += m.group(1)

    return(alignment)

def _get_results_Pfam_alignments(blast_results):

    # Initialize
    results_Pfam_alignments = {}

    # For each BLAST result...        
    for blast_result in blast_results:

        # Unwind
        DBDs = [v[0] for v in jaspar[blast_result[1]]["pfam"]]
        alignments = [v[1] for v in jaspar[blast_result[1]]["pfam"]]
        results_Pfam_alignments.setdefault(blast_result[1], [DBDs, alignments])

    return(results_Pfam_alignments)

def _filter_results(blast_results, pfam_results, DBDs, n=5):

    # Intialize
    filtered_results = []

    # For each filtered result...
    for filtered_result in _filter_results_by_Rost(blast_results, n=5):

        if filtered_result[1] in pfam_results:
            if pfam_results[filtered_result[1]][0] == DBDs:
                filtered_results.append(filtered_result)

    return(filtered_results)

def _filter_results_by_Rost(blast_results, n=5):

    # Initialize
    blast_homologs = []

    # For each result...
    for result in blast_results:

        # Initialize
        percent_identity = result[6]
        alignment_length = result[7]

        # If homologs...
        if _is_alignment_over_Rost_seq_id_curve(percent_identity, alignment_length, n):

            # Add homolog
            blast_homologs.append(result)

    return(blast_homologs)

def _is_alignment_over_Rost_seq_id_curve(percent_identity, L, n=5):
    """
    This function returns whether an alignment is over the Rost's pairwise
    sequence identity curve or not.
    """
    return(percent_identity >= _get_Rost_cutoff_percent_identity(L, n))

def _get_Rost_cutoff_percent_identity(L, n=5):
    """
    This function returns the Rost's cut-off percentage of identical residues
    for an alignment of length "L".
    """
    return(n + (480 * pow(L, -0.32 * (1 + pow(math.e, float(-L) / 1000)))))

def _is_alignment_over_Rost_seq_sim_curve(percent_similarity, L, n=12):
    """
    This function returns whether an alignment is over the Rost's pairwise
    sequence similarity curve or not.
    """
    return(percent_similarity >= _get_Rost_cutoff_percent_similarity(L, n))

def _get_Rost_cutoff_percent_similarity(L, n=12):
    """
    This function returns the Rost's cut-off percentage of "similar" residues
    for an alignment of length "L".
    """
    return(n + (420 * pow(L, -0.335 * (1 + pow(math.e, float(-L) / 2000)))))

# def _SeqRecord_profile_inference(seq_record, uniacc, files_dir):

#     # Initialize
#     inference_results = {}

#     # Load JSON files
#     global domains, jaspar
#     try:
#         domains, jaspar
#     except NameError:
#         domains, jaspar = _load_json_files(files_dir)

#     # If domains...
#     if uniacc in domains:
#         # For each domain...
#         for domain in domains[uniacc][0]:
#             # For each pairwise alignment...
#             for alignment in _pairwise_alignment(
#                 seq_record.seq, domain):
#                 # If alignment does not satisfy the threshold...
#                 identities = _get_alignment_identities(
#                     alignment[0], alignment[1]) / float(len(domain))
#                 if identities >= float(domains[uniacc][1]):
#                     # For each JASPAR matrix... #
#                     for matrix, gene_name in jaspar[uniacc]:
#                         # Infer matrix
#                         inference_results.setdefault((gene_name, matrix), identities)
#                         if identities > inference_results[(gene_name, matrix)]:
#                             inference_results[(gene_name, matrix)] = identities

#     return [[i[0], i[1], inference_results[i]] for i in inference_results]

# def _pairwise_alignment(A, B):
#     """
#     This function returns the alignments between two sequences "A" and "B" using
#     dynamic programming.
#     """
#     try:
#         # Parameters from EMBOSS needle
#         return pairwise2.align.globalds(A, B, blosum62, -10.0, -0.5)
#     except:
#         return []

# def _get_alignment_identities(A, B):
#     """
#     This function returns the number of identities between two aligned sequences
#     "A" and "B". If "A" and "B" have different lengths, returns None.
#     """
#     if len(A) == len(B):
#         return len([i for i in range(len(A)) if A[i] == B[i]])

#     return None

def _removeLowercase(s):

    return(s.translate(str.maketrans("", "", string.ascii_lowercase)))

def _fetchXs(seq1, seq2, similarity="identity"):
    """
    Called for each comparison to compare the strings.
    """

    # Fill Xs with zeroes
    scores = [0] * len(seq1)

    # Strings should have same length
    if len(seq1) == len(seq2):

        for n in range(len(seq1)):

            if similarity == "identity":
                scores[n] = _IDscoring(seq1[n], seq2[n])

            elif similarity == "blosum62":
                scores[n] = _BLOSUMscoring(seq1[n], seq2[n])

    else:
        # there is something wrong
        print("Strings have different length!\n\tA: %s\n\tB: %s" % (seq1, seq2))
        exit(0)

    return(scores)

def _IDscoring(aa1, aa2):

    if aa1 == aa2 and aa1 != "-":
        return(1)
    else:
        return(0)

def _BLOSUMscoring(aa1, aa2):

    if aa1 == "-" and aa2 == "-":
        return(1)
    elif aa1 == "-" or aa2 == "-":
        return(-4)
    else:
        if (aa1, aa2) in blosum62:
            return(blosum62[(aa1, aa2)])
        else:
            return(blosum62[(aa2, aa1)])

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()