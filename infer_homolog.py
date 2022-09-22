#!/usr/bin/env python

import argparse
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.SubsMat.MatrixInfo import blosum62
import copy
from functools import partial
import json
import math
from multiprocessing import Pool
import numpy as np
import os
import re
import shutil
import string
import subprocess
import sys
from tqdm import tqdm
bar_format="{percentage:3.0f}%|{bar:20}{r_bar}"
import warnings

# Defaults
root_dir = os.path.dirname(os.path.realpath(__file__))
files_dir = os.path.join(root_dir, "files")
# models_dir = os.path.join(root_dir, "models")

# Append JASPAR-profile-inference to path
sys.path.append(root_dir)

# Import globals
from __init__ import CisBP2Pfam, Jglobals, ReadSRModel, ScoreAlignmentResult
from infer_profile import __get_SeqRecord_Pfam_alignments

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns
    an {argparse} object.
    """

    # Initialize
    parser = argparse.ArgumentParser()

    # Mandatory args
    parser.add_argument("query",
        help="query sequence(s) in FASTA format")
    parser.add_argument("target",
        help="target sequence(s) in FASTA format")

    # Optional args
    parser.add_argument("--dummy-dir", default="/tmp/", metavar="DIR", 
        help="dummy directory (default = /tmp/)")
    parser.add_argument("--files-dir", default=files_dir, metavar="DIR",
        help="files directory from get_files.py (default = ./files/)")
    # # parser.add_argument("--models-dir", default=models_dir)
    parser.add_argument("--output-file", metavar="FILE",
        help="output file (default = STDOUT)")
    parser.add_argument("--threads", default=1, metavar="INT", type=int,
        help="number of threads to use (default = 1)")
    parser.add_argument("-w", "--warnings", action="store_true",
        help="issue warnings (default = False)")

    # Inference args
    inference_group = parser.add_argument_group("inference arguments")
    parser.add_argument("--no-blast", action="store_true",
        help="do not search for homologs using BLAST+ (default = False)")
    inference_group.add_argument("--rost", default=5, metavar="INT",
        help="\"n\" parameter for the Rost's curve (default = 5)")

    args = parser.parse_args()

    return(args)

def main():

    # Parse arguments
    args = parse_args()

    # Warnings
    if not args.warnings:
        warnings.filterwarnings("ignore")

    # Infer homologs
    infer_homologs(args.query, args.target, args.dummy_dir, args.files_dir,
        args.output_file, args.threads, args.no_blast, args.rost)

def infer_homologs(query_file, target_file, dummy_dir="/tmp/",
    files_dir=files_dir, output_file=None, threads=1, no_blast=False, n=5):

    # Initialize
    base_name = os.path.basename(__file__)
    pid = os.getpid()

    # Load data
    cisbp = __load_CisBP_models(files_dir)

    # Format BLAST+ database
    if not no_blast:
        __format_BLAST_database(target_file)

    # Create dummy dir
    dummy_dir = os.path.join(dummy_dir, "%s.%s" % (base_name, pid))
    dummy_file = os.path.join(dummy_dir, "inferred_homologs.tsv")
    if not os.path.exists(dummy_dir):
        os.makedirs(dummy_dir)

    # Get sequences as SeqRecords
    # Note: https://biopython.org/wiki/SeqRecord
    query_seq_records = []
    for seq_record in Jglobals.parse_fasta_file(query_file):
        query_seq_records.append(seq_record)
    target_seq_records = []
    for seq_record in Jglobals.parse_fasta_file(target_file):
        target_seq_records.append(seq_record)

    # Get SeqRecords Pfam alignments
    pfam_alignments = {}
    kwargs = {"total": len(query_seq_records), "bar_format": bar_format}
    pool = Pool(min([threads, len(query_seq_records)]))
    p = partial(__get_SeqRecord_Pfam_alignments, files_dir=files_dir,
        dummy_dir=dummy_dir)
    for id_alignments in tqdm(pool.imap(p, query_seq_records), **kwargs):
        if id_alignments[0] in pfam_alignments or len(id_alignments[1]) == 0:
            continue
        pfam_alignments.setdefault(id_alignments[0], {})
        pfam_alignments[id_alignments[0]].setdefault(id_alignments[1][0][0], [])
        pfam_alignments[id_alignments[0]][id_alignments[1][0][0]].append(id_alignments[1][0][1])
    kwargs = {"total": len(target_seq_records), "bar_format": bar_format}
    pool = Pool(min([threads, len(target_seq_records)]))
    for id_alignments in tqdm(pool.imap(p, target_seq_records), **kwargs):
        if id_alignments[0] in pfam_alignments or len(id_alignments[1]) == 0:
            continue
        pfam_alignments.setdefault(id_alignments[0], {})
        pfam_alignments[id_alignments[0]].setdefault(id_alignments[1][0][0], [])
        pfam_alignments[id_alignments[0]][id_alignments[1][0][0]].append(id_alignments[1][0][1])

    # Write
    columns = ["Query", "Target", "E-value", "Query Start-End",
        "Target Start-End", "DBD %ID"]
    Jglobals.write(dummy_file, "\t".join(columns))

    # Infer SeqRecords homologs
    kwargs = {"total": len(query_seq_records), "bar_format": bar_format}
    pool = Pool(min([threads, len(query_seq_records)]))
    p = partial(infer_SeqRecord_homologs, target_file=target_file,
        pfam_alignments=pfam_alignments, cisbp=cisbp, dummy_dir=dummy_dir,
        files_dir=files_dir, no_blast=no_blast, n=n)
    for inferences in tqdm(pool.imap(p, query_seq_records), **kwargs):
        for inference in inferences:
            Jglobals.write(dummy_file, "\t".join(map(str, inference)))
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

def __load_CisBP_models(files_dir="./files/"):

    # Initialize
    cisbp = {}

    for json_file in os.listdir(os.path.join(files_dir, "cisbp")):
        model = ReadSRModel(os.path.join(files_dir, "cisbp", json_file))
        cisbp.setdefault(CisBP2Pfam[model["Family_Name"]], model)

    return(cisbp)

def __format_BLAST_database(target_file):

    if not os.path.exists(f"{target_file}.phr"):

        # Make BLAST+ database
        cmd = "makeblastdb -in %s -dbtype prot" % target_file
        process = subprocess.run([cmd], shell=True, stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL)

def infer_SeqRecord_homologs(seq_record, target_file, pfam_alignments, cisbp,
    dummy_dir="/tmp/", files_dir="./files/", no_blast=False, n=5):

    # Initialize
    blast_results = {}
    inference_results = []

    # BLAST+ search
    if not no_blast:
        for r in blast(seq_record, target_file, n):
            blast_results.setdefault(r[1], r)

    for target_id in pfam_alignments:

        if not no_blast:
            if target_id not in blast_results:
                continue

        for DBD in pfam_alignments[seq_record.id]:

            if DBD not in pfam_alignments[target_id]:
                continue

            # Clean sequences
            seq1 = pfam_alignments[seq_record.id][DBD]
            seq1 = list([__remove_insertions(s) for s in seq1])
            seq2 = pfam_alignments[target_id][DBD]
            seq2 = list([__remove_insertions(s) for s in seq2])

            # Similarity regression alignment
            sr_alignment = {}
            ByPosPctID = __get_X(copy.copy(seq1), copy.copy(seq2), "identity")
            sr_alignment.setdefault("ByPos.PctID", ByPosPctID)
            ByPosAvgB62 = __get_X(copy.copy(seq1), copy.copy(seq2), "blosum62")
            sr_alignment.setdefault("ByPos.AvgB62", ByPosAvgB62)
            PctID_L = sum(ByPosPctID) / len(ByPosPctID)
            sr_alignment.setdefault("PctID_L", PctID_L)

            # Inference: Cis-BP
            if DBD in cisbp:
                model = cisbp[DBD]
            else:
                model = cisbp[None]
            _, Classification = ScoreAlignmentResult(sr_alignment, model)

            # i.e. inferred result
            if not Classification == "HSim":
                continue
            if not no_blast:
                inference_results.append([seq_record.id, target_id,
                    blast_results[target_id][4], blast_results[target_id][2],
                    blast_results[target_id][3],
                    round(sr_alignment["PctID_L"], 3)])
            else:
                inference_results.append([seq_record.id, target_id, "N/A",
                    "N/A", "N/A", round(sr_alignment["PctID_L"], 3)])
            break
    
    # Sort
    inference_results.sort(key=lambda x: (x[2], -x[-1]))

    return(inference_results)

def __make_seq_file(seq_record, file_name=".seq.fa"):

    # Remove seq file if exists...
    if os.path.exists(file_name):
        os.remove(file_name)

    # Write
    Jglobals.write(file_name, seq_record.format("fasta"))

def hmmscan(seq_file, hmm_file, dummy_dir="/tmp/",
    non_overlapping_domains=False):

    # Initialize
    out_file = os.path.join(dummy_dir, ".%s.out.txt" % os.getpid())

    # Scan
    cmd = "hmmscan --domtblout %s %s %s" % (out_file, hmm_file, seq_file)
    process = subprocess.run([cmd], shell=True, stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL)

    # Read domains
    domains = __read_domains(out_file)

    # Remove output file
    if os.path.exists(out_file):
        os.remove(out_file)

    # Filter overlapping domains
    if non_overlapping_domains:
        domains = __get_non_overlapping_domains(domains)

    # Yield domains one by one
    for pfam_ac, start, end, evalue in sorted(domains, key=lambda x: x[1]):

        yield(pfam_ac, start, end, evalue)

def __read_domains(file_name):
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

def __get_non_overlapping_domains(domains):
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

def hmmalign(seq_file, hmm_file):

    # Align
    cmd = "hmmalign --outformat PSIBLAST %s %s" % (hmm_file, seq_file)
    process = subprocess.check_output([cmd], shell=True, universal_newlines=True)

    return(__read_PSIBLAST_format(process))

def __read_PSIBLAST_format(psiblast_alignment):

    # Initialize
    alignment = ""

    # For each chunk...
    for chunk in psiblast_alignment.split("\n"):

        # If alignment substring...
        m = re.search("\s+(\S+)$", chunk)
        if m:
            alignment += m.group(1)

    return(alignment)

def blast(seq_record, target_file, n=5):

    # Initialize
    blast_results = set()
    outfmt = "sseqid pident length qstart qend sstart send evalue bitscore ppos qlen slen"

    # Run BLAST+
    cmd = "blastp -db %s -outfmt \"6 %s\"" % (target_file, outfmt)
    process = subprocess.Popen([cmd], shell=True, stdin=subprocess.PIPE,
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
        query_aligned_residues = int(blast_record[4]) - \
            int(blast_record[3]) + 1
        query_length = float(blast_record[10])
        target_aligned_residues = int(blast_record[6]) - \
            int(blast_record[5]) + 1
        target_length = float(blast_record[11])
        query_coverage = query_aligned_residues * 100 / query_length
        target_coverage = target_aligned_residues * 100 / target_length
        joint_coverage = math.sqrt(query_coverage * target_coverage)

        # Add BLAST+ record to search results
        blast_results.add((seq_record.id, target_id, query_start_end,
            target_start_end, e_value, score, percent_identities,
            alignment_length, percent_similarity, joint_coverage))

    # Return filtered results sorted by score
    blast_results = sorted(blast_results, key=lambda x: x[-1], reverse=True)

    return(__filter_blast_results_by_Rost(blast_results))

def __filter_blast_results_by_Rost(blast_results, n=5):

    # Initialize
    blast_homologs = []

    # For each result...
    for result in blast_results:

        # Initialize
        pid = result[6]
        L = result[7]

        # If homologs...
        if __is_alignment_over_Rost_seq_id_curve(pid, L, n):

            # Add homolog
            blast_homologs.append(result)

    return(blast_homologs)

def __is_alignment_over_Rost_seq_id_curve(pid, L, n=5):
    """
    This function returns whether an alignment is over the Rost's pairwise
    sequence identity curve or not.
    """
    return(pid >= __get_Rost_cutoff_percent_identity(L, n))

def __get_Rost_cutoff_percent_identity(L, n=5):
    """
    This function returns the Rost's cut-off percentage of identical residues
    for an alignment of length "L".
    """
    return(n + (480 * pow(L, -0.32 * (1 + pow(math.e, float(-L) / 1000)))))

# def __is_alignment_over_Rost_seq_sim_curve(psim, L, n=12):
#     """
#     This function returns whether an alignment is over the Rost's pairwise
#     sequence similarity curve or not.
#     """
#     return(psim >= _get_Rost_cutoff_percent_similarity(L, n))

# def __get_Rost_cutoff_percent_similarity(L, n=12):
#     """
#     This function returns the Rost's cut-off percentage of "similar" residues
#     for an alignment of length "L".
#     """
#     return(n + (420 * pow(L, -0.335 * (1 + pow(math.e, float(-L) / 2000)))))

def __get_CisBP_models(DBDs, cisbp):
    """
    Return Cis-BP highly-similar ("hsim") and dissimilar ("dis") thresholds
    based on DBD percentage of sequence identity ("pid") and similarity
    regression ("sr") models. For similarity regression models, it further
    returns the feature scaling "mean" and standard deviation ("sd"), and the
    features "intercept" and "weights".
    """

    # Initialize
    thresholds = {}

    # For each DBD...
    for DBD in DBDs:

        # Initialize
        thresholds.setdefault(DBD, {})

        # Get thresholds
        if DBD in cisbp:
            m = cisbp[DBD]
        else:
            m = cisbp[None]
        for what_threshold in ["dis", "hsim"]:
            t = m.get_model_threshold("pid", what_threshold)
            if t is None:
                t = cisbp[None].get_model_threshold("pid", what_threshold)
            thresholds[DBD].setdefault("pid", {})
            thresholds[DBD]["pid"].setdefault(what_threshold, t)
        thresholds[DBD].setdefault("sr", m.get_model("sr"))

    return(thresholds)

def __remove_insertions(s):
    """
    Remove insertions (i.e. lower case letters)
    """
    return(s.translate(str.maketrans("", "", string.ascii_lowercase)))

def __get_X(seq1, seq2, similarity="identity"):
    """
    Compare DBDs.
    """

    # Initialize
    scores = []

    # Reassign seq with more DBDs to seq1
    seq1, seq2 = __reassign(seq1, seq2)

    for i in range(len(seq1) - len(seq2) + 1):

        # Initialize
        arr = [0] * len(seq1[i])

        for j in range(len(seq2)):

            for k in range(len(seq1[i+j])):

                arr[k] += __score(seq1[i+j][k], seq2[j][k], similarity)

        # Append scores and rescale by the number of DBDs in seq2
        scores.append(np.array(arr))
        scores[-1] = scores[-1] / len(seq2)

    # Sort
    scores.sort(key=lambda x: sum(x), reverse=True)

    return(scores[0])

def __reassign(seq1, seq2):

    if len(seq1) < len(seq2):
        return(seq2, seq1)

    return(seq1, seq2)

def __score(aa1, aa2, similarity="identity"):

    if similarity == "identity":
        if aa1 == aa2 and aa1 != "-":
            return(1)
        else:
            return(0)
    elif similarity == "blosum62":
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