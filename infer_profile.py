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
from __init__ import Jglobals, CisBP

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
    parser.add_argument("sequences",
        help="input sequence(s) in FASTA format")

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
    inference_group.add_argument("-l", "--latest", action="store_true",
        help="return the latest version of each profile")
    inference_group.add_argument("--rost", default=5, metavar="INT",
        help="\"n\" parameter for the Rost's curve (default = 5)")
    inference_group.add_argument("--taxon", nargs="*", default=Jglobals.taxons,
        metavar="STR", help="return profiles from given taxon (default = all)")

    args = parser.parse_args()

    return(args)

def main():

    # Parse arguments
    args = parse_args()

    # Warnings
    if not args.warnings:
        warnings.filterwarnings("ignore")

    # Infer profiles
    infer_profiles(args.sequences, args.dummy_dir, args.files_dir,
        args.output_file, args.threads, args.latest, args.rost, args.taxon)

def infer_profiles(fasta_file, dummy_dir="/tmp/", files_dir=files_dir,
    output_file=None, threads=1, latest=False, n=5, taxons=Jglobals.taxons):

    # Initialize
    base_name = os.path.basename(__file__)
    pid = os.getpid()

    # Load data
    cisbp = __load_CisBP_models(files_dir)
    jaspar = __load_JASPAR_files(files_dir, taxons)

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

    # Write
    columns = ["Query", "TF Name", "TF Matrix", "E-value", "Query Start-End",
        "TF Start-End", "DBD %ID", "Similarity Regression"]
    Jglobals.write(dummy_file, "\t".join(columns))

    # Infer SeqRecord profiles
    kwargs = {"total": len(seq_records), bar_format: bar_format}
    pool = Pool(min([threads, len(seq_records)]))
    p = partial(infer_SeqRecord_profiles, cisbp=cisbp, dummy_dir=dummy_dir,
        files_dir=files_dir, jaspar=jaspar, latest=latest, n=n, taxons=taxons)
    for inferences in tqdm(pool.imap(p, seq_records), **kwargs):
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
        model = CisBP(os.path.join(files_dir, "cisbp", json_file))
        cisbp.setdefault(model.family, model)

    return(cisbp)

def __load_JASPAR_files(files_dir="./files/", taxons=["fungi", "insects",
    "nematodes", "plants", "vertebrates"]):

    # Initialize
    jaspar = {}
    pfams = {}
    profiles = {}
    uniprots = {}

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

    return(jaspar)

def infer_SeqRecord_profiles(seq_record, cisbp, jaspar, dummy_dir="/tmp/",
    files_dir="./files/", latest=False, n=5, taxons=["fungi", "insects",
    "nematodes", "plants", "vertebrates"]):

    # Initialize
    pfam_alignments = []
    inference_results = []

    # Get SeqRecord Pfam DBDs
    alignments = __get_SeqRecord_Pfam_alignments(seq_record,
        files_dir, dummy_dir)
    if len(alignments) == 0:
        return(inference_results)
    pfam_alignments.append({})
    for alignment in alignments:
        pfam_alignments[0].setdefault(alignment[0], [])
        pfam_alignments[0][alignment[0]].append(alignment[1])

    # BLAST+ search
    blast_results = blast(seq_record, files_dir, taxons, n)

    # Get Pfam DBDs of BLAST+ (filtered) results
    pfam_alignments.append(__get_blast_results_Pfam_alignments(blast_results,
        jaspar))

    # Get Cis-BP thresholds
    models = __get_CisBP_models(list(pfam_alignments[0].keys()), cisbp)

    for result in blast_results:
        scores = [None, None]
        for DBD in pfam_alignments[0]:

            if DBD not in pfam_alignments[1][result[1]]:
                continue

            # Clean sequences
            seq1 = pfam_alignments[0][DBD]
            seq1 = list([__remove_insertions(s) for s in seq1])
            seq2 = pfam_alignments[1][result[1]][DBD]
            seq2 = list([__remove_insertions(s) for s in seq2])

            # Inference: percentage of sequence identity
            X = __get_X(copy.copy(seq1), copy.copy(seq2), "identity")
            score = sum(X) / len(X)
            if score >= models[DBD]["pid"]["hsim"]:
                scores[0] = round(score, 3)

            # Inference: similarity regression
            if models[DBD]["sr"] is not None:
                if models[DBD]["sr"]["similarity"] == "blosum62":
                    X = __get_X(copy.copy(seq1), copy.copy(seq2), "blosum62")
                mean = models[DBD]["sr"]["mean"]
                sd = models[DBD]["sr"]["sd"]
                intercept = models[DBD]["sr"]["intercept"]
                weights = models[DBD]["sr"]["weights"]
                # i.e. length of Pfam domains might be different
                if len(X) == len(weights):
                    score = __get_similarity_regression_score(X, mean, sd,
                        intercept, weights)
                    if score >= models[DBD]["sr"]["hsim"]:
                        scores[1] = round(score, 3)

            # Skip
            if scores[0] is None and scores[1] is None:
                continue

            # Add inferred result
            for matrix, gene_name in jaspar[result[1]]["profiles"]:
                inference_results.append([result[0], gene_name, matrix,
                    result[4], result[2], result[3], scores[0], scores[1]])
            break
    
    # If use the lastest version of JASPAR...
    if latest:
        # Sort
        inference_results.sort(key=lambda x: (x[3], x[1], -float(x[2][2:])))
        # Remove profiles from older versions
        for i in sorted(frozenset(range(len(inference_results))), reverse=True):
            if inference_results[i][2][:6] == inference_results[i - 1][2][:6]:
                inference_results.pop(i)
    # ... Else, just sort...
    else:
        inference_results.sort(key=lambda x: (x[3], x[1], float(x[2][2:])))

    return(inference_results)

def __get_SeqRecord_Pfam_alignments(seq_record, files_dir="./files/",
    dummy_dir="/tmp/"):

    # Initialize
    pfam_alignments = []
    hmm_db = os.path.join(files_dir, "pfam", "All.hmm")

    # Make seq file
    seq_file = os.path.join(dummy_dir, ".%s.seq.fasta" % os.getpid())
    __make_seq_file(seq_record, seq_file)

    # For each DBD...
    for pfam_id_std, start, end, evalue in hmmscan(seq_file, hmm_db, dummy_dir,
        non_overlapping_domains=True):

        # Initialize
        hmm_file = os.path.join(files_dir, "pfam", "%s.hmm" % pfam_id_std)

        # Make seq file
        sub_seq_record = SeqRecord(seq_record.seq[start:end], id=seq_record.id,
            name=seq_record.name, description=seq_record.description)
        __make_seq_file(sub_seq_record, seq_file)

        # Add DBDs
        alignment = hmmalign(seq_file, hmm_file)
        pfam_alignments.append((pfam_id_std, alignment, start+1, end, evalue))

    return(pfam_alignments)

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

def blast(seq_record, files_dir="./files/", taxons=["fungi", "insects",
    "nematodes", "plants", "vertebrates"], n=5):

    # Initialize
    blast_results = set()
    outfmt = "sseqid pident length qstart qend sstart send evalue bitscore ppos qlen slen"

    # For each taxon...
    for taxon in taxons:

        # Taxon db
        taxon_db = os.path.join(files_dir, "%s.fa" % taxon)

        # Run BLAST+
        cmd = "blastp -db %s -outfmt \"6 %s\"" % (taxon_db, outfmt)
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

def __get_blast_results_Pfam_alignments(blast_results, jaspar):

    # Initialize
    pfam_alignments = {}

    # For each BLAST result...
    for blast_result in blast_results:
        pfam_alignments.setdefault(blast_result[1], {})
        for alignment in jaspar[blast_result[1]]["pfam"]:
            pfam_alignments[blast_result[1]].setdefault(alignment[0], [])
            pfam_alignments[blast_result[1]][alignment[0]].append(alignment[1])

    return(pfam_alignments)

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

def __get_similarity_regression_score(X, mean, sd, intercept, weights):
    """
    From "Evaluate Heldout Predictions.ipynb" (i.e. as in Lambert et al.)
    """

    # Z-scoring
    Z = (X - mean) / sd
    Z[np.isnan(Z)] = 0

    return(intercept + np.dot(weights, Z))

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()
