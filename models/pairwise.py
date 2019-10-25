#!/usr/bin/env python

import argparse
from Bio.SubsMat.MatrixInfo import blosum62
import json
import math
import numpy as np
from numpy import log10 as log
import os
import shutil
import string
import sys

# Defaults
out_dir = os.path.dirname(os.path.realpath(__file__))
root_dir = os.path.join(out_dir, os.pardir)
files_dir = os.path.join(root_dir, "files")

# Append JASPAR-profile-inference to path
sys.path.append(root_dir)

# Import globals
from __init__ import Jglobals
from infer_profile import (
    _is_alignment_over_Rost_seq_id_curve,
    _is_alignment_over_Rost_seq_sim_curve
)

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns
    an {argparse} object.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-f", default=files_dir, metavar="DIR",
        help="files directory (from get_files.py; default=../files/)")
    parser.add_argument("-o", default=out_dir, metavar="DIR",
        help="output directory (default = ./)")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Pairwise
    pairwise(os.path.abspath(args.f), os.path.abspath(args.o))

def pairwise(files_dir=files_dir, out_dir=out_dir):

    # Skip if pairwise JSON file already exists
    gzip_file = os.path.join(out_dir, "pairwise.json.gz")
    if not os.path.exists(gzip_file):

        # Initialize
        pairwise = {}

        # Create output dir
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # Get groups by DBD composition
        groups = _get_DBD_groups(files_dir)

        # Get Tomtom groups
        global tomtom
        tomtom = _get_Tomtom_groups(files_dir)

        # Get BLAST+ groups
        global blast
        blast = _get_BLAST_groups(files_dir)

        # For each key, values...
        for key, values in groups.items():

            if key != "zf-C2H2":
                continue

            # Initialize
            Xss = {}
            BLASTXss = {}
            Ys = []
            uniaccs = []

            # For each TF...
            for i in range(len(values) - 1):

                # For each next TF...
                for j in range(i + 1, len(values)):

                    if values[j][2] != "P49711":
                        continue

                    # For each sequence similarity representation...
                    for similarity in ["identity", "blosum62"]:

                        # # Initialize
                        # Xs = []

                        # # Inner most loop for examining EACH different component...
                        # for k in range(len(values[i][1])):

                        #     # Get Xs
                        #     seq1 = _removeLowercase(values[i][1][k])
                        #     seq2 = _removeLowercase(values[j][1][k])
                        #     Xs.extend(_fetchXs(seq1, seq2, similarity))

                        # Get Xs
                        Xs = _fetchXs(values[i][1], values[j][1], similarity)

                        # Get BLAST+ Xs
                        BLASTXs = _fetchBLASTXs(values[i][2], values[j][2], similarity)

                        # Append Xs
                        Xss.setdefault(similarity, [])
                        Xss[similarity].append(Xs)
                        BLASTXss.setdefault(similarity, [])
                        BLASTXss[similarity].append(BLASTXs)

                    # Get Y
                    y = _fetchY(values[i][0], values[j][0])

                    # Append, y, uniaccs
                    Ys.append(y)
                    uniaccs.append((values[i][2], values[j][2]))

            # Add to pairwise
            pairwise.setdefault(key, [Xss, BLASTXss, Ys, uniaccs])

        exit(0)

        # Write
        Jglobals.write(
            gzip_file[:-3],
            json.dumps(pairwise, sort_keys=True, indent=4, separators=(",", ": "))
        )
        fi = Jglobals._get_file_handle(gzip_file[:-3], "rb")
        fo = Jglobals._get_file_handle(gzip_file, "wb")
        shutil.copyfileobj(fi, fo)
        fi.close()
        fo.close()
        os.remove(gzip_file[:-3])

def _get_DBD_groups(files_dir=files_dir):

    # Load JSON file
    json_file = os.path.join(files_dir, "groups.DBDs.json.gz")
    handle = Jglobals._get_file_handle(json_file)
    groups = json.load(handle)
    handle.close()

    return(groups)

def _get_Tomtom_groups(files_dir=files_dir):

    # Initialize
    tomtom_filtered = {}

    # Load JSON file
    json_file = os.path.join(files_dir, "groups.tomtom.json.gz")
    handle = Jglobals._get_file_handle(json_file)
    tomtom_unfiltered = json.load(handle)
    handle.close()

    for matrix_id in tomtom_unfiltered:

        for m, e in tomtom_unfiltered[matrix_id]:

            tomtom_filtered.setdefault(matrix_id, {})
            tomtom_filtered[matrix_id].setdefault(m, e)

    return(tomtom_filtered)

def _get_BLAST_groups(files_dir=files_dir):

    # Initialize
    blast_filtered = {}

    # Load JSON file
    json_file = os.path.join(files_dir, "groups.blast.json.gz")
    handle = Jglobals._get_file_handle(json_file)
    blast_unfiltered = json.load(handle)
    handle.close()

    for uniacc in blast_unfiltered:

        for t, qse, tse, e, s, pid, al, psim, jc in blast_unfiltered[uniacc]:

            # (t) identifier of target sequence;
            # (qse, tse) start and end-position in query and in target;
            # (e) E-value;
            # (s) bit score;
            # (pid) percentage of identical matches;
            # (al) alignment length;
            # (psim) percentage of positive-scoring matches; and
            # (jc) joint coverage (i.e. square root of the coverage
            # on the query and the target).
            # s = s / al # transformation of bit score

            # Get Rost's verdict
            rost_seq_id = _is_alignment_over_Rost_seq_id_curve(pid, al)
            rost_seq_sim = _is_alignment_over_Rost_seq_sim_curve(psim, al)

            # Reformat percentages & coverage
            pid = pid / 100.0
            psim = psim / 100.0
            jc = jc / 100.0

            blast_filtered.setdefault(uniacc, {})
            blast_filtered[uniacc].setdefault(t, [[pid, rost_seq_id], [psim, rost_seq_sim], [jc]])

    return(blast_filtered)

def _removeLowercase(s):

    return(s.translate(str.maketrans("", "", string.ascii_lowercase)))

def _fetchXs(seqs1, seqs2, similarity="identity"):

    # Initialize
    means = []
    scores = []
    A = seqs1
    B = seqs2

    if len(A) < len(B):
        A = seqs2
        B = seqs1

    for a in range(len(A) - len(B) + 1):

        scores.append([])

        for b in range(len(B)):

            seqA = _removeLowercase(A[a+b])
            seqB = _removeLowercase(B[b])

            if similarity == "identity":
                scores[-1].append([0] * len(seqA))
            elif similarity == "blosum62":
                scores[-1].append([-4] * len(seqA))

            for n in range(len(seqA)):

                if similarity == "identity":
                    scores[-1][-1][n] = _IDscoring(seqA[n], seqB[n])
                elif similarity == "blosum62":
                    scores[-1][-1][n] = _BLOSUMscoring(seqA[n], seqB[n])

    # From Lambert et al.
    # For a multi-DBD alignment, the feature vector is generated by the
    # average score (identity or similarity) in each position of the DBD
    # alignment from all DBD arrays, normalizing by the DBD length of the
    # longest protein (Supplementary Fig. 2b).
    for s in scores:
        means.append([float(sum(z))/len(A) for z in zip(*s)])

    # Sort 
    for m in sorted(means, key=lambda x: sum(x), reverse=True):
        print(m, sum(m))
        # for seq2 in seqs2:

        #     seq2 = _removeLowercase(seq2)

        #     # Strings should have same length
        #     if len(seq1) == len(seq2):

        #         scores = [0] * len(seq1)

        #         for s in range(len(scores)):

        #             if similarity == "identity":
        #                 scores[s] = _IDscoring(seq1[s], seq2[s])

        #             elif similarity == "blosum62":
        #                 scores[s] = _BLOSUMscoring(seq1[s], seq2[s])

        #         pairs.append((scores, sum([seq1.count("-"), seq2.count("-")])))
    exit(0)

    # # Inner most loop for examining EACH different component...
    # for k in range(len(values[i][1])):

    #     # Get Xs
    #     seq1 = _removeLowercase(values[i][1][k])
    #     seq2 = _removeLowercase(values[j][1][k])
    #     Xs.extend(_fetchXs(seq1, seq2, similarity))
    pass

# def _fetchXs(seq1 , seq2, similarity="identity"):
#     """
#     Called for each comparison to compare the strings.
#     """

#     # Fill Xs with zeroes
#     scores = [0] * len(seq1)

#     # Strings should have same length
#     if len(seq1) == len(seq2):

#         for n in range(len(seq1)):

#             if similarity == "identity":
#                 scores[n] = _IDscoring(seq1[n], seq2[n])

#             elif similarity == "blosum62":
#                 scores[n] = _BLOSUMscoring(seq1[n], seq2[n])

#     else:
#         # there is something wrong
#         print("Strings have different length!\n\tA: %s\n\tB: %s" % (seq1, seq2))
#         exit(0)

#     return(scores)

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

def _fetchBLASTXs(uacc1, uacc2, similarity="identity"):

    # Initialize
    BLASTXs = [0.0, False, 0.0]

    if uacc1 in blast:
        if uacc2 in blast[uacc1]:
            if similarity == "identity":
                BLASTXs = blast[uacc1][uacc2][0] + blast[uacc1][uacc2][2]
            elif similarity == "blosum62":
                BLASTXs = blast[uacc1][uacc2][1] + blast[uacc1][uacc2][2]

    return(BLASTXs)

def _fetchY(maIDlist1, maIDlist2):
    """
    Returns "0 - log10(Tomtom e-value)" as the dependent value (i.e. Y).
    """

    # Initialize
    y = 0 - log(1963.0)

    # Iterate through profiles...
    for maID1 in maIDlist1:

        # Skip matrix ID
        if maID1 not in tomtom:
            continue

        for maID2 in maIDlist2:

            # Skip matrix ID
            if maID2 not in tomtom:
                continue

            # If profiles were clustered together...
            if maID2 in tomtom[maID1] and maID1 in tomtom[maID2]:

                joint_evalue = math.sqrt(tomtom[maID1][maID2] * tomtom[maID2][maID1])
                log_evalue = log(joint_evalue)

                # return(max(min(1.0, 0 - log_evalue / 10), 0.0))

                if 0 - log_evalue > y:
                    y = 0 - log_evalue

    return(y)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()