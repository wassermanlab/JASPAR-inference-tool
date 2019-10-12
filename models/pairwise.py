#!/usr/bin/env python

import argparse
from Bio.SubsMat.MatrixInfo import blosum62
from collections import Counter
import json
import math
from numpy import log10 as log
import os
import pickle
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

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns
    an {argparse} object.
    """

    parser = argparse.ArgumentParser()

    # parser.add_argument("-e", default=0.05, type=float, metavar="FLOAT",
    #     help="e-value threshold (default = 0.05)")
    parser.add_argument("-f", default=files_dir, metavar="DIR",
        help="files directory (from get_files.py; default=../files/)")
    parser.add_argument("-o", default=out_dir, metavar="DIR",
        help="output directory (default = ./)")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Pairwise
    # pairwise(args.e, os.path.abspath(args.f), os.path.abspath(args.o))
    pairwise(os.path.abspath(args.f), os.path.abspath(args.o))

# def pairwise(evalue=0.05, files_dir=files_dir, out_dir=out_dir):
def pairwise(files_dir=files_dir, out_dir=out_dir):

    # Skip if pickle file already exists
    pickle_file = os.path.join(out_dir, "pairwise.pickle")
    if not os.path.exists(pickle_file):

        # Initialize
        pairwise = {}

        # Create output dir
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # Get groups by DBD composition
        groups = _get_DBD_groups(files_dir)

        # Get Tomtom groups
        global tomtom
        # tomtom = _get_Tomtom_groups(evalue, files_dir)
        tomtom = _get_Tomtom_groups(files_dir)

        # Get BLAST+ groups
        global blast
        blast = _get_BLAST_groups(files_dir)

        # For each key, values...
        for key, values in groups.items():

            # Initialize
            Xss = {}
            BLASTXss = {}
            Ys = []
            uniaccs = []

            # # Skip if not enough TFs
            # if len(values) < 3:
            #     continue

            # For each TF...
            for i in range(len(values) - 1):

                # For each next TF...
                for j in range(i + 1, len(values)):

                    # For each sequence similarity representation...
                    for similarity in ["identity", "blosum62"]:

                        # Initialize
                        Xs = []
                        BLASTXs = []

                        # Inner most loop for examining EACH different component...
                        for k in range(len(values[i][1])):

                            # Get Xs
                            seq1 = _removeLowercase(values[i][1][k])
                            seq2 = _removeLowercase(values[j][1][k])
                            Xs.extend(_fetchXs(seq1, seq2, similarity))

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

            # # Skip if not enough classes in the data
            # if len(set([tuple(y) for y in Ys])) < 3:
            #     continue

            # Add to pairwise
            # pairwise.setdefault(key, [Xss, Ys, uniaccs])
            pairwise.setdefault(key, [Xss, BLASTXss, Ys, uniaccs])

        # Write pickle file
        with open(pickle_file, "wb") as f:
            pickle.dump(pairwise, f)

def _get_DBD_groups(files_dir=files_dir):

    # Load JSON file
    groups_json_file = os.path.join(files_dir, "groups.DBDs.json")
    with open(groups_json_file) as f:
        groups = json.load(f)

    return(groups)

# def _get_Tomtom_groups(evalue=0.05, files_dir=files_dir):
def _get_Tomtom_groups(files_dir=files_dir):

    # Initialize
    tomtom_filtered = {}

    # Load JSON file
    groups_json_file = os.path.join(files_dir, "groups.tomtom.json")
    with open(groups_json_file) as f:
        tomtom_unfiltered = json.load(f)

    for matrix_id in tomtom_unfiltered:

        # # Initialize
        # matrix_ids = set()

        for m, e in tomtom_unfiltered[matrix_id]:

        #     if e <= evalue:
        #         matrix_ids.add(m)

        # if len(matrix_ids) > 0:
        #     tomtom_filtered.setdefault(matrix_id, matrix_ids)
            tomtom_filtered.setdefault(matrix_id, {})
            tomtom_filtered[matrix_id].setdefault(m, e)

    return(tomtom_filtered)

def _get_BLAST_groups(files_dir=files_dir):

    # Initialize
    blast_filtered = {}

    # Load JSON file
    groups_json_file = os.path.join(files_dir, "groups.blast.json")
    with open(groups_json_file) as f:
        blast_unfiltered = json.load(f)

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
            pid = pid / 100.0
            psim = psim / 100.0
            jc = jc / 100.0

            blast_filtered.setdefault(uniacc, {})
            # blast_filtered[uniacc].setdefault(t, [s, pid, psim, jc])
            blast_filtered[uniacc].setdefault(t, [pid, psim, jc])

    return(blast_filtered)

def _removeLowercase(s):

    return(s.translate(str.maketrans("", "", string.ascii_lowercase)))

def _fetchXs(seq1 , seq2, similarity="identity"):
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

def _fetchBLASTXs(uacc1, uacc2, similarity="identity"):

    # Initialize
    BLASTXs = [0.0, 0.0]

    if uacc1 in blast:
        if uacc2 in blast[uacc1]:
            if similarity == "identity":
                BLASTXs[0] = blast[uacc1][uacc2][0]
                BLASTXs[1] = blast[uacc1][uacc2][2]
            elif similarity == "blosum62":
                BLASTXs[0] = blast[uacc1][uacc2][1]
                BLASTXs[1] = blast[uacc1][uacc2][2]

    return(BLASTXs)

def _fetchY(maIDlist1, maIDlist2):
    """
    get the list of maID from 1 and 2 to see if any of them is correlated.
    """

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
                return(0 - log_evalue)

    return(-2)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()