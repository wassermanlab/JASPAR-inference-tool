#!/usr/bin/env python

import argparse
from Bio.SubsMat.MatrixInfo import blosum62
import json
import os
import string
import sys

# Defaults
out_dir = os.path.dirname(os.path.realpath(__file__))
root_dir = os.path.join(out_dir, os.pardir, os.pardir)
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
    This function parses arguments provided via the command line and returns an {argparse} object.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-c", choices=["rsat", "tomtom"], default="tomtom", help="cluster profiles using \"rsat\" matrix-clustering or \"tomtom\" (i.e. default)")
    parser.add_argument("-f", default=files_dir, help="files directory (from get_files.py; default=../files/)", metavar="DIR")
    parser.add_argument("-o", default=out_dir, help="output directory (default = ./)", metavar="DIR")
    parser.add_argument("-r", choices=["id", "sim"], default="id", help="logistic regression based on sequence identity (i.e. \"id\"; default) or similarity (i.e. \"sim\")")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Pairwise
    pairwise(args.c, os.path.abspath(args.f), os.path.abspath(args.o), args.r)

def pairwise(cluster="tomtom", files_dir=files_dir, out_dir=out_dir, regression="id"):

    # Skip if pairwise JSON file already exists
    pairwise_json_file = os.path.join(out_dir, "pairwise.json")
    if not os.path.exists(pairwise_json_file):

        # Initialize
        pairwise = {}

        # Create output dir
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # Get clusters
        global clusters
        clusters = _get_clusters(cluster, files_dir)

        # Get domains
        groups = _get_groups(files_dir)

        # # For each key, values...
        for key, values in groups.items():
            # Skip if not enough values:
            # From PMID:1692833;
            # DNA-binding domains vary in length. Among the shortest the AT hook,
            # which recognizes sequences in the minor groove utilizing fewer than
            # a dozen amino acids.
            # 
            # 10 values gives you a training set of 45 
            # >>> from itertools import combinations
            # >>> digits = "".join([str(i) for i in range(10)])
            # >>> len(list(combinations(digits, 2)))
            # 45

            # Initialize
            Xs = []
            Ys = []
            uniaccs = []

            # Skip if not enough TFs
            if len(values) < 3:
                continue

            # For each TF...
            for i in range(len(values) - 1):

                # For next TF...
                for j in range(i + 1, len(values)):

                    # clean lists
                    x = []

                    # Inner most loop for examining EACH different component...
                    for k in range(len(values[i][1])):

                        # Get X
                        seq1 = _removeLowercase(values[i][1][k])
                        seq2 = _removeLowercase(values[j][1][k])
                        x.extend(_fetchXs(seq1, seq2, regression))

                    # Get Y
                    y = _fetchY(values[i][0], values[j][0])

                    # Get uniacc
                    uniacc = "{}*{}".format(values[i][2], values[j][2])

                    # Append Xs, Ys, uniaccs
                    Xs.append(x)
                    Ys.append(y)
                    uniaccs.append(uniacc)

            # Add to pairwise
            pairwise.setdefault(key, [Xs, Ys, uniaccs])

        # Write
        Jglobals.write(
            pairwise_json_file,
            json.dumps(pairwise, sort_keys=True, indent=4, separators=(",", ": "))
        )

def _get_clusters(cluster="tomtom", files_dir=files_dir):

    # Load JSON file
    clusters_json_file = os.path.join(files_dir, "clusters.%s.json" % cluster)
    with open(clusters_json_file) as f:
        clusters = json.load(f)

    return(clusters)

def _get_groups(files_dir=files_dir):

    # Load JSON file
    groups_json_file = os.path.join(files_dir, "groups.json")
    with open(groups_json_file) as f:
        groups = json.load(f)

    return(groups)

def _removeLowercase(strings):

    return(strings.translate(str.maketrans("", "", string.ascii_lowercase)))

def _fetchXs(seq1 , seq2, regression="id"):
    """
    Called for each comparison to compare the string.
    """

    # Fill Xs with zeroes
    scores = [0] * len(seq1)

    # Strings should have same length
    if len(seq1) == len(seq2):
        for n in range(len(seq1)):
            if regression == "id":
                scores[n] = _IDscoring(seq1[n], seq2[n])
            else:
                scores[n] = _BLOSUMScoring(seq1[n], seq2[n])
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

def _fetchY(maIDlist1, maIDlist2):
    """
    get the list of maID from 1 and 2 to see if any of them is correlated.
    """

    # Iterate through profiles...
    for maID1 in maIDlist1:

        # Skip matrix ID
        if maID1 not in clusters:
            continue

        for maID2 in maIDlist2:

            # Skip matrix ID
            if maID2 not in clusters:
                continue

            # If profiles were clustered together...
            if maID2 in clusters[maID1] and maID1 in clusters[maID2]:
                return(True)

    return(False)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()