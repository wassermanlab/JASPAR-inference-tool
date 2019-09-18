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

    parser.add_argument("-f", default=files_dir, help="files directory (from get_files.py; default=../files/)", metavar="DIR")
    parser.add_argument("-o", default=out_dir, help="output directory (default = ./)", metavar="DIR")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Pairwise
    pairwise(os.path.abspath(args.f), os.path.abspath(args.o))

def pairwise(files_dir=files_dir, out_dir=out_dir):

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
        clusters = _get_clusters(files_dir)

        # Get domains
        groups = _get_groups(files_dir)

        # For each key, values...
        for key, values in groups.items():
            # From PMID:1692833;
            # DNA-binding domains vary in length. Among the shortest the AT
            # hook, which recognizes sequences in the minor groove utilizing
            # fewer than a dozen amino acids.

            # Initialize
            Xss = {}
            Ys = []
            uniaccs = []

            # Skip if not enough TFs
            if len(values) < 3:
                continue

            # For each TF...
            for i in range(len(values) - 1):

                # For each next TF...
                for j in range(i + 1, len(values)):

                    # For each sequence similarity representation...
                    for similarity in ["identity", "blosum62"]:

                        # Initialize
                        Xs = []

                        # Inner most loop for examining EACH different component...
                        for k in range(len(values[i][1])):

                            # Get X
                            seq1 = _removeLowercase(values[i][1][k])
                            seq2 = _removeLowercase(values[j][1][k])
                            Xs.extend(_fetchXs(seq1, seq2, similarity))

                        # Append Xs
                        Xss.setdefault(similarity, [])
                        Xss[similarity].append(Xs)

                    # Get Y
                    y = _fetchY(values[i][0], values[j][0])

                    # Append y, uniaccs
                    Ys.append(y)
                    uniaccs.append("{}*{}".format(values[i][2], values[j][2]))

            # Skip if only one class in the data
            if len(set(Ys)) == 1:
                continue

            # Add to pairwise
            pairwise.setdefault(key, [Xss, Ys, uniaccs])

        # Write
        Jglobals.write(
            pairwise_json_file,
            json.dumps(pairwise, sort_keys=True, indent=4, separators=(",", ": "))
        )

def _get_clusters(files_dir=files_dir):

    # Load JSON file
    clusters_json_file = os.path.join(files_dir, "clusters.json")
    with open(clusters_json_file) as f:
        clusters = json.load(f)

    return(clusters)

def _get_groups(files_dir=files_dir):

    # Load JSON file
    groups_json_file = os.path.join(files_dir, "groups.json")
    with open(groups_json_file) as f:
        groups = json.load(f)

    return(groups)

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