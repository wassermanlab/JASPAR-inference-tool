#!/usr/bin/env python

import argparse
from Bio.SubsMat.MatrixInfo import blosum62
import json
import math
import numpy as np
from numpy import log10 as log
import os
import shutil
import sys

# Defaults
out_dir = os.path.dirname(os.path.realpath(__file__))
root_dir = os.path.join(out_dir, os.pardir)
files_dir = os.path.join(root_dir, "files")

# Append JASPAR-profile-inference to path
sys.path.append(root_dir)

# Import globals
from __init__ import Jglobals
from infer_profile import _fetchXs, _removeLowercase

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

        # For each key, values...
        for key, values in groups.items():

            # Initialize
            Xss = {}
            Ys = []
            TFpairs = []

            # For each TF...
            for i in range(len(values) - 1):

                # For each next TF...
                for j in range(i + 1, len(values)):

                    # For each sequence similarity representation...
                    for similarity in ["identity", "blosum62"]:

                        # Initialize
                        Xs = []

                        # Innermost loop examining each DBD sequence...
                        for k in range(len(values[i][1])):

                            # Get Xs
                            seq1 = _removeLowercase(values[i][1][k])
                            seq2 = _removeLowercase(values[j][1][k])
                            Xs.extend(_fetchXs(seq1, seq2, similarity))

                        # Append Xs
                        Xss.setdefault(similarity, [])
                        Xss[similarity].append(Xs)

                    # Get Y and append
                    y = _fetchY(values[i][0], values[j][0])
                    Ys.append(y)

                    # Get TF pair
                    TFpairs.append((values[i][2], values[j][2]))

            # Add to pairwise
            pairwise.setdefault(key, [Xss, Ys, TFpairs])

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

def _fetchY(maIDlist1, maIDlist2):
    """
    Returns "0 - log10(Tomtom e-value)" as the dependent value (i.e. Y).
    """

    # Initialize
    y = -3 # i.e. ~size of JASPAR

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

                if 0 - log_evalue > y:
                    y = 0 - log_evalue

    return(y)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()