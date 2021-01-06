#!/usr/bin/env python

import argparse
from Bio.SubsMat.MatrixInfo import blosum62
import copy
import json
import math
import numpy as np
from numpy import log10 as log
import os
import shutil
import sys

# Defaults
out_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "pairs")
root_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir)
groups_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "groups")

# Append JASPAR-profile-inference to path
sys.path.append(root_dir)

# Import globals
from __init__ import Jglobals
from infer_profile import __fetchXs, __removeLowercase

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns
    an {argparse} object.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("--groups-dir", default=groups_dir, metavar="DIR",
        help="files directory (from get_groups.py)")
    parser.add_argument("-o", default=out_dir, metavar="DIR",
        help="output directory (default = ./)")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Get pairs
    get_pairs(os.path.abspath(args.groups_dir), os.path.abspath(args.o))

def get_pairs(groups_dir=groups_dir, out_dir=out_dir):

    # Globals
    global cwd
    cwd = os.getcwd()

    # Create output dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Get pairs
    __pair(groups_dir, out_dir)

def __pair(groups_dir=groups_dir, out_dir=out_dir):

    # Skip if pairs JSON file already exists
    gzip_file = os.path.join(out_dir, "pairs.json.gz")
    if not os.path.exists(gzip_file):

        # Initialize
        pairs = {}

        # Move to output directory
        os.chdir(out_dir)

        # Load JSON files
        json_file = os.path.join(groups_dir, "groups.DBDs.json.gz")
        handle = Jglobals._get_file_handle(json_file)
        DBDs = json.load(handle)
        handle.close()
        json_file = os.path.join(groups_dir, "groups.clusters.json.gz")
        handle = Jglobals._get_file_handle(json_file)
        clusters = json.load(handle)
        handle.close()

        print(clusters)
        exit(0)

        # For each DBD, values...
        for DBD, values in DBDs.items():

            # Test C2H2 zinc fingers
            if DBD != "zf-C2H2":
                continue

            # Initialize
            features = {}
            labels = []
            pairings = []
            uniaccs = list(values.keys())

            # For each UniProt Accession...
            for i in range(len(uniaccs) - 1):

                # Initialize
                seqi = []

                # Innermost loop examining each DBD sequence...
                for k in range(len(values[uniaccs[i]])):
                    seqi.append(__removeLowercase(values[uniaccs[i]][k]))

                # For each other UniProt Accession...
                for j in range(i + 1, len(uniaccs)):

                    # Initialize
                    seqj = []

                    # Innermost loop examining each DBD sequence...
                    for k in range(len(values[uniaccs[j]])):
                        seqj.append(__removeLowercase(values[uniaccs[j]][k]))

                    # For each sequence similarity representation...
                    for similarity in ["identity", "blosum62"]:

                        # Get Xs
                        Xs = __fetchXs(seqi, seqj, similarity)
                        features.setdefault(similarity, [])
                        features[similarity].append(Xs)

                    # Get Y and append
                    y = __fetchY(uniaccs[i], uniaccs[j], clusters)
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

def _get_DBD_groups(groups_dir=groups_dir):

    # Load JSON file
    json_file = os.path.join(files_dir, "groups.DBDs.json.gz")
    handle = Jglobals._get_file_handle(json_file)
    groups = json.load(handle)
    handle.close()

    return(groups)

# def _get_Tomtom_groups(files_dir=files_dir):

#     # Initialize
#     tomtom_filtered = {}

#     # Load JSON file
#     json_file = os.path.join(files_dir, "groups.tomtom.json.gz")
#     handle = Jglobals._get_file_handle(json_file)
#     tomtom_unfiltered = json.load(handle)
#     handle.close()

#     for matrix_id in tomtom_unfiltered:

#         for m, e in tomtom_unfiltered[matrix_id]:

#             tomtom_filtered.setdefault(matrix_id, {})
#             tomtom_filtered[matrix_id].setdefault(m, e)

#     return(tomtom_filtered)

def __fetchY(uniacc1, uniacc2, clusters):
    """
    Returns whether or not the matrices of two TFs are clustered together as
    the dependent value (i.e. Y).
    """
    print(clusters[uniacc2])
    exit(0)

    return(len(a.intersection(b)) > 0)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()