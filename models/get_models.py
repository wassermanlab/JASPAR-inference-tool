#!/usr/bin/env python

import argparse
import json
import os
import sys

# Defaults
out_dir = os.path.dirname(os.path.realpath(__file__))
files_dir = os.path.join(out_dir, os.pardir, "files")

# Append JASPAR-profile-inference to path
sys.path.append(os.path.join(out_dir, os.pardir))

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

    # Get models
    get_models(os.path.abspath(args.f), os.path.abspath(args.o))

def get_models(files_dir=files_dir, out_dir=out_dir):

    # Create output dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    ##
    ## Step 1: group alignments by DBD composition
    ##
    _group_alignments(files_dir, out_dir)

    ##
    ## Step 2: get Y's for logistic regression
    ##
    _get_Ys(files_dir, out_dir)  

    ##
    ## Step 2: get one matrix per DBD composition
    ##

def _group_alignments(files_dir=files_dir, out_dir=out_dir):

    # Skip if alignments JSON file already exists
    alignments_json_file = os.path.join(out_dir, "alignments.json")
    if not os.path.exists(alignments_json_file):

        # Initialize
        alignments = {}

        # For each taxon...
        for taxon in Jglobals.taxons:

            # Ininitialize
            pfam_file_ext = ".pfam.json"
            pfam_json_file = os.path.join(files_dir, taxon + pfam_file_ext)

            # Load JSON file
            with open(pfam_json_file) as f:
                pfams = json.load(f)

            # For each uniacc...
            for uniacc, values in pfams.items():

                # Unwind
                domains = ";".join([v[0] for v in values])
                sequences = "".join([v[1] for v in values])

                # Add alignments
                alignments.setdefault(domains, [])
                alignments[domains].append([uniacc, sequences])

        # Write
        Jglobals.write(
            alignments_json_file,
            json.dumps(alignments, sort_keys=True, indent=4, separators=(",", ": "))
        )

def _get_Ys():

    # Skip if Y's JSON file already exists
    ys_json_file = os.path.join(out_dir, "y.json")
    if not os.path.exists(alignments_json_file):

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()