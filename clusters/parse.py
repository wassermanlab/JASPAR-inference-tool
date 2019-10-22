#!/usr/bin/env python

import argparse
import json
import os
import re
import sys

# Defaults
out_dir = os.path.dirname(os.path.realpath(__file__))
root_dir = os.path.join(out_dir, os.pardir)

# Append JASPAR-profile-inference to path
sys.path.append(root_dir)

# Import JASPAR-profile-inference functions
from __init__ import Jglobals

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns an {argparse} object.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-c", metavar="DIR",
        help="clusters directory (from downloads.py)")
    parser.add_argument("-o", default=out_dir, metavar="DIR",
        help="output directory (default = ./)")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    parse_clusters(os.path.abspath(args.c), os.path.abspath(args.o))

def parse_clusters(clusters_dir, output_dir=out_dir):

    # Initialize
    cwd = os.getcwd()

    # Create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # For each taxon...
    for taxon in Jglobals.taxons:

        # Create taxon directory
        taxon_dir = os.path.join(clusters_dir, taxon)

        # Move to taxon directory
        os.chdir(taxon_dir)

        # For each JSON file...
        for json_file in os.listdir(taxon_dir):

            # Initialize
            cluster = set()

            # For each line...
            for line in Jglobals.parse_file(json_file):

                # Get matrix ID
                m = re.search("\"label\": \".+(MA\d{4}_\d)\",", line)
                if m:
                    matrix_id = m.group(1).replace("_", ".")
                    cluster.add(matrix_id)
            print(cluster)

        exit(0)

        # Return to original directory
        os.chdir(cwd)

        # Write
        Jglobals.write(
            groups_json_file,
            json.dumps(overlaps, sort_keys=True, indent=4, separators=(",", ": "))
        )

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()