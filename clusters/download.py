#!/usr/bin/env python

import argparse
import os
from pathlib import Path
import re
import shutil
import sys

# Defaults
out_dir = os.path.dirname(os.path.realpath(__file__))
root_dir = os.path.join(out_dir, os.pardir)

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

    parser.add_argument("-o", default=out_dir, metavar="DIR",
        help="output directory (default = ./)")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Get clusters
    get_clusters(os.path.abspath(args.o))

def get_clusters(output_dir):

    # Initialize
    cwd = os.getcwd()

    # Create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # For each taxon...
    for taxon in Jglobals.taxons:

        # Initialize
        url = "http://folk.uio.no/jamondra/JASPAR_2020_clusters/%s/" % taxon
        clusters_file = "JASPAR_2020_matrix_clustering_%s_archive.zip" % taxon

        # Create taxon directory
        taxon_dir = os.path.join(out_dir, taxon)
        if not os.path.exists(taxon_dir):
            os.makedirs(taxon_dir)

        # Move to taxon directory
        os.chdir(taxon_dir)

        # If interactive trees directory does not exist...
        if not os.path.isdir("interactive_trees"):

            # Download clusters
            os.system("curl --silent -O %s" % os.path.join(url, "interactive_trees", clusters_file))

            # Unzip
            os.system("unzip -qq %s" % clusters_file)

            # Remove SQL files
            os.remove(clusters_file)

        # For each cluster...
        for cluster_file in Path(taxon_dir).glob("interactive_trees/*/parsed_tree_cluster_*.json"):

            # Initialize
            m = re.search("(parsed_tree_cluster_\d+.json)$", str(cluster_file))

            # If cluster does not exist...
            if not os.path.exists(m.group(1)):
                shutil.copy(str(cluster_file), m.group(1))

        # Return to original directory
        os.chdir(cwd)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()