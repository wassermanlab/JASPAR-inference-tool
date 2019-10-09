#!/usr/bin/env python

import argparse
import os

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns an {argparse} object.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-o", default="./", help="output directory (default = ./)", metavar="DIR")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Get Cis-BP
    get_cisbp(os.path.abspath(args.o))

def get_cisbp(output_dir):

    # Initialize
    cwd = os.getcwd()

    # Create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Move to output directory
    os.chdir(output_dir)

    # Skip if TFs file already exists
    if not os.path.exists("cisbp_1.02.tfs.sql"):

        # Download SQL files
        os.system("curl --silent -O http://cisbp.ccbr.utoronto.ca/data/1.02/DataFiles/SQLDumps/SQLArchive_cisbp_1.02.zip")

        # Unzip
        os.system("unzip -qq SQLArchive_cisbp_1.02.zip")

        # Remove SQL files
        os.remove("SQLArchive_cisbp_1.02.zip")

        # For each ZIP file...
        for zip_file in frozenset(os.listdir(os.getcwd())):

            # Skip non-zip files
            if not zip_file.endswith(".zip"):
                continue

            # Unzip
            os.system("unzip -qq %s" % zip_file)
            os.remove(zip_file)

    # Skip if E-scores file already exists
    if not os.path.exists("Escores.txt"):

        # Download E-scores file
        os.system("curl --silent -O http://cisbp.ccbr.utoronto.ca/data/1.02/DataFiles/Bulk_downloads/EntireDataset/Escores.txt.zip")

        # Unzip
        os.system("unzip -qq Escores.txt.zip")

        # Remove E-scores file
        os.remove("Escores.txt.zip")

    # Skip if PWMs dir already exists
    if not os.path.isdir("pwms"):

        # Download PWMs file
        os.system("curl --silent -O http://cisbp.ccbr.utoronto.ca/data/1.02/DataFiles/Bulk_downloads/EntireDataset/PWMs.zip")

        # Unzip
        os.system("unzip -qq PWMs.zip")

        # Remove PWMs file
        os.remove("PWMs.zip")

    # Return to original directory
    os.chdir(cwd)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()