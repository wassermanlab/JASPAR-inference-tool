#!/usr/bin/env python

import argparse
import os
import re
import shutil
import subprocess
import sys

# Defaults
out_dir = os.path.dirname(os.path.realpath(__file__))

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

    parser.add_argument("-d", "--devel", action="store_true", help="development mode (uses hfaistos; default = False)")
    parser.add_argument("-o", default=out_dir, help="output directory (default = ./)", metavar="DIR")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Make Tomtom files
    get_tomtom(args.devel, os.path.abspath(args.o))

def get_tomtom(devel=False, out_dir=out_dir):

    # Globals
    global cwd
    cwd = os.getcwd()

    # Create output dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Get JASPAR database
    _get_jaspar_database(devel, out_dir)

    # Get Tomtom e-values
    _get_tomtom_evalues(out_dir)

def _get_jaspar_database(devel=False, out_dir=out_dir):

    # Skip if JASPAR database already exists
    jaspar_database = os.path.join(out_dir, "jaspar.meme")
    if not os.path.exists(jaspar_database):

        # Get JASPAR URL
        jaspar_url = "http://jaspar.genereg.net/download/CORE/"
        if devel:
            jaspar_url = "http://hfaistos.uio.no:8002/download/CORE/"

        # For each taxon...
        for taxon in Jglobals.taxons:

            # Skip if taxon directory already exists
            taxon_dir = os.path.join(out_dir, taxon)
            if not os.path.exists(taxon_dir):

                # Initialize
                jaspar_file = "JASPAR2018_CORE_%s_redundant_pfms_meme.zip" % taxon
                if devel:
                    jaspar_file = "JASPAR2020_CORE_%s_redundant_pfms_meme.zip" % taxon

                # Create taxon directory
                os.makedirs(taxon_dir)

                # Move to taxon directory
                os.chdir(taxon_dir)

                # Download JASPAR profiles in MEME format
                _get_profiles(jaspar_url, jaspar_file)

                # For each MEME file...
                for meme_file in os.listdir("."):

                    # For each line...
                    for line in Jglobals.parse_file(meme_file):

                        # Write
                        Jglobals.write(jaspar_database, line)

                # Change dir
                os.chdir(cwd)

def _get_profiles(url, file_name):

        # Get JASPAR profiles
        os.system("curl --silent -O %s" % os.path.join(url, file_name))

        # Unzip
        os.system("unzip -qq %s" % file_name)

        # Remove zip files
        os.remove("%s" % file_name)

def _get_tomtom_evalues(out_dir=out_dir):
    """
    From http://meme-suite.org/doc/tomtom.html;
    In order to compute the scores, Tomtom needs to know the frequencies of the letters of the sequence alphabet in
    the database being searched (the 'background' letter frequencies). By default, the background letter frequencies
    included in the query motif file are used. The scores of columns that overlap for a given offset are summed.
    This summed score is then converted to a p-value. The reported p-value is the minimal p-value over all possible
    offsets. To compensate for multiple testing, each reported p-value is converted to an E-value by multiplying it
    by twice the number of target motifs. As a second type of multiple-testing correction, q-values for each match
    are computed from the set of p-values and reported.

    From PMID:17324271;
    [...]
    We show that Tomtom correctly assigns E values less than 0.01 to a large percentage of positive matches.
    [...]
    """

    # Initialize
    tomtom = {}
    jaspar_database = "jaspar.meme"
    tomtom_dir = "./tomtom_out/"

    # Skip if Tomtom JSON file already exists
    tomtom_json_file = os.path.join(out_dir, "tomtom.json")
    if not os.path.exists(tomtom_json_file):

        # Move to output directory
        os.chdir(out_dir)

        # Remove Tomtom directory
        if os.path.isdir(tomtom_dir):
            shutil.rmtree(tomtom_dir)

        # # For each taxon...
        # for taxon in Jglobals.taxons:

        #     # For each MEME file...
        #     for meme_file in os.listdir(taxon):

        #         # Initialize
        #         m = re.search("(MA\d{4}.\d).meme", meme_file)
        #         tomtom.setdefault(m.group(1), [])

        # Run Tomtom
        cmd = "tomtom -thresh 0.01 -evalue %s %s" % (jaspar_database, jaspar_database)
        process = subprocess.run([cmd], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # For each line...
        for line in Jglobals.parse_tsv_file(os.path.join(tomtom_dir, "tomtom.tsv")):

            # Skip header
            if line[0] == "Query_ID":
                continue

            # Skip self
            if line[0] == line[1]:
                continue

            # Add to Tomtom
            tomtom.setdefault(line[0], [])
            tomtom[line[0]].append(line[1])

        # Remove Tomtom directory
        shutil.rmtree(tomtom_dir)

        # Change dir
        os.chdir(cwd)

        # Write
        Jglobals.write(
            pfam_json_file,
            json.dumps(tomtom, sort_keys=True, indent=4, separators=(",", ": "))
        )

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()