#!/usr/bin/env python

import argparse
import json
import os
import sys

# Defaults
out_dir = os.path.dirname(os.path.realpath(__file__))
files_dir = os.path.join(out_dir, os.pardir, "files")
profiles_dir = os.path.join(out_dir, os.pardir, "profiles")

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
    parser.add_argument("-p", default=profiles_dir, help="profiles directory (from get_profiles.py; default=../profiles/)", metavar="DIR")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Get models
    get_models(os.path.abspath(args.f), os.path.abspath(args.o), os.path.abspath(args.p))

def get_models(files_dir=files_dir, out_dir=out_dir, profiles_dir=profiles_dir):

    ##
    ## Step 1: create one file per DBD composition with aligned sequences
    ##

    _get_alignments()

def _get_alignments(files_dir=files_dir, out_dir=out_dir):

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
    for k in alignments.keys():
        print(k, len(alignments[k]))
    exit(0)

    # store the aligned sequences into the list that is in a dictionary


    with open("/homed/home/nzhang/JASPAR-tools/jaspartools/files/pfam_DBD_noVer.json",'r') as DBDs:
        DBDdic = json.loads(DBDs.read())

    # new dic
    alignedD = {k:[] for k in DBDdic}

    # loop through the pfam hmm files
    animals = ['fungi','nematodes','vertebrates','plants','insects']
    for taxa in animals:
        with open('/homed/home/nzhang/JASPAR-tools/jaspartools/files/'+ taxa +'.pfam.hmm.json','r') as pfamF:
            pfamDic = json.loads(pfamF.read())

        # get the sequence of the pfam and store in dict
        for key,item in pfamDic.items():
            for lst in item:
                pfamID = lst[0]
                if pfamID in alignedD:
                    alignedD[pfamID].extend(lst[1:])

    # write the output as json (DBD hmm aligned)
    with open("/homed/home/nzhang/JASPAR-tools/jaspartools/files/pfamAlignVisual/pfam_DBDAlignment.json",'w+') as outfile:
        json.dump(alignedD , outfile , sort_keys = True , indent = 4, separators = (",",": "))

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()