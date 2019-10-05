#!/usr/bin/env python

import argparse
import json
import os
import pickle
import sys

# Defaults
out_dir = os.path.dirname(os.path.realpath(__file__))
root_dir = os.path.join(out_dir, os.pardir)
files_dir = os.path.join(root_dir, "files")
models_dir = os.path.join(root_dir, "models")

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
    parser.add_argument("-m", default=models_dir, help="models directory (from models.py; default=../models/)", metavar="DIR")
    parser.add_argument("-o", default=out_dir, help="output directory (default = ./)", metavar="DIR")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Make figures & tables
    figs_n_tabs(os.path.abspath(args.f), os.path.abspath(args.m), os.path.abspath(args.o))

def figs_n_tabs(files_dir=files_dir, models_dir=models_dir, out_dir=out_dir):

    # Skip if figure already exists
    figure_file = os.path.join(out_dir, "Fig2B.png")
    if not os.path.exists(figure_file):

        # Load Cis-BP file
        cisbp_file = os.path.join(out_dir, "Fig2B.json")
        with open(cisbp_file) as f:
            cisbp = json.load(f)

        # Load models file
        models_file = os.path.join(models_dir, "models.pickle")
        with open(pairwise_file, "rb") as f:
            models = pickle.load(f)

        print(cisbp)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()