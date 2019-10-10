#!/usr/bin/env python

import argparse
import json
import numpy
import operator
import os
import pickle
import re
import shutil
import socket
import subprocess
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

    parser.add_argument("-c", help="cis-bp directory (from downloads.py)", metavar="DIR")
    parser.add_argument("-o", default="./", help="output directory (default = ./)", metavar="DIR")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    parse_cisbp(os.path.abspath(args.c), os.path.abspath(args.o))

def parse_cisbp(cisbp_dir, output_dir="./"):

    # Initialize
    global cwd
    cwd = os.getcwd()

    # Get k-mers
    _get_kmers(cisbp_dir, output_dir)

    # Cis-BP to MEME format
    _reformat_to_meme(cisbp_dir, output_dir)

    # Group TFs by TF family
    _group_by_TF_family(cisbp_dir, output_dir)

def _get_kmers(cisbp_dir, output_dir="./"):

    # If k-mers directory does not exist...
    kmers_dir = os.path.join(output_dir, "kmers")
    if not os.path.exists(kmers_dir):

        # Create k-mers dir
        os.makedirs(kmers_dir)

        # Get escores
        escores = _get_escores(cisbp_dir)

        # Transpose matrix
        escores = list(zip(*escores))

        # Get k-mers
        kmers = escores.pop(0)

        # Change dir
        os.chdir(kmers_dir)

        # For each motif...
        for escore in escores:

            # Skip if pickle file already exists
            pickle_file = "%s.pickle" % escore[0]
            if not os.path.exists(pickle_file):

                # Initialize
                positive_kmers = set()

                # For each k-mer...
                for k in range(1, len(kmers)):

                    # If E-score < 0.45...
                    try:
                        if float(escore[k]) >= 0.45:
                            positive_kmers.add(kmers[k])
                    except:
                        pass

                # Write pickle file
                with open(pickle_file, "wb") as f:
                    pickle.dump(positive_kmers, f)

        # Change dir
        os.chdir(cwd)

def _get_escores(cisbp_dir):

    # Initialize
    escores = []

    # For each line...
    for line in Jglobals.parse_file(os.path.join(cisbp_dir, "Escores.txt")):

        # Initialize
        line = line.split("\t")

        # Get k-mer
        kmer = line.pop(0)

        # If 1st line...
        if len(escores) == 0:
            escores.append([kmer] + [i[:10] for i in line])

        # ... Else...
        else:

            escores.append([kmer])

            # For each E-score
            for E in line:
                try:
                    escores[-1].append('{0:.3g}'.format(float(E)))
                except:
                    escores[-1].append(None)

    return(escores)

def _reformat_to_meme(cisbp_dir, output_dir="./"):

    # Initialize
    kmers_dir = os.path.join(output_dir, "kmers")

    # If MEME directory does not exis...
    meme_dir = os.path.join(output_dir, "meme")
    if not os.path.exists(meme_dir):

        # Create TFs dir
        os.makedirs(meme_dir)

        # Change dir
        os.chdir(os.path.dirname(os.path.realpath(__file__)))

        # For each k-mers pickle file...
        for pickle_file in os.listdir(kmers_dir):

            # If valid file...
            m = re.search("^(M\d{4}_1.02).pickle$", pickle_file)
            if m:

                # Reformat PWM to MEME
                pwm_file = os.path.join(cisbp_dir, "pwms", "%s.txt" % m.group(1))
                meme_file = os.path.join(meme_dir, "%s.meme" % m.group(1))
                cmd = "python reformat2meme.py -i %s -m %s -o %s" % (pwm_file,
                    m.group(1), meme_file)
                process = subprocess.run([cmd], shell=True, stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL)

        # Change dir
        os.chdir(cwd)

def _group_by_TF_family(cisbp_dir, output_dir="./"):

    # Skip if groups JSON file already exists
    groups_json_file = os.path.join(out_dir, "groups.families.json")
    if not os.path.exists(groups_json_file):

        # Initialize
        groups = {}

        # Get TF motifs
        motifs = _get_motifs(cisbp_dir, output_dir)

        # Get TF families
        families = _get_families(cisbp_dir)

        # For each line...
        for line in Jglobals.parse_file(os.path.join(cisbp_dir, "cisbp_1.02.tfs.sql")):

            # If valid line...
            m = re.search("\('(.+)', '(.+)', '.+', '.+', '.+', '.+', '\w'\),*", line)
            if m:

                if m.group(1) in motifs:

                    # Add member
                    groups.setdefault(families[m.group(2)], [])
                    groups[families[m.group(2)]].append([motifs[m.group(1)], m.group(1)])

        # Write
        Jglobals.write(
            groups_json_file,
            json.dumps(groups, sort_keys=True, indent=4, separators=(",", ": "))
        )

        # Change dir
        os.chdir(cwd)

def _get_motifs(cisbp_dir, output_dir="./"):

    # Initialize
    motifs = {}
    kmers_dir = os.path.join(output_dir, "kmers")

    # For each line...
    for line in Jglobals.parse_file(os.path.join(cisbp_dir, "cisbp_1.02.motifs.sql")):

        # If valid line...
        m = re.search("\('(.+)', '(.+)', '.+', '.+', 'PBM', '.+', '.+', '.+'\),*", line)
        if m:

            # Ignore if no k-mers pickle file...
            if not os.path.exists(os.path.join(kmers_dir, "%s.pickle" % m.group(1))):
                continue

            motifs.setdefault(m.group(2), [])
            motifs[m.group(2)].append(m.group(1))

    return(motifs)

def _get_families(cisbp_dir):

    # Initialize
    families = {}

    # For each line...
    for line in Jglobals.parse_file(os.path.join(cisbp_dir, "cisbp_1.02.tf_families.sql")):

        # If valid line...
        m = re.search("\('(.+)', '(.+)', '.+', \d+, .+\),*", line)
        if m:
            families.setdefault(m.group(1), "+".join(list(sorted(m.group(2).split(",")))))

    return(families)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()