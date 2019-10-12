#!/usr/bin/env python

import argparse
from functools import partial
import json
from multiprocessing import Pool
import numpy
import operator
import os
from pathlib import Path
import pickle
import re
import shutil
import socket
import subprocess
import sys
from tqdm import tqdm

# Defaults
out_dir = os.path.dirname(os.path.realpath(__file__))
root_dir = os.path.join(out_dir, os.pardir)

# Append JASPAR-profile-inference to path
sys.path.append(root_dir)

# Import JASPAR-profile-inference functions
from __init__ import Jglobals
from files.get_files import Tomtom, _get_Tomtom_hits

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns an {argparse} object.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-c", help="cis-bp directory (from downloads.py)", metavar="DIR")
    parser.add_argument("-o", default=out_dir, metavar="DIR",
        help="output directory (default = ./)")
    parser.add_argument("--threads", default=1, type=int, metavar="INT",
        help="threads to use (default = 1)")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    parse_cisbp(os.path.abspath(args.c), os.path.abspath(args.o), args.threads)

def parse_cisbp(cisbp_dir, output_dir=out_dir, threads=1):

    # Initialize
    global cwd
    cwd = os.getcwd()

    # Get k-mers
    _get_kmers(cisbp_dir, output_dir)

    # Cis-BP to MEME format
    _reformat_to_meme(cisbp_dir, output_dir)

    # Group TFs by TF family
    _group_by_TF_family(cisbp_dir, output_dir)

    # Group matrices by Tomtom similarity
    _group_by_Tomtom(output_dir, threads)

    # Group matrices by k-mer overlap
    _group_by_kmer_overlap(output_dir)

def _get_kmers(cisbp_dir, output_dir=out_dir):

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

def _reformat_to_meme(cisbp_dir, output_dir=out_dir):

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

def _group_by_TF_family(cisbp_dir, output_dir=out_dir):

    # Skip if groups JSON file already exists
    groups_json_file = os.path.join(output_dir, "groups.families.json")
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

def _get_motifs(cisbp_dir, output_dir=out_dir):

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
            family = "+".join(list(sorted(m.group(2).split(","))))
            # Fix AP-2
            if family == "AP-2":
                family = "AP2"
            families.setdefault(m.group(1), family)

    return(families)

def _group_by_Tomtom(output_dir=out_dir, threads=1):

    # Skip if groups JSON file already exists
    groups_json_file = os.path.join(output_dir, "groups.tomtom.json")
    if not os.path.exists(groups_json_file):

        # Initialize
        tomtom = {}

        # Get all profiles
        profiles = [str(p) for p in Path(output_dir).glob("*/*.meme")]

        # Skip if MEME database already exists
        database = os.path.join(output_dir, "cisbp.meme")
        if not os.path.exists(database):

            # For each profile...
            for profile in profiles:

                # Cat to database
                os.system("cat %s >> %s" % (profile, database))

        # Parallelize
        pool = Pool(threads)
        parallelized = partial(Tomtom, database=database, out_dir=output_dir)
        for _ in tqdm(
            pool.imap(parallelized, profiles), desc="Tomtom", total=len(profiles)
        ):
            pass
        pool.close()
        pool.join()

        # Move to output directory
        os.chdir(output_dir)

        # For each profile...
        for profile in profiles:

            # Initialize
            m = re.search("(M\d{4}_1.02).meme$", profile)
            tomtom_dir = ".%s" % m.group(1)

            # Get hits
            tomtom.setdefault(m.group(1), _get_Tomtom_hits(tomtom_dir))

            # Remove Tomtom directory
            shutil.rmtree(tomtom_dir)

        # Write
        Jglobals.write(
            groups_json_file,
            json.dumps(tomtom, sort_keys=True, indent=4, separators=(",", ": "))
        )

        # For each taxon...
        for taxon in Jglobals.taxons:

            # Remove taxon directory
            if os.path.isdir(taxon):
                shutil.rmtree(taxon)

        # Change dir
        os.chdir(cwd)

def _group_by_kmer_overlap(output_dir=out_dir):

    # Initialize
    kmers_dir = os.path.join(output_dir, "kmers")

    # Skip if groups JSON file already exists
    groups_json_file = os.path.join(output_dir, "groups.overlap.json")
    if not os.path.exists(groups_json_file):

        # Initialize
        kmers = {}
        overlaps = {}

        # Load JSON file
        tomtom_json_file = os.path.join(output_dir, "groups.tomtom.json")
        with open(tomtom_json_file) as f:
            tomtom = json.load(f)

        # For each k-mers pickle file...
        for pickle_file in os.listdir(kmers_dir):

            # If valid file...
            m = re.search("^(M\d{4}_1.02).pickle$", pickle_file)
            if m:

                # Load pickle file
                with open(os.path.join(kmers_dir, pickle_file), "rb") as f:
                    kmers.setdefault(m.group(1), pickle.load(f))

        # For each matrix ID...
        for matrix in tomtom:

            overlaps.setdefault(matrix, [])

            # For each matrix ID...
            for next_matrix, evalue in tomtom[matrix]:
                
                intersection = len(kmers[matrix].intersection(kmers[next_matrix]))
                union = len(kmers[matrix].union(kmers[next_matrix]))
                overlaps[matrix].append([next_matrix, float(intersection)/union])

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