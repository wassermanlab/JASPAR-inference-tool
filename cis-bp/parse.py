#!/usr/bin/env python

import argparse
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
# Class       #
#-------------#

class TF(object):
    """
    This class defines an {TF} object.

    """

    def __init__(self, id=None, name=None, specie=None, family=None, motifs=None, sources=None, sequences=None, file_name=None):
        self._file = file_name
        self._id = id
        self._name = name
        self._specie = specie
        self._family = family
        self._motifs = motifs
        self._sources = sources
        self._sequences = sequences

        if self._file is not None:
            self._parse_file()

    def _parse_file(self):
        for line in functions.parse_file(self._file):
            if line.startswith("#"): continue
            line = line.split(";")
            self._id = line[0]
            self._name = line[1]
            self._specie = line[2]
            self._family = line[3].split(",")
            self._motifs = line[4].split(",")
            self._sources = line[5].split(",")
            self._sequences = line[6].split(",")

    def get_id(self):
        return self._id

    def get_name(self):
        return self._name

    def get_family(self, format=True):
        if self._family[0] == "AP-2": return "AP2" # Fixes weird AP-2 family annotation
        
        if format:
            family = self._family[0]
            family = family.replace("/", "_")
            family = family.replace(" ", "_")
            return family
            
        return self._family[0]

    def write(self, file_name):
        # Initialize #
        if os.path.exists(file_name): os.remove(file_name)
        functions.write(file_name, "#id;name;specie;family;motifs;sources;sequences")
        functions.write(file_name, "%s;%s;%s;%s;%s;%s;%s" % (self._id, self._name, self._specie, ",".join(self._family), ",".join(self._motifs), ",".join(self._sources), ",".join(self._sequences)))

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
    cwd = os.getcwd()

    # If k-mers dir does not exist...
    global kmers_dir
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

    # If TFs dir does not exist...
    tfs_dir = os.path.join(output_dir, "tfs")
    if not os.path.exists(tfs_dir):

        # Create TFs dir
        os.makedirs(tfs_dir)

        # Change dir
        os.chdir(tfs_dir)

        # Get TF motifs
        motifs = _get_motifs(cisbp_dir)

        # Get TF families
        families = _get_families(cisbp_dir)

        # For each line...
        for line in Jglobals.parse_file(os.path.join(cisbp_dir, "cisbp_1.02.tfs.sql")):

            # If valid line...
            m = re.search("\('(.+)', '(.+)', '.+', '.+', '(.+)', '(.+)', 'D'\),*", line)
            if m:

                print(line)
                exit(0)

                # For each motif...
                for motif in motifs:
                    if motifs[motif][0] == m.group(1):
                        tf_motifs.append(motif)
                        tf_sources.append(sources[motifs[motif][1]])
                        tf_sequences.append(motifs[motif][2])
                if len(tf_motifs) > 0:
                    tf_obj = TF(m.group(1), m.group(3), re.sub("_", " ", m.group(4)), families[m.group(2)], tf_motifs, tf_sources, tf_sequences)
                    # Skip if TF file already exists #
                    tf_file = os.path.join(os.path.abspath(options.output_dir), "tfs", "%s.txt" % m.group(1))
                    if not os.path.exists(tf_file):
                        tf_obj.write(tf_file)


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

def _get_motifs(cisbp_dir):

    # Initialize
    motifs = {}

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
            families.setdefault(m.group(1), set(m.group(2).split(",")))

    return(families)

# ###############################
# # 2. Get associated PBM files #
# ###############################
# if options.verbose: sys.stdout.write("\nGet associated PBM files...\n")
# # Skip if starts later #
# if options.start_step <= 2:
#     # Verbose mode #
#     if options.verbose: sys.stdout.write("\n")
#     # For each TF... #
#     for tf_file in sorted(os.listdir(os.path.join(os.path.abspath(options.output_dir), "tfs"))):
#         # Get TF object #
#         tf_obj = TF(file_name=os.path.join(os.path.abspath(options.output_dir), "tfs", tf_file))
#         # Verbose mode... #
#         if options.verbose: sys.stdout.write("\t%s...\n" % tf_obj.get_id())
#         # For each sequence... #
#         for i in range(len(tf_obj._sequences)):
#             ##############################
#             # 2.1 Align positive k-mers  #
#             ##############################
#             if options.verbose: sys.stdout.write("\t\t-- aligning positive k-mers...\n")
#             # Skip if k-mers file already exists #
#             kmers_file = os.path.join(os.path.abspath(options.output_dir), "kmers", "%s.%s.txt" % (tf_obj.get_id(), str(i)))
#             if not os.path.exists(kmers_file):
#                 # Initialize #
#                 positive_kmers = []
#                 matches = set()
#                 # For each line... #
#                 for line in functions.parse_file(os.path.join(os.path.abspath(options.output_dir), "motifs", "%s.txt" % tf_obj._motifs[i])):
#                     if line.startswith("#"): continue
#                     line = line.split(";")
#                     try:
#                         # If positive k-mer... #
#                         if float(line[2]) >= float(config.get("Parameters", "min_escore_positives")):
#                             positive_kmers.append([line[0], float(line[2])])
#                             positive_kmers.append([line[1], float(line[2])])
#                     except: pass
#                 # If positive k-mers... #
#                 if len(positive_kmers) > 0:
#                     # Initialize #
#                     pwm = []
#                     done = set()
#                     pwm_file = os.path.join(os.path.abspath(options.cisbp_dir), "pwms", "%s.txt" % tf_obj._motifs[i])
#                     functions.write(kmers_file, "#kmer;position")
#                     # For each line... #
#                     for line in functions.parse_file(pwm_file):
#                         line = line.split("\t")
#                         position = line.pop(0)
#                         if position != "Pos":
#                             pwm.append(map(str, [round(float(j) * 1000, 0) for j in line]))
#                     # For each sliding window... #
#                     for j in range(len(pwm) - 8 + 1):
#                         # Initialize #
#                         dummy_file = os.path.join(os.path.abspath(options.dummy_dir), "%s.%s.pfm" % (os.path.basename(__file__), os.getpid()))
#                         # For each row... #
#                         sub_pwm = pwm[j:j + 8]
#                         for k in zip(*sub_pwm):
#                             functions.write(dummy_file, "\t".join(k))
#                         # Get motif #
#                         f = open(dummy_file)
#                         motif = motifs.read(f, "pfm")
#                         f.close()
#                         # Remove dummy file #
#                         os.remove(dummy_file)
#                         # Get pseudocounts #
#                         motif.pseudocounts = motifs.jaspar.calculate_pseudocounts(motif)
#                         # Get motif max. score, min. score and score threshold #
#                         max_score = motif.pssm.max
#                         min_score = motif.pssm.min
#                         motif_score_threshold = (max_score - min_score) * 0.8 + min_score
#                         # For each k-mer... #
#                         for kmer in sorted(positive_kmers, key=lambda x: x[1], reverse=True):
#                             # Format sequence #
#                             sequence = Seq("".join(kmer[0]), IUPAC.unambiguous_dna)
#                             # Get PWM matches #
#                             pwm_matches = [[position, score] for position, score in motif.pssm.search(sequence, threshold=motif_score_threshold)]
#                             # For each matching position, score... #
#                             for position, score in sorted(pwm_matches, key=lambda x: x[-1], reverse=True):
#                                 if position == 0:
#                                     matches.add((kmer[0], j, kmer[1], score))
#                     # For each hit... #
#                     for match in sorted(matches, key=operator.itemgetter(2, 3), reverse=True):
#                         # Get best match #
#                         if match[0] in done: continue
#                         functions.write(kmers_file, "%s;%s" % (match[0], match[1]))
#                         done.add(match[0])
#             ##############################
#             # 2.2 Extract TF sequences   #
#             ##############################
#             if options.verbose: sys.stdout.write("\t\t-- extracting TF sequences...\n")
#             # Skip if k-mers file does not exists #
#             if not os.path.exists(kmers_file): continue
#             # Skip if sequence file already exists #
#             sequence_file = os.path.join(os.path.abspath(options.output_dir), "sequences", "%s.%s.fa" % (tf_obj.get_id(), str(i)))
#             if not os.path.exists(sequence_file):
#                 functions.write(sequence_file, ">%s\n%s" % (tf_obj.get_id(), tf_obj._sequences[i]))

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()