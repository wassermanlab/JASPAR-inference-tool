#!/usr/bin/env python

import argparse
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
    This function parses arguments provided via the command line and returns
    an {argparse} object.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", help="input pwm (in any format)", metavar="FILE")
    parser.add_argument("-m", metavar="STR",
        help="motif id (i.e. what to place in the \"MOTIF\" field)")
    parser.add_argument("-o", help="output pwm (in MEME format)", metavar="FILE")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Reformat PWM to MEME
    pwm_to_meme(os.path.abspath(args.i), os.path.abspath(args.o), args.m)

def pwm_to_meme(input_file, output_file, motif_id):
    """
    This function reformats (in theory) any PWM to be used by MEME.

    @input:
    input_file {FILE}
    output_file {FILE}
    motif_id {STR} what to place in the "MOTIF" field
    """

    # Read matrix
    matrix = _read_matrix(input_file)

    # Transpose matrix
    if len(matrix) > len(matrix[0]):
        matrix = list(zip(*matrix))

    # Remove extra columns
    while len(matrix) > 4:
        matrix.pop(0)

    # Get number of sites
    nsites = _get_nsites(matrix)

    # Reformat matrix
    matrix = _reformat_matrix(matrix, nsites)

    # If not count matrix, set sites to 100
    if nsites == 1:
        nsites = 100

    # Write
    Jglobals.write(output_file, "MEME version 4")
    Jglobals.write(output_file, "")
    Jglobals.write(output_file, "ALPHABET= ACGT")
    Jglobals.write(output_file, "")
    Jglobals.write(output_file, "strands: + -")
    Jglobals.write(output_file, "")
    Jglobals.write(output_file, "Background letter frequencies")
    Jglobals.write(output_file, "A 0.25 C 0.25 G 0.25 T 0.25")
    Jglobals.write(output_file, "")
    Jglobals.write(output_file, "MOTIF %s" % motif_id)

    Jglobals.write(output_file, "letter-probability matrix: alength= 4 w= %s nsites= %s E= 0" % (len(matrix), nsites))

    # For each row...
    for row in matrix:
        Jglobals.write(output_file, "%s" % " ".join(["{:9.6f}".format(c) for c in row]))

def _read_matrix(matrix_file):

    # Initialize
    matrix = []

    # For each line
    for line in Jglobals.parse_file(matrix_file):

        # Skip comments
        if line.startswith("#") or line.startswith(">"):
            continue
    
        # Get data
        data = re.findall(r'[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?', line)

        # If data...
        if len(data) != 0:

            # Initialize
            matrix.append([])

            # For each value...
            for i in data:
                matrix[-1].append(i)

    return(matrix)

def _get_nsites(matrix):

    return(int(round(sum([float(r[0]) for r in matrix]))))

def _reformat_matrix(matrix, nsites):

    # Initialize
    format_matrix = []

    # For each row...
    for row in list(zip(*matrix)):
        format_matrix.append([float(i)/nsites for i in row])

    return(format_matrix)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()