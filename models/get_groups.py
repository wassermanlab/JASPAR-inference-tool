#!/usr/bin/env python

import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import json
import os
import shutil
import sys
from tqdm import tqdm

# Defaults
out_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)))
root_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir)
files_dir = os.path.join(root_dir, "files")

# Append JASPAR-profile-inference to path
sys.path.append(root_dir)

# Import globals
from __init__ import Jglobals
from infer_profile import hmmAlign, hmmScan, __makeSeqFile

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns
    an {argparse} object.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("--files-dir", default=files_dir, metavar="DIR",
        help="files directory (from get_files.py)")
    parser.add_argument("-o", default=out_dir, metavar="DIR",
        help="output directory (default = ./)")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Get groups
    get_groups(os.path.abspath(args.files_dir), os.path.abspath(args.o))

def get_groups(files_dir=files_dir, out_dir=out_dir):

    # Globals
    global cwd
    cwd = os.getcwd()

    # Create output dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Get groups
    __get_groups(files_dir, out_dir)

def __get_groups(files_dir=files_dir, out_dir=out_dir):

    # Skip if JSON file already exists
    gzip_file = os.path.join(out_dir, "groups.json.gz")
    if not os.path.exists(gzip_file):

        # Initialize
        groups = {}

        # Move to output directory
        os.chdir(out_dir)

        # For each taxon...
        for taxon in Jglobals.taxons:

            # Get Pfam alignments
            pfams = __get_Pfam_alignments(taxon, files_dir, out_dir)

            # For each UniProt Accession...
            for uniacc, values in pfams.items():

                for value in values:

                    # Add DBD
                    groups.setdefault(value[0], {})
                    groups[value[0]].setdefault(uniacc, [])
                    groups[value[0]][uniacc].append(value[1])

        # Write
        Jglobals.write(
            gzip_file[:-3],
            json.dumps(groups, sort_keys=True, indent=4)
        )
        fi = Jglobals._get_file_handle(gzip_file[:-3], "rb")
        fo = Jglobals._get_file_handle(gzip_file, "wb")
        shutil.copyfileobj(fi, fo)
        fi.close()
        fo.close()
        os.remove(gzip_file[:-3])

        # Change dir
        os.chdir(cwd)

def __get_Pfam_alignments(taxon, files_dir=files_dir, out_dir=out_dir):

    # Change dir
    os.chdir(out_dir)

    # Initialize
    pfams = {}
    seq_file = ".seq.fasta"
    hmm_db = os.path.join(files_dir, "pfam-DBDs", "all_DBDs.hmm")

    # Load JSON files
    json_file = os.path.join(files_dir, "%s.uniprot.json" % taxon)
    with open(json_file) as f:
        uniaccs = json.load(f)

    # For each uniacc...
    for uniacc in uniaccs:

        # Initialize
        pfams.setdefault(uniacc, [])

        # Make seq file
        seq = Seq(uniaccs[uniacc][1], IUPAC.protein)
        seq_record = SeqRecord(seq, id=uniacc, name=uniacc,
            description=uniacc)
        __makeSeqFile(seq_record, seq_file)

        # For each DBD...
        for pfam_ac, start, end, evalue in hmmScan(seq_file, hmm_db,
            non_overlapping_domains=True):

            # Initialize
            hmm_file = os.path.join(files_dir, "pfam-DBDs", "%s.hmm" % pfam_ac)

            # Make seq file
            sub_seq = seq[start:end]
            seq_record = SeqRecord(sub_seq, id=uniacc, name=uniacc,
                description=uniacc)
            __makeSeqFile(seq_record, seq_file)

            # Add DBDs
            alignment = hmmAlign(seq_file, hmm_file)
            pfams[uniacc].append((pfam_ac, alignment, start+1, end, evalue))

    # Remove seq file
    if os.path.exists(seq_file):
        os.remove(seq_file)

    # Change dir
    os.chdir(cwd)

    return(pfams)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()