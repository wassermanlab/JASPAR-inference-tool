#!/usr/bin/env python

import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import json
import os
import re
import shutil
import subprocess
import sys

# Defaults
out_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)))
root_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir)
files_dir = os.path.join(root_dir, "files")
lib_dir = os.path.join(root_dir, "lib")

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

    # Globals
    global cwd
    cwd = os.getcwd()

    # Create output dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Get groups
    get_groups(os.path.abspath(args.files_dir), os.path.abspath(args.o))

def get_groups(files_dir=files_dir, out_dir=out_dir):

    # Get groups
    __group_by_DBD(files_dir, out_dir)
    __group_by_matrix_cluster(files_dir, out_dir)
    # __group_by_sequence(files_dir, out_dir)

def __group_by_DBD(files_dir=files_dir, out_dir=out_dir):

    # Skip if JSON file already exists
    gzip_file = os.path.join(out_dir, ".groups.DBDs.json.gz")
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

def __group_by_matrix_cluster(files_dir=files_dir, out_dir=out_dir):

    # Skip if JSON file already exists
    gzip_file = os.path.join(out_dir, ".groups.matrix-clusters.json.gz")
    if not os.path.exists(gzip_file):

        # Initialize
        groups = {}

        # Create groups dir
        if not os.path.isdir("groups"):
            os.makedirs("groups")

        # Move to groups directory
        os.chdir("groups")

        # Skip if already done
        jaspar_profiles = "JASPAR2020_CORE_profiles.jaspar"
        if not os.path.exists(jaspar_profiles):

            # For each taxon...
            for taxon in Jglobals.taxons:

                # For each JASPAR profile...
                taxon_dir = os.path.join(files_dir, taxon)
                for jaspar_profile in os.listdir(taxon_dir):

                    # Skip
                    if not jaspar_profile.endswith(".jaspar"):
                        continue

                    # Write
                    jaspar_profile = os.path.join(taxon_dir, jaspar_profile)
                    for line in Jglobals.parse_file(jaspar_profile):
                        Jglobals.write(jaspar_profiles, line)

        # For stringency criterion...
        for criterion in ["stringent", "lenient"]:

            # Initialize
            prefix = "%s_clustering" % criterion

            # Skip if already done
            leafs_file = os.path.join(
                "%s_tables" % prefix, "leaf_to_cluster.tab"
            )
            if not os.path.exists(leafs_file):

                # Initialize
                if criterion ==  "stringent":
                    cor = 0.8
                    Ncor = 0.65
                else:
                    cor = 0.6
                    Ncor = 0.4

                # From RSAT matrix-clustering (PMID: 28591841)
                # Based on this study, we defined the default parameters:
                # the motif-to-motif similarity matrix is computed using Ncor
                # with a minimal alignment width of 5 columns, the motif tree
                # is built with the average linkage rule, and the partitioning
                # criterion combines thresholds on two metrics: cor ≥ 0.6 and
                # Ncor ≥ 0.4. [...] To obtain non-redundant motifs whilst
                # preserving specificity, we used more stringent partitioning
                # criteria than the default (cor ≥ 0.8 and Ncor ≥ 0.65).
                cmd = """
                    RSAT matrix-clustering -v 1 \
                    -matrix %s %s jaspar \
                    -hclust_method average \
                    -calc sum \
                    -title '%s' \
                    -metric_build_tree Ncor \
                    -lth w 5 -lth cor %s -lth Ncor %s \
                    -label_in_tree name \
                    -return json \
                    -quick \
                    -radial_tree_only \
                    -o %s
                """ % (prefix, jaspar_profiles, prefix, cor, Ncor, prefix)
                process = subprocess.run(
                    [cmd], shell=True, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, timeout=(3600 * 5)
                )

            # For each matrix...
            for matrix_id, cluster_id in Jglobals.parse_tsv_file(leafs_file):

                # Get matrix id
                m = re.search("(MA\d{4}\.\d)$", matrix_id)

                # Add cluster
                groups.setdefault(criterion, {})
                groups[criterion].setdefault(int(cluster_id), [])
                groups[criterion][int(cluster_id)].append(m.group(1))

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

def __group_by_sequence(files_dir=files_dir, out_dir=out_dir):

    # Skip if JSON file already exists
    gzip_file = os.path.join(out_dir, ".groups.sequence.json.gz")
    if not os.path.exists(gzip_file):

        # Initialize
        groups = {}
        prefix = "JASPAR2020_CORE"

        # Create groups dir
        if not os.path.isdir("groups"):
            os.makedirs("groups")

        # Move to groups directory
        os.chdir("groups")   

        # Skip if already done
        clusters_file = "%s_cluster.tsv" % prefix
        if not os.path.exists(clusters_file):

            # Format db
            cmd = """
                mmseqs easy-cluster %s/*.fa %s /tmp \
                -c 0.5 --cov-mode 1 --min-seq-id 0.2
            """ % (files_dir, prefix)
            process = subprocess.run([cmd], shell=True, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)

        # For each cluster...
        for uniacc, next_uniacc in Jglobals.parse_tsv_file(clusters_file):

            # Add cluster
            groups.setdefault(uniacc, [])
            groups[uniacc].append(next_uniacc)

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

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()