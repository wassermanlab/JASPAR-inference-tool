#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SubsMat.MatrixInfo import blosum62
import copy
import json
import numpy as np
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
from infer_profile import __get_X, __removeLowercase

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
    files_dir = os.path.abspath(args.files_dir)
    out_dir = os.path.abspath(args.o)

    # Create output dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Get data splits
    get_data_splits(files_dir, out_dir)

def get_data_splits(files_dir=files_dir, out_dir=out_dir):

    # Skip if JSON file already exists
    gzip_file = os.path.join(out_dir, "data_splits.json.gz")
    if not os.path.exists(gzip_file):

        # Initialize
        data_splits = {}

        # Create data splits dir
        if not os.path.isdir(os.path.join(out_dir, "data_splits")):
            os.makedirs(os.path.join(out_dir, "data_splits"))

        # Move to data splits directory
        os.chdir(os.path.join(out_dir, "data_splits"))

        # Load JSON file
        uniacc2matrix = {}
        uniacc2sequence = {}
        for taxon in Jglobals.taxons:
            json_file = os.path.join(files_dir, "%s.uniprot.json" % taxon)
            handle = Jglobals._get_file_handle(json_file)
            for uniacc, values in json.load(handle).items():
                uniacc2matrix.setdefault(uniacc, values[0])
                uniacc2sequence.setdefault(uniacc, values[1])

        # Group uniaccs by their Pfam DBD (i.e. Xs)
        dbd2uniaccs = __group_uniaccs_by_DBD(uniacc2sequence, files_dir)

        # Group matrices by their similarity (i.e. ys)
        matrix2cluster = __group_matrices_by_similarity(files_dir)

        # For each DBD, values...
        for DBD, values in dbd2uniaccs.items():

            # # Test C2H2 zinc fingers
            # if DBD != "zf-C2H2":
            #     continue

            # Initialize
            uniaccs = list(values.keys())

            # Get training / test sets
            training, test = __get_training_test_sets(
                DBD, uniaccs, uniacc2sequence
            )

            # For each UniProt Accession...
            for i in range(len(uniaccs) - 1):

                # Initialize
                seqi = []

                # For each DBD sequence...
                for k in range(len(values[uniaccs[i]])):
                    seqi.append(__removeLowercase(values[uniaccs[i]][k]))

                # For each other UniProt Accession...
                for j in range(i + 1, len(uniaccs)):

                    # Initialize
                    seqj = []

                    # For each other DBD sequence...
                    for k in range(len(values[uniaccs[j]])):
                        seqj.append(__removeLowercase(values[uniaccs[j]][k]))

                    # Train or test?
                    if uniaccs[i] in training and uniaccs[j] in training:
                        data_split = "training"
                    else:
                        data_split = "test"

                    # Get y
                    y = __get_y(
                        uniaccs[i], uniaccs[j], uniacc2matrix, matrix2cluster
                    )

                    # For each sequence similarity representation...
                    for similarity in ["identity", "blosum62"]:

                        # Get X
                        X = __get_X(seqi, seqj, similarity)

                        # Add data
                        data_splits.setdefault(DBD, {})
                        data_splits[DBD].setdefault(data_split, [])
                        data_splits[DBD][data_split].append(
                            [uniaccs[i], uniaccs[j], list(X), y, similarity]
                        )

            # Remove DBDs with:
            if DBD in data_splits:
                # 1) lack of data
                if len(data_splits[DBD]) != 2:
                    data_splits.pop(DBD, None)
                # 2) lack of positive or negative examples
                else:
                    for k in frozenset(data_splits[DBD]):
                        ys = set([v[3][0][0] for v in data_splits[DBD][k]])
                        if len(ys) != 2:
                            data_splits.pop(DBD, None)
                            break

        # Write
        Jglobals.write(
            gzip_file[:-3],
            json.dumps(data_splits, sort_keys=True, indent=4)
        )
        fi = Jglobals._get_file_handle(gzip_file[:-3], "rb")
        fo = Jglobals._get_file_handle(gzip_file, "wb")
        shutil.copyfileobj(fi, fo)
        fi.close()
        fo.close()
        os.remove(gzip_file[:-3])

def __group_uniaccs_by_DBD(uniacc2sequence, files_dir=files_dir):

    # Skip if JSON file already exists
    gzip_file = "groups.DBD.json.gz"
    if not os.path.exists(gzip_file):

        # Initialize
        dbd2uniaccs = {}

       # For each uniacc, sequence...
        for uniacc, sequence in uniacc2sequence.items():

            # Get Pfam alignments
            alignments = __get_Pfam_alignments(uniacc, sequence, files_dir)

            # For each alignment...
            for alignment in alignments:

                # Add DBD
                dbd2uniaccs.setdefault(alignment[0], {})
                dbd2uniaccs[alignment[0]].setdefault(uniacc, [])
                dbd2uniaccs[alignment[0]][uniacc].append(alignment[1])

        # Write
        Jglobals.write(
            gzip_file[:-3],
            json.dumps(dbd2uniaccs, sort_keys=True, indent=4)
        )
        fi = Jglobals._get_file_handle(gzip_file[:-3], "rb")
        fo = Jglobals._get_file_handle(gzip_file, "wb")
        shutil.copyfileobj(fi, fo)
        fi.close()
        fo.close()
        os.remove(gzip_file[:-3])

    # Return groupings
    handle = Jglobals._get_file_handle(gzip_file)
    return(json.load(handle))

def __get_Pfam_alignments(uniacc, sequence, files_dir=files_dir):

    # Initialize
    alignments = []
    seq_file = ".seq.fasta"
    hmm_db = os.path.join(files_dir, "pfam-DBDs", "all_DBDs.hmm")

    # Make seq file
    seq = Seq(sequence, IUPAC.protein)
    record = SeqRecord(seq, id=uniacc, name=uniacc, description=uniacc)
    __makeSeqFile(record, seq_file)

    # For each DBD...
    for pfam_ac, start, end, evalue in hmmScan(seq_file, hmm_db,
        non_overlapping_domains=True):

        # Initialize
        hmm_file = os.path.join(files_dir, "pfam-DBDs", "%s.hmm" % pfam_ac)

        # Make seq file
        sub_seq = seq[start:end]
        record = SeqRecord(sub_seq, id=uniacc, name=uniacc, description=uniacc)
        __makeSeqFile(record, seq_file)

        # Add DBDs
        alignment = hmmAlign(seq_file, hmm_file)
        alignments.append((pfam_ac, alignment, start+1, end, evalue))

    # Remove seq file
    if os.path.exists(seq_file):
        os.remove(seq_file)

    return(alignments)

def __group_matrices_by_similarity(files_dir=files_dir):

    # Skip if JSON file already exists
    gzip_file = "groups.matrix.json.gz"
    if not os.path.exists(gzip_file):

        # Initialize
        cluster2matrix = {}

        # Skip if JASPAR profiles already exist
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
                    stderr=subprocess.PIPE
                )

            # For each matrix...
            for matrix_id, cluster_id in Jglobals.parse_tsv_file(leafs_file):

                # Get matrix id
                m = re.search("(MA\d{4}\.\d)$", matrix_id)

                # Add cluster
                cluster2matrix.setdefault(criterion, {})
                cluster2matrix[criterion].setdefault(int(cluster_id), [])
                cluster2matrix[criterion][int(cluster_id)].append(m.group(1))

        # Remove JASPAR profiles
        if os.path.exists(jaspar_profiles):
            os.remove(jaspar_profiles)

        # Write
        Jglobals.write(
            gzip_file[:-3],
            json.dumps(cluster2matrix, sort_keys=True, indent=4)
        )
        fi = Jglobals._get_file_handle(gzip_file[:-3], "rb")
        fo = Jglobals._get_file_handle(gzip_file, "wb")
        shutil.copyfileobj(fi, fo)
        fi.close()
        fo.close()
        os.remove(gzip_file[:-3])

    # Return groupings
    matrix2cluster = {}
    handle = Jglobals._get_file_handle(gzip_file)
    for criterion, clusters in json.load(handle).items():
        matrix2cluster.setdefault(criterion, {})
        for cluster, matrix_ids in clusters.items():
            for matrix_id in matrix_ids:
                matrix2cluster[criterion].setdefault(matrix_id, int(cluster))
    handle.close()
    return(matrix2cluster)

def __get_training_test_sets(DBD, uniaccs, uniacc2sequence):

    # Skip if TXT file already exists
    txt_file = os.path.join(DBD, "repset.txt")
    if not os.path.exists(txt_file):

        # Initialize
        records = []
        seq_file = ".seq.fasta"

        # Create DBD dir
        if not os.path.exists(DBD):
            os.makedirs(DBD)

        # For each uniacc...
        for uniacc in uniaccs:

            # Make record
            seq = Seq(uniacc2sequence[uniacc], IUPAC.protein)
            record = SeqRecord(seq, id=uniacc, name=uniacc, description=uniacc)
            records.append(record)

        # Write sequences
        SeqIO.write(records, seq_file, "fasta")

        # Get representative sequences
        cmd = "%s/repset.py --outdir %s --seqs %s" % (
            lib_dir, os.path.abspath(DBD), seq_file
        )
        process = subprocess.run(
            [cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )

        # Remove seq file
        if os.path.exists(seq_file):
            os.remove(seq_file)

    # Get training (90%) and test (10%) sets
    seqs = [l for l in Jglobals.parse_file(txt_file)]
    training = seqs[:len(seqs) * 90 // 100]
    test = seqs[len(training):]

    return(set(training), set(test))

def __train_or_test_split(uniacc1, uniacc2, non_redundant):

    if uniacc1 in non_redundant:
        if len(non_redundant[uniacc1]) > 2:
            return("test")
    if uniacc2 in non_redundant:
        if len(non_redundant[uniacc2]) > 2:
            return("test")

    return("train")

def __get_y(uniacc1, uniacc2, uniacc2matrix, matrix2cluster):
    """
    Returns whether or not the matrices of two TFs belong to the same cluster
    as the dependent value (i.e. Y).
    """

    # Initialize
    y = [[0.], [0.]]
    criteria = ["stringent", "lenient"]

    for i in range(len(criteria)):
        for matrix_id1 in uniacc2matrix[uniacc1]:
            cluster1 = matrix2cluster[criteria[i]][matrix_id1]
            for matrix_id2 in uniacc2matrix[uniacc2]:
                if cluster1 == matrix2cluster[criteria[i]][matrix_id2]:
                    y[i][0] = 1.
        break

    return(y)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()