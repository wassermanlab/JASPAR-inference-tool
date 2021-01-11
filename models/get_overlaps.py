#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from functools import partial
import json
import math
from multiprocessing import Pool
import numpy as np
import os
import pandas as pd
import random
import re
import shutil
import subprocess
import sys
from tqdm import tqdm


# Defaults
out_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "overlaps")
root_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir)
files_dir = os.path.join(root_dir, "files")
jaspar_ucsc_tracks = os.path.join(root_dir, "submodules", "JASPAR-UCSC-tracks")

# Append JASPAR-profile-inference to path
sys.path.append(root_dir)

# Import globals
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

    parser.add_argument("--files-dir", default=files_dir, metavar="DIR",
        help="files directory (from get_files.py)")
    parser.add_argument("-o", default=out_dir, metavar="DIR",
        help="output directory (default = ./)")
    parser.add_argument("--threads", default=1, metavar="INT",
        help="threads to use (default = 1)")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Get overlaps
    get_overlaps(
        os.path.abspath(args.files_dir), os.path.abspath(args.o),
        int(args.threads)
    )

def get_overlaps(files_dir=files_dir, out_dir=out_dir, threads=1):

    # Globals
    global cwd
    cwd = os.getcwd()

    # Create output dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Get overlaps
    __get_overlaps(files_dir, out_dir, threads)

def __get_overlaps(files_dir=files_dir, out_dir=out_dir, threads=1):

    # Skip if JSON file already exists
    gzip_file = os.path.join(out_dir, os.pardir, "overlaps.json.gz")
    if not os.path.exists(gzip_file):

        # Initialize
        overlaps = {}

        # Move to output directory
        os.chdir(out_dir)

        # Get TFBS predictions
        global predictions
        predictions = __get_TFBS_predictions(
            __get_random_sequences(out_dir, threads), out_dir, threads
        )

        # Get TFBS overlaps
        overlaps = __get_TFBS_overlaps(out_dir, threads)

        # Write
        Jglobals.write(
            gzip_file[:-3],
            json.dumps(overlaps, sort_keys=True, indent=4)
        )
        fi = Jglobals._get_file_handle(gzip_file[:-3], "rb")
        fo = Jglobals._get_file_handle(gzip_file, "wb")
        shutil.copyfileobj(fi, fo)
        fi.close()
        fo.close()
        os.remove(gzip_file[:-3])

        # Change dir
        os.chdir(cwd)

def __get_random_sequences(out_dir=out_dir, threads=1):

    # Skip if FASTA file already exists
    gzip_file = ".random_seqs.fa.gz"
    if not os.path.exists(gzip_file):

        # Initialize
        records = []
        random_sequences = set()
        kwargs = {
            "total": 101,
            "ncols": 100
        }

        # STDOUT
        Jglobals.write(None, "\n*** Create random sequences")

        # Get 1M random sequences
        pool = Pool(threads)
        for sequences in tqdm(
            pool.imap(__create_random_sequences, range(101)), **kwargs
        ):
            total_sequences = len(random_sequences)
            while len(random_sequences) < total_sequences + 9901:
                random_sequences.add(sequences.pop(0))
        pool.close()
        pool.join()

        # For each sequence...
        for i, s in enumerate(random_sequences):
            record = SeqRecord(Seq(s), id=str(i + 1), name="", description="")
            records.append(record)

        # Write
        handle = Jglobals._get_file_handle(gzip_file[:-3], "wt")
        SeqIO.write(sequences=records, handle=handle, format="fasta")
        handle.close()
        process = subprocess.run(["gzip %s" % gzip_file[:-3]], shell=True)

    return(gzip_file)

def __create_random_sequences(gc_content=0, length=30, total_sequences=4**7):

    # Initialize
    random.seed(gc_content)
    records = []
    sequences = set()

    # Compute all possible cumulative weights
    cum_weights = []
    for i in range(100 - gc_content + 1):
        for j in range(gc_content + 1):
            A=i; T=(100-gc_content-A); C=j; G=(gc_content-j)
            cum_weights.append([A, A+C, A+C+G, A+C+G+T])

    # While not enough sequences...
    while len(sequences) < total_sequences:

        # Get sequence
        sequence = random.choices(
            list("ACGT"), cum_weights=random.choice(cum_weights), k=length
        )

        # Add sequence
        sequences.add("".join(sequence))

    return(list(sequences))

    # # For each sequence...
    # for i, sequence in enumerate(sequences):
    #     record = SeqRecord(
    #         Seq(sequence),
    #         id=str(i + (total_sequences * gc_content)),
    #         name="",
    #         description="")
    #     records.append(record)

    # return(records)

def __get_TFBS_predictions(fasta_file, out_dir=out_dir, threads=1):

    # Initialize
    predictions = {}

    # Get predictions
    d = "./predictions/"
    if not os.path.isdir(d):

        # STDOUT
        Jglobals.write(None, "\n*** Scan JASPAR profiles")

        # Create FASTA file
        tmp_file = "tmp.fa"
        cmd = "zless %s > %s" % (fasta_file, tmp_file)
        process = subprocess.run([cmd], shell=True)

        # Scan sequences
        cmd = """
            %s/scan_sequence.py \
            --fasta-file %s \
            --profiles-dir %s/profiles/ \
            --output-dir %s \
            --threads %s
        """ % (jaspar_ucsc_tracks, tmp_file, jaspar_ucsc_tracks, d, threads)
        process = subprocess.run([cmd], shell=True)
        os.remove(tmp_file)

    # STDOUT
    Jglobals.write(None, "\n*** Load TFBS predictions")

    # Load predictions files
    files = []
    for predictions_file in os.listdir(d):
        files.append(os.path.join(d, predictions_file))
    kwargs = {"total": len(files), "ncols": 100}
    pool = Pool(threads)
    for matrix_id, p in tqdm(
        pool.imap(__load_TFBS_predictions, files), **kwargs
    ):
        predictions.setdefault(matrix_id, p)
    pool.close()
    pool.join()

    return(predictions)

def __load_TFBS_predictions(scans_file):

    # Initialize
    predictions = []
    matrix_id = re.search("(MA\d{4}\.\d)", scans_file)

    # Load TFBS predictions
    try:
        df = pd.read_csv(scans_file, sep="\t", header=None)
        df = df.sort_values(by=[4], ascending=False)
        return(
            matrix_id.group(1), np.array(pd.Series(df[0].to_list()).unique())
        )
    except:
        return(matrix_id.group(1), None)

def __get_TFBS_overlaps(out_dir=out_dir, threads=1):

    # Initialize
    overlaps = {}
    matrix_ids = list(predictions.keys())
    kwargs = {
        "total": len(matrix_ids),
        "ncols": 100
    }

    # Create directory
    d = "./overlaps/"
    if not os.path.isdir(d):
        os.makedirs(d)

    # STDOUT
    Jglobals.write(None, "\n*** Load TFBS overlaps")

    # Get overlaps
    pool = Pool(threads)
    p = partial(__compute_TFBS_overlaps, overlaps_dir=d)
    for matrix_id, o in tqdm(pool.imap(p, matrix_ids), **kwargs):
        overlaps.setdefault(matrix_id, o)
    pool.close()
    pool.join()

    return(overlaps)

def __compute_TFBS_overlaps(matrix_id, overlaps_dir):

    # Skip if JSON file already exists
    gzip_file = os.path.join("./overlaps/", "%s.json.gz" % matrix_id)
    if not os.path.exists(gzip_file):

        # Initialize
        overlaps = {}
        ar1 = predictions[matrix_id]

        # For each other matrix...
        for next_matrix_id in predictions:

            # Compute overlap
            ar2 = predictions[next_matrix_id]
            if ar1 is None or ar2 is None:
                overlaps.setdefault(next_matrix_id, 0.)
            else:
                intersect = np.intersect1d(ar1[:len(ar2)], ar2[:len(ar1)])
                union = np.union1d(ar1[:len(ar2)], ar2[:len(ar1)])
                overlap = round(intersect.size / float(union.size), 3)
                overlaps.setdefault(next_matrix_id, overlap)

        # Write
        Jglobals.write(
            gzip_file[:-3],
            json.dumps(overlaps, sort_keys=True, indent=4)
        )
        fi = Jglobals._get_file_handle(gzip_file[:-3], "rb")
        fo = Jglobals._get_file_handle(gzip_file, "wb")
        shutil.copyfileobj(fi, fo)
        fi.close()
        fo.close()
        os.remove(gzip_file[:-3])

    # Load JSON file
    handle = Jglobals._get_file_handle(gzip_file)
    overlaps = json.load(handle)

    return(matrix_id, overlaps)

def __TFBS_in_sequence(motif, seq, threshold):

    # Search
    for position, score in motif.pssm.search(seq, threshold):
        return(True)

    return(False)

def __get_cosine_similarities(vectors, out_dir=out_dir):

    # Skip if JSON file already exists
    gzip_file = ".cosine_similarities.json.gz"
    if not os.path.exists(gzip_file):

        # Initialize
        data = []
        matrix_ids = sorted(vectors.keys())

        # For each matrix...
        for i in range(len(matrix_ids) - 1):

            # Test CTCF
            if matrix_ids[i] != "MA0139.1":
                continue

            # Initialize
            u = vectors[matrix_ids[i]]
            print(sum(u))

            # For each other matrix...
            for j in range(i + 1, len(matrix_ids)):

                # Add cosine distance
                v = vectors[matrix_ids[j]]
                print(sum(v))
                cosine_distance = round(1 - distance.cosine(u, v), 3)
                data.append([matrix_ids[i], matrix_ids[j], cosine_distance])
                data.append([matrix_ids[j], matrix_ids[i], cosine_distance])

            data.sort(key=lambda x: x[-1], reverse=True)
            print(data[0])
            exit(0)

def __group_by_cluster(files_dir=files_dir, out_dir=out_dir):

    # Skip if groups JSON file already exists
    gzip_file = os.path.join(out_dir, "groups.clusters.json.gz")
    if not os.path.exists(gzip_file):

        # Initialize
        clusters = {}
        motifs = {}

        # Move to output directory
        os.chdir(out_dir)

        # Load JSON files
        json_file = os.path.join(out_dir, "groups.DBDs.json.gz")
        handle = Jglobals._get_file_handle(json_file)
        DBDs = json.load(handle)
        handle.close()

        # For each taxon...
        for taxon in Jglobals.taxons:

            # Load JSON files
            json_file = os.path.join(files_dir, "%s.uniprot.json" % taxon)
            with open(json_file) as f:
                uniaccs = json.load(f)

            # For each UniProt Accession...
            for uniacc, values in uniaccs.items():
                motifs.setdefault(uniacc, [taxon, values[0]])

        # For each DBD...
        for DBD in sorted(DBDs):

            # # Test Forkhead
            # if DBD != "Forkhead":
            #     continue

            # Get clusters
            uniaccs = list(DBDs[DBD].keys())
            clusters = __get_DBD_clusters(DBD, uniaccs, motifs, files_dir)
        exit(0)

        #         # For each matrix ID...
        #         for matrix_id, values in clusters.items():

        #             for value in values:

        #                 groups.setdefault(matrix_id, [])
        #                 groups[matrix_id].append((taxon, value))

        # # Write
        # Jglobals.write(
        #     gzip_file[:-3],
        #     json.dumps(groups, sort_keys=True, indent=4)
        # )
        # fi = Jglobals._get_file_handle(gzip_file[:-3], "rb")
        # fo = Jglobals._get_file_handle(gzip_file, "wb")
        # shutil.copyfileobj(fi, fo)
        # fi.close()
        # fo.close()
        # os.remove(gzip_file[:-3])

        # # Change dir
        # os.chdir(cwd)

def __get_DBD_clusters(DBD, uniaccs, motifs, files_dir=files_dir):

    # Initialize
    matrix_ids = set()

    # Create output dir
    if not os.path.isdir(DBD):
        os.makedirs(DBD)

    # Skip if already done
    jaspar_profiles = os.path.abspath(os.path.join(DBD, "%s.jaspar" % DBD))
    if not os.path.exists(jaspar_profiles):

        # Initialize
        fo = Jglobals._get_file_handle(jaspar_profiles, "a")

        # For each UniProt Accession...
        for uniacc in uniaccs:

            # For each motif...
            for matrix_id in motifs[uniacc][1]:

                # Skip
                if matrix_id in matrix_ids:
                    continue

                # Write
                file_name = os.path.join(
                    files_dir, motifs[uniacc][0], "%s.jaspar" % matrix_id
                )
                fi = Jglobals._get_file_handle(file_name, "r")
                fo.write(fi.read())
                fi.close()

                # Done
                matrix_ids.add(matrix_id)

    # Skip if already done
    leafs_file = os.path.join(DBD, "%s_tables" % DBD, "leaf_to_cluster.tab")
    if not os.path.exists(leafs_file):

        # RSAT matrix-clustering
        cmd = """
            RSAT matrix-clustering -v 2 \
            -matrix %s %s jaspar \
            -hclust_method average \
            -calc sum \
            -title %s \
            -metric_build_tree Ncor \
            -lth w 5 -lth cor 0.6 -lth Ncor 0.4 \
            -label_in_tree name \
            -return json \
            -quick \
            -radial_tree_only \
            -o %s
        """ % (DBD, jaspar_profiles, DBD, DBD)
        process = subprocess.run([cmd], shell=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, cwd=DBD)

        # except:
        #     # i.e. to avoid crushing due to error
        #     pass

# def _group_by_Tomtom(out_dir=out_dir, threads=1):

#     # Skip if groups JSON file already exists
#     gzip_file = os.path.join(out_dir, "groups.tomtom.json.gz")
#     if not os.path.exists(gzip_file):

#         # Initialize
#         tomtom = {}

#         # Get all JASPAR profiles
#         jaspar_profiles = _get_profiles_from_latest_version(Path(out_dir).glob("*/*.meme"))
#         # jaspar_profiles = [str(f) for f in Path(out_dir).glob("*/*.meme")]

#         # Skip if JASPAR MEME database already exists
#         database = os.path.join(out_dir, "jaspar.meme")
#         if not os.path.exists(database):

#             # For each JASPAR profile...
#             for jaspar_profile in jaspar_profiles:

#                 # Cat to database
#                 os.system("cat %s >> %s" % (jaspar_profile, database))

#         # Parallelize
#         pool = Pool(threads)
#         parallelized = partial(Tomtom, database=database, out_dir=out_dir)
#         for _ in tqdm(
#             pool.imap(parallelized, jaspar_profiles), desc="Tomtom", total=len(jaspar_profiles)
#         ):
#             pass
#         pool.close()
#         pool.join()

#         # Move to output directory
#         os.chdir(out_dir)

#         # For each JASPAR profile...
#         for jaspar_profile in jaspar_profiles:

#             # Initialize
#             m = re.search("(MA\d{4}.\d).meme$", jaspar_profile)
#             tomtom_dir = ".%s" % m.group(1)

#             # Get hits
#             tomtom.setdefault(m.group(1), _get_Tomtom_hits(tomtom_dir))

#             # Remove Tomtom directory
#             shutil.rmtree(tomtom_dir)

#         # Write
#         Jglobals.write(
#             gzip_file[:-3],
#             json.dumps(tomtom, sort_keys=True, indent=4)
#         )
#         fi = Jglobals._get_file_handle(gzip_file[:-3], "rb")
#         fo = Jglobals._get_file_handle(gzip_file, "wb")
#         shutil.copyfileobj(fi, fo)
#         fi.close()
#         fo.close()
#         os.remove(gzip_file[:-3])

#         # For each taxon...
#         for taxon in Jglobals.taxons:

#             # Remove taxon directory
#             if os.path.isdir(taxon):
#                 shutil.rmtree(taxon)

#         # Change dir
#         os.chdir(cwd)

# def Tomtom(meme_file, database, out_dir=out_dir):
#     """
#     From http://meme-suite.org/doc/tomtom.html;
#     In order to compute the scores, Tomtom needs to know the frequencies of
#     the letters of the sequence alphabet in the database being searched (the
#     "background" letter frequencies). By default, the background letter fre-
#     quencies included in the query motif file are used. The scores of columns
#     that overlap for a given offset are summed. This summed score is then con-
#     verted to a p-value. The reported p-value is the minimal p-value over all
#     possible offsets. To compensate for multiple testing, each reported p-value
#     is converted to an E-value by multiplying it by twice the number of target
#     motifs. As a second type of multiple-testing correction, q-values for each
#     match arecomputed from the set of p-values and reported.
#     """

#     # Skip if output directory already exists
#     m = re.search("(MA\d{4}.\d|M\d{4}_1.02).meme$", meme_file)
#     output_dir = os.path.join(out_dir, ".%s" % m.group(1))
#     if not os.path.isdir(output_dir):

#         # Run Tomtom
#         cmd = "tomtom -o %s -thresh 10000 -evalue %s %s" % (output_dir, meme_file, database)
#         process = subprocess.run([cmd], shell=True, stdout=subprocess.DEVNULL,
#             stderr=subprocess.DEVNULL)

# def _get_profiles_from_latest_version(jaspar_profiles):

#     # Initialize
#     done = set()
#     latest_version_profiles = []

#     # For each profile...
#     for jaspar_profile in sorted(jaspar_profiles, reverse=True):

#         # Initialize
#         m = re.search("(MA\d{4}).\d.\S+$", str(jaspar_profile))
#         matrix_id = m.group(1)

#         # Skip if done
#         if matrix_id in done:
#             continue

#         # i.e. a profile from the latest version
#         latest_version_profiles.append(str(jaspar_profile))

#         # Done
#         done.add(matrix_id)

#     return(latest_version_profiles)

# def _get_Tomtom_hits(tomtom_dir):

#     # Intialize
#     hits = []

#     # For each line...
#     for line in Jglobals.parse_tsv_file(os.path.join(tomtom_dir, "tomtom.tsv")):

#         # Skip comments
#         if line[0].startswith("#"):
#             continue

#         # Skip header
#         if line[0] == "Query_ID":
#             continue

#         # Skip self
#         if line[0][:6] == line[1][:6]:
#             continue

#         # Add to hits
#         hits.append([line[1], float(line[4])])

#     return(hits)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()