#!/usr/bin/env python

import argparse
from Bio.SubsMat.MatrixInfo import blosum62
import json
import math
import numpy as np
from numpy import log10 as log
import os
import shutil
import string
import sys

# Defaults
out_dir = os.path.dirname(os.path.realpath(__file__))
root_dir = os.path.join(out_dir, os.pardir)
files_dir = os.path.join(root_dir, "files")

# Append JASPAR-profile-inference to path
sys.path.append(root_dir)

# Import globals
from __init__ import Jglobals
from infer_profile import (
    _is_alignment_over_Rost_seq_id_curve,
    _is_alignment_over_Rost_seq_sim_curve
)

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns
    an {argparse} object.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-f", default=files_dir, metavar="DIR",
        help="files directory (from get_files.py; default=../files/)")
    parser.add_argument("-o", default=out_dir, metavar="DIR",
        help="output directory (default = ./)")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Pairwise
    pairwise(os.path.abspath(args.f), os.path.abspath(args.o))

def pairwise(files_dir=files_dir, out_dir=out_dir):

    # Skip if pairwise JSON file already exists
    gzip_file = os.path.join(out_dir, "pairwise.json.gz")
    if not os.path.exists(gzip_file):

        # Initialize
        pairwise = {}

        # Create output dir
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # Get groups by DBD composition
        groups = _get_DBD_groups(files_dir)

        # Get Tomtom groups
        global tomtom
        tomtom = _get_Tomtom_groups(files_dir)

        # # Get BLAST+ groups
        # global blast
        # blast = _get_BLAST_groups(files_dir)

        # For each key, values...
        for key, values in groups.items():

            # Initialize
            Xss = {}
            # BLASTXss = {}
            Ys = []
            TFpairs = []

            # For each TF...
            for i in range(len(values) - 1):

                # For each next TF...
                for j in range(i + 1, len(values)):

                    # Get Xs
                    identityXs, blosum62Xs = _fetchXs(values[i][1], values[j][1])
                    Xss.setdefault("identity", [])
                    Xss.setdefault("blosum62", [])
                    if identityXs is not None and blosum62Xs is not None:
                        Xss["identity"].append(identityXs)
                        Xss["blosum62"].append(blosum62Xs)
                    else:
                        # Skip this pair!
                        continue

                    # # Get BLAST+ Xs
                    # identityBLASTXs, blosum62BLASTXs = _fetchBLASTXs(values[i][2], values[j][2])
                    # BLASTXss.setdefault("identity", [])
                    # BLASTXss.setdefault("blosum62", [])
                    # BLASTXss["identity"].append(identityBLASTXs)
                    # BLASTXss["blosum62"].append(blosum62BLASTXs)

                    # Get Y
                    y = _fetchY(values[i][0], values[j][0])
                    Ys.append(y)

                    # Get TF pair
                    TFpairs.append((values[i][2], values[j][2]))

            # Add to pairwise
            # pairwise.setdefault(key, [Xss, BLASTXss, Ys, TFpairs])
            pairwise.setdefault(key, [Xss, Ys, TFpairs])

        # Write
        Jglobals.write(
            gzip_file[:-3],
            json.dumps(pairwise, sort_keys=True, indent=4, separators=(",", ": "))
        )
        fi = Jglobals._get_file_handle(gzip_file[:-3], "rb")
        fo = Jglobals._get_file_handle(gzip_file, "wb")
        shutil.copyfileobj(fi, fo)
        fi.close()
        fo.close()
        os.remove(gzip_file[:-3])

def _get_DBD_groups(files_dir=files_dir):

    # Load JSON file
    json_file = os.path.join(files_dir, "groups.DBDs.json.gz")
    handle = Jglobals._get_file_handle(json_file)
    groups = json.load(handle)
    handle.close()

    return(groups)

def _get_Tomtom_groups(files_dir=files_dir):

    # Initialize
    tomtom_filtered = {}

    # Load JSON file
    json_file = os.path.join(files_dir, "groups.tomtom.json.gz")
    handle = Jglobals._get_file_handle(json_file)
    tomtom_unfiltered = json.load(handle)
    handle.close()

    for matrix_id in tomtom_unfiltered:

        for m, e in tomtom_unfiltered[matrix_id]:

            tomtom_filtered.setdefault(matrix_id, {})
            tomtom_filtered[matrix_id].setdefault(m, e)

    return(tomtom_filtered)

# def _get_BLAST_groups(files_dir=files_dir):

#     # Initialize
#     blast_filtered = {}

#     # Load JSON file
#     json_file = os.path.join(files_dir, "groups.blast.json.gz")
#     handle = Jglobals._get_file_handle(json_file)
#     blast_unfiltered = json.load(handle)
#     handle.close()

#     for uniacc in blast_unfiltered:

#         for t, qse, tse, e, s, pid, al, psim, jc in blast_unfiltered[uniacc]:

#             # (t) identifier of target sequence;
#             # (qse, tse) start and end-position in query and in target;
#             # (e) E-value;
#             # (s) bit score;
#             # (pid) percentage of identical matches;
#             # (al) alignment length;
#             # (psim) percentage of positive-scoring matches; and
#             # (jc) joint coverage (i.e. square root of the coverage
#             # on the query and the target).
#             # s = s / al # transformation of bit score

#             # Get Rost's verdict
#             rost_seq_id = _is_alignment_over_Rost_seq_id_curve(pid, al)
#             rost_seq_sim = _is_alignment_over_Rost_seq_sim_curve(psim, al)

#             # Reformat percentages & coverage
#             pid = pid / 100.0
#             psim = psim / 100.0
#             jc = jc / 100.0

#             blast_filtered.setdefault(uniacc, {})
#             blast_filtered[uniacc].setdefault(t, [[pid, rost_seq_id], [psim, rost_seq_sim], [jc]])

#     return(blast_filtered)

def _fetchXs(seqs1, seqs2):

    # Initialize
    permA = []
    permB = []
    scores = []

    # Find the max. consecutive domains (i.e. not None)
    max1 = _get_max_consecutive_not_Nones(seqs1)
    max2 = _get_max_consecutive_not_Nones(seqs2)

    if max1 < max2:
        A = seqs2
        B = seqs1
        step = max1
    else:
        A = seqs1
        B = seqs2
        step = max2

    # Avoid overhangs
    for a in range(len(A) - step + 1):
        permA.append(list(range(a, a+step)))
    for b in range(len(B) - step + 1):
        permB.append(list(range(b, b+step)))


    # For each sequence similarity representation...
    for similarity in ["identity", "blosum62"]:

        scores.append([])

        # Avoid overhangs
        for a in permA:

            for b in permB:

                # i.e. gap; skip
                if None in A[a[0]:a[-1]+1] or None in B[b[0]:b[-1]+1]:
                    continue

                scores[-1].append([])

                # For each DBD...
                for d in range(len(a)):

                    seqA = _removeLowercase(A[a[d]])
                    seqB = _removeLowercase(B[b[d]])

                    if similarity == "identity":
                        scores[-1][-1].append([0] * len(seqA))
                    else:
                        scores[-1][-1].append([-4] * len(seqA))

                    for n in range(len(seqA)):

                        if similarity == "identity":
                            scores[-1][-1][-1][n] = _IDscoring(seqA[n], seqB[n])
                        else:
                            scores[-1][-1][-1][n] = _BLOSUMscoring(seqA[n], seqB[n])

    # From Lambert et al.
    # For TF families that have DBDs present in arrays [...] the best ungapped
    # and overlapping pairwise alignment of DBD arrays is found by selecting
    # the alignment offset with the maximum amino acid identity.
    ix = None
    max_identities = 0
    for s in range(len(scores[0])):
        if sum([sum(a) for a in scores[0][s]]) > max_identities:
            ix = s
            max_identities = sum([sum(a) for a in scores[0][s]])

    if ix is not None:

        # From Lambert et al.
        # For a multi-DBD alignment, the feature vector is generated by the
        # average score (identity or similarity) in each position of the DBD
        # alignment from all DBD arrays, normalizing by the DBD length of the
        # longest protein.
        identityXs = [float(sum(z))/len(A) for z in zip(*scores[0][ix])]
        blosum62Xs = [float(sum(z))/len(A) for z in zip(*scores[1][ix])]
        
        return(identityXs, blosum62Xs)

    return(None, None)

def _get_max_consecutive_not_Nones(seq):

    # Intialize
    notNones = [0]    

    # For each domain...
    for d in seq:
        if d is None:
            notNones.append(0)
        else:
            notNones[-1] += 1

    return(max(notNones))

def _removeLowercase(s):

    return(s.translate(str.maketrans("", "", string.ascii_lowercase)))

def _IDscoring(aa1, aa2):

    if aa1 == aa2 and aa1 != "-":
        return(1)
    else:
        return(0)

def _BLOSUMscoring(aa1, aa2):

    if aa1 == "-" and aa2 == "-":
        return(1)
    elif aa1 == "-" or aa2 == "-":
        return(-4)
    else:
        if (aa1, aa2) in blosum62:
            return(blosum62[(aa1, aa2)])
        else:
            return(blosum62[(aa2, aa1)])

# def _fetchBLASTXs(uacc1, uacc2):

#     if uacc1 in blast:
#         if uacc2 in blast[uacc1]:
#             return(blast[uacc1][uacc2][0] + blast[uacc1][uacc2][2], blast[uacc1][uacc2][1] + blast[uacc1][uacc2][2])

#     return([0.0, False, 0.0], [0.0, False, 0.0])

def _fetchY(maIDlist1, maIDlist2):
    """
    Returns "0 - log10(Tomtom e-value)" as the dependent value (i.e. Y).
    """

    # Initialize
    y = 0 - log(1645.0)

    # Iterate through profiles...
    for maID1 in maIDlist1:

        # Skip matrix ID
        if maID1 not in tomtom:
            continue

        for maID2 in maIDlist2:

            # Skip matrix ID
            if maID2 not in tomtom:
                continue

            # If profiles were clustered together...
            if maID2 in tomtom[maID1] and maID1 in tomtom[maID2]:

                joint_evalue = math.sqrt(tomtom[maID1][maID2] * tomtom[maID2][maID1])
                log_evalue = log(joint_evalue)

                # return(max(min(1.0, 0 - log_evalue / 10), 0.0))

                if 0 - log_evalue > y:
                    y = 0 - log_evalue

    return(y)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()