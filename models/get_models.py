#!/usr/bin/env python

import argparse
from collections import Counter
import copy
from glmnet import LogitNet
from itertools import chain
import json
import numpy as np
import os
import pickle
import shutil
from sklearn.metrics import precision_recall_curve
import sys

# Defaults
out_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)))
root_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir)
files_dir = os.path.join(root_dir, "files")
data_splits_file = os.path.join(out_dir, "data_splits.json.gz")

# Append JASPAR-profile-inference to path
sys.path.append(root_dir)

# Import globals
from __init__ import Jglobals

#-------------#
# Classes     #
#-------------#

class LeaveOneTfOut:

    def __init__(self, n_splits=3, shuffle=False, random_state=None):
        self.n_splits = n_splits

    def split(self, X=None, y=None, groups=None):

        for train, test in iterator:
            yield(train, test)

    def get_n_splits(self, X=None, y=None, groups=None):
        return self.n_splits

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns
    an {argparse} object.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("--data-splits", default=data_splits_file,
        metavar="FILE", help="compressed json file (from get_data_splits.py)")
    parser.add_argument("--files-dir", default=files_dir, metavar="DIR",
        help="files directory (from get_files.py)")
    parser.add_argument("-o", default=out_dir, metavar="DIR",
        help="output directory (default = ./)")
    parser.add_argument("-v", "--verbose", action="store_true",
        help="verbose mode (default = False)")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Globals
    global cwd
    cwd = os.getcwd()
    out_dir = os.path.abspath(args.o)
    global data_splits
    handle = Jglobals._get_file_handle(os.path.abspath(args.data_splits))
    data_splits = json.load(handle)
    handle.close()
    global lambdas
    lambdas = 10**np.linspace(6, -6, 100)
    global verbose
    verbose = args.verbose
    global cisbp
    cisbp = {}
    cisbp_dir = os.path.join(files_dir, "cisbp")
    for json_file in os.listdir(cisbp_dir):
        handle = Jglobals._get_file_handle(os.path.join(cisbp_dir, json_file))
        dictionary = json.load(handle)
        cisbp.setdefault(dictionary["Family_Name"], dictionary)
        handle.close()

    # Create output dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Get models
    get_models(out_dir)

def get_models(out_dir=out_dir):

    # Skip if JSON file already exists
    json_file = os.path.join(out_dir, "models.json")
    if not os.path.exists(json_file):

        # Initialize
        models = {
            "Keys": "DBD",
            "Values": {
                "Keys": "Similarity",
                "Values": [
                    "Coefficients",
                    "Best lambda",
                    "Threshold @ 75% Precision",
                    "Recall @ 75% Precision",
                    "Precision",
                    "Recall",
                    "Thresholds"
                ]
            }
        }

        # Create data splits dir
        if not os.path.isdir(os.path.join(out_dir, "models")):
            os.makedirs(os.path.join(out_dir, "models"))

        # Move to data splits directory
        os.chdir(os.path.join(out_dir, "models"))

        # For each DBD composition...
        for DBD in data_splits:

            # Verbose mode
            if verbose:
                Jglobals.write(None, "\n%s..." % DBD)

            # # Skip if pickle file already exists
            # pickle_file = "%s.pickle" % DBD
            # if not os.path.exists(pickle_file):

            #     # Initialize
            #     pickles = {}
            #     models.setdefault(DBD, {})

            # For each sequence similarity representation...
            for similarity in ["identity", "blosum62"]:

                # Train model
                m, coefficients, lamdabest = __train_model(DBD, similarity)
                # pickles.setdefault(similarity, m)

                # Compute statistics
                statistics = __compute_statistics(DBD, similarity, m)
                # models[DBD].setdefault(similarity, [])
                # models[DBD][similarity].append([
                #     coefficients,
                #     lamdabest,
                #     statistics["Threshold @ 75% Precision"],
                #     statistics["Recall @ 75% Precision"],
                #     statistics["Precision"],
                #     statistics["Recall"],
                #     statistics["Thresholds"]
                # ])

            #     # Write
            #     handle = Jglobals._get_file_handle(pickle_file, "wb")
            #     pickle.dump(pickles, handle)
            #     handle.close()

            # Compute Cis-BP statistics
            # statistics = __compute_CisBP_statistics(DBD)
        #     models[DBD].setdefault("cisbp", [])
        #     models[DBD][similarity].append([
        #         coefficients,
        #         lamdabest,
        #         statistics["Cis-BP Threshold"],
        #         statistics["Precision @ Cis-BP Threshold"],
        #         statistics["Recall @ Cis-BP Threshold"],
        #         statistics["Precision"],
        #         statistics["Recall"]
        #     ])

        #     exit(0)

        # # Write
        # Jglobals.write(
        #     json_file,
        #     json.dumps(models, sort_keys=True, indent=4)
        # )

def __train_model(DBD, similarity):

    # Verbose mode
    if verbose:
        Jglobals.write(None, "\t*** similarity: %s" % similarity)

    # Get pairs, Xs, ys
    pairs, Xs, ys = __get_pairs_Xs_ys(DBD, similarity, "training")

    try:

        # Get iterator for cross-validation
        global iterator
        iterator = __leave_one_TF_out(pairs)

        # Get Xs, ys, weights
        Xs_train = np.asarray(Xs)
        ys_train = np.asarray(ys["stringent"])
        weights = np.asarray(__get_weights(ys["stringent"]))

        # Initialize model
        m = LogitNet(
            alpha=0, n_splits=len(set(list(chain.from_iterable(pairs)))),
            lambda_path=lambdas, lower_limits=0
        )

        # Set custom cross-validation
        m.CV = LeaveOneTfOut

        # Fit
        m_fit = m.fit(Xs_train, ys_train.ravel(), sample_weight=weights)

        # Verbose mode
        if verbose:
            Jglobals.write(None, "\t*** LogitNet training: SUCCESS!")

        return(m_fit, list(m_fit.coef_[0]), m_fit.lambda_max_)

    except:

        # Verbose mode
        if verbose:
            Jglobals.write(None, "\t*** LogitNet training: FAIL!")

        return(None, None, None) 

def __get_pairs_Xs_ys(DBD, similarity, data_split="train"):

    # Initialize
    pairs = []
    Xs = []
    ys = {}

    # For each criterion...
    for criterion in ["stringent", "lenient"]:
        ys.setdefault(criterion, [])

    # For each pair, X, y...
    for TFa, TFb, X, y, s in data_splits[DBD][data_split]:

        # Skip
        if s != similarity:
            continue

        # Add pair, X, y
        pairs.append([TFa, TFb])
        Xs.append(X)
        ys["stringent"].append(y[0])
        ys["lenient"].append(y[1])

    return(pairs, Xs, ys)

def __leave_one_TF_out(pairs):
    """
    Leave one TF out.
    """

    # Initialize
    d = {}
    iterator = []

    # 1. Store TFs and their indices into a dict
    for pair in pairs:
        for TF in pair:
            d.setdefault(TF, [])
    for TF, idxs in d.items():
        for idx in range(len(pairs)):
            if TF in pairs[idx]:
                idxs.append(idx)

    # 2. For each TF, figure out which indices to include
    idxs = set(list(range(len(pairs))))
    for TF in d:
        iterator.append((np.array(list(idxs - set(d[TF]))), np.array(d[TF])))

    return(iterator)

def __get_weights(ys):
    """
    Weight samples 1/freq.
    """

    # Initialize
    weights = {}

    for k, v in Counter(list(chain.from_iterable(ys))).items():
        weights.setdefault(k, 1/(float(v)/len(ys)))

    return([weights[y[0]] for y in ys])

def __compute_statistics(DBD, similarity, m=None):

    # Initialize
    statistics = {
        "Threshold @ 75% Precision": None,
        "Recall @ 75% Precision": None,
        "Precision": None,
        "Recall": None,
        "Thresholds": None,
    }

    if m is not None:

        # Get pairs, Xs, ys
        pairs, Xs, ys = __get_pairs_Xs_ys(DBD, similarity, "test")

        # Get best lambda (i.e. max)
        lambdabest = m.lambda_max_
        if verbose:
            Jglobals.write(None, "\t*** best lambda: %s" % round(lambdabest, 8))

        # Get Xs, ys
        Xs_test = np.asarray(Xs)
        ys_test = np.asarray(ys["stringent"])

        # Predict
        p = m.predict_proba(Xs_test, lamb=lambdabest)
        p_transposed = p.transpose()

        # Statistics
        precision, recall, thresholds = precision_recall_curve(
            ys_test.ravel(), p_transposed[1]
        )
        statistics["Precision"] = list(precision)
        statistics["Recall"] = list(recall)
        statistics["Thresholds"] = list(thresholds)
        for idx in range(len(thresholds)):
            if precision[idx] < 0.75:
                continue
            statistics["Threshold @ 75% Precision"] = thresholds[idx]
            statistics["Recall @ 75% Precision"] = recall[idx]
            break

        # Verbose
        if verbose:
            if statistics["Threshold @ 75% Precision"] is not None:
                recall = round(statistics["Threshold @ 75% Precision"] * 100, 2)
                Jglobals.write(
                    None,
                    "\t*** Recall @ 75% Precision: {0:.2f}%".format(recall)
                )
            else:
                Jglobals.write(None, "\t*** Recall @ 75% Precision: FAIL!")

    return(statistics)

def __compute_CisBP_statistics(DBD):

    # Initialize
    statistics = {
        "sequence identity": {
            "Threshold": None,
            "Precision @ Threshold": None,
            "Recall @ Threshold": None,
            "Precision": None,
            "Recall": None,
            "Thresholds": None
        },
        "similarity regression": {
            "Threshold": None,
            "Precision @ Threshold": None,
            "Recall @ Threshold": None,
            "Precision": None,
            "Recall": None,
            "Thresholds": None
        }
    }

    # Get pairs, Xs, ys
    pairs, Xs, ys = __get_pairs_Xs_ys(DBD, "identity", "test")

    # Get Xs, ys
    Xs_test = np.asarray(Xs)
    ys_test = np.asarray(ys["stringent"])
    ys_true = float(sum(ys_test.ravel()))

    for k in statistics.keys():

        # Verbose mode
        if verbose:
            Jglobals.write(None, "\t*** Cis-BP: %s" % k)

        # Initialize
        precision = []
        recall = []

        if k == "sequence identity":

            # Initialize
            threshold = __get_CisBP_sequence_identity_threshold(DBD)
            statistics[k]["Threshold"] = threshold

            # Get Xs
            Xs = np.array([float(sum(X)) / len(X) for X in Xs_test])

        else:

            continue

            # Initialize
            threshold, weights = __get_CisBP_similarity_regression_weights(DBD)
            statistics[k]["Threshold"] = threshold

            # Get Xs
            Xs = np.array([sum(X * weights) for X in Xs_test])

        # Get thresholds
        thresholds = sorted(set(Xs))

        # For each threshold...
        for t in thresholds:

            # Get statistics
            idxs = np.where(Xs >= t)
            p = (Xs[idxs] > statistics[k]["Threshold"]).astype(int)
            ys = ys_test[idxs].ravel().astype(int)
            precision.append(sum(p/len(p)))
            recall.append(sum(p/ys_true))

        # Statistics
        statistics[k]["Precision"] = precision
        statistics[k]["Precision"].append(1.)
        statistics[k]["Recall"] = recall
        statistics[k]["Recall"].append(0.)
        statistics[k]["Thresholds"] = thresholds
        for idx in range(len(thresholds)):
            if thresholds[idx] < statistics[k]["Threshold"]:
                continue
            statistics[k]["Precision @ Threshold"] = precision[idx]
            statistics[k]["Recall @ Threshold"] = recall[idx]
            break

        # Verbose
        if verbose:
            Jglobals.write(
                None, "\t*** Threshold: {0:.2f}".format(
                    statistics[k]["Threshold"]
                )
            )
            if statistics[k]["Precision @ Threshold"] is not None:
                P = round(statistics[k]["Precision @ Threshold"] * 100, 2)
                Jglobals.write(
                    None, "\t*** Precision @ Threshold: {0:.2f}%".format(P)
                )
                R = round(statistics[k]["Recall @ Threshold"] * 100, 2)
                Jglobals.write(
                    None, "\t*** Recall @ Threshold: {0:.2f}%".format(R)
                )
            else:
                Jglobals.write(None, "\t*** Precision @ Threshold: FAIL!")
                Jglobals.write(None, "\t*** Recall @ Threshold: FAIL!")

    return(statistics)

def __get_CisBP_sequence_identity_threshold(DBD):

    # Initialize
    domains = {
        "B3": None,
        "CP2": None,
        "CUT": None,
        "E2F_TDP": "E2F",
        "HLH": "bHLH",
        "HMG_box": "Sox",
        "HSF_DNA-bind": "HSF",
        "IRF": None,
        "MADF_DNA_bdg": "MADF",
        "MH1": None,
        "Myb_DNA-binding": "Myb/SANT",
        "NAM": "NAC/NAM",
        "PAX": None,
        "Pou": "Homeodomain,POU",
        "RFX_DNA_binding": "RFX",
        "RHD_DNA_bind": None,
        "SRF-TF": None,
        "STAT_bind": None,
        "TCR": "TCR/CxC",
        "TEA": None,
        "Zn_clus": "Zinc cluster",
        "bZIP_1": "bZIP",
        "zf-C2H2": "C2H2 ZF",
        "zf-C4": None,
        "zf-Dof": "Dof"
    }

    if DBD in cisbp:
        domain = DBD
    elif DBD in domains:
        if domains[DBD] in cisbp:
            domain = domains[DBD]
        else:
            domain = "NO_THRESHOLD"

    if "Baseline" in cisbp[domain]:
        return(cisbp[domain]["Baseline"]["Threshold.HSim"])

    return(cisbp[domain]["Threshold.HSim"])

def __get_CisBP_similarity_regression_weights(DBD):

    # Initialize
    domains = {
        "B3": None,
        "CP2": None,
        "CUT": None,
        "E2F_TDP": "E2F",
        "HLH": "bHLH",
        "HMG_box": "Sox",
        "HSF_DNA-bind": "HSF",
        "IRF": None,
        "MADF_DNA_bdg": "MADF",
        "MH1": None,
        "Myb_DNA-binding": "Myb/SANT",
        "NAM": "NAC/NAM",
        "PAX": None,
        "Pou": "Homeodomain,POU",
        "RFX_DNA_binding": "RFX",
        "RHD_DNA_bind": None,
        "SRF-TF": None,
        "STAT_bind": None,
        "TCR": "TCR/CxC",
        "TEA": None,
        "Zn_clus": "Zinc cluster",
        "bZIP_1": "bZIP",
        "zf-C2H2": "C2H2 ZF",
        "zf-C4": None,
        "zf-Dof": "Dof"
    }

    if DBD in cisbp:
        domain = DBD
    elif DBD in domains:
        if domains[DBD] in cisbp:
            domain = domains[DBD]
        else:
            domain = "NO_THRESHOLD"

    if "SR.Weights" in cisbp[domain]:
        return(
            cisbp[domain]["Threshold.HSim"],
            np.array(cisbp[domain]["SR.Weights"])
        )

    return(None, None)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()