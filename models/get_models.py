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

    parser.add_argument("--data-splits-file", default=data_splits_file,
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
    handle = Jglobals._get_file_handle(os.path.abspath(args.data_splits_file))
    data_splits = json.load(handle)
    handle.close()
    global lambdas
    lambdas = 10**np.linspace(6, -6, 100)
    global verbose
    verbose = args.verbose
    global domains
    domains = {}
    domains_file = os.path.abspath(os.path.join(files_dir, "pfam-DBDs.json"))
    handle = Jglobals._get_file_handle(domains_file)
    for key, values in json.load(handle).items():
        domains.setdefault(values[0], float(values[1]))
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
                    "Recall"
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

            # Skip if pickle file already exists
            pickle_file = "%s.pickle" % DBD
            if not os.path.exists(pickle_file):

                # Initialize
                pickles = {}
                models.setdefault(DBD, {})

                # For each sequence similarity representation...
                for similarity in ["identity", "blosum62"]:

                    # Train model
                    m, coefficients, lamdabest = __train_model(DBD, similarity)
                    pickles.setdefault(similarity, m)

                    # Compute statistics
                    statistics = __compute_statistics(DBD, similarity, m)
                    models[DBD].setdefault(similarity, [])
                    models[DBD][similarity].append([
                        coefficients,
                        lamdabest,
                        statistics["Threshold @ 75% Precision"],
                        statistics["Recall @ 75% Precision"],
                        statistics["Precision"],
                        statistics["Recall"]
                    ])

                # Write
                handle = Jglobals._get_file_handle(pickle_file, "wb")
                pickle.dump(pickles, handle)
                handle.close()

            # Compute Cis-BP statistics
            statistics = __compute_CisBP_statistics(DBD)
            models[DBD].setdefault("cisbp", [])
            models[DBD][similarity].append([
                coefficients,
                lamdabest,
                statistics["Threshold @ 75% Precision"],
                statistics["Recall @ 75% Precision"],
                statistics["Precision"],
                statistics["Recall"]
            ])

            exit(0)

        # Write
        Jglobals.write(
            json_file,
            json.dumps(models, sort_keys=True, indent=4)
        )

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
            Jglobals.write(None, "\t*** LogitNet training: success!")

        return(m_fit, list(m_fit.coef_[0]), m_fit.lambda_max_)

    except:

        # Verbose mode
        if verbose:
            Jglobals.write(None, "\t*** LogitNet training: fail!")

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
                Jglobals.write(None, "\t*** Recall @ 75% Precision: fail!")

    return(statistics)

def __compute_CisBP_statistics(DBD):

    # Initialize
    precision = []
    recall = []
    statistics = {
        "Cis-BP Threshold": None,
        "Precision @ Cis-BP Threshold": None,
        "Recall @ Cis-BP Threshold": None,
        "Precision": None,
        "Recall": None,
    }

    # Get pairs, Xs, ys
    pairs, Xs, ys = __get_pairs_Xs_ys(DBD, "identity", "test")

    # Get Xs, ys
    Xs_test = np.asarray(Xs)
    ys_test = np.asarray(ys["stringent"])

    # Get thresholds
    thresholds = sorted(set([float(sum(X)) / len(X) for X in Xs_test]))

    # For each threshold...
    for t in thresholds:

        # Initialize
        tp = 0; fp = 0

        # For each X, y...
        for X, y in zip(Xs_test, ys_test):

            # Initialize
            pid = float(sum(X)) / len(X)

            # If %ID is smaller than threshold = negative
            if pid < t:
                continue

            # If y is one = true
            if y[0] == 1: tp += 1
            else: fp += 1

        # Precision, recall
        precision.append(float(tp) / (tp + fp))
        recall.append(float(tp) / sum(ys_test.ravel()))

    # Statistics
    precision.append(1.)
    recall.append(0.)
    for idx in range(len(thresholds)):
        if thresholds[idx] < domains[DBD]:
            continue
        statistics["Cis-BP Threshold"] = domains[DBD]
        statistics["Precision @ Cis-BP Threshold"] = precision[idx]
        statistics["Recall @ Cis-BP Threshold"] = recall[idx]
        break

    # Verbose
    if verbose:
        threshold = statistics["Cis-BP Threshold"]
        Jglobals.write(
            None,
            "\t*** Cis-BP Threshold: {0:.2f}%".format(round(threshold * 100, 2))
        )
        P = round(statistics["Precision @ Cis-BP Threshold"] * 100, 2)
        Jglobals.write(
            None, "\t*** Precision @ Cis-BP Threshold: %s" % P
        )
        R = round(statistics["Recall @ Cis-BP Threshold"] * 100, 2)
        Jglobals.write(
            None, "\t*** Recall @ Cis-BP Threshold: %s" % R
        )

    print(thresholds)
    print(precision)
    print(recall)
    exit(0)

    return(statistics)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()