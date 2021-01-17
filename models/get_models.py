#!/usr/bin/env python

import argparse
from collections import Counter
import copy
from glmnet import LogitNet
from itertools import chain
import json
import numpy as np
import os
import pandas as pd
import pickle
import shutil
from sklearn.metrics import precision_recall_curve
import sys

# Defaults
out_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)))
root_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir)
data_splits_dir = copy.copy(out_dir)

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

    parser.add_argument(
        "--data-splits-dir", default=data_splits_dir, metavar="DIR",
        help="data splits directory (\"-o\" from get_data_splits.py)")
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
    global data_splits
    json_file = os.path.join(args.data_splits_dir, ".data_splits.json.gz")
    handle = Jglobals._get_file_handle(json_file)
    data_splits = json.load(handle)
    handle.close()
    global lambdas
    lambdas = 10**np.linspace(6, -6, 100)
    global verbose
    verbose = args.verbose

    # Create output dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Get models
    get_models(os.path.abspath(args.o))

def get_models(out_dir=out_dir):

    # Skip if JSON file already exists
    json_file = os.path.join(out_dir, "models.json.gz")
    if not os.path.exists(json_file):

        # Initialize
        models = {
            "Keys": "DBD",
            "Values": [
                "similarity",
                "coefficients",
                "Y cut-off",
                "coverage @ 75% precision"
            ]
        }

        # Create groups dir
        if not os.path.isdir("models"):
            os.makedirs("models")

        # Move to groups directory
        os.chdir("models")

        # For each DBD composition...
        for DBD in data_splits:

            # Verbose mode
            if verbose:
                Jglobals.write(None, "\n%s..." % DBD)

            # Skip if pickle file already exists
            pickle_file = "%s.pickle" % DBD
            if not os.path.exists(pickle_file):

                # Initialize
                models = {}

                # For each sequence similarity representation...
                for similarity in ["identity", "blosum62"]:

                    # Train model
                    m = __train_model(DBD, similarity)
                    models.setdefault(similarity, m)

                    # Compute statistics
                    statistics = __compute_statistics(DBD, similarity, m)

                # Write
                handle = Jglobals._get_file_handle(pickle_file, "wb")
                pickle.dump(models, handle)
                handle.close()

            # # Skip if JSON file already exists
            # if DBD == "AP2":
            #     print(models["identity"].coef_)
            #     exit(0)

            # # Train models
            # # similarity, model, lambdabest, y = _train_LinReg_models(values)
            # similarity, coeffs, Y, coverage = _train_LinReg_models(values)

            # # Add model
            # # models.setdefault(domain, [similarity, model, lambdabest, y])
            # models.setdefault(domain, [similarity, coeffs, Y, coverage])

        # # Write pickle file
        # with open(gzip_file, "wb") as f:
        #     pickle.dump(models, f)
        # # Write
        # Jglobals.write(
        #     json_file[:-3],
        #     json.dumps(models, sort_keys=True, indent=4)
        # )
        # fi = Jglobals._get_file_handle(json_file[:-3], "rb")
        # fo = Jglobals._get_file_handle(json_file, "wb")
        # shutil.copyfileobj(fi, fo)
        # fi.close()
        # fo.close()
        # os.remove(json_file[:-3])

        # Verbose mode
        if verbose:
            Jglobals.write(None, "\n")

# def __get_lambda_path(min_lambda=1e-3, max_lambda=1e+3, reg_step = 0.01):

#     lambdas = np.arange(log(min_lambda), log(max_lambda), reg_step)
#     lambdas = list(np.power(10, lambdas))
#     lambdas.append(max_lambda)
#     lambdas.sort(reverse=True)

#     return(lambdas)

def __train_model(DBD, similarity):

    # Verbose mode
    if verbose:
        Jglobals.write(None, "\t*** LogitNet: %s" % similarity)

    # Get pairs, Xs, ys
    pairs, Xs, ys = __get_pairs_Xs_ys(DBD, similarity, "train")

    # Skip training
    if len(list(chain.from_iterable(ys["stringent"]))) < 2:

        # Verbose mode
        if verbose:
            Jglobals.write(None, "\t*** not enough data: skip training!")

        return(None)

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
            Jglobals.write(None, "\t*** model training: successful!")

        return(m_fit)

    except:

        # Verbose mode
        if verbose:
            Jglobals.write(None, "\t*** model training: failed!")

        return(None) 

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

def __compute_statistics(DBD, similarity, m):

    # Skip
    if m is None:

        # Verbose mode
        if verbose:
            Jglobals.write(None, "\t*** not enough data: skip statistics!")

        return(None, None, None, None)

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
    if verbose:
        success = False
        for idx in range(len(thresholds)):
            if precision[idx] < 0.75:
                continue
            r = round(recall[idx] * 100, 2)
            Jglobals.write(
                None, "\t*** model recall @ 75% precision: {0:.2f}%".format(r)
            )
            success = True
            break
        if not success:
            Jglobals.write(None, "\t*** model does not achieve 75% precision!")

    return(lambdabest, idx, precision, recall, thresholds)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()