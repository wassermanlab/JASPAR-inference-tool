#!/usr/bin/env python

# Reference:
# https://towardsdatascience.com/
# building-a-logistic-regression-in-python-step-by-step-becd4d56c9c8
#
# Notes:
# For similarity representation "%ID", we don't really apply a regression.
# Instead, we set a threshold on the Y at a precision >= 75%.

import argparse
import json
import numpy as np
import os
import pickle
import re
from sklearn.linear_model import LogisticRegressionCV, RidgeCV
from sklearn.metrics import precision_recall_curve
import sys

# Defaults
out_dir = os.path.dirname(os.path.realpath(__file__))
root_dir = os.path.join(out_dir, os.pardir)

# Append JASPAR-profile-inference to path
sys.path.append(root_dir)

# Import globals
from __init__ import Jglobals

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns an {argparse} object.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-o", default=out_dir, help="output directory (default = ./)", metavar="DIR")
    parser.add_argument("-p", help="pickle file from pairwise.py", metavar="PICKLE")
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose mode (default = False)")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Make Pfam files
    train_models(os.path.abspath(args.p), os.path.abspath(args.o), args.verbose)

def train_models(pairwise_file, out_dir=out_dir, verbose=False):

    # Skip if JSON/pickle files already exists
    json_file = os.path.join(out_dir, "results.json")
    pickle_file = os.path.join(out_dir, "models.pickle")
    if not os.path.exists(json_file) or not os.path.exists(pickle_file):

        # Initialize
        models = {
            "Keys": "DBD composition",
            "Values": {
                (
                    "regression approach",
                    "similarity representation"
                ) : ("Y at 75% precision", "model")
            }
        }
        results = {
            "Keys": "DBD composition",
            "Values": {
                (
                    "regression approach",
                    "similarity representation"
                ) : (["precisions"], ["recalls"], ["Ys"], ["weights"])
            }
        }

        # Load pickle file
        with open(pairwise_file, "rb") as f:
            pairwise = pickle.load(f)

        # For each DBD composition...
        for domains, values in pairwise.items():

            # Initialize
            Xs = {}
            models.setdefault(domains, {})
            results.setdefault(domains, {})
            Ys = np.array(values[1])
            Ys_int = Ys * 1
            tfPairs = values[2]

            # Verbose mode
            if verbose:
                Jglobals.write(None, "\nRegressing %s..." % domains)

            # Leave one TF out:
            # 1. Get the TF names into a dict
            # 2. For each TF, figure out which index to include
            tfIdxs = {}
            for tfPair in tfPairs:
                for tf in tfPair.split("*"):
                    tfIdxs.setdefault(tf, [])
            for tf, idxs in tfIdxs.items():
                for i in range(len(tfPairs)):
                    tfs = set(tfPairs[i].split("*"))
                    if tf in tfs:
                        idxs.append(i)

            # Get iterator for cross validation
            myCViterator = _leaveOneTFOut(tfIdxs, len(tfPairs))

            # For each sequence similarity representation...
            for similarity in ["identity", "blosum62"]:

                # Add Xs
                Xs.setdefault(similarity, np.asarray(values[0][similarity]))

                # Verbose mode
                if verbose:
                    a, b = Xs[similarity].shape
                    Jglobals.write(None, "\t*** Xs (%s): %s / %s" % (similarity, a, b))

            # Verbose mode
            if verbose:
                Jglobals.write(None, "\t*** Ys: %s" % Ys_int.shape)
                Jglobals.write(None, "\t*** TF pairs: %s" % len(tfPairs))

            # For each regression approach...
            for regression in ["linear", "logistic"]:

                # If linear regression...
                if regression == "linear":
                    alphas = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0] # i.e. reg. coeff. from Andrew Ng
                    regModel = RidgeCV(alphas=alphas, cv=myCViterator)

                # ... Else...
                else:
                    regModel = LogisticRegressionCV(Cs=10, cv=myCViterator, max_iter=50000)

                # For each sequence similarity representation...
                for similarity in ["identity", "blosum62", "%ID"]:

                    # Initialize
                    if similarity == "%ID":
                        if regression == "logistic":
                            continue
                        myXs = []
                        for pairwise in Xs["identity"]:
                            myXs.append([float(sum(pairwise)) / len(pairwise)])
                        myXs = np.array(myXs)
                    else:
                        myXs = Xs[similarity]

                    # Fit model...
                    fitRegModel = regModel.fit(myXs, Ys_int)

                    # If linear regression...
                    if regression == "linear":
                        predictions = fitRegModel.predict(myXs)

                    # ... Else...
                    else:
                        predictions = fitRegModel.predict_proba(myXs)[:,1]

                    # Get precision-recall curve
                    Prec, Rec, Ys = precision_recall_curve(Ys_int, predictions)
                    recall, y = _get_recall_and_y_at_precision_threshold(Prec, Rec, Ys, threshold=0.75)

                    # Verbose mode
                    if verbose:
                        Jglobals.write(None, "\t*** Recall at 75% Precision threshold ({} + {}): {}".format(regression, similarity, recall))

                    # Add fitRegModel
                    models[domains].setdefault((regression, similarity), (y, fitRegModel))
                    results[domains].setdefault((regression, similarity), (Prec, Rec, Ys, fitRegModel.coef_.tolist()[0]))

        # Write JSON
        Jglobals.write(
            json_file,
            json.dumps(results, sort_keys=True, indent=4, separators=(",", ": "))
        )

        # Write pickle file
        with open(pickle_file, "wb") as f:
            pickle.dump(models, f)

def _leaveOneTFOut(tfIdxs, l):

    myCViterator = []

    # For each TF...
    for tf, idxs in tfIdxs.items():

        s = [i for i in range(l)]
        testIdx = idxs
        trainIdx = list(set(s) - set(testIdx))

        myCViterator.append((trainIdx, testIdx))

    return(myCViterator)

def _get_recall_and_y_at_precision_threshold(Prec, Rec, Ys, threshold=0.75):

    # For each y...
    for i in range(len(Ys)):

        if Prec[i] >= threshold:

            return(Rec[i], Ys[i])

    # Return worst case scenario:
    # i.e. recall = 1 and y = 1
    return(0.0, 1.0)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()