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

    # Skip if pickle file already exists
    models_file = os.path.join(out_dir, "models.pickle")
    # results_file = os.path.join(out_dir, "results.pickle")
    # if not os.path.exists(models_file) or not os.path.exists(results_file):
    if not os.path.exists(models_file):

        # Initialize
        models = {
            "Keys": "DBD composition",
            "Values": {
                (
                    "regression approach",
                    "similarity representation"
                ) : (["precisions"], ["recalls"], ["Ys"], ["TF recalls"], "recall at 75% precision", "Y at 75% precision", "TF recall at 75% precision", "model")
            }
        }
        # results = {
        #     "Keys": "DBD composition",
        #     "Values": {
        #         (
        #             "regression approach",
        #             "similarity representation"
        #         ) : (["precisions"], ["recalls"], ["Ys"], ["weights"])
        #     }
        # }

        # Load pickle file
        with open(pairwise_file, "rb") as f:
            pairwise = pickle.load(f)

        # For each DBD composition...
        for domains, values in pairwise.items():

            # Initialize
            Xs = {}
            # BLASTXs = {}
            models.setdefault(domains, {})
            # results.setdefault(domains, {})
            # Ys = np.array(values[2])
            Ys = np.array(values[1])
            Ys_int = Ys * 1
            # tfPairs = values[3]
            tfPairs = values[2]

            # Verbose mode
            if verbose:
                Jglobals.write(None, "\nRegressing %s..." % domains)

            # Leave one TF out:
            # 1. Get the TF names into a dict
            # 2. For each TF, figure out which index to include
            tfIdxs = {}
            for tfPair in tfPairs:
                for tf in tfPair:
                    tfIdxs.setdefault(tf, [])
            for tf, idxs in tfIdxs.items():
                for i in range(len(tfPairs)):
                    if tf in tfPairs[i]:
                        idxs.append(i)

            # Get iterator for cross validation
            myCViterator = _leaveOneTFOut(tfIdxs, len(tfPairs))

            # For each sequence similarity representation...
            for similarity in ["identity", "blosum62"]:

                # Add Xs
                Xs.setdefault(similarity, np.asarray(values[0][similarity]))
                # BLASTXs.setdefault(similarity, np.array(values[1][similarity]))

                # Verbose mode
                if verbose:
                    a, b = Xs[similarity].shape
                    Jglobals.write(None, "\t*** Xs (%s): %s / %s" % (similarity, a, b))
                    # a, b = BLASTXs[similarity].shape
                    # Jglobals.write(None, "\t*** Xs (BLAST+): %s / %s" % (a, b))

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
                # for similarity in ["identity", "blosum62", "%ID"]:
                for similarity in ["identity", "blosum62"]:


                    # # For use BLAST+ Xs...
                    # for use_blast_Xs in [False, True]:

                    # Initialize
                    # if similarity == "%ID":
                    #     if regression == "logistic":
                    #         continue
                    #     # if use_blast_Xs:
                    #     #     continue
                    #     myXs = []
                    #     for pairwise in Xs["identity"]:
                    #         myXs.append([float(sum(pairwise)) / len(pairwise)])
                    #     myXs = np.array(myXs)
                    # # elif use_blast_Xs:
                    # #     myXs = np.concatenate((Xs[similarity], BLASTXs[similarity]), axis=1)
                    # else:
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
                    tfRec = _get_tf_recall_curve(tfPairs, predictions, Ys)
                    recall = _get_value_at_precision_threshold(Prec, Rec, threshold=0.75)
                    y = _get_value_at_precision_threshold(Prec, Ys, threshold=0.75)
                    tf_recall = _get_value_at_precision_threshold(Prec, tfRec, threshold=0.75)

                    # Verbose mode
                    if verbose:
                        # Jglobals.write(None, "\t*** Recall at 75% Precision threshold ({} + {} + BLAST+ = {}): {}".format(regression, similarity, use_blast_Xs, recall))
                        Jglobals.write(None, "\t*** Recall at 75% Precision threshold ({} + {}): {}".format(regression, similarity, recall))
                        Jglobals.write(None, "\t*** Recalled TFs at 75% Precision threshold ({} + {}): {}".format(regression, similarity, tf_recall))

                    # Add fitRegModel
                    # models[domains].setdefault((regression, similarity), (recall, tf_recall, y, fitRegModel))
                    # results[domains].setdefault((regression, similarity), (Prec, Rec, Ys, tfRec, fitRegModel.coef_.tolist()[0]))
                    models[domains].setdefault((regression, similarity), (Prec, Rec, Ys, tfRec, recall, y, tf_recall, fitRegModel))

        # Write pickle file
        with open(models_file, "wb") as f:
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

def _get_tf_recall_curve(tfPairs, predictions, Ys):
    """
    From sklearn.metrics.precision_recall_curve

    The last precision and recall values are "1." and "0." respectively and do
    not have a corresponding threshold. This ensures that the graph starts on
    the y axis.
    """

    # Initialize
    tf_recall = []
    tfs = set(tf for tfPair in tfPairs for tf in tfPair)

    # For each y...
    for y in Ys:

        # Initialize
        tfs_recalled = set()

        # For each predictions...
        for i in range(len(predictions)):
            if predictions[i] >= y:
                tfs_recalled.add(tfPairs[i][0])
                tfs_recalled.add(tfPairs[i][1])

        tf_recall.append(float(len(tfs_recalled)) / len(tfs))

    # Append a last recall value of 0.0
    tf_recall.append(0.0)

    return(np.array(tf_recall))

def _get_value_at_precision_threshold(Prec, values, threshold=0.75):

    # For each y...
    for i in range(len(Prec)):

        if Prec[i] >= threshold or i +1 == len(values):
            break

    return(values[i])

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()