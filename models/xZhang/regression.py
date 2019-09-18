#!/usr/bin/env python

# Reference:
# https://towardsdatascience.com/
# building-a-logistic-regression-in-python-step-by-step-becd4d56c9c8

import argparse
import json
import matplotlib.pyplot as plt 
import numpy as np
import os
import pandas as pd
import pickle
import re
import seaborn as sns
import sklearn
from sklearn.linear_model import LogisticRegressionCV, RidgeCV
from sklearn.metrics import confusion_matrix
import sys

# Defaults
out_dir = os.path.dirname(os.path.realpath(__file__))
root_dir = os.path.join(out_dir, os.pardir, os.pardir)

# Defaults for plotting
plt.rc("font", size=14)
sns.set(style="white")
sns.set(style="whitegrid",color_codes=True)

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

    parser.add_argument("-j", help="json file (from pairwise.py)", metavar="JSON")
    parser.add_argument("-o", default=out_dir, help="output directory (default = ./)", metavar="DIR")
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose mode (default = False)")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Make Pfam files
    train_models(os.path.abspath(args.j), os.path.abspath(args.o), args.verbose)

def train_models(json_file, out_dir=out_dir, verbose=False):

    # Initialize
    models = {"Regressions": ["Logistic", "Ridge"]}
    results = {"Regressions": ["Logistic", "Ridge"]}

    # Get cluster and regression
    m = re.search("pairwise.(rsat|tomtom)\+(id|sim).json", json_file)
    cluster = m.group(1)
    regression = m.group(2)

    # Skip if JSON/pickle files already exists
    results_json_file = os.path.join(out_dir, "results.%s+%s.json" % (cluster, regression))
    models_pickle_file = os.path.join(out_dir, "models.%s+%s.pickle" % (cluster, regression))
    if not os.path.exists(results_json_file) or not os.path.exists(models_pickle_file):

        # Load JSON file
        with open(json_file) as f:
            pairwise = json.load(f)

        # For each DBD composition...
        for domains, values in pairwise.items():

            # Verbose mode
            if verbose:
                Jglobals.write(None, "\nRegressing %s..." % domains)

            # Initialize
            models.setdefault(domains, [])
            results.setdefault(domains, [])
            Xs = np.asarray(values[0])
            Ys = np.array(values[1])
            tfPairs = values[2]

            # Verbose mode
            if verbose:
                Jglobals.write(None, "\t*** Xs: %s / %s" % Xs.shape)
                Jglobals.write(None, "\t*** Ys: %s" % Ys.shape)
                Jglobals.write(None, "\t*** TF pairs: %s" % len(tfPairs))

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

            # 1) Logistic regression
            logReg = LogisticRegressionCV(Cs=10, cv=myCViterator, max_iter=1000)
            logRegModel = logReg.fit(Xs, Ys)
            predictions = logRegModel.predict(Xs)
            ##
            ## From sklearn.metrics.confusion_matrix:
            ##
            ## The count of true negatives is [0][0], false negatives is [1][0],
            ## true positives is [1][1] and false positives is [0][1].
            ##
            matrix = confusion_matrix(Ys, predictions)
            tn = float(matrix[0][0])
            fn = float(matrix[1][0])
            tp = float(matrix[1][1])
            fp = float(matrix[0][1])
            accuracy = 0
            if (tp + tn + fp + fn) > 0:
                # Fix zeroDivisionError
                accuracy =  (tp + tn) / (tp + tn + fp + fn)
            precission = 0
            if (tp + fp) > 0:
                # Fix zeroDivisionError
                precission = tp / (tp + fp)
            recall = 0
            if (tp + fn) > 0:
                # Fix zeroDivisionError
                recall = tp / (tp + fn)
            results[domains].append(([accuracy, precission, recall], logRegModel.coef_.tolist()[0]))
            models[domains].append(logRegModel)

            # Verbose mode
            if verbose:
                Jglobals.write(None, "\t*** Accuracy (LogisticRegressionCV): %s" % accuracy)
                Jglobals.write(None, "\t*** Precission (LogisticRegressionCV): %s " % precission)
                Jglobals.write(None, "\t*** Recall (LogisticRegressionCV): %s" % recall)

            # 2) Ridge regression (i.e. linear)
            # alphas = [0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
            alphas = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0] # i.e. reg. coeff. from Andrew Ng
            linReg = RidgeCV(alphas=alphas, cv=myCViterator)
            linRegModel = linReg.fit(Xs, Ys)
            predictions = linRegModel.predict(Xs)
            matrix = confusion_matrix(Ys, predictions >= 0.5)
            tn = float(matrix[0][0])
            fn = float(matrix[1][0])
            tp = float(matrix[1][1])
            fp = float(matrix[0][1])
            accuracy = 0
            if (tp + tn + fp + fn) > 0:
                # Fix zeroDivisionError
                accuracy =  (tp + tn) / (tp + tn + fp + fn)
            precission = 0
            if (tp + fp) > 0:
                # Fix zeroDivisionError
                precission = tp / (tp + fp)
            recall = 0
            if (tp + fn) > 0:
                # Fix zeroDivisionError
                recall = tp / (tp + fn)
            results[domains].append(([accuracy, precission, recall], logRegModel.coef_.tolist()[0]))
            models[domains].append(logRegModel)

            # Verbose mode
            if verbose:
                Jglobals.write(None, "\t*** Accuracy (RidgeCV): %s" % accuracy)
                Jglobals.write(None, "\t*** Precission (RidgeCV): %s " % precission)
                Jglobals.write(None, "\t*** Recall (RidgeCV): %s" % recall)

        # Write JSON
        Jglobals.write(
            results_json_file,
            json.dumps(results, sort_keys=True, indent=4, separators=(",", ": "))
        )

        # Write pickle file
        with open(models_pickle_file, "wb") as f:
            pickle.dump(models, f)

def _leaveOneTFOut(tfIdxs, l): # dict

    myCViterator = []

    # For each TF...
    for tf, idxs in tfIdxs.items():

        s = [i for i in range(l)]
        testIdx = idxs
        trainIdx = list(set(s) - set(testIdx))

        myCViterator.append((trainIdx, testIdx))

    return(myCViterator)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()