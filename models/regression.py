#!/usr/bin/env python

# Reference:
# https://towardsdatascience.com/
# building-a-logistic-regression-in-python-step-by-step-becd4d56c9c8

import argparse
import json
import numpy as np
import os
import pickle
import re
from sklearn.linear_model import LogisticRegressionCV, RidgeCV
from sklearn.metrics import average_precision_score
import sys

# Defaults
out_dir = os.path.dirname(os.path.realpath(__file__))
root_dir = os.path.join(out_dir, os.pardir)

# # Defaults for plotting
# import matplotlib.pyplot as plt 
# import pandas as pd
# import seaborn as sns
# plt.rc("font", size=14)
# sns.set(style="white")
# sns.set(style="whitegrid",color_codes=True)

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

def train_models(pairwise_file, out_dir=out_dir, verbose=False):

    # Skip if JSON/pickle files already exists
    json_file = os.path.join(out_dir, "results.json")
    pickle_file = os.path.join(out_dir, "models.pickle")
    if not os.path.exists(json_file) or not os.path.exists(pickle_file):

        # Initialize
        models = {}
        results = {}

        # Load JSON file
        with open(pairwise_file) as f:
            pairwise = json.load(f)

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
                Jglobals.write(None, "\t*** Ys: %s" % Ys.shape)
                Jglobals.write(None, "\t*** TF pairs: %s" % len(tfPairs))

            # For each regression approach...
            for regression in ["linear", "logistic"]:

                # If linear regression...
                if regression == "linear":
                    alphas = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0] # i.e. reg. coeff. from Andrew Ng
                    regModel = RidgeCV(alphas=alphas, cv=myCViterator)

                # ... Else...
                else:
                    regModel = LogisticRegressionCV(Cs=10, cv=myCViterator, max_iter=1000)

                # For each sequence similarity representation...
                for similarity in ["identity", "blosum62"]:          

                    # Fit model...
                    fitRegModel = regModel.fit(Xs[similarity], Ys_int)

                    # If linear regression...
                    if regression == "linear":
                        predictions = fitRegModel.predict(Xs[similarity])

                    # ... Else...
                    else:
                        predictions = fitRegModel.predict_proba(Xs[similarity])[:,1]

                    # i.e. area under the precision-recall curve
                    avg_precision = average_precision_score(Ys_int, predictions)
                    
                    # Verbose mode
                    if verbose:
                        Jglobals.write(None, "\t*** Avg. Precision (%s + %s): %s" % (regression, similarity, avg_precision))

                    # Add fitRegModel
                    models[domains].setdefault((regression, similarity), fitRegModel)
                    results[domains].setdefault((regression, similarity), (avg_precision, fitRegModel.coef_.tolist()[0]))

            exit(0)

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