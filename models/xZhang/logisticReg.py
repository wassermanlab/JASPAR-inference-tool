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
from sklearn import linear_model
from sklearn.linear_model import LogisticRegressionCV as logistic
from sklearn.linear_model import RidgeCV as linear
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
    """
    Model information:
    # logistic and linear regression for each regression
    # Only TOMTOM result vs TOMTOM result & Homology as Y values
    # L2 regularization strengthen w cross validation: Leave one out
    """

    # Initialize
    models = {"Regressions": ["Logistic", "Ridge"]}
    regScores = {"Regressions": ["Logistic", "Ridge"]}

    # Get cluster and regression
    m = re.search("pairwise.(rsat|tomtom)\+(id|sim).json", json_file)
    cluster = m.group(1)
    regression = m.group(2)

    # Load JSON file
    with open(json_file) as f:
        pairwise = json.load(f)

    # For each DBD composition...
    for domains, values in pairwise.items():

        if domains != "Homeodomain":
            continue

        # Verbose mode
        if verbose:
            Jglobals.write(None, "\nRegressing %s..." % domains)

        # Initialize
        models.setdefault(domains, [])
        regScores.setdefault(domains, [])
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
        logreg = logistic(Cs = 10, cv = myCViterator, scoring = "accuracy", max_iter = 10000)
        result = logreg.fit(Xs, Ys)
        score = result.score(Xs, Ys)
        regScores[domains].append((score, result.coef_.tolist()[0]))
        models[domains].append(result)

        # Verbose mode
        if verbose:
            Jglobals.write(None, "\t*** Logistic score: %s" % regScores[domains][0][0])
            Jglobals.write(None, "\t*** Weights:")
            Jglobals.write(None, ", ".join(map(str, regScores[domains][0][1])))

        # 2) Linear regression (i.e. Ridge)
        alphas = [0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
        linreg = linear(alphas=alphas, cv=myCViterator)
        result = linreg.fit(Xs, Ys)
        score = result.score(Xs, Ys)
        regScores[domains].append((score, result.coef_.tolist()[0]))
        models[domains].append(result)

        # Verbose mode
        if verbose:
            Jglobals.write(None, "\t*** Ridge (linear) score: %s" % regScores[domains][1][0])

        exit(0)

    # # use other score than accuracy?? since linear would be >0.5 to be 1 else .


    # # output the regression/logistic score (accuracy)
    # with open("regression_result.json",'w+') as f_ou:
    #     json.dump(regScore,f_ou, sort_keys = True , indent = 4, separators = (",",": "))
        
    # # store models
    # with open("trained_models.pickle",'wb') as wf:
    #     pickle._dump(model,wf)
    # # to load = just use pickle.load(f,..)
    # # each family has four models

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