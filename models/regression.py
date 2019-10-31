#!/usr/bin/env python

# Reference:
# https://towardsdatascience.com/
# building-a-logistic-regression-in-python-step-by-step-becd4d56c9c8

import argparse
from collections import Counter
from glmnet import ElasticNet
import json
import numpy as np
from numpy import log10 as log
from operator import itemgetter 
import os
import pickle
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
# Classes     #
#-------------#

class LeaveOneTfOut:

    def __init__(self, n_splits=3, shuffle=False, random_state=None):
        self.n_splits = n_splits

    def split(self, X=None, y=None, groups=None):

        for train, test in CViterator:
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

    parser.add_argument("-j", metavar="JSON",
        help="compressed JSON file from pairwise.py")
    parser.add_argument("-o", default=out_dir, metavar="DIR",
        help="output directory (default = ./)")
    parser.add_argument("--threads", type=int, default=1, metavar="INT",
        help="threads to use (default = 1)")
    parser.add_argument("-v", "--verbose", action="store_true",
        help="verbose mode (default = False)")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Make Pfam files
    train_models(os.path.abspath(args.j), os.path.abspath(args.o),
        args.threads, args.verbose)

def train_models(pairwise_file, out_dir=out_dir, threads=1, verbose=False):

    # Initialize
    global lambdas
    lambdas = _get_lambda_path()
    global threshPos
    threshPos = 8.

    # Skip if pickle file already exists
    models_file = os.path.join(out_dir, "models.pickle")
    if not os.path.exists(models_file):

        # Initialize
        models = {
            "Keys": "DBD composition",
            "Values": ["similarity", "model", "lambdabest", "y cut-off"]
        }

        # Load JSON file
        global pairwise
        handle = Jglobals._get_file_handle(pairwise_file)
        pairwise = json.load(handle)
        handle.close()

        # For each DBD composition...
        for domain, values in pairwise.items():

            # Verbose mode
            if verbose:
                Jglobals.write(None, "\nRegressing %s..." % domain)

            # Train models
            similarity, model, lambdabest, y = _train_LinReg_models(values, threads, verbose)

            # Add model
            models.setdefault(domain, [similarity, model, lambdabest, y])

        # Write pickle file
        with open(models_file, "wb") as f:
            pickle.dump(models, f)

        # Verbose mode
        if verbose:
            Jglobals.write(None, "\n")

def _get_lambda_path(min_lambda=1e-3, max_lambda=1e+3, reg_step = 0.01):

    lambdas = np.arange(log(min_lambda), log(max_lambda), reg_step)
    lambdas = list(np.power(10, lambdas))
    lambdas.append(max_lambda)
    lambdas.sort(reverse=True)

    return(lambdas)

def _train_LinReg_models(values, threads=1, verbose=False):

    # Initialize
    global TFpairs
    results = []
    Xs = values[0]
    Ys = np.array(values[1])
    TFpairs = values[2]

    # Get weights
    weights = _get_weights(Ys)

    # Get unique TFs
    TFs = set()
    for TFpair in TFpairs:
        TFs.add(TFpair[0])
        TFs.add(TFpair[1])

    # Get CV iterator
    global CViterator
    CViterator = _leaveOneTfOut(TFpairs)

    # For each sequence similarity representation...
    for similarity in ["identity", "blosum62"]:

        # Verbose mode
        if verbose:
            Jglobals.write(None, "\t*** ElasticNet: %s" % similarity)

        try:

            # Initialize
            Xss = np.asarray(Xs[similarity])
            m = ElasticNet(alpha=0, n_splits=len(TFs), lambda_path=lambdas,
                lower_limits=np.zeros(len(Xss[0])), n_jobs=threads)

            # Set custom cross-validation
            m.CV = LeaveOneTfOut

            # Fit
            mFit = m.fit(Xss, Ys, sample_weight=weights)

            # Get best lambda (i.e. max)
            lambdabest = mFit.lambda_max_
            if verbose:
                Jglobals.write(None, "\t\t- lambdabest: %s" % lambdabest)

            # Predict
            p = mFit.predict(Xss, lamb=lambdabest)

            # PPV
            prec, rec, thresh = precision_recall_curve(Ys >= threshPos, p)
            for x in range(len(thresh)):
                if prec[x] >= 0.75:
                    results.append([similarity, mFit, lambdabest, thresh, rec[x]])
                    break

            # If results...
            if len(results) > 0:
                if verbose:
                    Jglobals.write(None, "\t\t- recall @ 75% precision: {0:.2f}%".\
                        format(round(results[-1][-1] * 100, 2)))
            # Else...
            else:
                if verbose:
                    Jglobals.write(None, "\t\t- could not reach 75% precision!")
        except:
            pass

    # No results...
    if len(results) == 0:
        if verbose:
            Jglobals.write(None, "\t*** similarity regression models could not be trained!")
        return(None, None, None, None)

    # Sort by recall
    results.sort(key=lambda x: x[-1], reverse=True)

    return(results[0][:4])

def _leaveOneTfOut(TFpairs):
    """
    Leave one TF out.
    """

    # Initialize
    CViterator = []

    # 1. Get the TF names into a dict
    TFidxs = {}
    for tfPair in TFpairs:
        for tf in tfPair:
            TFidxs.setdefault(tf, [])
    for tf, idxs in TFidxs.items():
        for i in range(len(TFpairs)):
            if tf in TFpairs[i]:
                idxs.append(i)

    # 2. For each TF, figure out which index to include
    for tf, idxs in TFidxs.items():
        s = [i for i in range(len(TFpairs))]
        idxTest = idxs
        idxTrain = list(set(s) - set(idxTest))
        CViterator.append((np.array(idxTrain), np.array(idxTest)))

    return(CViterator)

def _get_weights(Ys):
    """
    Weight samples 1/freq.
    """

    labels = Ys >= threshPos
    c = Counter(labels)
    w = {}
    for l in c:
        w.setdefault(l, 1/(float(c[l])/len(labels)))

    return[w[l] for l in labels]

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()