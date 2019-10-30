#!/usr/bin/env python

# Reference:
# https://towardsdatascience.com/
# building-a-logistic-regression-in-python-step-by-step-becd4d56c9c8

import argparse
from collections import Counter
from glmnet import ElasticNet, LogitNet
import json
import numpy as np
from numpy import log10 as log
from operator import itemgetter 
import os
import pickle
from sklearn.metrics import mean_squared_error, precision_recall_curve
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

    parser.add_argument("--domain", metavar="STR",
        help="regress models for given domain (default = all)")
    parser.add_argument("-o", default=out_dir, metavar="DIR",
        help="output directory (default = ./)")
    parser.add_argument("-p", metavar="JSON",
        help="compressed JSON file from pairwise.py")
    parser.add_argument("--threads", type=int, default=1, metavar="INT",
        help="threads to use (default = 1)")
    parser.add_argument("--thresh-neg", type=float, default=1., metavar="FLT",
        help="threads to use (default = 1.0)")
    parser.add_argument("--thresh-pos", type=float, default=6., metavar="FLT",
        help="threads to use (default = 6.0)")
    parser.add_argument("-v", "--verbose", action="store_true",
        help="verbose mode (default = False)")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Make Pfam files
    train_models(os.path.abspath(args.p), args.domain, os.path.abspath(args.o),
        args.threads, args.thresh_neg, args.thresh_pos, args.verbose)

def train_models(pairwise_file, dbd=None, out_dir=out_dir, threads=1, thresh_neg=1., 
    thresh_pos=6., verbose=False):

    # Initialize
    global lambdas
    lambdas = _get_lambda_path()

    # Skip if pickle file already exists
    models_file = os.path.join(out_dir, "models.pickle")
    results_file = os.path.join(out_dir, "results.json")
    if not os.path.exists(models_file) or not os.path.exists(results_file):

        # Initialize
        models = {
            "Keys": "DBD composition",
            "Values": {
                (
                    "regression approach",
                    "similarity representation"
                ) : "model"
            }
        }
        results = {
            "Keys": "DBD composition",
            "Values": {
                (
                    "regression approach",
                    "similarity representation"
                ) : ("Ys", "predictions")
            }
        }

        # Load JSON file
        global pairwise
        handle = Jglobals._get_file_handle(pairwise_file)
        pairwise = json.load(handle)
        handle.close()

        # E-value thresholds from cisbp.ipynb:
        # For (+) 6.702705320605285 and 6.34620565834759
        global threshPos
        threshPos = thresh_pos
        # For (-) 0.9352972238139634 and -0.6065759798508961
        global threshNeg
        threshNeg = thresh_neg

        # For each DBD composition...
        for domain, values in pairwise.items():
            # Zn_clus
            # PAX
            # Forkhead
            # HMG_box
            # zf-C2H2+zf-C2H2+zf-C2H2
            # zf-C4

            if dbd is not None:
                if domain != dbd:
                    continue

            # Verbose mode
            if verbose:
                Jglobals.write(None, "\nRegressing %s..." % domain)

            # Train models
            _train_LinReg_models(values, threads, verbose)
            _train_LogReg_models(values, threads, verbose)
            #_train_BLAST_models(values, threads, verbose)

    # # Write pickle file
    # with open(models_file, "wb") as f:
    #     pickle.dump(models, f)

def _train_LinReg_models(values, threads=1, verbose=False):

    # Initialize
    global TFpairs
    Xs = values[0]
    Ys = np.array(values[2])
    TFpairs = values[3]

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

        # Initialize
        Xss = np.asarray(Xs[similarity])
        m = ElasticNet(alpha=0, n_splits=len(TFs), lambda_path=lambdas,
            lower_limits=np.zeros(len(Xss[0])))

        # Set custom cross-validation
        m.CV = LeaveOneTfOut

        # Fit
        mFit = m.fit(Xss, Ys, sample_weight=weights)

        # Get best lambda (i.e. max)
        lambdabest = mFit.lambda_max_
        print("\t\t- lambdabest: %s" % lambdabest)

        # Predict
        p = mFit.predict(Xss, lamb=lambdabest)
        prec, rec, thresh = precision_recall_curve(Ys >= threshPos, p)
        for x in range(len(thresh)):
            if prec[x] >= 0.75:
                print(prec[x], rec[x], thresh[x])
                break

def _train_LogReg_models(values, threads=1, verbose=False):

    # Initialize
    global TFpairs
    Xs = values[0]
    Ys = np.array(values[2])
    TFpairs = values[3]

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
            Jglobals.write(None, "\t*** LogitNet: %s" % similarity)

        # Initialize
        Xss = np.asarray(Xs[similarity])
        m = LogitNet(alpha=0, n_splits=len(TFs), lambda_path=lambdas,
            lower_limits=np.zeros(len(Xss[0])))

        # Set custom cross-validation
        m.CV = LeaveOneTfOut

        # Fit
        mFit = m.fit(Xss, Ys >= threshPos, sample_weight=weights)

        # Get best lambda (i.e. max)
        lambdabest = mFit.lambda_max_
        print("\t\t- lambdabest: %s" % lambdabest)

        # Predict
        p = mFit.predict_proba(Xss, lamb=lambdabest)[:,1]
        prec, rec, thresh = precision_recall_curve(Ys >= threshPos, p)
        for x in range(len(thresh)):
            if prec[x] >= 0.75:
                print(prec[x], rec[x], thresh[x])
                break

# def _transform_Ys(Ys):

#     # Initialize
#     Ys_transformed = []

#     for y in Ys:
#         if y < 0:
#             y = 0
#         elif y >= 10:
#             y = 10
#         Ys_transformed.append(y)

#     return(Ys_transformed)

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
    Weight positive samples 1/freq (or 10x whichever is higher) higher 
    than negatives because they are usually at a much lower frequency.
    """

    labels = Ys >= threshPos
    c = Counter(labels)
    w = {}
    for l in c:
        w.setdefault(l, 1/(float(c[l])/len(labels)))
    # if w[True] < 10:
    #     w[True] = 10.
    # w[False] = 1.

    return[w[l] for l in labels]

def _get_lambda_path(min_lambda=1e-3, max_lambda=1e+3, reg_step = 0.01):

    lambdas = np.arange(log(min_lambda), log(max_lambda), reg_step)
    lambdas = list(np.power(10, lambdas))
    lambdas.append(max_lambda)
    lambdas.sort(reverse=True)

    return(lambdas)

# def _get_best_lambda(CVs, logistic=False):

#     # Initialize
#     MSEs = []
#     lambdabest = max(lambdas)
#     MSEbest = None

#     # For each cross-validation...
#     for cv in range(len(CVs)):

#         MSEs.append([])

#         # 0: xTrain
#         # 1: yTrain
#         # 2: wTrain
#         # 3: xTest
#         # 4: yTest
#         # 5: wTest
#         # 6: mFit
#         xTest = CVs[cv][3]
#         yTrue = list(CVs[cv][4])
#         mFit = CVs[cv][6]
#         if logistic:
#             yPred = list(map(list, zip(*(yPred * 1))))
#             print(yPred)
#         else:
#             yPred = mFit.predict(xTest, lamb=mFit.lambda_path_)
#             yPred = list(map(list, zip(*yPred)))


#         for l in range(len(lambdas)):

#             if logistic:
#                 MSEs[-1].append(mean_squared_error(yTrue, yPred[l]))
#             else:
#                 MSEs[-1].append(mean_squared_error(yTrue, yPred[l]))

    # Transpose
    MSEs = list(map(list, zip(*MSEs)))

    for l in range(len(lambdas)):

        MSEavg = np.mean(MSEs[l])

        if MSEbest is None or MSEavg < MSEbest:
            lambdabest = lambdas[l]
            MSEbest = np.mean(MSEs[l])

    return(lambdabest, MSEbest)

def _get_prc(CVs, lambdabest, logistic=False, npv=False):
    """
    From sklearn.metrics.precision_recall_curve

    The last precision and recall values are "1." and "0." respectively and do
    not have a corresponding threshold. This ensures that the graph starts on
    the y axis.
    """

    # Initialize
    yTrue = []
    yPred = []

    # For each cross-validation...
    for cv in range(len(CVs)):

        # 0: xTrain
        # 1: yTrain
        # 2: wTrain
        # 3: xTest
        # 4: yTest
        # 5: wTest
        # 6: mFit
        xTest = CVs[cv][3]
        yTrue += list(CVs[cv][4])
        mFit = CVs[cv][6]
        if logistic:
            yPred += list(mFit.predict_proba(xTest, lamb=lambdabest))
        else:
            yPred += list(mFit.predict(xTest, lamb=lambdabest))

    if logistic:
        if npv:
            yTrue = np.array(yTrue) * -1
        else:
            yTrue = np.array(yTrue) * 1
    if not logistic:
        if npv:
            yTrue = np.array(yTrue) < threshNeg
        else:
            yTrue = np.array(yTrue) >= threshPos

    return(precision_recall_curve(yTrue, np.array(yPred)))

def _get_tf_recall_curve(tfPairs, labels, predictions, Ys):
    """
    From sklearn.metrics.precision_recall_curve

    The last precision and recall values are "1." and "0." respectively and do
    not have a corresponding threshold. This ensures that the graph starts on
    the y axis.
    """

    # Initialize
    tf_recall = []
    tfs = set(tf for i in range(len(tfPairs)) if labels[i] == 1 for j in tfPair)

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