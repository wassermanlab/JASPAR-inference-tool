#!/usr/bin/env python

# Reference:
# https://towardsdatascience.com/
# building-a-logistic-regression-in-python-step-by-step-becd4d56c9c8

import argparse
from collections import Counter
from functools import partial
from glmnet import ElasticNet
from multiprocessing import Pool
import numpy as np
from operator import itemgetter 
import os
import pickle
from sklearn.metrics import mean_squared_error, precision_recall_curve
import sys
from tqdm import tqdm

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
    This function parses arguments provided via the command line and returns
    an {argparse} object.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-o", default=out_dir, metavar="DIR",
        help="output directory (default = ./)")
    parser.add_argument("-p", metavar="PICKLE",
        help="pickle file from pairwise.py")
    parser.add_argument("--threads", type=int, default=1, metavar="INT",
        help="threads to use (default = 1)")
    parser.add_argument("-v", "--verbose", action="store_true",
        help="verbose mode (default = False)")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Make Pfam files
    train_models(os.path.abspath(args.p), os.path.abspath(args.o), args.threads, args.verbose)

def train_models(pairwise_file, out_dir=out_dir, threads=1, verbose=False):

    # Skip if pickle file already exists
    models_file = os.path.join(out_dir, "models.pickle")
    results_file = os.path.join(out_dir, "results.pickle")
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

        # Load pickle file
        global pairwise
        with open(pairwise_file, "rb") as f:
            pairwise = pickle.load(f)

        # For each DBD composition...
        for domains, values in pairwise.items():

            if domains != "GATA":
                continue

            # Train models
            _train_SR_models(domains, values, threads, verbose)
            _train_BLAST_models(domains, values, threads, verbose)

    # # Write pickle file
    # with open(models_file, "wb") as f:
    #     pickle.dump(models, f)

def _train_SR_models(domains, values, threads=1, verbose=False):

    # Initialize
    Xs = values[0]
    Ys = np.array(values[2])
    # Ys = np.array(_transform_Ys(Ys))
    TFpairs = values[3]

    # E-value thresholds from cisbp.ipynb:
    # For (+) 6.702705320605285 and 6.34620565834759
    # For (-) 0.9352972238139634 and -0.6065759798508961
    threshPos = 6.0
    threshNeg = 1.0

    # Verbose mode
    if verbose:
        Jglobals.write(None, "\nRegressing %s..." % domains)

    # Get weights
    weights = _get_weights(Ys, threshPos)

    # # For each sequence similarity representation...
    # for similarity in ["identity", "blosum62"]:

    #     # Add Xs
    #     Xs.setdefault(similarity, np.asarray(values[0][similarity]))
    #     BLASTXs.setdefault(similarity, np.array(values[1][similarity]))

    #     # Verbose mode
    #     if verbose:
    #         a, b = Xs[similarity].shape
    #         Jglobals.write(None, "\t*** Xs (%s): %s / %s" % (similarity, a, b))
    #         a, b = BLASTXs[similarity].shape
    #         Jglobals.write(None, "\t*** Xs (BLAST+): %s / %s" % (a, b))

    # # Verbose mode
    # if verbose:
    #     Jglobals.write(None, "\t*** Ys: %s" % Ys.shape)
    #     Jglobals.write(None, "\t*** TF pairs: %s" % len(tfPairs))

    # For each sequence similarity representation...
    for similarity in ["identity", "blosum62"]:

        # Initialize
        Xss = np.asarray(Xs[similarity])
        limits = np.zeros(len(Xss[0]))

        # Get lambdas for cross-validation
        model = ElasticNet(alpha=0, lower_limits=limits, standardize=False)
        modelFit = model.fit(Xss, Ys, sample_weight=weights)
        lambdas = modelFit.lambda_path_

        # Initialize cross-validation model 
        modelCV = ElasticNet(alpha=0, lambda_path=lambdas, lower_limits=limits,
            standardize=False)

        # Get cross-validations
        CVs = _get_cross_validations(Xss, Ys, weights, modelCV, TFpairs)

        # Get best lambda
        idx, lambdabest, mse = _get_best_lambda(CVs, lambdas)

        # Get Precision, Recall, thresholds
        prec, rec, threshs = _get_prc(CVs, idx, threshPos, npv=False)
        # nprec, nrec, nthreshs = _get_prc(CVs, idx, threshNeg, npv=True)

        # For each profile...
        for i in range(len(prec)):
            print(prec[i], rec[i], threshs[i])
        continue


        # # Verbose mode
        # if verbose:
        #     # Jglobals.write(None, "\t*** Recall at 75% Precision threshold ({} + {} + BLAST+ = {}): {}".format(regression, similarity, use_blast_Xs, recall))
        #     Jglobals.write(None, "\t*** Recall at 75% Precision threshold ({} + {}): {}".format(regression, similarity, recall))
        #     Jglobals.write(None, "\t*** Recalled TFs at 75% Precision threshold ({} + {}): {}".format(regression, similarity, tf_recall))

        # # Add fitRegModel
        # # models[domains].setdefault((regression, similarity), (recall, tf_recall, y, fitRegModel))
        # # results[domains].setdefault((regression, similarity), (Prec, Rec, Ys, tfRec, fitRegModel.coef_.tolist()[0]))
        # models[domains].setdefault((regression, similarity), (Prec, Rec, Ys, tfRec, recall, y, tf_recall, fitRegModel))

def _train_BLAST_models(domains, values, threads=1, verbose=False):

    # Initialize
    Xs = values[0]
    Ys = np.array(values[2])
    # Ys = np.array(_transform_Ys(Ys))
    TFpairs = values[3]

    # E-value thresholds from cisbp.ipynb:
    # For (+) 6.702705320605285 and 6.34620565834759
    # For (-) 0.9352972238139634 and -0.6065759798508961
    threshPos = 6.0
    threshNeg = 1.0

    # Verbose mode
    if verbose:
        Jglobals.write(None, "\nRegressing %s..." % domains)

    # Get weights
    weights = _get_weights(Ys, threshPos)

    # # For each sequence similarity representation...
    # for similarity in ["identity", "blosum62"]:

    #     # Add Xs
    #     Xs.setdefault(similarity, np.asarray(values[0][similarity]))
    #     BLASTXs.setdefault(similarity, np.array(values[1][similarity]))

    #     # Verbose mode
    #     if verbose:
    #         a, b = Xs[similarity].shape
    #         Jglobals.write(None, "\t*** Xs (%s): %s / %s" % (similarity, a, b))
    #         a, b = BLASTXs[similarity].shape
    #         Jglobals.write(None, "\t*** Xs (BLAST+): %s / %s" % (a, b))

    # # Verbose mode
    # if verbose:
    #     Jglobals.write(None, "\t*** Ys: %s" % Ys.shape)
    #     Jglobals.write(None, "\t*** TF pairs: %s" % len(tfPairs))

    # For each sequence similarity representation...
    for similarity in ["identity", "blosum62"]:

        if similarity == "blosum62":
            continue

        # Initialize
        Xss = []
        limits = np.zeros(1)

        # For each X...
        for X in Xs[similarity]:
            Xss.append(sum(X)/float(len(X)))
        Xss = np.asarray(Xss)
        Xss = Xss.reshape(-1, 1)

        # Get lambdas for cross-validation
        model = ElasticNet(alpha=0, lower_limits=limits, standardize=False)
        modelFit = model.fit(Xss, Ys, sample_weight=weights)
        lambdas = modelFit.lambda_path_

        # Initialize cross-validation model 
        modelCV = ElasticNet(alpha=0, lambda_path=lambdas, lower_limits=limits,
            standardize=False)

        # Get cross-validations
        CVs = _get_cross_validations(Xss, Ys, weights, modelCV, TFpairs)

        # Get best lambda
        idx, lambdabest, mse = _get_best_lambda(CVs, lambdas)

        # Get Precision, Recall, thresholds
        prec, rec, threshs = _get_prc(CVs, idx, threshPos, npv=False)
        # nprec, nrec, nthreshs = _get_prc(CVs, idx, threshNeg, npv=True)

        # For each profile...
        for i in range(len(prec)):
            print(prec[i], rec[i], threshs[i])
        exit(0)


        # # Verbose mode
        # if verbose:
        #     # Jglobals.write(None, "\t*** Recall at 75% Precision threshold ({} + {} + BLAST+ = {}): {}".format(regression, similarity, use_blast_Xs, recall))
        #     Jglobals.write(None, "\t*** Recall at 75% Precision threshold ({} + {}): {}".format(regression, similarity, recall))
        #     Jglobals.write(None, "\t*** Recalled TFs at 75% Precision threshold ({} + {}): {}".format(regression, similarity, tf_recall))

        # # Add fitRegModel
        # # models[domains].setdefault((regression, similarity), (recall, tf_recall, y, fitRegModel))
        # # results[domains].setdefault((regression, similarity), (Prec, Rec, Ys, tfRec, fitRegModel.coef_.tolist()[0]))
        # models[domains].setdefault((regression, similarity), (Prec, Rec, Ys, tfRec, recall, y, tf_recall, fitRegModel))

def _transform_Ys(Ys):

    # Initialize
    Ys_transformed = []

    for y in Ys:
        if y < 0:
            y = 0
        elif y >= 10:
            y = 10
        Ys_transformed.append(y)

    return(Ys_transformed)

def _get_weights(Ys, threshold=6.0):
    """
    Weight positive samples 1/freq (or 10x whichever is higher) higher 
    than negatives because they are usually at a much lower frequency.
    """

    labels = Ys >= threshold
    c = Counter(labels)
    w = {}
    for l in c:
        w.setdefault(l, 1/(float(c[l])/len(labels)))
    if w[True] < 10:
        w[True] = 10.
    w[False] = 1.

    return[w[l] for l in labels]

def _get_cross_validations(x, y, w, m, pairs):

    # Initialize
    CVs = []

    # Get CV iterator
    CViterator = _get_CV_iterator(pairs)

    # For each cross-validation...
    for cv in range(len(CViterator)):

        # Get train/test Xs/Ys
        xTrain = np.asarray(itemgetter(*CViterator[cv][0])(x))
        yTrain = np.asarray(itemgetter(*CViterator[cv][0])(y))
        wTrain = np.asarray(itemgetter(*CViterator[cv][0])(w))
        xTest = np.asarray(itemgetter(*CViterator[cv][1])(x))
        yTest = np.asarray(itemgetter(*CViterator[cv][1])(y))
        wTest = np.asarray(itemgetter(*CViterator[cv][1])(w))

        # Fit cross-validation model
        mFit = m.fit(xTrain, yTrain, sample_weight=wTrain)

        # Predict
        CVs.append((xTrain, yTrain, wTrain, xTest, yTest, wTest, mFit, mFit.predict(xTest, lamb=mFit.lambda_path_)))

    return(CVs)

def _get_CV_iterator(TFpairs):
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
        CViterator.append((idxTrain, idxTest))

    return(CViterator)

def _get_best_lambda(CVs, lambdas):

    # Initialize
    regStrength = []

    # Get Ys
    yTrue, yPred = _get_Ys_CV(CVs)

    # For each lambda...
    for l in range(len(lambdas)):
        regStrength.append((l, lambdas[l], mean_squared_error(yTrue, yPred[l])))

    # Sort by mean squared error
    regStrength.sort(key=lambda x: x[-1])

    return(regStrength[0])

def _get_prc(CVs, idx, threshold, npv=False):

    # Get Ys
    yTrue, yPred = _get_Ys_CV(CVs)

    if not npv:
        labels = (np.array(yTrue) >= threshold) * 1
    else:
        labels = (np.array(yTrue) < threshold) * 1

    prec, rec, threshs = precision_recall_curve(labels, yPred[idx])
    threshs = np.append(np.array([min(yPred[idx])]), threshs)

    return(prec, rec, threshs)

def _get_Ys_CV(CVs):

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
        # 7: mFit.predict(xTest, lamb=mFit.lambda_path_)))
        yTrue += list(CVs[cv][4])
        yPred += list(CVs[cv][7])

    # Transpose
    yPred = list(map(list, zip(*yPred)))

    return(yTrue, yPred)

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