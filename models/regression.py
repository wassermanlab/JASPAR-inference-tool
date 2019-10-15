#!/usr/bin/env python

# Reference:
# https://towardsdatascience.com/
# building-a-logistic-regression-in-python-step-by-step-becd4d56c9c8

import argparse
from functools import partial
from glmnet import ElasticNet
from multiprocessing import Pool
import numpy as np
from operator import itemgetter 
import os
import pickle
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
        evalue_threshold = 7.0 # from cisbp.ipynb: 6.7027 and 6.3462
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
        with open(pairwise_file, "rb") as f:
            pairwise = pickle.load(f)

        # For each DBD composition...
        for domains, values in pairwise.items():

            if domains != "WRKY":
                continue

            # Initialize
            Xs = {}
            BLASTXs = {}
            models.setdefault(domains, {})
            results.setdefault(domains, {})
            Ys = np.array(values[2])
            tfPairs = values[3]

            # Verbose mode
            if verbose:
                Jglobals.write(None, "\nRegressing %s..." % domains)

            # Get CV iterator
            myCViterator = _get_CV_iterator(tfPairs)

            # For each sequence similarity representation...
            for similarity in ["identity", "blosum62"]:

                # Add Xs
                Xs.setdefault(similarity, np.asarray(values[0][similarity]))
                BLASTXs.setdefault(similarity, np.array(values[1][similarity]))

                # Verbose mode
                if verbose:
                    a, b = Xs[similarity].shape
                    Jglobals.write(None, "\t*** Xs (%s): %s / %s" % (similarity, a, b))
                    a, b = BLASTXs[similarity].shape
                    Jglobals.write(None, "\t*** Xs (BLAST+): %s / %s" % (a, b))

            # Verbose mode
            if verbose:
                Jglobals.write(None, "\t*** Ys: %s" % Ys.shape)
                Jglobals.write(None, "\t*** TF pairs: %s" % len(tfPairs))

            # For each sequence similarity representation...
            for similarity in ["identity", "blosum62"]:

                # Initialize
                myXs = Xs[similarity]
                myBLASTXs = BLASTXs[similarity]
                lower_limits = np.zeros(len(myXs[0]))

                # For each regression approach...
                for regression in ["linear", "logistic"]:

                    if regression == "linear":
                        m = ElasticNet(alpha=0, lower_limits=lower_limits, standardize=False)

                    elif regression == "logistic":
                        continue

                    # For each cross-validation...
                    for i in range(len(myCViterator)):

                        # Initialize
                        x_train = np.asarray(itemgetter(*myCViterator[i][0])(myXs))
                        y_train = np.asarray(itemgetter(*myCViterator[i][0])(Ys))
                        x_test = np.asarray(itemgetter(*myCViterator[i][1])(myXs))
                        y_test = np.asarray(itemgetter(*myCViterator[i][1])(Ys))

                        m = m.fit(x_train, y_train)
                        p = m.predict(x_test)
                        print(p)

                        # if regression == "logistic":
                        #     y_train = y_train >= evalue_threshold
                        #     y_test = y_test >= evalue_threshold

                        # Fit model...
                        # fitRegModel = OneVsRestClassifier(regModel).fit(myXs, Ys_transform)
                        #m_cv = m.fit(x_train, y_train)

                        # Predict
                        # if regression == "linear":
                        # elif regression == "logistic":
                        #     p = m_cv.predict_proba(x_test)[:,1]

                        # print(x_test)
                        # print(y_test)
                        # print(p)
                        # for i in range(len(y_test)):
                        #     print(y_test[i], p[i])
                    #     print(m_cv.predict([[0 for i in range(57)]]))
                    #     print(m_cv.predict([[1 for i in range(57)]]))
                    #     exit(0)

                    #     # Add CV predictions
                    #     predictions.append((p, y_test))

                    # exit(0)

                    # # # Predict
                    # # if regression == "linear":
                    # #     predictions = fitRegModel.predict(myXs)

                    # # # ... Else...
                    # # else:
                    # #     predictions = fitRegModel.predict_proba(myXs)

                    # for i in range(100):
                    #     print(Ys[i], predictions[i], abs(Ys[i] - predictions[i]))
                    # exit(0)

                    # # Get precision-recall curve
                    # Prec, Rec, Ys = precision_recall_curve(Ys_int, predictions)
                    
                    # tfRec = _get_tf_recall_curve(tfPairs, Ys_int, predictions, Ys)
                    # recall = _get_value_at_precision_threshold(Prec, Rec, threshold=0.75)
                    # y = _get_value_at_precision_threshold(Prec, Ys, threshold=0.75)
                    # tf_recall = _get_value_at_precision_threshold(Prec, tfRec, threshold=0.75)

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

def _get_CV_iterator(tfPairs):
    """
    Leave one TF out.
    """

    # Initialize
    myCViterator = []

    # 1. Get the TF names into a dict
    tfIdxs = {}
    for tfPair in tfPairs:
        for tf in tfPair:
            tfIdxs.setdefault(tf, [])
    for tf, idxs in tfIdxs.items():
        for i in range(len(tfPairs)):
            if tf in tfPairs[i]:
                idxs.append(i)

    # 2. For each TF, figure out which index to include
    for tf, idxs in tfIdxs.items():

        s = [i for i in range(len(tfPairs))]
        testIdx = idxs
        trainIdx = list(set(s) - set(testIdx))

        myCViterator.append((trainIdx, testIdx))

    return(myCViterator)

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