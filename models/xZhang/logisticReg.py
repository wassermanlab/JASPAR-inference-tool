#!/usr/bin/env python

# this python script will import the data from the 
# json (pairwise_output.json) and apply logistic regression 
# on it

# Reference:
# https://towardsdatascience.com/building-a-logistic-regression-in-python-step-by-step-becd4d56c9c8

import argparse
import json
import matplotlib.pyplot as plt 
import numpy as np
import os
import pandas as pd
import pickle
import seaborn as sns
import sklearn
from sklearn import linear_model
from sklearn.linear_model import LogisticRegressionCV
from sklearn.linear_model import RidgeCV
import sys

# Defaults
out_dir = os.path.dirname(os.path.realpath(__file__))
root_dir = os.path.join(out_dir, os.pardir, os.pardir)
files_dir = os.path.join(root_dir, "files")

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

    parser.add_argument("-c", choices=["rsat", "tomtom"], default="tomtom", help="cluster profiles using \"rsat\" matrix-clustering or \"tomtom\" (i.e. default)")
    parser.add_argument("-f", default=files_dir, help="files directory (from get_files.py; default=../files/)", metavar="DIR")
    parser.add_argument("-o", default=out_dir, help="output directory (default = ./)", metavar="DIR")

    return(parser.parse_args())


def main():

    # Parse arguments
    args = parse_args()

    # Make Pfam files
    train_models(args.c, os.path.abspath(args.f), os.path.abspath(args.o))

def train_models(cluster="tomtom", files_dir=files_dir, out_dir=out_dir):
    """
    Model information:
    # logistic and linear regression for each regression
    # Only TOMTOM result vs TOMTOM result & Homology as Y values
    # L2 regularization strengthen w cross validation: Leave one out
    """

# read in files
with open("pairwise_output.IDscore.V3.json",'r') as pairF:
    pairD = json.loads(pairF.read())

# generator def
def cv_leaveOneTFOut(tfNames,totalLen): # dict
    # generator yields TFs
    myCViterator = []
    for name,value in tfNames.items():
        s = [i for i in range(totalLen)]
        testindex = value
        trainindex = list(set(s)-set(testindex))
        myCViterator.append((trainindex,testindex))
    return(myCViterator)


# Def for ridge
# def ridgeRegression(x,y,alpha):
#     ridgereg = Ridge(alpha = alpha)
#     # did not normalize

#     ridgereg.fit(x,y)
#     # TBD 

# the score (accuracy) is stored in the dictionary as a list
regScore = dict()
regScore["Header_columns"] = ['logistic with homolog as x','Ridge with homolog as x','logistic with merged homolog and tomtom','Ridge with merged homolog and tomtom']

# model
model = dict()

# output score and the weight matrix in the form of tuples

# for each family
for key, value in pairD.items():
    print("*******Family ",key," is regressed*********")
    regScore[key] = list()
    model[key] = list()

    # read in x y
    # key = protein DBD name sequence
    positionMat = np.asarray(value[0]) # list of list
    print("position matrix dim: ",positionMat.shape)
    tomtom = np.array(value[1])
    print("tomtom dim: ",tomtom.shape)
    rost = np.array(value[2])
    print("Rost homolo dim: ",rost.shape)
    tfPairs = value[3] # just as a list
    print("tfPairs dim: ",len(tfPairs))

    # z which is combination of tomtom and rost
    z = list()
    for i in range(len(value[1])):
        if value[1][i] and value[2][i]:
            z.append(True)
        else:
            z.append(False)

    zArr = np.array(z)

    # leave one out as in leave one TF out!!!
    # 1. get the names of the TF into a dict
    # 2. (during iteration) if contained that then continue to next loop
    tfNames = dict()
    for pair in tfPairs:
        names = pair.split('*')
        if not names[0] in tfNames:
            tfNames[names[0]]=list()
        if not names[1] in tfNames:
            tfNames[names[1]]=list() 

    # Figure out for each entry of the dict, which column to include and which to not
    indexPresent = list() # index to keep track which ones does not contain TF
    for name,lst in tfNames.items():
        # go through list of tfPairs
        for i in range(len(tfPairs)):
            if name in tfPairs[i]:
                lst.append(i) # lst consists of test set


    # part 1: treat rost as an x value

    print("Logistic regression .///")
    
    #iterator for cross validation
    myCViterator = cv_leaveOneTFOut(tfNames,len(tfPairs))

    # logistic regression
    logreg = LogisticRegressionCV(Cs = 10,cv = myCViterator, scoring = "accuracy", max_iter = 10000)

    # add the rost score as an additional column
    #X_data = positionMat
    #X_data[:,:-1] = rost.T
    #X_data = np.hstack((positionMat,rost.T))
    X_data = np.concatenate((positionMat,rost[:,None]),axis=1) # Value error: Transpose the adding one
    #X_data = np.hstack()
    #print("X data dimension: ",X_data.shape)

    # fit
    try:
        resultclf = logreg.fit(X_data,tomtom)
    except:
        print("Not enough class: ",X_data.shape)
        continue
    resultclf = logreg.fit(X_data,tomtom)
    score = resultclf.score(X_data,tomtom)
    print("Logistic score: ",score)
    regScore[key].append((score,resultclf.coef_.tolist()[0]))
    print("Weights of things: ",resultclf.coef_.tolist()[0])
    
    model[key].append(resultclf)
    #exit(0)

    ################ Commented out: no cross valid on alpha ##########################
    # # train and test split data
    # positionMat_train, positionMat_test = positionMat.take(train_index, axis=0), positionMat.take(test_index, axis=0)
    # rost_train, rost_test = rost[train_index], rost[test_index]
    # tomtom_train, tomtom_test = tomtom[train_index], tomtom[test_index]
    
    # # treat rost just as another column...

    # X_train = np.append(positionMat_train,rost_train, axis=1)
    # print("Dimention of training: ",X_train.shape)
    # logreg.fit(X_train, tomtom_train)

    # # Assess accuracy
    # X_test = np.append(positionMat_test,rost_test, axis=1)
    # print("Dimension of testing: ", X_test.shape)
    # Y_pred = logreg.predict(X_test)

    # scoring = logreg.score(X_test,tomtom_test)
    # print("Accuracy score: ",scoring)
    # # original paper needs at least 75% precision threshold

    # # confusion matrix
    # confusion_matrix = confusion_matrix(tomtom_test,Y_pred)

    # # Matthews correlation coefficient
    # mattcoef = sklearn.metrics.matthews_corrcoef(tomtom_test,Y_pred)
    # print("Matthew's Correlation coefficient: ",mattcoef)

        # select the best scoring one
    

    # part 1 b: Ridge regularization
    
    # CV splits are done above at myCViterator
    ridge = RidgeCV(alphas = [0.01,0.05,0.1,0.3,0.5,0.7,0.9,1.0],cv = myCViterator)
    ridgefit = ridge.fit(X_data,tomtom)

    rgscore = ridgefit.score(X_data,tomtom)
    print("Ridge Score: ",rgscore)
    regScore[key].append((rgscore,ridgefit.coef_.tolist()[0]))
    model[key].append(ridgefit)
    # # train and test split data
    # positionMat_train, positionMat_test = positionMat.take(train_index, axis=0), positionMat.take(test_index, axis=0)
    # rost_train, rost_test = rost[train_index], rost[test_index]
    # tomtom_train, tomtom_test = tomtom[train_index], tomtom[test_index]
    
    # # treat rost just as another column...

    # X_train = np.append(positionMat_train,rost_train, axis=1)
    # print("Dimention of training: ",X_train.shape)
    # X_test = np.append(positionMat_test,rost_test, axis=1)
    # print("Dimension of testing: ", X_test.shape)

    # ridgeRegression(X_train, tomtom_train,alpha)

    # # Assess accuracy
    
    # Y_pred = logreg.predict(X_test)

    # scoring = logreg.score(X_test,tomtom_test)
    # print("Accuracy score: ",scoring)
    # # original paper needs at least 75% precision threshold
    # # define higher and lower threshold from each model
    # # from the model:
    # # >95% => highly similar
    # # 25% ID dissimilar 

    # # confusion matrix
    # confusion_matrix = confusion_matrix(tomtom_test,Y_pred)

    # # Matthews correlation coefficient
    # mattcoef = sklearn.metrics.matthews_corrcoef(tomtom_test,Y_pred)
    # print("Matthew's Correlation coefficient: ",mattcoef)

    # part 2: treat zArr as the y value, positionMat as the x value

    # logistic
    logreg = LogisticRegressionCV(Cs = 10,cv = myCViterator, scoring = "accuracy", max_iter = 10000)

    # add the rost score as an additional column
    #X_data = positionMat
    #X_data[:,:-1] = rost.T
    #X_data = np.hstack((positionMat,rost.T))
    X_data = positionMat
    print("X data dimension: ",X_data.shape)

    # fit
    resultclf = logreg.fit(X_data,zArr)
    score = resultclf.score(X_data,zArr)
    print("zArr as y value , logistic score: ",score)
    regScore[key].append((score,resultclf.coef_.tolist()[0]))
    print("Weights of things: ",resultclf.coef_.tolist()[0])
    model[key].append(resultclf)

    # ridge z array
    ridge = RidgeCV(alphas = [0.01,0.05,0.1,0.3,0.5,0.7,0.9,1.0],cv = myCViterator)
    ridgefit = ridge.fit(X_data,zArr)

    rgscore = ridgefit.score(X_data,zArr)
    print("ZArr as value, Ridge Score: ",rgscore)
    regScore[key].append((score,ridgefit.coef_.tolist()[0]))
    model[key].append(ridgefit)
    #exit(0)
    # visualize the data (as in )

# use other score than accuracy?? since linear would be >0.5 to be 1 else .


# output the regression/logistic score (accuracy)
with open("regression_result.json",'w+') as f_ou:
    json.dump(regScore,f_ou, sort_keys = True , indent = 4, separators = (",",": "))
    
# store models
with open("trained_models.pickle",'wb') as wf:
    pickle._dump(model,wf)
# to load = just use pickle.load(f,..)
# each family has four models

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()