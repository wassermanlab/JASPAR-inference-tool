# this file will turn the pfam_DBDAlignment.json file to a readable number (Identity)
import pandas as pd
import numpy as np

import json
import subprocess
import string
from array import *
import os
from Bio.SubsMat.MatrixInfo import blosum62

file_header = os.path.dirname(os.path.abspath(__file__))
# global var

#pfamAhmm = "/sp4work/nicole/Pfam-A.hmm" #change if this become part of the tool

hmmf = "tmp.hmm"

# Already got consensus sequence from getConsensus.py script
# def getPfamSeq(pfamid): #get consensus sequence (for each pfam hmm)
#     global hmmf
#     global pfamAhmm
#     # grep the real name (with versions)
#     hmm_name = subprocess.check_output(['grep',pfamid,pfamAhmm],universal_newlines=True).split()[1]
    
#     # run hmmfetch and save as tmp.hmm
#     hmm = subprocess.run(['hmmfetch','-o',hmmf,pfamAhmm,hmm_name])

#     # get and parse output from hmmemit
#     seq = subprocess.check_output(['hmmemit','-c',hmmf])


def getDiffChar(string1,string2): # return the positions where the string differs
    diff = [0]*len(string1)
    # string1 should have same length as string 2
    if (len(string1) == len(string2)):
        # proceed to for loop
        for i in range(len(string1)):
            diff[i] = BLOSUMScoring(string1[i],string2[i])
    else:
        # there is something wrong
        print(string1,string2," is not equal length")
        exit(0)
    print(string1)
    print(string2)
    print("difference is:", diff)
    return(diff)

def BLOSUMScoring(aa1,aa2):
    if aa1 == '-' and aa2 =='-':
        return(1)
    elif aa1 == '-' or aa2 =='-':
        return(-4)
    else:
        try:
            return(blosum62[(aa1,aa2)])
        except:
            return(blosum62[(aa2,aa1)])

def scoringID(testseq1,testseq2): # called for each comparison to compare the string
    # compare and score
    positionScore = getDiffChar(testseq1,testseq2)
    return(positionScore) # returns an list of 1 or 0 (on the length)

def recursionSeq(startseq,listSeq): # the recursion that calls itself..

    #base case
    if (len(listSeq) == 1):
        return([scoringID(startseq,listSeq[0])])

    curScore = scoringID(startseq,listSeq[0])
    fulScore = [curScore]+recursionSeq(startseq,listSeq[1:])

    return(fulScore)

def permutationSequences(listOfseq): # 
    permutationScore = list()
    for i in range(len(listOfseq)-1):
        if (i == 0):
            print(listOfseq[i],"* and *",listOfseq[i+1:])
            permutationScore = recursionSeq(listOfseq[i],listOfseq[i+1:])
            #exit(0)
        else:
            permutationScore = permutationScore + recursionSeq(listOfseq[i],listOfseq[i+1:])

    return(permutationScore)


# Pairwise comparison so don't need the majority consens sequences
# # read in consensus json
# with open(file_header + "consensusDBDseq.json",'r') as consensusF:
#     consensusD = json.loads(consensusF.read())

# read in alignment json

# Change fileName comment when testing
fileName = file_header + "/pfam_DBDAlignment.json"
#fileName = file_header + "/test.json"

with open(fileName,'r') as alignF:
    alignD = json.loads(alignF.read())

# dictionary for the matrices
matrixD = {k:[] for k in alignD}

# process the file :
# 1. remove the lowercase
# 2. convert them to a list of 1 and 0
# 3. store it as a list of json

for key,item in alignD.items(): # key does not have version .xx number

    # get the real pfam sequence
    #realseq = consensusD[key]
    
    #remove lower case letters in the list
    alignments = [x.translate(str.maketrans('','',string.ascii_lowercase)) for x in item]
    print(alignments)
    permScore = permutationSequences(alignments)
    matrixD[key] = permScore
    #print(permScore)
    #print(permScore)

    #exit(0)
# storing the matrix...
with open(file_header+"/BLOSUMScoreMatrix.json",'w+') as wf:
    json.dump(matrixD, wf , sort_keys = True , indent = 4, separators = (",",": "))
