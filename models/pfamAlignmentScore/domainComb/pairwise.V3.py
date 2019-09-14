# this is a new script that does pairwise alignments of the sequences ()
# the empty string set is removed ... no pfam ID
# this pairwise comparison is different... - It is by tomtom output

###### V2 update #######
# This version includes Rost curves as well (read in as a dict) and 
# incorporates the rating into a three component matrix
########################


###### V3 #############
# This will concatenate the lists (of multiple domain) of ID to one


import json
import string # needed for string remove lower case
import os # needed for string remove lower case
from Bio.SubsMat.MatrixInfo import blosum62

file_header = os.path.dirname(os.path.abspath(__file__))

with open('all_against_all.json','r') as yOut: # reciprocal hits result
    yDic = json.loads(yOut.read())

newYDic = dict()
values = []

##maidDNE = ["MA0549","MA0506","MA0322"] # matrix ID that does not exist

# need to process y dic...
for key,val in yDic.items():
    newKey = key.split(".")[0]
    for associate in val:
        newA = associate.split(".")[0]
        values.append(newA)
        print(newA)
    if newKey in newYDic:
        # error: duplicative key
        print(newKey," in dict")
        exit(0)

    newYDic[newKey] = values
    values = []

print(newYDic)
#exit(0)

with open('domainHMMalign_withMAID.V2.json','r') as domainF:
    domainD = json.loads(domainF.read())




####################### V2 block #########################
# Read in rost curve dict
with open('homologs.json','r') as rostF:
    rostD = json.loads(rostF.read())

def crossRost(unip1,unip2):
    global rostD
    if unip1 in rostD[unip2] and unip2 in rostD[unip1]:
        return(True)
    else:
        return(False)

##########################################################


# helper scoring functions
def getDiffChar(string1,string2): # return the positions where the string differs
    diff = [0]*len(string1)
    # string1 should have same length as string 2
    if (len(string1) == len(string2)):
        # proceed to for loop
        for i in range(len(string1)):
            #diff[i] = BLOSUMScoring(string1[i],string2[i]) # if want to BLOSUM score
            diff[i] = IDscoring(string1[i],string2[i]) # if want ID score
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

def IDscoring(aa1,aa2):
    if aa1==aa2 and aa1!='-':
        return(1)
    else:
        return(0)

def scoringID(testseq1,testseq2): # called for each comparison to compare the string
    # compare and score
    positionScore = getDiffChar(testseq1,testseq2)
    return(positionScore) # returns an list of 1 or 0 (on the length)

# the order of calling the functions
# 1. scoring ID (give sequences) -> 2. getDiffChar -> 3. BLOSUM -> -> 
# return list of score (each position is a score)


def eachMAcomp(maID1,maID2):
    global newYDic
    if maID2 in newYDic[maID1] and maID1 in newYDic[maID2]: # both in other's cluster
        return(True)
    else:
        return(False)

def fetchY(maIDlist1,maIDlist2):
    # global maidDNE
    # get the list of maID from 1 and 2 to see if any of them is correlated
    for maID in maIDlist1:
        maIDn = maID.split('.')[0]

        for maID2 in maIDlist2:
            
            print("MAID ",maIDn)
            maID2n = maID2.split('.')[0]
            print("MAID ",maID2n)
            #if maID in maidDNE: # if it is in the maidDNE (maid that does not exist)
            #    return(False)
            if eachMAcomp(maIDn,maID2n):
                return(True)
    return(False)

def removeLower(stringsss):
    # basically remove lowercase...
    return(stringsss.translate(str.maketrans('','',string.ascii_lowercase)))

dataD = dict() # output dictionary
matrixScoring = list() # will be lists of lists
yscoring = list()
rostScoring = list()
# added unipID
unipIDs = list()

# V3 variables ##################################
multReg = list() # V3L for linking region directly
ysingleScore = list() # record all singles (k segments)
rostSingleScore = list() 
stringent = True # if True, all must be True then True; This is on the TOMTOM score
# the Rost result would not be affected by the stringent parameter


for key,value in domainD.items():
    # domainD organization
    # key = pfam ID sequences ex. GATA,GATA
    # value[x] = [[MA],[sequence],unip]
    # value[x][0] = [list of MA]
    # value[x][1] = [list of sequence]
    # V2: value[x][2] = unip

    # aimed output:
    # key = pfam ID sequence
    # value = [[matrix score],[y score]] # the numbers should correspond

    if len(value) < 3:
        continue # skip things that are too small

    # loop for comparing elements
    for i in range(len(value)):
        # nested while to keep comparing sequences 
        for j in range(i+1,len(value)):

            # clean lists
            multReg = list()
            rostPass = False
            yPass = False
            ysingleScore = list()
            rostSingleScore = list()

            # inner most loop for examining EACH different components
            print(value[i])
            for k in range(len(value[i][1])):
                # for each sequence in the y...

                #print(value[i][1][k],' and ',value[j][1][k])
                scores = scoringID(removeLower(value[i][1][k]),removeLower(value[j][1][k]))
                #print(scores)

                print(value[i][0],' and ',value[j][0])
                yscore = fetchY(value[i][0],value[j][0])
                print(yscore)

                ########### V2 ############
                print("unip ",value[i][2],' and ',value[j][2])
                unipComp = value[i][2]+"*"+value[j][2]
                rostScore = crossRost(value[i][2],value[j][2])
                print(rostScore)

                ############# V3 ##############
                ysingleScore.append(yscore)
                rostSingleScore.append(rostScore)
                multReg.extend(scores) 
                #matrixScoring.append(scores)
            
            if stringent:
                yPass = all(x==True for x in ysingleScore)
                rostPass = all(x==True for x in rostSingleScore)
            else:
                # not stringent
                yPass = any(x==True for x in ysingleScore)
                rostPass = any(x==True for x in rostSingleScore)

            # append must be made at a higher level
            yscoring.append(yPass)
            rostScoring.append(rostPass)
            unipIDs.append(unipComp) # don't need to be changed bc would be the same
            matrixScoring.append(multReg)

            #print("length of yscore, rost score, matrix Score",yscoring,' ',rostScoring,' ',matrixScoring)

            

    dataD[key] = [matrixScoring,yscoring,rostScoring,unipIDs]
    matrixScoring = list()
    yscoring = list()
    rostScoring = list()
    unipIDs = list()

wantedID = "Homeodomain"
for stuff in dataD[wantedID][0]:
    if len(stuff)!=len(dataD[wantedID][0][0]):
        print(len(stuff))

with open(file_header+'/pairwise_output.IDscore.V3.json','w+') as wf:
    json.dump(dataD,wf,sort_keys = True , indent = 4, separators = (",",": "))

# key error MA0549