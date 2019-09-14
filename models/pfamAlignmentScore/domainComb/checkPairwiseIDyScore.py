# Since everything is true. This py will pull from
# domainHMMalign_withMAID.json and all_against_all.json 
# to figure out if all is as 1

import json

# y hits
with open("all_against_all.json",'r') as yF:
    yD = json.loads(yF.read())

# domain hits
with open("domainHMMalign_withMAID.json",'r') as domainF:
    domainD = json.loads(domainF.read())

def eachMAcomp(maID1,maID2):
    global yD
    print(maID1,", ",maID2)
    if maID2 in yD[maID1] and maID1 in yD[maID2]: # both in other's cluster
        return(True)
    else:
        return(False)

def fetchY(maIDlist1,maIDlist2):
    # global maidDNE
    # get the list of maID from 1 and 2 to see if any of them is correlated
    for maID in maIDlist1:
        for maID2 in maIDlist2:
            #maID = maID.split('.')[0]
            #maID2 = maID.split('.')[0]
            #if maID in maidDNE: # if it is in the maidDNE (maid that does not exist)
            #    return(False)
            if eachMAcomp(maID,maID2):
                return(True)
    return(False)

newYDic = dict()
newVal = list()
for key,val in yD.items():
    newKey = key.split(".")[0]
    for associate in val:
        newA = associate.split(".")[0]
        newVal.append(newA)
        print(newA)
    if newKey in newYDic:
        # error: duplicative key
        print(newKey," in dict")
        exit(0)

    newYDic[newKey] = newVal
    newVal= []

print(newYDic)


# check for every 
for key,value in domainD.items():
    # domainD organization
    # key = pfam ID sequences ex. GATA,GATA
    # value[x] = [[MA],[sequence]]
    # value[x][0] = [list of MA]
    # value[x][1] = [list of sequence]

    # aimed output:
    # key = pfam ID sequence
    # value = [[matrix score],[y score]] # the numbers should correspond

    if len(value) < 3:
        continue # skip things that are too small

    # loop for comparing elements
    for i in range(len(value)):
        # nested while to keep comparing sequences 
        for j in range(i+1,len(value)):

            # inner most loop for examining EACH different components
            print(value[i])
            for k in range(len(value[i][1])):
                # for each sequence in the y...

                # print(value[i][1][k],' and ',value[j][1][k])
                # scores = scoringID(removeLower(value[i][1][k]),removeLower(value[j][1][k]))
                # print(scores)

                print(value[i][0],' and ',value[j][0])
                yscore = fetchY(value[i][0],value[j][0])
                print(yscore)

                # print(yscoring)

    # dataD[key] = [matrixScoring,yscoring]
    # matrixScoring = list()
    yscoring = list()
