
# this file will visualize the pfam alignments (spots of insertion or gap .)
import sys
import json

# store the aligned sequences into the list that is in a dictionary


with open("/homed/home/nzhang/JASPAR-tools/jaspartools/files/pfam_DBD_noVer.json",'r') as DBDs:
    DBDdic = json.loads(DBDs.read())

# new dic
alignedD = {k:[] for k in DBDdic}

# loop through the pfam hmm files
animals = ['fungi','nematodes','vertebrates','plants','insects']
for taxa in animals:
    with open('/homed/home/nzhang/JASPAR-tools/jaspartools/files/'+ taxa +'.pfam.hmm.json','r') as pfamF:
        pfamDic = json.loads(pfamF.read())

    # get the sequence of the pfam and store in dict
    for key,item in pfamDic.items():
        for lst in item:
            pfamID = lst[0]
            if pfamID in alignedD:
                alignedD[pfamID].extend(lst[1:])

# write the output as json (DBD hmm aligned)
with open("/homed/home/nzhang/JASPAR-tools/jaspartools/files/pfamAlignVisual/pfam_DBDAlignment.json",'w+') as outfile:
    json.dump(alignedD , outfile , sort_keys = True , indent = 4, separators = (",",": "))


    