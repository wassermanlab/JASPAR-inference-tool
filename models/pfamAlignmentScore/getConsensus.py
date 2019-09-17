
# this file gets all the DBD Pfam and gets the majority consensus sequence from them

import subprocess
import json
from pathlib import Path
import os

file_header = os.path.dirname(os.path.abspath(__file__))
pfamAhmm = "/sp4work/nicole/Pfam-A.hmm" #change if this become part of the tool
hmmf = file_header + "tmp.hmm"

def parseFASTA(fasta):
    #print("input parse: "+ fasta)
    content = fasta.split('\n')
    #print(content)
    seq = ''
    for line in content:
        if not '>' in line:
            seq += line
            #print(seq)
    #print("parsed fasta" + line)
    return(seq)

def hmmemit(pfamid):

    global pfamAhmm
    global hmmf
    
    #remove pre-existing file
    config = Path(hmmf)
    if config.is_file():
        subprocess.run(['rm','-f',hmmf])

    hmm_name = subprocess.check_output(['grep',pfamid,pfamAhmm],universal_newlines=True).split()[1]
    
    # run hmmfetch and save as tmp.hmm
    hmm = subprocess.run(['hmmfetch','-o',hmmf,pfamAhmm,hmm_name])

    # get and parse output from hmmemit
    seq = parseFASTA(subprocess.check_output(['hmmemit','-c',hmmf],universal_newlines=True))
    #print("emit: ",seq)
    return(seq)


with open("/homed/home/nzhang/JASPAR-tools/jaspartools/files/pfam_DBD_noVer.json",'r') as pfamIDf:
    pfamIDdic = json.loads(pfamIDf.read())

consensusD = {k:"" for k in pfamIDdic}

for key,item in pfamIDdic.items():
    seq = hmmemit(key)
    consensusD[key] = seq

with open("/homed/home/nzhang/JASPAR-tools/jaspartools/files/pfamAlignVisual/consensusDBDseq.json",'w+') as outf:
    json.dump(consensusD , outf , sort_keys = True , indent = 4, separators = (",",": "))


