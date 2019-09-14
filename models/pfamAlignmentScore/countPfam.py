# used to count dictionary number of entry
import json
import os

fileHeader = os.path.dirname(os.path.abspath(__file__))
#
filename = "/pfam_DBDAlignment.json"

with open(fileHeader+filename,'r') as fin:
    jD = json.loads(fin.read())

print(len(jD))