#!/usr/bin/env python

import os
import sys
import re

# Defaults
out_dir = os.path.dirname(os.path.realpath(__file__))
root_dir = os.path.join(out_dir, os.pardir)

# Append JASPAR-profile-inference to path
sys.path.append(root_dir)

# Import globals
from __init__ import Jglobals, CisBP2Pfam, ReadSRModel

cisbp = {}

for json_file in os.listdir("cisbp"):
    model = ReadSRModel(os.path.join("cisbp", json_file))
    cisbp.setdefault(CisBP2Pfam[model["Family_Name"]], model)

for pfam in cisbp:
    if cisbp[pfam]["Model.Class"] == "SimilarityRegression":
        hmm_file = os.path.join("pfam", f"{pfam}.hmm")
        for line in Jglobals.parse_file(hmm_file):
            m = re.search("^LENG\s+(\d+)", line)
            if m:
                print(
                    pfam,
                    int(m.group(1)) == len(cisbp[pfam]["SR.FeatureScales.mean"])
                )