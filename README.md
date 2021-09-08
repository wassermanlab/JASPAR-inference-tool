# JASPAR profile inference tool
This repository contains the data and code used by the JASPAR profile inference tool. For more information please refer to the supplementary data from JASPAR [2016](https://academic.oup.com/nar/article/44/D1/D110/2502663) and [2020](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz1001/5614568).

## News
01/09/2021 We have updated the profile inference tool with new profiles for the 2022 release of JASPAR.
31/01/2021 We have updated the profile inference tool as described in the similarity regression [manuscript](https://www.nature.com/articles/s41588-019-0411-1).
~~01/09/2019 We have improved the profile inference tool by implementing our own [similarity regression](https://www.nature.com/articles/s41588-019-0411-1) method.~~

## Content
* The `conda` folder contains contains the [`environment.yml`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/conda/environment.yml) file used to develop the profile inference tool for JASPAR 2020 (see installation)
* The `examples` folder contains the sequences of two transcription factors (TFs) and one protein that is not a transcription factor, such as the human serine/threonine-protein kinase [mTOR](https://www.uniprot.org/uniprot/P42345)
* The `files` folder contains the output of the script [`get-files.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/files/get-files.py), which downloads TF sequences from [UniProt](https://www.uniprot.org/), DNA-binding domains (DBDs) from [Pfam](https://pfam.xfam.org/), retrieves infernece models from [Cis-BP](http://cisbp.ccbr.utoronto.ca/), etc.
* ~~The `models` folder contains the similarity regression models created by calling the script [`pairwise.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/models/pairwise.py) followed by [`regression.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/models/regression.py)~~
* The script [`infer-profile.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/infer-profile.py) takes as input ~~the folders `files` and `models`, plus~~ one or more proteic sequences in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format) (_e.g._ a proteome), and infers DNA-binding profiles from JASPAR 

The original scripts used for the publication of [JASPAR 2016](https://doi.org/10.1093/nar/gkv1176) have been placed in the folder [`version-1.0`](https://github.com/wassermanlab/JASPAR-profile-inference/tree/master/version-1.0).

## Dependencies
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [HMMER](http://hmmer.org/) (version ≥3.0)
* [Python 3](https://www.python.org/download/releases/3/) with the following libraries: [Biopython](http://biopython.org), [CoreAPI](http://www.coreapi.org), [GitPython](https://gitpython.readthedocs.io/en/stable/), [joblib](https://joblib.readthedocs.io/en/latest/), ~~[glmnet](https://github.com/civisanalytics/python-glmnet)~~, [NumPy](https://numpy.org/), [pandas](https://pandas.pydata.org/), [ProDy](http://prody.csb.pitt.edu/), ~~[SciPy](https://www.scipy.org/)~~, ~~[scikit-learn](https://scikit-learn.org/stable/)~~ and [tqdm](https://tqdm.github.io) 
* ~~The [RSAT matrix-clustering](http://pedagogix-tagc.univ-mrs.fr/rsat/matrix-clustering_form.cgi) tool~~
* ~~[Tomtom](http://meme-suite.org/doc/tomtom.html) as distributed in the [MEME](http://meme-suite.org/index.html) suite (version ≥5.0)~~

Note that for running `infer_profile.py`, the CoreAPI, GitPython, ~~glmnet, SciPy and scikit-learn,~~ and ProDy python packages are not required.

## Installation
All dependencies can be installed through the [conda](https://docs.conda.io/en/latest/) package manager:
```
conda env create -f ./conda/environment.yml
```

## Usage
To illustrate how the profile inference tool can be used, we provide an example for the [zebra fish](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?&id=7955) TF [egr1](https://www.uniprot.org/uniprot/P26632), and the [fission yeast](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?&id=4896) TF [tbp](https://www.uniprot.org/uniprot/P17871):
```
$ ./infer-profile.py --files-dir ./files/ --latest ./examples/egr1+tbp1.fa 
100%|████████████████████| 2/2 [00:08<00:00,  4.28s/it]
Query   TF Name TF Matrix       E-value Query Start-End TF Start-End    DBD %ID
sp|P26632|EGR1_DANRE    EGR1    MA0162.2        0.0     1-511   1-543   0.971
sp|P26632|EGR1_DANRE    EGR3    MA0732.1        6.81e-89        57-410  38-374  0.899
sp|P26632|EGR1_DANRE    Egr2    MA0472.1        5.95e-72        55-398  38-424  0.942
sp|P26632|EGR1_DANRE    EGR4    MA0733.1        9.11e-51        306-401 478-573 0.783
sp|P17871|TBP_SCHPO     SPT15   MA0386.1        8.18e-126       17-230  29-239  0.912
sp|P17871|TBP_SCHPO     TBP     MA0108.2        3.66e-109       8-230   114-337 0.771
```
The tool infers that the motif of `sp|P26632|EGR1_DANRE` should be similar to [EGR1](http://jaspar.genereg.net/matrix/MA0162.4/), [EGR2](http://jaspar.genereg.net/matrix/MA0472.1/), [EGR3](http://jaspar.genereg.net/matrix/MA0732.1/) and [EGR4](http://jaspar.genereg.net/matrix/MA0733.1/), and that the motif of `sp|P17871|TBP_SCHPO` should be similar to [SPT15](http://jaspar.genereg.net/matrix/MA0386.1/) and [TBP](http://jaspar.genereg.net/matrix/MA0108.2/).

#### As a Python module
```
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import importlib
infer_profile = importlib.import_module("infer-profile")

# Transcription factor Sox-3-B of Xenopus laevis
# https://www.uniprot.org/uniprot/Q5FWM3.fasta
seq = [
    "MYSMLDTDMKSPVQQSNALSGGPGTPGGKGNTSTPDQDRVKRPMNAFMVWSRGQRRKMAQ",
    "ENPKMHNSEISKRLGADWKLLSDSEKRPFIDEAKRLRAVHMKDYPDYKYRPRRKTKTLLK",
    "KDKYSLPGNLLAPGINPVSGGVGQRIDTYPHMNGWTNGAYSLMQEQLGYGQHPAMNSSQM",
    "QQIQHRYDMGGLQYSPMMSSAQTYMNAAASTYSMSPAYNQQSSTVMSLASMGSVVKSEPS",
    "SPPPAITSHTQRACLGDLRDMISMYLPPGGDAGDHSSLQNSRLHSVHQHYQSAGGPGVNG",
    "TVPLTHI"
]

# Load data
cisbp = infer_profile.__load_CisBP_models()
jaspar = infer_profile.__load_JASPAR_files_n_models()

# Infer profiles
seq_record = SeqRecord(Seq("".join(seq)), id="Sox-3-B")
inferred_profiles = infer_profile.infer_SeqRecord_profiles(
    seq_record, cisbp, jaspar, latest=True)

# Print
rows = [["Query", "TF Name", "TF Matrix", "E-value", "Query Start-End",
        "TF Start-End", "DBD %ID"]]
for inferred_profile in inferred_profiles:
    rows.append(inferred_profile)
for row in rows:
    print("\t".join(map(str, row)))

Query   TF Name TF Matrix       E-value Query Start-End TF Start-End    DBD %ID
Sox-3-B Sox3    MA0514.1        3.91e-129       1-307   1-375   0.942
Sox-3-B POU2F1::SOX2    MA1962.1        9.56e-115       1-307   1-317   0.913
Sox-3-B SOX2    MA0143.4        9.56e-115       1-307   1-317   0.913
Sox-3-B Pou5f1::Sox2    MA0142.1        7.27e-112       1-307   1-319   0.913
Sox-3-B Sox2    MA0143.1        7.27e-112       1-307   1-319   0.913
Sox-3-B Sox1    MA0870.1        6.37e-81        1-307   1-391   0.884
Sox-3-B SOX21   MA0866.1        6.71e-54        38-127  6-95    0.899
Sox-3-B SOX14   MA1562.1        2.24e-53        38-127  6-95    0.884
Sox-3-B D       MA0445.1        7.12e-45        32-145  134-239 0.826
Sox-3-B SOX15   MA1152.1        3.64e-44        22-117  34-126  0.812
Sox-3-B SRY     MA0084.1        1.11e-41        27-187  51-198  0.667
Sox-3-B Sox11   MA0869.1        4.58e-35        40-117  49-126  0.696
Sox-3-B SOX18   MA1563.1        3.38e-34        25-117  70-162  0.551
Sox-3-B SOX4    MA0867.1        1.36e-33        40-117  59-136  0.696
Sox-3-B Sox17   MA0078.1        1.6e-33 18-143  50-168  0.58
Sox-3-B SOX12   MA1561.1        1.11e-32        40-116  40-116  0.681
Sox-3-B SOX9    MA0077.1        9.43e-32        31-114  96-179  0.638
Sox-3-B SOX8    MA0868.1        1.56e-31        40-114  102-176 0.652
Sox-3-B SOX10   MA0442.1        1.74e-31        31-128  95-192  0.638
Sox-3-B Sox6    MA0515.1        6.99e-26        40-126  620-706 0.551
Sox-3-B Sox5    MA0087.1        1.03e-25        40-126  556-642 0.551
Sox-3-B SOX13   MA1120.1        4.97e-25        40-120  424-504 0.551
```