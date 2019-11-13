# JASPAR profile inference tool
This repository contains the data and code used by the JASPAR profile inference tool. For more information please refer to the supplementary data from JASPAR [2016](https://academic.oup.com/nar/article/44/D1/D110/2502663) and [2020](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz1001/5614568).

## News
01/09/2019 We have improved the profile inference tool by implementing our own [similarity regression](https://www.nature.com/articles/s41588-019-0411-1) method.

## Content
* The `examples` folder contains the sequences of two transcription factors (TFs) and one protein that is not a transcription factor, such as the human serine/threonine-protein kinase [mTOR](https://www.uniprot.org/uniprot/P42345)
* The `files` folder contains the output of the script [`get_files.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/files/get_files.py), which downloads TF sequences from [UniProt](https://www.uniprot.org/), DNA-binding domains (DBDs) from [Pfam](https://pfam.xfam.org/), and retrieves cut-offs on the DBD percentage of sequence identity from [Cis-BP](http://cisbp.ccbr.utoronto.ca/), etc.
* The `models` folder contains the similarity regression models created by calling the script [`pairwise.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/models/pairwise.py) followed by [`regression.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/models/regression.py)
* The script [`infer_profile.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/infer_profile.py) takes as input the folders `files` and `models`, plus one or more proteic sequences in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format) (_e.g._ a proteome), and infers DNA-binding profiles from JASPAR 
* The file [`environment.yml`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/environment.yml) contains the conda environment used to develop the profile inference tool for JASPAR 2020 (see installation)

The original scripts used for the publication of [JASPAR 2016](https://doi.org/10.1093/nar/gkv1176) have been placed in the folder [`version-1.0`](https://github.com/wassermanlab/JASPAR-profile-inference/tree/master/version-1.0).

## Dependencies
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [HMMER](http://hmmer.org/) (version ≥3.0)
* [Python 3.7](https://www.python.org/download/releases/3.7/) with the following libraries: [Biopython](http://biopython.org) (<1.74), [CoreAPI](http://www.coreapi.org), [glmnet](https://github.com/civisanalytics/python-glmnet), [NumPy](https://numpy.org/), [pandas](https://pandas.pydata.org/), [ProDy](http://prody.csb.pitt.edu/), [SciPy](https://www.scipy.org/), [scikit-learn](https://scikit-learn.org/stable/) and [tqdm](https://tqdm.github.io) 
* [Tomtom](http://meme-suite.org/doc/tomtom.html) as distributed in the [MEME](http://meme-suite.org/index.html) suite (version ≥5.0)

Note that for running `infer_profile.py`, neither Tomtom nor the CoreAPI, glmnet, ProDy and scikit-learn python packages are required.

## Installation
All dependencies can be installed through the [conda](https://docs.conda.io/en/latest/) package manager:
```
conda env create -f ./environment.yml
```

## Usage
To illustrate how the profile inference tool can be used, we provide an example for the [brown rat](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=10116&lvl=3&lin=f&keep=1&srchmode=1&unlock) TF [Egr1](https://www.uniprot.org/uniprot/P08154):
```
./infer_profile.py --fasta-file ./examples/Egr1.fa --files-dir ./files/ --models-dir ./models/ --latest
100%|█████████████████████████████████████████████████████████████████████| 1/1 [00:03<00:00,  3.54s/it]
Query   TF Name   TF Matrix   E-value   Query Start-End   TF Start-End   DBD %ID   Similarity Regression
Egr1    EGR1      MA0162.4    0.0     	1-508       	  29-543    	 1.0       20.153
Egr1    EGR3      MA0732.1    4.53e-90  69-422       	  46-385    	 0.884     18.268
Egr1    EGR2      MA0472.1    2.27e-74  62-398       	  45-424    	 0.957     19.287
Egr1    EGR4      MA0733.1    1.12e-51  306-401      	  478-573    	 0.812     None
```
The tool infers that the DNA-binding preferences of the `Query` (_i.e._ Egr1) are similar to those from the JASPAR TFs [EGR1](http://jaspar.genereg.net/matrix/MA0162.4/), [EGR2](http://jaspar.genereg.net/matrix/MA0472.1/), [EGR3](http://jaspar.genereg.net/matrix/MA0732.1/) and [EGR4](http://jaspar.genereg.net/matrix/MA0733.1/). The inference is based on the percentage of identical residues between the DBDs of their TFs (_i.e._ `DBD %ID`) and, except for EGR4, a linear regression model trained on the DBD pairwise sequence identities of JASPAR TFs with 3x `zf-C2H2` domains (_i.e._ `Similarity Regression`).
