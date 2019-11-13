# JASPAR profile inference tool
This repository contains the data and code used by the JASPAR profile inference tool. For more information refer to the Supplementary Data of the JASPAR [2016](https://academic.oup.com/nar/article/44/D1/D110/2502663) and [2020](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz1001/5614568) manuscripts.

## News
11/11/2019 We improved the profile inference tool using our own implementation of the recently described [similarity regression](https://www.nature.com/articles/s41588-019-0411-1) method.

## Content
* The folder `examples` contains the sequences of two transcription factors (TFs) and that of a negative example such as the serine/threonine-protein kinase [mTOR](https://www.uniprot.org/uniprot/P42345)
* The folder `files` contains the output from the script [`get_files.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/files/get_files.py), which creates [JSON files](https://en.wikipedia.org/wiki/JSON) for inference and downloads profile inference cut-offs on the percentage of sequence identity from [Cis-BP](http://cisbp.ccbr.utoronto.ca/), TF sequences and DNA-binding domains (DBDs) from [UniProt](https://www.uniprot.org/) and [Pfam](https://pfam.xfam.org/), respectively, etc.
* The folder `models` contains the outputs from the scripts [`pairwise.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/models/pairwise.py), which creates the pairwise alignment of DBDs, and [`regression.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/models/regression.py), which builts the similarity regression models
* The script [`infer_profile.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/infer_profile.py) takes as input the `files` and `models` folders and a proteic sequence, in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format), and provides profile inferences
* The file [`environment.yml`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/environment.yml) contains the conda environment (see Installation) used to develop the profile inference tool for JASPAR 2020

The original scripts used for the publication of [JASPAR 2016](https://doi.org/10.1093/nar/gkv1176) have been placed in the folder [`version-1.0`](https://github.com/wassermanlab/JASPAR-profile-inference/tree/master/version-1.0).

## Dependencies
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [HMMER](http://hmmer.org/) (version ≥3.0)
* [Python 3.7](https://www.python.org/download/releases/3.7/) with the [Biopython](http://biopython.org) (<1.74), [CoreAPI](http://www.coreapi.org), [glmnet](https://github.com/civisanalytics/python-glmnet), [NumPy](https://numpy.org/), [pandas](https://pandas.pydata.org/), [ProDy](http://prody.csb.pitt.edu/), [SciPy](https://www.scipy.org/), [scikit-learn](https://scikit-learn.org/stable/) and [tqdm](https://tqdm.github.io) libraries
* [Tomtom](http://meme-suite.org/doc/tomtom.html) as distributed in the [MEME](http://meme-suite.org/index.html) suite (version ≥5.0)

Note that for running `infer_profile.py`, neither Tomtom nor any of the following Python dependencies, CoreAPI, glmnet, ProDy and scikit-learn, are required.

## Installation
All dependencies can be installed through the [conda](https://docs.conda.io/en/latest/) package manager:
```
conda env create -f ./environment.yml
```

## Usage
To illustrate the use of the profile inference tool, we provide an example for the [brown rat](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=10116&lvl=3&lin=f&keep=1&srchmode=1&unlock) TF [Egr1](https://www.uniprot.org/uniprot/P08154):
```
./infer_profile.py --fasta-file ./examples/Egr1.fa --files-dir ./files/ --models-dir ./models/ --latest
100%|█████████████████████████████████████████████████████████████████████| 1/1 [00:03<00:00,  3.54s/it]
Query   TF Name TF Matrix       E-value Query Start-End TF Start-End    DBD %ID Similarity Regression
EGR1_RAT        EGR1    MA0162.4        0.0     1-508   29-543  1.0     20.153
EGR1_RAT        EGR3    MA0732.1        4.53e-90        69-422  46-385  0.884   18.268
EGR1_RAT        EGR2    MA0472.1        2.27e-74        62-398  45-424  0.957   19.287
EGR1_RAT        EGR4    MA0733.1        1.12e-51        306-401 478-573 0.812   None
```
Based on the sequence of Egr1, the tool inferred the JASPAR profiles of [EGR1](http://jaspar.genereg.net/matrix/MA0162.4/), [EGR2](http://jaspar.genereg.net/matrix/MA0472.1/) and [EGR3](http://jaspar.genereg.net/matrix/MA0732.1/) based on:
* The percentage of identical residues between the DBDs of Egr1 and those of the JASPAR TFs (_i.e._ `DBD %ID`); and
* A linear regression model trained on the pairwise sequence identities of JASPAR TFs with the DBD composition `3x zf-C2H2` (_i.e._ `Similarity Regression`)

Note that the profile [EGR4](http://jaspar.genereg.net/matrix/MA0733.1/) could not be inferred by Similarity Regression.