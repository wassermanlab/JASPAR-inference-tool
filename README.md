# JASPAR profile inference tool
This repository contains the data and code used by the JASPAR profile inference tool.
</br>
For more information refer to the Supplementary Data of the JASPAR [2016](https://academic.oup.com/nar/article/44/D1/D110/2502663) and [2020](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz1001/5614568) manuscripts.

## News
11/11/2019 We improved the profile inference tool using our own implementation of the recently described [similarity regression](https://www.nature.com/articles/s41588-019-0411-1) method.

## Content
* The folder `examples` contains the sequences of two transcription factors and that of a negative example (_i.e._ [MTOR](https://www.uniprot.org/uniprot/P42345))
* The folder `files` contains the output from the script [`get_files.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/files/get_files.py), which creates [JSON files](https://en.wikipedia.org/wiki/JSON) for inference (_i.e._ *.json.gz) and downloads profile inference cut-offs on the percentage of sequence identity from [Cis-BP](http://cisbp.ccbr.utoronto.ca/), transcription factor sequences and DNA-binding domains (DBDs) from UniProt and [Pfam](https://pfam.xfam.org/), respectively, etc.
* The folder `models` contains the outputs from the scripts [`pairwise.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/models/pairwise.py), which creates the pairwise alignment of DBDs, and [`regression.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/models/regression.py), which builts the similarity (linear) regression models
* The script [`infer_profile.py`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/infer_profile.py) takes as input the `files` and `models` folders and a proteic sequence, in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format), and provides profile inferences
* The file [`environment.yml`](https://github.com/wassermanlab/JASPAR-profile-inference/blob/master/environment.yml) contains the conda environment (see Installation) to run the profile inference tool as of JASPAR 2020

The original scripts used for the publication of [JASPAR 2016](https://doi.org/10.1093/nar/gkv1176) have been placed in the folder [`version-1.0`](https://github.com/wassermanlab/JASPAR-profile-inference/tree/master/version-1.0).

## Dependencies
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [HMMER](http://hmmer.org/) (version 3+)
* [Python 3.7](https://www.python.org/download/releases/3.7/) with the [Biopython](http://biopython.org) (<1.74), [CoreAPI](http://www.coreapi.org), [glmnet](https://github.com/civisanalytics/python-glmnet), [NumPy](https://numpy.org/), [pandas](https://pandas.pydata.org/), [ProDy](http://prody.csb.pitt.edu/), [SciPy](https://www.scipy.org/), [scikit-learn](https://scikit-learn.org/stable/) and [tqdm](https://tqdm.github.io) libraries
* [`Tomtom`](http://meme-suite.org/doc/tomtom.html) as distributed in the [`MEME`](http://meme-suite.org/index.html) suite

Note that for running `infer_profile.py`, the Python dependencies CoreAPI, glmnet, ProDy and scikit-learn, and Tomtom are not required.

## Installation
All dependencies can be installed through the [conda](https://docs.conda.io/en/latest/) package manager:
```
conda env create -f ./environment.yml
```

## Usage
To illustrate the use of the profile inference tool, we provide an example for the [brown rat](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=10116&lvl=3&lin=f&keep=1&srchmode=1&unlock) transcription factor [Egr1](https://www.uniprot.org/uniprot/P08154):
* Create pairwise 
```
./infer_profile.py --fasta-file ./examples/Egr1.fa --files-dir ./files/ --models-dir ./models/ --latest --taxon vertebrates 
100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 1/1 [00:01<00:00,  1.66s/it]
Query   TF Name TF Matrix       E-value Query Start-End TF Start-End    DBD %ID Similarity Regression
EGR1_RAT        EGR1    MA0162.4        0.0     1-508   29-543  1.0     20.152809255049007
EGR1_RAT        EGR3    MA0732.1        4.53e-90        69-422  46-385  0.8840579710144927      18.268131142177957
EGR1_RAT        EGR2    MA0472.1        2.27e-74        62-398  45-424  0.9565217391304349      19.2874303040493
EGR1_RAT        EGR4    MA0733.1        1.12e-51        306-401 478-573 0.8115942028985507      None
```
For this example, the profiles .
* Create the genomic track
```
./scans2bigBed -c ./genomes/sacCer3/sacCer3.chrom.sizes -i ./tracks/sacCer3/ -o ./tracks/sacCer3.bb -t 4
```
TFBS predictions from the previous step are merged into a [bigBed track file](https://genome.ucsc.edu/goldenPath/help/bigBed.html). As scores (column 5), we use <i>p</i>-values from PWMScan (scaled between 0-1000, where 0 corresponds to <i>p</i>-value = 1 and 1000 to <i>p</i>-value ≤ 10-10). This allows for comparison of prediction confidence across TFBSs. Again, for this example, this step should be completed within a few minutes, while for larger genomes it can take a few hours.

**Important note:** both disk space and memory requirements for large genomes (*i.e.* danRer11, hg19, hg38 and mm10) are substantial. In these cases, we highly recommend allocating at least 1Tb of disk space and 512Gb of ram.

## Similarity regression (results)
```
Regressing AFT...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing AFT+AFT...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing AP2...
	*** ElasticNet: identity
		- lambdabest: 0.001
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 0.6456542290345665
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing AP2+AP2...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing AP2+B3...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing ARID...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing AT_hook+AT_hook+AT_hook+AT_hook...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing B3...
	*** ElasticNet: identity
		- lambdabest: 6.4565422903453475
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 7.2443596007485285
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing B3+B3...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing BAF1_ABF1...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing BrkDBD...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing CBFB_NFYA...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing CENP-B_N...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing CG-1...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing COE1_DBD...
	*** ElasticNet: identity
		- lambdabest: 114.81536214965975
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 588.8436553554222
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing CP2...
	*** ElasticNet: identity
		- lambdabest: 999.9999999997058
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 999.9999999997058
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing CUT+CUT+CUT+Homeodomain...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing CUT+Homeodomain...
	*** ElasticNet: identity
		- lambdabest: 0.37153522909712566
		- recall @ 75% precision: 50.00%
	*** ElasticNet: blosum62
		- lambdabest: 1.9952623149685564
		- recall @ 75% precision: 50.00%

Regressing Copper-fist...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing DM...
	*** ElasticNet: identity
		- lambdabest: 999.9999999997049
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 999.9999999997049
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing DM+DM...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing E2F_TDP...
	*** ElasticNet: identity
		- lambdabest: 3.235936569295725
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 5.370317963701544
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing E2F_TDP+E2F_TDP...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing Ets...
	*** ElasticNet: identity
		- lambdabest: 2.951209226665883
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 0.12302687708122553
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing FAR1...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing Forkhead...
	*** ElasticNet: identity
		- lambdabest: 2.511886431509161
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 3.548133892335137
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing GAGA_bind...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing GATA...
	*** ElasticNet: identity
		- lambdabest: 0.001
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 0.10715193052374997
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing GATA+GATA...
	*** ElasticNet: identity
		- lambdabest: 3.467368504524714
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 4.265795188015167
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing GCM...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing GCR1_C...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing HLH...
	*** ElasticNet: identity
		- lambdabest: 5.370317963701544
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 7.585775750290393
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing HLH+HLH...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing HMG_box...
	*** ElasticNet: identity
		- lambdabest: 13.489628825913801
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 22.908676527672828
		- recall @ 75% precision: 33.33%

Regressing HMG_box+HMG_box...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing HPD...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing HSF_DNA-bind...
	*** ElasticNet: identity
		- lambdabest: 0.004168693834703227
		- recall @ 75% precision: 80.00%
	*** ElasticNet: blosum62
		- lambdabest: 23.98832919018975
		- recall @ 75% precision: 80.00%

Regressing Homeodomain...
	*** ElasticNet: identity
		- lambdabest: 5.495408738575237
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 3.235936569295725
		- recall @ 75% precision: 1.23%

Regressing Homeodomain+CUT+CUT+CUT+Homeodomain...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing Homeodomain+Homeodomain...
	*** ElasticNet: identity
		- lambdabest: 0.001
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 5.370317963701544
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing IRF...
	*** ElasticNet: identity
		- lambdabest: 4.265795188015167
		- recall @ 75% precision: 100.00%
	*** ElasticNet: blosum62
		- lambdabest: 8.511380382022123
		- recall @ 75% precision: 100.00%

Regressing KilA-N...
	*** ElasticNet: identity
		- lambdabest: 0.001
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 10.964781961429676
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing LAG1-DNAbind...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing LAG1-DNAbind+LAG1-DNAbind...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing LOB...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing MADF_DNA_bdg...
	*** ElasticNet: identity
		- lambdabest: 50.11872336271567
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 75.85775750290021
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing MADF_DNA_bdg+MADF_DNA_bdg...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing MADF_DNA_bdg+MADF_DNA_bdg+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing MH1...
	*** ElasticNet: identity
		- lambdabest: 0.001
		- recall @ 75% precision: 50.00%
	*** ElasticNet: blosum62
		- lambdabest: 0.001
		- recall @ 75% precision: 50.00%

Regressing Myb_DNA-binding...
	*** ElasticNet: identity
		- lambdabest: 10.232929922805527
		- recall @ 75% precision: 4.00%
	*** ElasticNet: blosum62
		- lambdabest: 16.982436524613917
		- recall @ 75% precision: 8.00%

Regressing Myb_DNA-binding+Myb_DNA-binding...
	*** ElasticNet: identity
		- lambdabest: 3.162277660167836
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 2.3988329190190925
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing Myb_DNA-binding+Myb_DNA-binding+Myb_DNA-binding...
	*** ElasticNet: identity
		- lambdabest: 50.11872336271567
		- recall @ 75% precision: 100.00%
	*** ElasticNet: blosum62
		- lambdabest: 173.78008287489288
		- recall @ 75% precision: 100.00%

Regressing Myb_DNA-binding+Myb_DNA-binding+Myb_DNA-binding+Myb_DNA-binding...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing Myb_DNA-binding+Myb_DNA-binding+Myb_DNA-binding+Myb_DNA-binding+Myb_DNA-binding+Myb_DNA-binding...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing Myb_DNA-binding+Rap1-DNA-bind...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing NAM...
	*** ElasticNet: identity
		- lambdabest: 1.2302687708121949
		- recall @ 75% precision: 100.00%
	*** ElasticNet: blosum62
		- lambdabest: 8.511380382022123
		- recall @ 75% precision: 75.00%

Regressing NDT80_PhoG...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing P53...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing PAX...
	*** ElasticNet: identity
		- lambdabest: 0.001
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 0.001
		- recall @ 75% precision: 100.00%

Regressing PAX+Homeodomain...
	*** ElasticNet: identity
		- lambdabest: 3.630780547700379
		- recall @ 75% precision: 100.00%
	*** ElasticNet: blosum62
		- lambdabest: 12.88249551692874
		- recall @ 75% precision: 100.00%

Regressing Pou+Homeodomain...
	*** ElasticNet: identity
		- lambdabest: 0.05248074602497283
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 0.69183097091884
		- recall @ 75% precision: 20.00%

Regressing Pou+Pou+Homeodomain...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing RFX_DNA_binding...
	*** ElasticNet: identity
		- lambdabest: 999.9999999997058
		- recall @ 75% precision: 25.00%
	*** ElasticNet: blosum62
		- lambdabest: 999.9999999997058
		- recall @ 75% precision: 25.00%

Regressing RHD_DNA_bind...
	*** ElasticNet: identity
		- lambdabest: 0.6165950018613977
		- recall @ 75% precision: 100.00%
	*** ElasticNet: blosum62
		- lambdabest: 7.943282347241294
		- recall @ 75% precision: 100.00%

Regressing RRM_1...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing Runt...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing SAM_LFY...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing SAND...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing SBP...
	*** ElasticNet: identity
		- lambdabest: 1.1220184543017955
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 2.2908676527673952
		- recall @ 75% precision: 50.00%

Regressing SRF-TF...
	*** ElasticNet: identity
		- lambdabest: 0.5011872336272059
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 0.33884415613916047
		- recall @ 75% precision: 4.76%

Regressing STAT_bind...
	*** ElasticNet: identity
		- lambdabest: 16.982436524613917
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 12.589254117939138
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing STAT_bind+STAT_bind...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing STE...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing T-box...
	*** ElasticNet: identity
		- lambdabest: 21.877616239490866
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 23.442288153194195
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing T-box+HLH...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing TBP+TBP...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing TCP...
	*** ElasticNet: identity
		- lambdabest: 74.13102413007404
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 58.88436553554511
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing TCR+TCR...
	*** ElasticNet: identity
		- lambdabest: 8.511380382022123
		- recall @ 75% precision: 100.00%
	*** ElasticNet: blosum62
		- lambdabest: 7.762471166285434
		- recall @ 75% precision: 100.00%

Regressing TEA...
	*** ElasticNet: identity
		- lambdabest: 1.3182567385562052
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 6.165950018613675
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing TF_AP-2...
	*** ElasticNet: identity
		- lambdabest: 831.7637711024294
		- recall @ 75% precision: 100.00%
	*** ElasticNet: blosum62
		- lambdabest: 999.9999999997049
		- recall @ 75% precision: 100.00%

Regressing THAP...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing WRKY...
	*** ElasticNet: identity
		- lambdabest: 1.9054607179629404
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 1.9498445997577305
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing WRKY+WRKY...
	*** ElasticNet: identity
		- lambdabest: 0.30199517204016485
		- recall @ 75% precision: 100.00%
	*** ElasticNet: blosum62
		- lambdabest: 4.466835921508831
		- recall @ 75% precision: 100.00%

Regressing Zn_clus...
	*** ElasticNet: identity
		- lambdabest: 999.9999999997049
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 0.9772372209556672
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing Zn_clus+Zn_clus...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing bZIP_1...
	*** ElasticNet: identity
		- lambdabest: 4.67735141287114
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 4.570881896147929
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing bZIP_1+bZIP_1...
	*** ElasticNet: identity
		- lambdabest: 0.6606934480075045
		- recall @ 75% precision: 100.00%
	*** ElasticNet: blosum62
		- lambdabest: 0.15848931924609425
		- recall @ 75% precision: 100.00%

Regressing bZIP_1+bZIP_1+bZIP_1...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing zf-BED...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing zf-C2H2+bZIP_1...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing zf-C2H2+zf-C2H2...
	*** ElasticNet: identity
		- lambdabest: 0.0616595001861428
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 0.001
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing zf-C2H2+zf-C2H2+zf-C2H2...
	*** ElasticNet: identity
		- lambdabest: 32.35936569295566
		- recall @ 75% precision: 65.91%
	*** ElasticNet: blosum62
		- lambdabest: 64.56542290345031
		- recall @ 75% precision: 65.91%

Regressing zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2...
	*** ElasticNet: identity
		- lambdabest: 0.8511380382022541
		- recall @ 75% precision: 100.00%
	*** ElasticNet: blosum62
		- lambdabest: 1.4791083881679772
		- recall @ 75% precision: 100.00%

Regressing zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+Homeodomain+zf-C2H2+zf-C2H2+zf-C2H2...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2...
	*** ElasticNet: identity
		- lambdabest: 9.120108393557325
		- recall @ 75% precision: 11.11%
	*** ElasticNet: blosum62
		- lambdabest: 19.05460717962847
		- recall @ 75% precision: 11.11%

Regressing zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2...
	*** ElasticNet: identity
		- lambdabest: 112.20184543016853
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 81.28305161639034
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2...
	*** ElasticNet: identity
		- lambdabest: 8.912509381335727
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 2.511886431509161
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2...
	*** ElasticNet: identity
		- lambdabest: 999.9999999997049
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 999.9999999997049
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+Homeodomain+Homeodomain+zf-C2H2+Homeodomain+zf-C2H2...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2...
	*** ElasticNet: identity
		- lambdabest: 0.001
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 42.657951880149575
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2...
	*** ElasticNet: identity
		- lambdabest: 26.915348039263304
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 30.902954325129095
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2...
	*** ElasticNet: identity
		- lambdabest: 999.9999999997049
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 999.9999999997049
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2+zf-C2H2...
	*** ElasticNet: identity
	*** ElasticNet: blosum62
	*** could not train similarity regression!

Regressing zf-C4...
	*** ElasticNet: identity
		- lambdabest: 2.3442288153195343
		- could not reach 75% precision!
	*** ElasticNet: blosum62
		- lambdabest: 5.128613839912715
		- could not reach 75% precision!
	*** could not train similarity regression!

Regressing zf-Dof...
	*** ElasticNet: identity
		- lambdabest: 1.2589254117939757
		- recall @ 75% precision: 6.52%
	*** ElasticNet: blosum62
		- lambdabest: 2.7542287033377013
		- recall @ 75% precision: 2.17%
```