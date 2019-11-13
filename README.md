# JASPAR profile inference tool
For the [2016 release](https://doi.org/10.1093/nar/gkv1176) of [JASPAR](http://jaspar.genereg.net/), we have incorporated the capacity of inferring a JASPAR TF binding profile recognized by a DNA binding domain. Following a similar approach than the [Cis-BP database](http://cisbp.ccbr.utoronto.ca) (please refer to the original [Cell paper](https://doi.org/10.1016/j.cell.2014.08.009) for more details), for a given TF, the profile inference tool compares the DBD sequence of that TF to those of homologous TFs stored in JASPAR and, wherever possible, infers the binding profile(s) of that TF from the best compared JASPAR homologous TFs. Please refer to the JASPAR 2016 manuscript for more details.

## News
**TO BE UPDATED**

## Content
The repository is organized as follows:
* The `examples` folder contains a TF (*i.e.* `MAX.fa`) and a non-TF proteic sequence (*i.e.* `MTOR.fa`) in FASTA format
* The `files` folder contains the output from `make_files.py`: *i.e.* `domains.json`, `jaspar.json` and several BLAST-formatted databases and JSONs
* The scripts `functions.py`, `make_files.py` and `profile_inferrer.py`

## Dependencies
The scripts for running the profile inference tool require the following dependencies:
* [`BLAST+`](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [`Python 2.7 or 3.x`](https://www.python.org) with the [`Biopython`](http://biopython.org), [`bioservices`](https://bioservices.readthedocs.io), [`CoreAPI`](http://www.coreapi.org) and [`tqdm`](https://pypi.org/project/tqdm/) libraries

## Usage
The script `profile_inferrer.py` infers one or more JASPAR TF binding profiles recognized by a sequence of interest. It requires the following inputs:
* A file containing one or more sequences in FASTA format
* The path to the `files` folder

Non-mandatory options include:
* The path to a "dummy" directory where to create temporary files (option `--dummy`, by default is set to the global temporary directory `/tmp`)
* The N parameter for the [Rost's sequence identity curve](https://doi.org/10.1093/protein/12.2.85) (option `-n`; by default is set to `5` to ensure ~99% of correctly assigned homologs)
* The name of a file to output the results (option `-o`; by default is set to the standard output stream)
* The number of threads to use (*i.e.* to speed-up the inference; option `--threads`)

JASPAR database options (for a customized profile inference):
* Options `--fungi`, `--insects`, `--nematodes`, `--plants` and `--vertebrates` restrict the inference of profiles to the specified taxons (by default, the inference tool uses profiles regardless of taxon)
* The option `-l` limits the inference of profiles to the latest JASPAR version

As a usage example, the inferred JASPAR profiles for the `MAX` TF can be obtained as follows: `./profile_inferrer.py ./examples/MAX.fa ./files/`.

The script returns (if any) the inferred JASPAR profiles for the specified TF(s) along with the details regarding the inference, including:
* The BLAST alignment between the query and JASPAR TF, including the the start and end amino acid positions, and the Expect value (E); and
* The % of identical residues between the query and the JASPAR TF DBDs

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