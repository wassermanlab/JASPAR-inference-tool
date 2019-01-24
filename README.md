# JASPAR profile inference tool
For the [2016 release](https://doi.org/10.1093/nar/gkv1176) of [JASPAR](http://jaspar.genereg.net/), we have incorporated the capacity of inferring a JASPAR TF binding profile recognized by a DNA binding domain. Following a similar approach than the [Cis-BP database](http://cisbp.ccbr.utoronto.ca) (please refer to the original [Cell paper](https://doi.org/10.1016/j.cell.2014.08.009) for more details), for a given TF, the profile inference tool compares the DBD sequence of that TF to those of homologous TFs stored in JASPAR and, wherever possible, infers the binding profile(s) of that TF from the best compared JASPAR homologous TFs. Please refer to the JASPAR 2016 manuscript for more details.

## News
**TO BE UPDATED**

## Content
The repository is organized as follows:
* The `examples` folder contains a TF (*i.e.* `MAX.fa`) and a non-TF proteic sequence (*i.e.* `MTOR.fa`) in FASTA format
* The `files` folder contains the output from `make_files.py`: *i.e.* `domains.json`, `jaspar.json` and several BLAST formatted databases
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
* The `blastp` alignment between the query and JASPAR TF, including the the start and end amino acid positions, and the Expect value (E); and
* The % of identical residues between the query and the JASPAR TF DBDs