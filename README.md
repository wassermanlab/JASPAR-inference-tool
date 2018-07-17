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
* [`Python 2.7`](https://www.python.org/download/releases/2.7/) with the [`Biopython`](http://biopython.org), [`CoreAPI`](http://www.coreapi.org), [`tqdm`](https://pypi.org/project/tqdm/) and [`UniProt`](https://github.com/boscoh/uniprot) libraries

## Usage
The script `profile_inferrer.py` infers one or more JASPAR TF binding profiles recognized by a sequence of interest. It requires the following inputs:
* The path to the `BLAST+` bin directory where `blastp` is located (option `-b`)
* The path to the `files` folder (option `-f`)
* A file containing one or more sequences in FASTA format (option `-i`)

Non-mandatory options include:
* The path to a "dummy" directory where to create temporary files (option `--dummy`, by default is set to the global temporary directory `/tmp`)
* The N parameter for the [Rost's sequence identity curve](https://doi.org/10.1093/protein/12.2.85) (option `-n`; by default is set to `0`)
* The name of a file to output the results (option `-o`; by default is set to the standard output stream)
* A taxonomic group (*i.e.* `fungi`, `insects`, `nematodes`, `plants`, or `vertebrates`) to which to restrict the profile inference (option `-t`; by default is set to ignore taxon restrictions)

Inference modes (for a more restricted profile inference):
* The option `-l` limits the inference of profiles to the latest JASPAR version
* The option `-s` ignores inferred profiles representing different TFs (*e.g.* hetetodimers)

As a usage example, the inferred JASPAR profiles for the `MAX` TF can be obtained as follows: `./profile_inferrer.py -b $BLAST_PATH -f ./files/ -i ./examples/MAX.fa`.

The script returns all inferred JASPAR TF binding profiles along with details regarding  the inference details, including:
* The `BLAST` alignment between the query and the JASPAR TF, including the aligned sequences, the start and end amino acid positions, and the Expect value (E); and
* The % of identical residues between the query and the JASPAR TF DBD

Note that for the `profile_inferrer.py` script to work with future versions of JASPAR other than 2018, first users would have to create a new `files` folder with the script `make_files.py`.