# JASPAR profile inference tool
For the [2016 release](https://doi.org/10.1093/nar/gkv1176) of [JASPAR](http://jaspar.genereg.net/), we have incorporated the capacity of inferring a JASPAR TF binding profile recognized by a DNA binding domain. Following a similar approach than the [Cis-BP database](http://cisbp.ccbr.utoronto.ca) (please refer to the original [Cell paper](https://doi.org/10.1016/j.cell.2014.08.009) for more details), for a given TF, the profile inference tool compares the DBD sequence that TF to those of homologous TFs stored in JASPAR, and, wherever possible, infers that TF binding profile(s) from the best compared JASPAR homologous TFs. Please refer to the JASPAR 2016 manuscript for more details.

## Content
The repository is organized as follows:
* The `examples` folder contains a positive (*i.e.* `MAX.fa`) and a negative example (*i.e.* `MTOR.fa`) for testing
* The `files` folder contains the output from `make_files.py`:
  - `domains.json` contains UniProt Accessions, their associated DBD sequences and thresholds on the % of DBD sequence identity for profile inference
  - `jaspar.json` contains UniProt Accessions, their associated JASPAR matrix IDs and TF Names
  - `fungi.fa`, `insects.fa`, `nematodes.fa`, `plants.fa` and `vertebrates.fa` are taxon-specific TF databases formatted with `makeblastdb` for `blastp` searches (the `sequence.fa` is a general database that groups TFs from all taxons) 
* The scripts `functions.py`, `make_files.py` and `profile_inferrer.py` (*i.e.* the profile inference tool)

## Dependencies
The scripts for running the profile inference tool require the following dependencies:
* [`BLAST+`](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [`Python 2.7`](https://www.python.org/download/releases/2.7/) with the [`Biopython`](http://biopython.org), [`CoreAPI`](http://www.coreapi.org) and [`UniProt`](https://github.com/boscoh/uniprot) libraries

## Usage
The inference tool requires the output from `make_files.py` as prerequisites:
* 

## Configuration
```
explain that users must exec make_files.py and create symbolic links to both blastp and makeblastdb
```

## Usage
```
some cool description with two examples:

python motif_inferrer.py -i examples/MAX.fa
```

If we don't want to infer JASPAR profiles representing different TFs (e.g. hetetodimers):
```
python motif_inferrer.py -b $BLAST_PATH -f $FILES_DIR -i examples/MAX.fa -s
```





