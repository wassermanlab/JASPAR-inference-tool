# JASPAR motif inference tool
```
some cool description
```

## Dependencies
The inference tool requires the following dependencies:
* [`BLAST+`](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [`Python 2.7`](https://www.python.org/download/releases/2.7/) with the [`Biopython`](http://biopython.org), [`CoreAPI`](http://www.coreapi.org) and [`UniProt`](https://github.com/boscoh/uniprot) libraries

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
python motif_inferrer.py -f $FILES_DIR -i examples/MAX.fa -s
```





