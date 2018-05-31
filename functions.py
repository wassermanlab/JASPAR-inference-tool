import os, sys, re

def parse_file(file_name):
    """
    This function parses any file and yields lines one by one.
    
    @input:
    file_name {string}
    @return:
    line {string}

    """
 
    if os.path.exists(file_name):
        # Open file handle #
        try: f = open(file_name, "rt")
        except: raise ValueError("Could not open file %s" % file_name)
        # For each line... #
        for line in f:
            yield line.strip("\n")
        f.close()
    else:
        raise ValueError("File %s does not exist!" % file_name)

def parse_fasta_file(file_name, clean=True, proteinogenize=True):
    """
    This function parses any FASTA file and yields sequences one by one
    in the form header, sequence.

    @input:
    file_name {string}
    @return:
    line {list} header, sequence

    """

    # Initialize #
    header = ""
    sequence = ""
    # For each line... #
    for line in parse_file(file_name):
        if len(line) == 0: continue
        if line.startswith("#"): continue
        if line.startswith(">"):
            if sequence != "":
                if clean:
                    sequence = re.sub("\W|\d", "X", sequence)
                yield header, sequence
            m = re.search("^>(.+)", line)
            header = m.group(1)
            sequence = ""
        else:
            sequence += line.upper()
    if clean:
        # Convert non-amino acids to Xs #
        sequence = re.sub("[^ACDEFGHIKLMNPQRSTUVWY]", "X", sequence)
    if proteinogenize:
        # Convert selenocysteines to Cs #
        sequence = re.sub("U", "C", sequence)

    yield header, sequence

def write(file_name=None, line=None):
    """
    This function adds a {line} to file (or to stdout if no file
    was provided).

    @input:
    file_name {string}
    line {string}

    """

    if file_name is None: sys.stdout.write("%s\n" % line)
    else:
        with open(file_name, "a") as out_file:
            out_file.write("%s\n" % line)
