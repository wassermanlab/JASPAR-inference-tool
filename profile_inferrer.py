#!/usr/bin/env python2.7
import os, re
from Bio.Blast import NCBIXML
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
import json
import math
from multiprocessing import Pool
import optparse
import shutil
import subprocess

# Import my functions #
import functions

#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("./%prog -f <files_dir> -i <input_file> [-b <blast_dir> -m <mmseqs_dir>] [--dummy=<dummy_dir> -n <n_parameter> -o <output_file> -t <taxon> --threads=<threads>] [-l -s]")

    parser.add_option("-f", action="store", type="string", dest="files_dir", help="Files directory (output directory from make_files.py)", metavar="<files_dir>")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (i.e. one or more sequences in FASTA format)", metavar="<input_file>")

    group = optparse.OptionGroup(parser, "Homology searching options")
    group.add_option("-b", action="store", type="string", dest="blast_dir", help="Full path to BLAST+ bin directory (i.e. where \"makeblastdb\" is located; e.g. $BLAST_PATH/bin)", metavar="<blast_dir>")
    group.add_option("-m", action="store", type="string", dest="mmseqs_dir", help="Full path to MMseqs2 bin directory (i.e. where \"mmseqs\" is located; e.g. $MMSEQS_PATH/bin)", metavar="<mmseqs_dir>")
    parser.add_option_group(group)
    
    group = optparse.OptionGroup(parser, "Non-mandatory options")
    group.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="<dummy_dir>")
    group.add_option("-n", default=0, action="store", type="int", dest="n_parameter", help="N parameter for the Rost's curve (e.g. n=5 ensures 99% of correctly assigned homologs; default = 0)", metavar="<n_parameter>")
    group.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="<output_file>")
    group.add_option("-t", action="store", dest="taxon", help="Taxonomic group (i.e. \"fungi\", \"insects\", \"nematodes\", \"plants\", or \"vertebrates\"; default = None)", metavar="<taxon>")
    group.add_option("--threads", default=1, action="store", type="int", dest="threads", help="Total number of cores to be used for the computation (default = 1)", metavar="<threads>")
    parser.add_option_group(group)
    
    group = optparse.OptionGroup(parser, "Inference modes")
    group.add_option("-l", "--latest", default=False, action="store_true", dest="latest", help="Latest mode (return the latest version of a profile; default = False)")
    group.add_option("-s", "--single", default=False, action="store_true", dest="single", help="Singleton mode (return profiles from a single TF; default = False)")

    parser.add_option_group(group)

    (options, args) = parser.parse_args()

    if options.files_dir is None or options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    if options.blast_dir is None and options.mmseqs_dir is None:
        parser.error("missing an homology searching method:\n\tspecify full path to either BLAST+ (option \"-b\") or MMseqs2 (option \"-m\") bin directory")
    
    if options.blast_dir is not None and options.mmseqs_dir is not None:
        parser.error("specify ONLY one homology searching method:\n\ti.e. full path to either BLAST+ (option \"-b\") or MMseqs2 (option \"-m\") bin directory")

    if options.taxon is not None:
        if options.taxon not in ["fungi", "insects", "nematodes", "plants", "vertebrates"]:
            parser.error("invalid taxon: %s\n\tvalid taxons include \"fungi\", \"insects\", \"nematodes\", \"plants\", and \"vertebrates\"" % options.taxon)

    return options

def is_alignment_over_Rost_sequence_identity_curve(identities, align_length, parameter=0):
    """
    This function evaluates whether an alignment is over {True} or 
    below {False} the Rost's sequence identity curve.
    
    @input:
    identities {int}
    align_length {int}
    parameter {int} N parameter in the curve (if > 0 more strict)
    @return: {boolean}
    
    """
    return identities >= get_Rost_ID_threshold(align_length, n=parameter)

def get_Rost_ID_threshold(L, n=0):
    """
    This function returns the Rost sequence identity threshold for a
    given alignment of length "L".

    @input:
    L {int} alignment length
    parameter {int} N parameter in the curve (if > 0 more strict)
    @return: {Decimal}
        
    """
    

    return n + (480*pow(L,float('-0.32')*(1+pow(float(repr(math.e)),float(repr(float(-L)/1000))))))

def get_alignment_identities(A, B):
    """
    This function returns the number of identities between a pair
    of aligned sequences {A} and {B}. If {A} and {B} have different
    lengths, returns None.

    @input:
    A {string} aligned sequence A (with residues and gaps)
    B {string} aligned sequence B (with residues and gaps)
    @return: {int} or None

    """
    if len(A) == len(B):
        return len([i for i in range(len(A)) if A[i] == B[i]])

    return None

def parallelize_inference(i):
    
    # Initialize #
    inferences = []
    
    # For each domain... #
    for domain in domains[homologs[i][1]][0]:
        # For each pairwise alignment... #
        for alignment in pairwise_alignment(sequences[alignments[i][0].strip()], domain):
            # If at least one DBD alignment passes threshold... #
            identities = get_alignment_identities(alignment[0], alignment[1])/float(len(domain))
            if identities >= float(domains[homologs[i][1]][1]):
                # For each JASPAR matrix... #
                for matrix, genename in jaspar[homologs[i][1]]:
                    # Infer matrix #
                    inferences.append([alignments[i][0], genename, matrix, homologs[i][6], "%s-%s" % (homologs[i][2], homologs[i][3]), "%s-%s" % (homologs[i][4], homologs[i][5]), identities])

                return inferences
    
    return inferences

def pairwise_alignment(A, B):
    """
    This function returns the alignments between a pair of sequences
    {A} and {B} using a dynamic programming algorithm.

    @input:
    A {string} sequence A
    B {string} sequence B
    @return: {Alignments}
        
    """
    
    try:
        # Gap open/extend parameters from EMBOSS needle #
        alignments = pairwise2.align.globalds(A, B, blosum62, -10.0, -0.5)
    except:
        raise ValueError("Needleman-Wunsch pairwise alignment failed:\n\tA: %s\n\tB: %s" % (A, B))

    return alignments

def clean_alignment(A, B):
    """
    Thiss function returns the alignment between a pair of sequences
    {A} and {B} ungapped (i.e. clean) with respect to {A}.

    @input:
    A {string} sequence A
    B {string} sequence B
    @return: cleaned {A} and {B}
        
    """

    matrix = [list(A), list(B)]
    # Transpose matrix #
    matrix_t = zip(*matrix)
    for i in sorted(range(len(matrix_t)), reverse=True):
        # If position is a gap... #
        if matrix_t[i][0] == "-": matrix_t.pop(i)
    # Transpose matrix #
    matrix = zip(*matrix_t)

    return "".join(matrix[0]), "".join(matrix[1])
    
#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Globals #
    global alignments, domains, homologs, sequences
    # Initialize #
    alignments = []
    homologs = []
    inferred_profiles = set()
    sequences = {}
    dummy_dir = os.path.join(os.path.abspath(options.dummy_dir), "%s.%s" % (os.path.basename(__file__), os.getpid()))
    if not os.path.exists(dummy_dir): os.makedirs(dummy_dir)
    alignments_tsv = os.path.join(dummy_dir, "alignments.tsv")
    alignments_msa = os.path.join(dummy_dir, "alignments.msa")
    alignments_xml = os.path.join(dummy_dir, "alignments.xml")
    domains_json = os.path.join(os.path.abspath(options.files_dir), "domains.json")
    domains = json.loads("\n".join([line for line in functions.parse_file(domains_json)]))
    jaspar_json = os.path.join(os.path.abspath(options.files_dir), "jaspar.json")
    jaspar = json.loads("\n".join([line for line in functions.parse_file(jaspar_json)]))
    jaspar_mmseqs_db = os.path.join(os.path.abspath(options.files_dir), "sequences.fa.db")
    jaspar_blast_db = os.path.join(os.path.abspath(options.files_dir), "sequences.fa")
    if options.taxon is not None:
        jaspar_mmseqs_db = os.path.join(os.path.abspath(options.files_dir), "%s.fa.db" % options.taxon)
        jaspar_blast_db = os.path.join(os.path.abspath(options.files_dir), "%s.fa" % options.taxon)
    query_mmseqs_db = os.path.join(dummy_dir, "query.db")
    query_mmseqs_ali = os.path.join(dummy_dir, "query.ali")
    results_file = os.path.join(dummy_dir, "results.tsv")

    # For each header, sequence... #
    for header, sequence in functions.parse_fasta_file(os.path.abspath(options.input_file)):
        sequences.setdefault(header, sequence)

    # If MMseqs2... #
    if options.mmseqs_dir is not None:
        # Create db #
        try:
            process = subprocess.check_output([os.path.join(os.path.abspath(options.mmseqs_dir), "mmseqs"), "createdb", os.path.abspath(options.input_file), query_mmseqs_db], stderr=subprocess.STDOUT)
        except:
            raise ValueError("Could not create MMseqs2 db: %s" % query_mmseqs_db)
        # Index db #
        try:
            process = subprocess.check_output([os.path.join(os.path.abspath(options.mmseqs_dir), "mmseqs"), "createindex", query_mmseqs_db, os.path.abspath(options.dummy_dir), "--threads", str(options.threads)], stderr=subprocess.STDOUT)
        except:
            raise ValueError("Could not index MMseqs2 db: %s" % query_db)
        # Search JASPAR db #
        try:
            process = subprocess.check_output([os.path.join(os.path.abspath(options.mmseqs_dir), "mmseqs"), "search", jaspar_mmseqs_db, query_mmseqs_db, query_mmseqs_ali, dummy_dir, "-s", "7.5", "--threads", str(options.threads)], stderr=subprocess.STDOUT)
        except:
            raise ValueError("Could not search JASPAR db!")
        # Reformat alignments #
        try:
            process = subprocess.check_output([os.path.join(os.path.abspath(options.mmseqs_dir), "mmseqs"), "convertalis", jaspar_mmseqs_db, query_mmseqs_db, query_mmseqs_ali, alignments_tsv, "--threads", str(options.threads)], stderr=subprocess.STDOUT)
        except:
            raise ValueError("Could not reformat alignments!")
        # Generate MSAs #
        try:
            process = subprocess.check_output([os.path.join(os.path.abspath(options.mmseqs_dir), "mmseqs"), "result2msa", jaspar_mmseqs_db, query_mmseqs_db, query_mmseqs_ali, alignments_msa, "--summarize", "--threads", str(options.threads)], stderr=subprocess.STDOUT)
        except:
            raise ValueError("Could not reformat alignments!")
        # For each line... #
        for line in functions.parse_file(alignments_msa):
            # Skip empty lines #
            if len(line) > 1:
                # If new target... #
                if line[0] == "#" or line[1] == "#":
                    target = None
                    target_seq = None
                # If header... #
                elif line.startswith(">"):
                    if target is None:
                        target = line[1:]
                    else: query = line[1:]
                # ... Else... #
                else:
                    if target_seq is None: target_seq = line
                    else:
                        query_seq = line
                        start = re.search("\W?\w", query_seq)
                        end = re.search("\w\W*$", query_seq)
                        alignments.append((query, query_seq[start.end() - 1:end.start() + 1], target, target_seq[start.end() - 1:end.start() + 1]))
        # For each line... #
        for line in functions.parse_file(alignments_tsv):
            # The file is formatted as a tab-separated list with 12 columns: (1,2) identifiers for query and target sequences/profiles, (3) sequence identity, (4) alignment length, (5) number of mismatches, (6) number of gap openings, (7-8, 9-10) domain start and end-position in query and in target, (11) E-value, and (12) bit score.
            target, query, perc_identity, alignment_length, mismatches, gaps, target_start, target_end, query_start, query_end, e_value, score = line.split("\t")
            # Skip if target does not have assigned domains #
            if target not in domains: continue
            # If structural homologs... #
            if is_alignment_over_Rost_sequence_identity_curve(int(int(alignment_length) - (int(mismatches) + int(gaps))), int(alignment_length), parameter=int(options.n_parameter)):
                homologs.append((query, target, query_start, query_end, target_start, target_end, float(e_value)))
            else :
                alignments.pop(len(homologs))
    # ... Else, use BLAST+... #
    else:
        # Exec blastp #
        try:
            process = subprocess.check_output([os.path.join(os.path.abspath(options.blast_dir), "blastp"), "-query", os.path.abspath(options.input_file), "-db", jaspar_blast_db, "-out", alignments_xml, "-outfmt", "5", "-num_threads", str(options.threads)], stderr=subprocess.STDOUT)
        except:
            raise ValueError("Could not search JASPAR db!")
        # Parse BLAST results #
        blast_records = NCBIXML.parse(open(alignments_xml))
        # For each BLAST record... #
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                # Skip if target does not have assigned domains #
                if alignment.hit_def not in domains: continue
                for hsp in alignment.hsps:
                    # If structural homologs... #
                    if is_alignment_over_Rost_sequence_identity_curve(get_alignment_identities(hsp.query, hsp.sbjct), len(hsp.query), parameter=int(options.n_parameter)):
                        homologs.append((blast_record.query, alignment.hit_def, hsp.query_start, hsp.query_end, hsp.sbjct_start, hsp.sbjct_end, float(hsp.expect)))
                        alignments.append((blast_record.query, hsp.query, alignment.hit_def, hsp.sbjct))
                        break

    # Parallelize inference #
    pool = Pool(options.threads)
    inferences = pool.map(parallelize_inference, [i for i in range(len(homologs))])

    # Write output #
    functions.write(results_file, "Query\tTF Name\tTF Matrix\tE-value\tQuery Start-End\tTF Start-End\tDBD %ID")
    # For each homolog... #
    for homolog in inferences:
        # For each inference... #
        for inference in sorted(homolog, key=lambda x: (x[3], -int(x[2][-1]))):
            # If latest mode... #
            if options.latest:
                if (inference[0], inference[2][:6]) in inferred_profiles: continue
            # If single mode... #
            if options.single:
                if "::" in inference[1]: continue
            # Write output #
            functions.write(results_file, "%s" % "\t".join(map(str, inference)))
            inferred_profiles.add((inference[0], inference[2][:6]))

    # Output #
    if options.output_file is not None:
        # Write output #
        shutil.copy(results_file, os.path.abspath(options.output_file))
    else:
        # For each line... #
        for line in functions.parse_file(results_file):
            # Write output #
            functions.write(None, line)
    # Remove files #
    shutil.rmtree(dummy_dir)
