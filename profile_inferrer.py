#!/usr/bin/env python2.7
import os, re
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
from itertools import chain
import json
import math
import numpy
from multiprocessing import Pool
import optparse
import shutil
import subprocess
from tqdm import tqdm

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

    parser = optparse.OptionParser("./%prog -b <blast_dir> -f <files_dir> -i <input_file> [--dummy=<dummy_dir> -n <n_parameter> -o <output_file> -t <taxon> --threads=<threads>] [-l -s]")

    parser.add_option("-b", action="store", type="string", dest="blast_dir", help="Full path to BLAST+ bin directory (i.e. where \"makeblastdb\" is located; e.g. $BLAST_PATH/bin)", metavar="<blast_dir>")
    parser.add_option("-f", action="store", type="string", dest="files_dir", help="Files directory (output directory from make_files.py)", metavar="<files_dir>")
    parser.add_option("-i", action="store", type="string", dest="input_file", help="Input file (i.e. one or more sequences in FASTA format)", metavar="<input_file>")

    group = optparse.OptionGroup(parser, "Non-mandatory options")
    group.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="<dummy_dir>")
    group.add_option("-n", default=0, action="store", type="int", dest="n_parameter", help="N parameter for the Rost's curve (e.g. n=5 ensures 99% of correctly assigned homologs; default = 0)", metavar="<n_parameter>")
    group.add_option("-o", action="store", type="string", dest="output_file", help="Output file (default = stdout)", metavar="<output_file>")
    group.add_option("-t", action="store", dest="taxon", help="Taxonomic group (i.e. \"fungi\", \"insects\", \"nematodes\", \"plants\", or \"vertebrates\"; default = None)", metavar="<taxon>")
    group.add_option("--threads", default=1, action="store", type="int", dest="threads", help="Total number of cores to be used for the computation (default = 1)", metavar="<threads>")
    parser.add_option_group(group)
    
    group = optparse.OptionGroup(parser, "Inference modes")
    group.add_option("-l", "--latest", default=False, action="store_true", dest="latest", help="Latest mode (return the latest version of a profile; default = False)")
    # Need to need of a mode to filter "::" heterodimers
#    group.add_option("-m", "--sensitive", default=False, action="store_true", dest="sensitive", help="Sensitive mode (infer a profile only if the query satisfies the sequence identity requirement with ALL DBDs of the JASPAR TF; default = False)")
#    group.add_option("-s", "--sensitive", default=False, action="store_true", dest="sensitive", help="Sensitive mode (infer a profile only if the query satisfies the sequence identity requirement with ALL DBDs of the JASPAR TF; default = False)")

    parser.add_option_group(group)

    (options, args) = parser.parse_args()

    if options.blast_dir is None or options.files_dir is None or options.input_file is None:
        parser.error("missing arguments: type option \"-h\" for help")

    if options.taxon is not None:
        if options.taxon not in ["fungi", "insects", "nematodes", "plants", "vertebrates"]:
            parser.error("invalid taxon: %s\n\tvalid taxons include \"fungi\", \"insects\", \"nematodes\", \"plants\", and \"vertebrates\"" % options.taxon)

    return options

def parallelize_homology_search(i):
    
    # Initialize #
    homology_search = []
    targets = set()
    (header, sequence) = sequences[i]
    jaspar_db = os.path.join(os.path.abspath(options.files_dir), "sequences.fa")
    if options.taxon is not None:
        jaspar_db = os.path.join(os.path.abspath(options.files_dir), "%s.fa" % options.taxon)

    # Exec blastp #
    try:
        process = subprocess.Popen([os.path.join(os.path.abspath(options.blast_dir), "blastp"), "-db", jaspar_db, "-outfmt", "6"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        process.stdin.write(">%s\n%s" % (header, sequence))
        (blast_records, blast_errors) = process.communicate()
        # For each BLAST+ record... #
        for blast_record in blast_records.split("\n"):
            # Skip if not a BLAST+ record #
            blast_record = blast_record.split("\t")
            if len(blast_record) != 12: continue
            # Each BLAST+ record is formatted as a tab-separated list with 12 columns: (1,2) identifiers for query and target sequences, (3) percentage sequence identity, (4) alignment length, (5) number of mismatches, (6) number of gap openings, (7-8, 9-10) start and end-position in query and in target, (11) E-value, and (12) bit score.
            (query, target, perc_identity, alignment_length, mismatches, gaps, query_start, query_end, target_start, target_end, e_value, score) = blast_record
            # Skip if target exists #
            if target in targets: continue
            # If structural homologs... #
            if is_alignment_over_Rost_sequence_identity_curve(round(int(alignment_length) * float(perc_identity)/100), int(alignment_length), options.n_parameter):
                homology_search.append([i, target, "%s-%s" % (query_start, query_end), "%s-%s" % (target_start, target_end), float(e_value)])
    except:
#        raise ValueError("Could not search JASPAR db!")
        return []
    
    return homology_search

def is_alignment_over_Rost_sequence_identity_curve(identities, align_length, n_parameter=0):
    """
    This function evaluates whether an alignment is over {True} or 
    below {False} the Rost's sequence identity curve.
    
    @input:
    identities {int}
    align_length {int}
    parameter {int} N parameter in the curve (if > 0 more strict)
    @return: {boolean}
    
    """

    return identities >= get_Rost_ID_threshold(align_length, n=n_parameter)

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

def parallelize_profile_inference(i):
    
    # Initialize #
    profile_inference = []
    (header, sequence) = sequences[homologs[i][0]]

    # If domains... #
    if homologs[i][1] in domains:
        # For each domain... #
        for domain in domains[homologs[i][1]][0]:
            # Skip if inferred profile #
            if len(profile_inference) > 0: break
            # For each pairwise alignment... #
            for alignment in pairwise_alignment(sequence, domain):
                # If DBD alignment does not satisfy the threshold... #
                identities = get_alignment_identities(alignment[0], alignment[1])/float(len(domain))
                if identities < float(domains[homologs[i][1]][1]):
                    pass
#                    # If sensitive mode... #
#                    if options.sensitive:
#                        return []
                else:
                    # For each JASPAR matrix... #
                    for matrix, genename in jaspar[homologs[i][1]]:
                        # Infer matrix #
                        profile_inference.append([header, genename, matrix, homologs[i][4], homologs[i][2], homologs[i][3], identities])

    return profile_inference

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
        # Parameters from EMBOSS needle #
        alignments = pairwise2.align.globalds(A, B, blosum62, -10.0, -0.5)
    except:
#        raise ValueError("Needleman-Wunsch pairwise alignment failed:\n\tA: %s\n\tB: %s" % (A, B))
        return []

    return alignments

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
    
#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Globals #
    global domains, homologs, jaspar, options, sequences

    # Arguments & Options #
    options = parse_options()

    # Initialize #
    inferred_profiles = set()
    dummy_dir = os.path.join(os.path.abspath(options.dummy_dir), "%s.%s" % (os.path.basename(__file__), os.getpid()))
    if not os.path.exists(dummy_dir): os.makedirs(dummy_dir)
    domains_json = os.path.join(os.path.abspath(options.files_dir), "domains.json")
    domains = json.loads("\n".join([line for line in functions.parse_file(domains_json)]))
    jaspar_json = os.path.join(os.path.abspath(options.files_dir), "jaspar.json")
    jaspar = json.loads("\n".join([line for line in functions.parse_file(jaspar_json)]))
    results_file = os.path.join(dummy_dir, "results.tsv")

    # For each header, sequence... #
    sequences = [(header, sequence) for header, sequence in functions.parse_fasta_file(os.path.abspath(options.input_file))]

    # Parallelize homology search #
    pool = Pool(options.threads)
    homologs = list(set([tuple(i) for i in chain.from_iterable(tqdm(pool.imap(parallelize_homology_search, range(len(sequences))), desc="BLAST+ search", total=len(sequences)))]))
    pool.close()
    pool.join()
    
    # Parallelize DBD inference #
    pool = Pool(options.threads)
    inferences = list(set([tuple (i) for i in chain.from_iterable(tqdm(pool.imap(parallelize_profile_inference, range(len(homologs))), desc="DBD inference", total=len(homologs)))]))
    pool.close()
    pool.join()

    # For each inference... #
    functions.write(results_file, "Query\tTF Name\tTF Matrix\tE-value\tQuery Start-End\tTF Start-End\tDBD %ID")
    for inference in sorted(inferences, key=lambda x: (x[0], x[3], x[1], -float(x[2][2:]))):
        # If latest mode... #
        if options.latest:
            if (inference[0], inference[2][:6]) in inferred_profiles: continue
#        # If single mode... #
#        if options.single:
#            if "::" in inference[1]: continue
        # Write output #
        functions.write(results_file, "%s" % "\t".join(map(str, inference)))
        inferred_profiles.add((inference[0], inference[2][:6]))

    # Write output #
    if options.output_file is not None:
        shutil.copy(results_file, os.path.abspath(options.output_file))
    else:
        # For each line... #
        for line in functions.parse_file(results_file):
            functions.write(None, line)

    # Remove dummy dir #
    shutil.rmtree(dummy_dir)
