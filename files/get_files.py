#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import coreapi
from git import Repo
import io
import json
import os
import pickle
import re
import requests
import shutil
import subprocess
import sys
from urllib.request import urlretrieve

# Defaults
out_dir = os.path.dirname(os.path.realpath(__file__))
root_dir = os.path.join(out_dir, os.pardir)

# Append JASPAR-profile-inference to path
sys.path.append(root_dir)

# Import globals
from __init__ import Jglobals
from infer_profile import __make_seq_file, hmmscan, hmmalign

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns
    an {argparse} object.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--out-dir",
        default="./",
        help="output directory (default = ./)",
    )
    parser.add_argument(
        "-u", "--update",
        default=False,
        help="update files (default = False)",
        action="store_true"
    )

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Globals
    global client
    client = coreapi.Client()
    global codec
    codec = coreapi.codecs.CoreJSONCodec()
    global cwd
    cwd = os.getcwd()
    global devnull
    devnull = subprocess.DEVNULL
    global jaspar_url
    jaspar_url = "https://jaspar.elixir.no/"
    global clusters_file_ext
    clusters_file_ext = ".clusters.json"
    global pfam_file_ext
    pfam_file_ext = ".pfam.json"
    global profiles_file_ext
    profiles_file_ext = ".profiles.json"
    global uniprot_file_ext
    uniprot_file_ext = ".uniprot.json"
    out_dir = os.path.abspath(args.out_dir)

    # Prepare for update
    if args.update:
        os.chdir(out_dir)
        os.system("./rm_files.sh")
        os.chdir(cwd)

    # Get files
    get_files(out_dir)

def get_files(out_dir=out_dir):

    # Create output dir
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    # Download Pfam DBD hidden Markov models
    __download_Pfam_DBD_HMMs(out_dir)

    # Download Cis-BP similarity regression models
    __download_CisBP_models(out_dir)

    # For each taxon...
    for taxon in Jglobals.taxons:

        # Download JASPAR files
        __download_JASPAR_profiles(taxon, out_dir)
        __get_profile_info(taxon, out_dir)

        # Download TF sequences from UniProt
        __download_UniProt_sequences(taxon, out_dir)

        # Format BLAST+ database
        __format_BLAST_database(taxon, out_dir)

        # Get Pfam alignments
        __get_Pfam_alignments(taxon, out_dir)

def __download_Pfam_DBD_HMMs(out_dir=out_dir):

    # Skip if Pfam file already exists
    json_file = os.path.join(out_dir, "pfam.json")
    if not os.path.exists(json_file):

        # Initialize
        pfams = {}
        pfam_ids = set()

        # Create Pfam dir
        pfam_dir = os.path.join(out_dir, "pfam")
        if not os.path.isdir(pfam_dir):
            os.makedirs(pfam_dir)

        # Change dir
        os.chdir(pfam_dir)

        # Skip if Cis-BP file already exists
        cisbp_file = "TF_Information_all_motifs.txt.zip"
        if not os.path.exists(cisbp_file):
            url = "http://cisbp.ccbr.utoronto.ca/data/2.00/" + \
                "DataFiles/Bulk_downloads/EntireDataset/"
            urlretrieve(os.path.join(url, cisbp_file), cisbp_file)

        # Get DBD/cut-off pairs
        cmd = "unzip -p %s | cut -f 11 | sort | uniq | grep -v DBDs" % \
            cisbp_file
        process = subprocess.run([cmd], shell=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)

        # For each line...
        for line in process.stdout.decode("utf-8").split("\n"):

            # For each Pfam ID...
            for pfam_id in line.split(","):

                # Skip if not Pfam ID
                if pfam_id == "UNKNOWN" or pfam_id == "":
                    continue

                # Add Pfam ID
                pfam_ids.add(pfam_id)

        # Skip if Pfam file already exists
        pfam_file = "Pfam-A.seed.gz"
        if not os.path.exists(pfam_file):
            url = "https://ftp.ebi.ac.uk/pub/databases/" + \
                "Pfam/releases/Pfam32.0/"
            urlretrieve(os.path.join(url, pfam_file), pfam_file)

        # Parse seed files
        pfam_id = None
        pfam_ac = None
        msa = []
        cmd = "zless %s" % pfam_file
        process = subprocess.run([cmd], shell=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)

        # For each line...
        for line in process.stdout.decode("latin-1").split("\n"):

            # Add line
            msa.append(line)

            # Get Pfam ID
            m = re.search("^#=GF\sID\s+(\S+)$", line)
            if m:
                pfam_id = m.group(1)

            # Get Pfam accession
            m = re.search("^#=GF\sAC\s+(PF\d{5}).\d+$", line)
            if m:
                pfam_ac = m.group(1)

            if line.startswith("//"):

                # If Pfam ID in Cis-BP...
                if pfam_id in pfam_ids:

                    # Write MSA
                    msa_file = "%s.msa" % pfam_id
                    for m in msa:
                        Jglobals.write(msa_file, m)

                    # HMM build
                    hmm_file = "%s.hmm" % pfam_id
                    cmd = "hmmbuild %s %s" % (hmm_file, msa_file)
                    process = subprocess.run([cmd], shell=True, stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL)

                    # HMM press
                    cmd = "hmmpress -f %s" % hmm_file
                    process = subprocess.run([cmd], shell=True, stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL)

                    # Add Pfam
                    pfams.setdefault(pfam_ac, pfam_id)

                    # Remove MSA file
                    os.remove(msa_file)

                pfam_id = None
                pfam_ac = None
                msa = []

        # Skip if HMM database of all DBDs already exists
        hmm_db = "All.hmm"
        if not os.path.exists(hmm_db):

            # For each HMM file...
            for hmm_file in os.listdir("."):

                # Skip if not HMM file
                if not hmm_file.endswith(".hmm"): continue

                # Add HMM to database
                for line in Jglobals.parse_file(hmm_file):
                    Jglobals.write(hmm_db, line)

            # HMM press
            cmd = "hmmpress -f %s" % hmm_db
            process = subprocess.run([cmd], shell=True,
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Remove Cis-BP file
        if os.path.exists(cisbp_file):
            os.remove(cisbp_file)

        # Remove Pfam file
        if os.path.exists(pfam_file):
            os.remove(pfam_file)

        # Write
        Jglobals.write(
            json_file, json.dumps(pfams, sort_keys=True, indent=4)
        )

        # Change dir
        os.chdir(cwd)

def __download_CisBP_models(out_dir=out_dir):

    # Skip if Cis-BP directory already exists
    cisbp_dir = os.path.join(out_dir, "cisbp")
    if not os.path.isdir(cisbp_dir):

        # Create Cis-BP dir
        if not os.path.isdir(cisbp_dir):
            os.makedirs(cisbp_dir)

        # Change dir
        os.chdir(cisbp_dir)

        # Clone repo to tmp dir
        url = "https://github.com/smlmbrt/SimilarityRegression.git"
        Repo.clone_from(url, "tmp")

        # For each JSON file...
        for json_file in os.listdir(os.path.join("tmp", "SRModels")):

            # Skip
            if not json_file.endswith(".json"):
                continue

            # Copy
            shutil.copy(os.path.join("tmp", "SRModels", json_file), json_file)

        # Remove tmp dir
        shutil.rmtree("tmp")

    # Change dir
    os.chdir(cwd)

def __download_JASPAR_profiles(taxon, out_dir=out_dir):
        
    # Skip if taxon directory already exists
    taxon_dir = os.path.join(out_dir, taxon)
    if not os.path.exists(taxon_dir):

        # Create taxon directory
        if os.path.isdir(taxon_dir):
            shutil.rmtree(taxon_dir)            
        os.makedirs(taxon_dir)

        # Move to taxon directory
        os.chdir(taxon_dir)

        # Initialize
        jaspar_file = "JASPAR%s_CORE_%s_redundant_pfms_jaspar.zip" % \
            (Jglobals.version, taxon)

        # Get JASPAR profiles
        if not os.path.exists(jaspar_file):
            urlretrieve(
                os.path.join(
                    jaspar_url,
                    "download",
                    "data",
                    str(Jglobals.version),
                    "CORE",
                    jaspar_file
                ), jaspar_file
            )

        # Unzip
        os.system("unzip -qq %s" % jaspar_file)

        # Remove zip files
        os.remove("%s" % jaspar_file)

        # Change dir
        os.chdir(cwd)

def __get_profile_info(taxon, out_dir=out_dir):

    # Skip if taxon profiles JSON file already exists
    profiles_json_file = os.path.join(out_dir, taxon + profiles_file_ext)
    if not os.path.exists(profiles_json_file):

        # Initialize
        profiles = {}
        url = os.path.join(jaspar_url, "api", "v1", "taxon", taxon)
        response = client.get(url)
        json_obj = json.loads(codec.encode(response))

        # While there are more pages...
        while json_obj["next"] is not None:

            # For each profile...
            for profile in json_obj["results"]:

                # Add profiles from the CORE collection...
                if profile["collection"] == "CORE":
                    profiles.setdefault(profile["matrix_id"], profile["name"])

            # Go to next page
            response = client.get(json_obj["next"])
            json_obj = json.loads(codec.encode(response))

        # Do last page
        for profile in json_obj["results"]:

            # Add profiles from the CORE collection...
                if profile["collection"] == "CORE":
                    profiles.setdefault(profile["matrix_id"], profile["name"])

        # Write
        Jglobals.write(
            profiles_json_file, json.dumps(profiles, sort_keys=True, indent=4)
        )

def __download_UniProt_sequences(taxon, out_dir=out_dir):

    # Change dir
    os.chdir(out_dir)

    # Skip if pickle file already exists
    pickle_file = ".%s.uniaccs.pickle" % taxon
    if not os.path.exists(pickle_file):

        # Initialize
        uniaccs = {}

        # Load JSON file
        profiles_json_file = taxon + profiles_file_ext
        with open(profiles_json_file) as f:
            profiles = json.load(f)

        # For each profile...
        for profile in sorted(profiles):

            # Get profile detailed info
            url = os.path.join(jaspar_url, "api", "v1", "matrix", profile)
            response = client.get(url)
            json_obj = json.loads(codec.encode(response))

            # For each UniProt Accession...
            for uniacc in json_obj["uniprot_ids"]:

                # Skip
                if uniacc == "":
                    continue

                # Initialize
                uniacc = uniacc.strip(" ")
                uniaccs.setdefault(uniacc, [[], None])

                # Add uniacc
                if profile not in uniaccs[uniacc][0]:
                    uniaccs[uniacc][0].append(profile)

        # Write pickle file
        with open(pickle_file, "wb") as f:
            pickle.dump(uniaccs, f)

    # Skip if taxon uniprot JSON file already exists
    uniprot_json_file = taxon + uniprot_file_ext
    if not os.path.exists(uniprot_json_file):

        # Load pickle file
        with open(pickle_file, "rb") as f:
            uniaccs = pickle.load(f)

        # For each sublist...
        for uniaccs_100 in __split_in_chunks(list(uniaccs.keys())):
        
            # Get UniProt sequences
            query = "+OR+".join(f"accession:{acc}" for acc in uniaccs_100)
            url = "https://rest.uniprot.org/uniprotkb/" + \
                f"stream?&format=fasta&query={query}"
            request = requests.get(url)
            handle = io.StringIO(request.text)
            records = list(SeqIO.parse(handle, "fasta"))

            # For each UniProt Accession...
            for record in records:
                uniacc = record.id.split("|")[1]
                if uniacc in uniaccs:
                    uniaccs[uniacc][1] = str(record.seq)

        # Fix faulty sequences
        for uniacc in uniaccs:
            if not uniaccs[uniacc][1]:
                print("faulty UniProt accession", taxon, uniacc)
                url = f"https://www.uniprot.org/uniprot/{uniacc}.fasta"
                response = requests.get(url)
                uniaccs[uniacc][1] = "".join(response.text.split("\n")[1:])

        # Write
        Jglobals.write(
            uniprot_json_file, json.dumps(uniaccs, sort_keys=True, indent=4)
        )

    # Change dir
    os.chdir(cwd)

def __split_in_chunks(l, n=100): 
    for i in range(0, len(l), n):  
        yield l[i:i + n] 

def __format_BLAST_database(taxon, out_dir=out_dir):

    # Skip if taxon FASTA file already exists
    fasta_file = os.path.join(out_dir, "%s.fa" % taxon)
    if not os.path.exists(fasta_file):

        # Load JSON file
        uniprot_json_file = taxon + uniprot_file_ext
        with open(uniprot_json_file) as f:
            uniaccs = json.load(f)

        # For each UniProt Accession...
        for uniacc in sorted(uniaccs):
            seq = uniaccs[uniacc][1]
            Jglobals.write(fasta_file, ">%s\n%s" % (uniacc, seq))

        # Make BLAST+ database
        cmd = "makeblastdb -in %s -dbtype prot" % fasta_file
        process = subprocess.run(
            [cmd], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def __get_Pfam_alignments(taxon, out_dir=out_dir):

    # Skip if Pfam JSON file already exists
    pfam_json_file = os.path.join(out_dir, taxon + pfam_file_ext)
    if not os.path.exists(pfam_json_file):

        # Change dir
        os.chdir(out_dir)

        # Initialize
        pfams = {}
        seq_file = ".seq.fasta"
        hmm_db = os.path.join("pfam", "All.hmm")
        uniprot_json_file = taxon + uniprot_file_ext

        # Load JSON file
        with open(uniprot_json_file) as f:
            uniaccs = json.load(f)

        # For each uniacc...
        for u in uniaccs:

            # Initialize
            pfams.setdefault(u, [])

            # Make seq file
            seq = Seq(uniaccs[u][1], IUPAC.protein)
            record = SeqRecord(seq, id=u, name=u, description=u)
            __make_seq_file(record, seq_file)

            # For each DBD...
            for pfam_id_std, start, end, evalue in hmmscan(seq_file, hmm_db,
                non_overlapping_domains=True):

                # Initialize
                hmm_file = os.path.join("pfam", "%s.hmm" % pfam_id_std)

                # Make seq file
                sub_seq = seq[start:end]
                record = SeqRecord(sub_seq, id=u, name=u, description=u)
                __make_seq_file(record, seq_file)

                # Add DBDs
                alignment = hmmalign(seq_file, hmm_file)
                pfams[u].append((pfam_id_std, alignment, start+1, end, evalue))

        # Write
        Jglobals.write(
            pfam_json_file, json.dumps(pfams, sort_keys=True, indent=4)
        )

        # Remove seq file
        if os.path.exists(seq_file):
            os.remove(seq_file)

        # Change dir
        os.chdir(cwd)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()