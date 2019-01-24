#!/usr/bin/env python

import os, re
import argparse
from bioservices import UniProt
import coreapi
import hashlib
import json
import shutil
import subprocess

# Import my functions
import functions

# Globals
taxons = ["fungi", "insects", "nematodes", "plants", "vertebrates"]

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line using argparse.
    """

    parser = argparse.ArgumentParser(description="creates files for profile inference tool.")

    # Optional args
    files_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "files")
    parser.add_argument("-o", metavar="OUTDIR", default=files_dir,
        help="Output directory (default = %s)" % files_dir)

    return parser.parse_args()

def make_files(out_dir=os.path.dirname(os.path.realpath(__file__))):

    # Initialize
    cwd = os.getcwd()
    matrix_ids = set()
    codec = coreapi.codecs.CoreJSONCodec()
    uniprot = UniProt(verbose=False, cache=False)

    # Create output dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # For each taxon...
    for taxon in taxons:
        # Skip if taxon profiles JSON file already exists
        profiles_json_file = os.path.join(
            out_dir, "%s.profiles.json" % taxon)
        if not os.path.exists(profiles_json_file):
            try:
                # Initialize
                profiles = {}
                client = coreapi.Client()        
                response = client.get(
                    "http://jaspar.genereg.net/api/v1/taxon/%s/" % taxon)
                json_obj = json.loads(codec.encode(response))
                # While there are more pages...
                while json_obj["next"] is not None:
                    # For each profile...
                    for profile in json_obj["results"]:
                        # If CORE collection profile...
                        if profile["collection"] == "CORE":
                            # Add profile
                            profiles.setdefault(profile["matrix_id"], profile["name"])
                    # Go to next page
                    response = client.get(json_obj["next"])
                    json_obj = json.loads(codec.encode(response))
                # Do last page
                for profile in json_obj["results"]:
                    # If CORE collection profile...
                    if profile["collection"] == "CORE":
                        # Add profile
                        profiles.setdefault(profile["matrix_id"], profile["name"])
                # Write
                functions.write(profiles_json_file, json.dumps(
                    profiles, sort_keys=True, indent=4, separators=(",", ": ")))
            except:
                raise ValueError("Could not fetch %s profiles from JASPAR" % taxon)

        # Skip if taxon uniprot JSON file already exists
        uniprot_json_file = os.path.join(out_dir, "%s.uniprot.json" % taxon)
        if not os.path.exists(uniprot_json_file):
            try:
                # Initialize
                uniaccs = {}
                client = coreapi.Client()
                # Load JSON file
                with open(profiles_json_file) as f:
                    profiles = json.load(f)
                # For each profile...
                for profile in sorted(profiles):
                    # Get profile detailed info
                    response = client.get(
                        "http://jaspar.genereg.net/api/v1/matrix/%s/" % profile)
                    json_obj = json.loads(codec.encode(response))
                    # Fix bugged cases
                    if json_obj["matrix_id"] == "MA0328.1":
                        json_obj["uniprot_ids"] = ["P0CY08"]
                    if json_obj["matrix_id"] == "MA0110.1":
                        json_obj["uniprot_ids"] = ["P46667"]
                    if json_obj["matrix_id"] == "MA0058.1":
                        json_obj["uniprot_ids"] = ["P61244"]
                    if json_obj["matrix_id"] == "MA0046.1":
                        json_obj["uniprot_ids"] = ["P20823"]
                    if json_obj["matrix_id"] == "MA0098.1":
                        json_obj["uniprot_ids"] = ["P14921"]
                    if json_obj["matrix_id"] == "MA0052.1":
                        json_obj["uniprot_ids"] = ["Q02078"]
                    if json_obj["matrix_id"] == "MA0024.1":
                        json_obj["uniprot_ids"] = ["Q01094"]
                    if json_obj["matrix_id"] == "MA0138.1":
                        json_obj["uniprot_ids"] = ["Q13127"]
                    # For each UniProt Accession...
                    for uniacc in json_obj["uniprot_ids"]:
                        # Strip
                        uniacc = uniacc.strip(" ")
                        # Add uniacc
                        uniaccs.setdefault(uniacc, [[], None])
                        if profile not in uniaccs[uniacc][0]:
                            uniaccs[uniacc][0].append(profile)
                # For each UniProt Accession...
                for uniacc in uniaccs:
                    # Get UniProt sequence
                    uniaccs[uniacc][1] = uniprot.get_fasta_sequence(uniacc)
                # Write
                functions.write(uniprot_json_file, json.dumps(
                    uniaccs, sort_keys=True, indent=4, separators=(",", ": ")))
            except:
                raise ValueError("Could not fetch %s sequences from UniProt" % taxon)

        # Skip if taxon FASTA file already exists
        fasta_file = os.path.join(out_dir, "%s.fa" % taxon)
        if not os.path.exists(fasta_file):
            # Load JSON file
            with open(uniprot_json_file) as f:
                uniaccs = json.load(f)
            # For each UniProt Accession...
            for uniacc in sorted(uniaccs):
                # Write
                functions.write(fasta_file,
                    ">%s\n%s" % (uniacc, uniaccs[uniacc][1]))
            # Create BLAST+ db
            try:
                process = subprocess.check_output([
                    "makeblastdb",
                    "-in", fasta_file,
                    "-dbtype", "prot"],
                    stderr=subprocess.STDOUT)
            except:
                raise ValueError("Could not create BLAST+ database: %s" % fasta_file)

    # Skip if Cis-BP JSON file already exists
    cisbp_json_file = os.path.join(out_dir, "cisbp.json")
    if not os.path.exists(cisbp_json_file):
        # Initialize
        proteins = {}
        prot_features = {}
        tfs = {}
        tf_families = {}
        # Create Cis-BP dir
        cisbp_dir = os.path.join(out_dir, "cisbp")
        if not os.path.exists(cisbp_dir):
            os.makedirs(cisbp_dir)
        # Change dir
        os.chdir(cisbp_dir)
        # Skip if TFs file already exists
        if not os.path.exists("cisbp_1.02.tfs.sql"):
            # Download SQL files
            os.system("curl --silent -O http://cisbp.ccbr.utoronto.ca/data/1.02/DataFiles/SQLDumps/SQLArchive_cisbp_1.02.zip")
            # Unzip
            os.system("unzip -qq SQLArchive_cisbp_1.02.zip")
            # Remove SQL files
            os.remove("SQLArchive_cisbp_1.02.zip")
            # For each ZIP file...
            for zip_file in frozenset(os.listdir(os.getcwd())):
                # Skip non-zip files
                if not zip_file.endswith(".zip"): continue
                # Unzip
                os.system("unzip -qq %s" % zip_file)
                os.remove(zip_file)
        # Return to original dir
        os.chdir(cwd)
        # Get protein features
        with open(os.path.join(cisbp_dir, "cisbp_1.02.prot_features.sql")) as f:
            # For each line...
            for line in f:
                m = re.search("\('.+', '(.+)', '.+', \d+, \d+, '(.+)'\)", line)
                if m:
                    prot_features.setdefault(m.group(1), set())
                    prot_features[m.group(1)].add(m.group(2))
        # Get TFs
        with open(os.path.join(cisbp_dir, "cisbp_1.02.tfs.sql")) as f:
            # For each line...
            for line in f:
                m = re.search("\('(.+)', '(.+)', '.+', '.+', '.+', '.+', '.+'\)", line)
                if m:
                    tfs.setdefault(m.group(1), m.group(2))
        # Get TF families
        with open(os.path.join(cisbp_dir, "cisbp_1.02.tf_families.sql")) as f:
            # For each line...
            for line in f:
                m = re.search("\('(.+)', '.+', '.+', \d+, (.+)\)", line)
                if m:
                    tf_families.setdefault(m.group(1), m.group(2))
        # Get proteins
        with open(os.path.join(cisbp_dir, "cisbp_1.02.proteins.sql")) as f:
            # For each line...
            for line in f:
                m = re.search("\('(.+)', '(.+)', '.+', '.+', '([A-Z]+)\W*'\)", line)
                if m:
                    if m.group(1) not in prot_features: continue
                    # Digest to MD5
                    h = hashlib.new("md5")
                    h.update(m.group(3).encode("utf-8"))
                    md5 = h.hexdigest() + m.group(3)[:4] + m.group(3)[-4:]
                    proteins.setdefault(md5, [tf_families[tfs[m.group(2)]], []])
                    # For each domain...
                    for domain in prot_features[m.group(1)]:
                        if domain not in proteins[md5][1]:
                            proteins[md5][1].append(domain)
        # Write
        functions.write(cisbp_json_file, json.dumps(
            proteins, sort_keys=True, indent=4, separators=(",", ": ")))
        # Remove Cis-BP dir
        shutil.rmtree(cisbp_dir)

    # Skip if JSON files already exist
    domains_json_file = os.path.join(out_dir, "domains.json")
    jaspar_json_file = os.path.join(out_dir, "jaspar.json")
    if not os.path.exists(domains_json_file) or not os.path.exists(jaspar_json_file):
        # Initialize
        domains = {}
        jaspar = {}
        # Load JSON file
        with open(cisbp_json_file) as f:
            cisbp = json.load(f)
        # Remove JSON files
        if os.path.exists(domains_json_file): domains_json_file
        if os.path.exists(jaspar_json_file): jaspar_json_file
        # For each taxon...
        for taxon in taxons:
            # Load JSON files
            profiles_json_file = os.path.join(
                out_dir, "%s.profiles.json" % taxon)
            with open(profiles_json_file) as f:
                profiles = json.load(f)
            uniprot_json_file = os.path.join(
                out_dir, "%s.uniprot.json" % taxon)
            with open(uniprot_json_file) as f:
                uniaccs = json.load(f)
            # For each UniProt Accession...
            for uniacc in sorted(uniaccs):
                # Skip if no sequence
                if uniaccs[uniacc][1] is None: continue
                # Digest to MD5
                h = hashlib.new("md5")
                h.update(uniaccs[uniacc][1].encode("utf-8"))
                md5 = h.hexdigest() + uniaccs[uniacc][1][:4] + uniaccs[uniacc][1][-4:]
                # If sequence in Cis-BP...
                if md5 in cisbp:
                    # Add to domains
                    domains.setdefault(uniacc, [cisbp[md5][1], cisbp[md5][0]])
                    # For each profile...
                    for profile in uniaccs[uniacc][0]:
                        # Add to JASPAR
                        jaspar.setdefault(uniacc, [])
                        jaspar[uniacc].append([profile, profiles[profile]])
        # Write
        functions.write(domains_json_file, json.dumps(
            domains, sort_keys=True, indent=4, separators=(",", ": ")))
        functions.write(jaspar_json_file, json.dumps(
            jaspar, sort_keys=True, indent=4, separators=(",", ": ")))

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Parse arguments
    args = parse_args()

    # Make files
    make_files(args.o)