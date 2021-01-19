#!/usr/bin/env python

import argparse
import coreapi
from git import Repo
import json
import os
import pickle
# Download of Pfam/UniProt via RESTFUL API
from prody.database import pfam, uniprot
import re
import shutil
import subprocess
import sys
import time
from urllib.request import urlretrieve

# Defaults
out_dir = os.path.dirname(os.path.realpath(__file__))
root_dir = os.path.join(out_dir, os.pardir)

# Append JASPAR-profile-inference to path
sys.path.append(root_dir)

# Import globals
from __init__ import Jglobals

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns
    an {argparse} object.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--devel", action="store_true",
        help="development mode (uses hfaistos; default = False)")
    parser.add_argument("-o", default=out_dir, metavar="DIR",
        help="output directory (default = ./)")

    return(parser.parse_args())

def main():

    # Parse arguments
    args = parse_args()

    # Get files
    get_files(args.devel, os.path.abspath(args.o))

def get_files(devel=False, out_dir=out_dir):

    # Globals
    global client
    global codec
    global cwd
    global profiles_file_ext
    global clusters_file_ext
    global uniprot_file_ext
    client = coreapi.Client()
    codec = coreapi.codecs.CoreJSONCodec()
    cwd = os.getcwd()
    profiles_file_ext = ".profiles.json"
    clusters_file_ext = ".clusters.json"
    uniprot_file_ext = ".uniprot.json"

    # Create output dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Get JASPAR URL
    global jaspar_url
    jaspar_url = "http://jaspar.genereg.net/"
    if devel:
        jaspar_url = "http://hfaistos.uio.no:8002/"

    # Download Pfam DBD hidden Markov models
    __download_Pfam_DBD_HMMs(out_dir)

    # Download Cis-BP similarity regression models
    __download_CisBP_models(out_dir)

    # For each taxon...
    for taxon in Jglobals.taxons:

        # Download JASPAR files
        __download_JASPAR_profiles(taxon, out_dir)
        __get_profile_info(taxon, out_dir)
        # __get_motif_clusters(taxon, out_dir)

        # Download TF sequences from UniProt
        __download_UniProt_sequences(taxon, out_dir)

        # Format BLAST+ database
        __format_BLAST_database(taxon, out_dir)

def __download_Pfam_DBD_HMMs(out_dir=out_dir):

    # Skip if Pfam DBD file already exists
    pfam_DBD_file = os.path.join(out_dir, "pfam-DBDs.json")
    if not os.path.exists(pfam_DBD_file):

        # Initialize
        cutoffs = {}
        pfam_DBDs = {}
        url = "http://cisbp.ccbr.utoronto.ca/data/2.00/DataFiles/Bulk_downloads/EntireDataset/"
        cisbp_file = "TF_Information_all_motifs.txt.zip"

        # Change dir
        os.chdir(out_dir)

        # Skip if Cis-BP file already exists
        if not os.path.exists(cisbp_file):
            urlretrieve(os.path.join(url, cisbp_file), cisbp_file)

        # Get DBD/cut-off pairs
        cmd = "unzip -p %s | cut -f 11,13 | sort | uniq | grep -v DBDs" % \
            cisbp_file
        process = subprocess.run([cmd], shell=True, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)

        # For each line...
        for line in process.stdout.decode("utf-8").split("\n"):

            # Initialize
            line = line.split("\t")
            if len(line) != 2:
                continue
            pfam_ids, cutoff = line[0], line[1]

            # For each Pfam ID...
            for pfam_id in pfam_ids.split(","):

                # Skip if not Pfam ID
                if pfam_id == "UNKNOWN":
                    continue

                # Add Pfam ID cut-off
                cutoffs.setdefault(pfam_id, float(cutoff))

        # Create Pfam dir
        pfam_dir = "pfam-DBDs"
        if not os.path.exists(pfam_dir):
            os.makedirs(pfam_dir)

        # Change dir
        os.chdir(pfam_dir)

        # For each Pfam ID...
        for pfam_id in sorted(cutoffs):

            # Fetch MSA from Pfam
            attempts = 0
            while attempts < 5:
                try:
                    msa_file = pfam.fetchPfamMSA(pfam_id, alignment="seed")
                    break
                except:
                    # i.e. try again in 5 seconds
                    attempts += 1
                    time.sleep(5)

            # For each line...
            for line in Jglobals.parse_file(msa_file):

                m = re.search("^#=GF\sID\s+(\S+)$", line)
                if m:
                    pfam_id_std = m.group(1)

                m = re.search("^#=GF\sAC\s+(PF\d{5}).\d+$", line)
                if m:
                    pfam_ac = m.group(1)
                    break

            # HMM build
            hmm_file = "%s.hmm" % pfam_id_std
            cmd = "hmmbuild %s %s" % (hmm_file, msa_file)
            process = subprocess.run([cmd], shell=True, stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL)

            # HMM press
            cmd = "hmmpress -f %s" % hmm_file
            process = subprocess.run([cmd], shell=True, stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL)

            # Add Pfam
            pfam_DBDs.setdefault(pfam_ac, [pfam_id_std, cutoffs[pfam_id]])

            # Remove MSA file
            os.remove(msa_file)

        # Skip if HMM database of all DBDs already exists
        hmm_db = "all_DBDs.hmm"
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

        # Write
        Jglobals.write(
            pfam_DBD_file,
            json.dumps(pfam_DBDs, sort_keys=True, indent=4)
        )

        # Change dir
        os.chdir(out_dir)

        # Remove Cis-BP file
        if os.path.exists(cisbp_file):
            os.remove(cisbp_file)

    # Change dir
    os.chdir(cwd)

def __download_CisBP_models(out_dir=out_dir):

    # Skip if Cis-BP directory already exists
    cisbp_dir = os.path.join(out_dir, "cisbp")
    if not os.path.isdir(cisbp_dir):

        # Create Cis-BP dir
        os.makedirs(cisbp_dir)

        # Change dir
        os.chdir(cisbp_dir)

        # Clone repo to tmp dir
        url = "https://github.com/smlmbrt/SimilarityRegression.git"
        repo = Repo.clone_from(url, "tmp")

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
        os.makedirs(taxon_dir)

        # Move to taxon directory
        os.chdir(taxon_dir)

        # Initialize
        jaspar_file = "JASPAR2020_CORE_%s_redundant_pfms_jaspar.zip" % taxon
        if "hfaistos.uio.no:8002" in jaspar_url:
            jaspar_file = "JASPAR2020_CORE_%s_redundant_pfms_jaspar.zip" % taxon

        # Get JASPAR profiles
        if not os.path.exists(jaspar_file):
            urlretrieve(
                os.path.join(jaspar_url, "download", "CORE", jaspar_file),
                jaspar_file
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
            profiles_json_file,
            json.dumps(profiles, sort_keys=True, indent=4)
        )

# def __get_motif_clusters(taxon, out_dir=out_dir):

#     # Skip if taxon clusters JSON file already exists
#     clusters_json_file = os.path.join(out_dir, taxon + clusters_file_ext)
#     if not os.path.exists(clusters_json_file):

#         # Create taxon directory
#         clusters_dir = os.path.join(out_dir, taxon, "clusters")
#         os.makedirs(clusters_dir)

#         # Move to clusters directory
#         os.chdir(clusters_dir)

#         # Initialize
#         clusters = {}
#         leafs_file = os.path.join(
#             "JASPAR_2020_matrix_clustering_%s_tables" % taxon,
#             "leaf_to_cluster.tab"
#         )
#         clusters_file = "JASPAR_2020_matrix_clustering_%s_archive.zip" % taxon
#         url = os.path.join(
#             jaspar_url, "static", "clustering", "JASPAR_2020_clusters", taxon,
#             "radial_trees"
#         )

#         # Skip if Cis-BP file already exists
#         if not os.path.exists(clusters_file):
#             urlretrieve(os.path.join(url, clusters_file), clusters_file)

#         # Unzip
#         os.system("unzip -qq %s" % clusters_file)

#         # For each line...
#         for line in Jglobals.parse_file(leafs_file):
#             m = re.search("JASPAR_2020_\w+_m\d+_(MA\d{4})_(\d)\t(\d+)", line)
#             matrix_id = "%s.%s" % (m.group(1), m.group(2))
#             clusters.setdefault(matrix_id, [])
#             clusters[matrix_id].append(int(m.group(3)))

#         # Remove clusters directory
#         shutil.rmtree(clusters_dir)

#         # Write
#         Jglobals.write(
#             clusters_json_file,
#             json.dumps(clusters, sort_keys=True, indent=4)
#         )

#         # Change dir
#         os.chdir(cwd)

def __download_UniProt_sequences(taxon, out_dir=out_dir):

    # # Initialize
    # faulty_profiles = {
    #     "MA0024.1": ["Q01094"],
    #     "MA0046.1": ["P20823"],
    #     "MA0052.1": ["Q02078"],
    #     "MA0058.1": ["P61244"],
    #     "MA0098.1": ["P14921"],
    #     "MA0110.1": ["P46667"],
    #     "MA0138.1": ["Q13127"],
    #     "MA0328.1": ["P0CY08"],
    #     "MA0529.1": ["Q94513"],
    #     "MA0529.2": ["Q94513"]
    # }
    faulty_sequences = {
        "B9GPL8": [
            "MEEVGAQVAAPIFIHEALSSRYCDMTSMAKKHDLSYQSPNSQLQQHQFLQASREKNWNSK",
            "AWDWDSVDDDGLGLNLGGSLTSVEEPVSRPNKRVRSGSPGNGSYPMCQVDNCKEDLSKAK",
            "DYHRRHKVCQVHSKATKALVGKQMQRFCQQCSRFHPLTEFDEGKRSCRRRLAGHNRRRRK",
            "TQPEDVTSRLLLPGNPDMNNNGNLDIVNLLTALARSQGKTYLPMIDFYVPPFVLTNCPTV",
            "PDKDQLIQILNKINSLPLPMDLAAKLSNIASLNVKNPNQPYLGHQNRLNGTASSPSTNDL",
            "LAVLSTTLAASAPDALAILSQRSSQSSDNDKSKLPGPNQVTVPHLQKRSNVEFPAVGVER",
            "ISRCYESPAEDSDYQIQESRPNLPLQLFSSSPENESRQKPASSGKYFSSDSSNPIEERSP",
            "SSSPPVVQKLFPLQSTAETMKSEKMSVSREVNANVEGDRSHGCVLPLELFRGPNREPDHS",
            "SFQSFPYRGGYTSSSGSDHSPSSQNSDPQDRTGRIIFKLFDKDPSHFPGTLRTKIYNWLS",
            "NSPSEMESYIRPGCVVLSVYLSMPSASWEQLERNLLQLVDSLVQDSDSDLWRSGRFLLNT",
            "GRQLASHKDGKVRLCKSWRTWSSPELILVSPVAVIGGQETSLQLKGRNLTGPGTKIHCTY",
            "MGGYTSKEVTDSSSPGSMYDEINVGGFKIHGPSPSILGRCFIEVENGFKGNSFPVIIADA",
            "SICKELRLLESEFDENAVVSNIVSEEQTRDLGRPRSREEVMHFLNELGWLFQRKSMPSMH",
            "EAPDYSLNRFKFLLIFSVERDYCVLVKTILDMLVERNTCRDELSKEHLEMLYEIQLLNRS",
            "VKRRCRKMADLLIHYSIIGGDNSSRTYIFPPNVGGPGGITPLHLAACASGSDGLVDALTN",
            "DPHEIGLSCWNSVLDANGLSPYAYAVMTKNHSYNLLVARKLADKRNGQISVAIGNEIEQA",
            "ALEQEHVTISQFQRERKSCAKCASVAAKMHGRFLGSQGLLQRPYVHSMLAIAAVCVCVCL",
            "FFRGAPDIGLVAPFKWENLNYGTI"
        ]
    }

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

            # # Fix faulty profiles
            # if json_obj["matrix_id"] in faulty_profiles:
            #     json_obj["uniprot_ids"] = faulty_profiles[json_obj["matrix_id"]]

            # For each UniProt Accession...
            for uniacc in json_obj["uniprot_ids"]:

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

        # For each UniProt Accession...
        for uniacc in uniaccs:

            # Fix faulty sequences
            if uniacc in faulty_sequences:
                uniaccs[uniacc][1] = "".join(faulty_sequences[uniacc]) 
                continue

            # Get UniProt sequence
            u = uniprot.queryUniprot(uniacc)
            uniaccs[uniacc][1] = "".join(u["sequence   0"].split("\n"))

        # Write
        Jglobals.write(
            uniprot_json_file,
            json.dumps(uniaccs, sort_keys=True, indent=4)
        )

    # Change dir
    os.chdir(cwd)

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
            Jglobals.write(fasta_file, ">%s\n%s" % (uniacc, uniaccs[uniacc][1]))

        # Make BLAST+ database
        cmd = "makeblastdb -in %s -dbtype prot" % fasta_file
        process = subprocess.run([cmd], shell=True, stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()