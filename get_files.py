#!/usr/bin/env python

import os, re
import argparse
from Bio import SearchIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import coreapi
import gzip
import json
from pathlib import Path
# Download of Pfam/UniProt via RESTFUL API
from prody.database import pfam, uniprot
import subprocess
import shutil
import wget
import zipfile # for unzip 

# Import globals
from __init__ import Jglobals

# Defaults
out_dir = os.path.dirname(os.path.realpath(__file__))

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns an {argparse} object.
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--devel", action="store_true", help="development mode (uses hfaistos; default = False)")
    parser.add_argument("-o", default=out_dir, help="output directory (default = ./)", metavar="DIR")

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
    global pfam_DBD_file
    global pfam_file_ext
    global profiles_file_ext
    global uniprot_file_ext
    client = coreapi.Client()
    codec = coreapi.codecs.CoreJSONCodec()
    cwd = os.getcwd()
    pfam_DBD_file = "pfam-DBDs.json"
    pfam_file_ext = ".pfam.json"
    profiles_file_ext = ".profiles.json"
    uniprot_file_ext = ".uniprot.json"

    # Create output dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Get JASPAR URL
    global jaspar_url
    jaspar_url = "http://jaspar.genereg.net/"
    if devel:
        jaspar_url = "http://hfaistos.uio.no:8002/"

    # Download Pfam DBDs
    _download_Pfam_DBDs(out_dir)
    exit(0)

    # For each taxon...
    for taxon in Jglobals.taxons:

        # Download JASPAR profiles
        _download_JASPAR_profiles(taxon, out_dir)

        # Download UniProt sequences
        _download_UniProt_sequences(taxon, out_dir)

        # Get Pfam alignments
        _get_Pfam_alignments(taxon, out_dir)
        exit(0)

def _download_Pfam_DBDs(out_dir=out_dir):

    # Initialize
    pfam_DBDs = {}
    pfam_ids = set()
    url = "http://cisbp.ccbr.utoronto.ca/data/2.00/DataFiles/Bulk_downloads/EntireDataset/"
    cisbp_file = "TF_Information_all_motifs.txt.zip"

    # Create Pfam dir
    pfam_dir = os.path.join(out_dir, "pfam-DBDs")
    if not os.path.exists(pfam_dir):
        os.makedirs(pfam_dir)

    # Skip if Pfam DBD file already exists
    pfam_DBD_file = os.path.join(out_dir, "pfam-DBDs.json")
    if not os.path.exists(pfam_DBD_file):

        # Change dir
        os.chdir(out_dir)

        # Skip if Cis-BP file already exists
        if not os.path.exists(cisbp_file):

            # Download TF info
            os.system("curl --silent -O %s%s" % (url, cisbp_file))

        # Get DBDs
        cmd = "unzip -p %s | cut -f 11 | sort | uniq | grep -v DBDs" % cisbp_file
        process = subprocess.run([cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # For each output line...
        for line in process.stdout.decode("utf-8").split("\n"):

            # Add Pfam IDs
            for pfam_id in line.split(","):
                    pfam_ids.add(pfam_id)
        print(pfam_ids)
        exit(0)

        # Change dir
        os.chdir(pfam_dir)

        # For each Pfam ID...
        for pfam_id in sorted(pfam_ids):

            try:

                # Fetch MSA from Pfam
                msa_file = pfam.fetchPfamMSA(pfam_id, alignment="seed")

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
                cmd = "hmmbuild %s %s &> /dev/null" % (hmm_file, msa_file)
                process = subprocess.run([cmd], shell=True)

                # HMM press
                cmd = "hmmpress -f %s &> /dev/null" % hmm_file
                process = subprocess.run([cmd], shell=True)

                # Add Pfam
                pfam_DBDs.setdefault(pfam_ac, pfam_id_std)

                # Remove MSA file
                os.remove(msa_file)

            except:
                print("\nCould not fetch MSA for id: %s\n" % pfam_id)

        # Change dir
        os.chdir(out_dir)

        # Write
        Jglobals.write(
            pfam_DBD_file,
            json.dumps(pfam_DBDs, sort_keys=True, indent=4, separators=(",", ": "))
        )

        # Remove Cis-BP file
        os.remove(cisbp_file)

    # Change dir
    os.chdir(cwd)

def _download_JASPAR_profiles(taxon, out_dir=out_dir):

    # Initialize
    url = os.path.join(jaspar_url, "api", "v1", "taxon", taxon)

    # Skip if taxon profiles JSON file already exists
    profiles_json_file = os.path.join(out_dir, taxon + profiles_file_ext)
    if not os.path.exists(profiles_json_file):

        # Initialize
        profiles = {}
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
            json.dumps(profiles, sort_keys=True, indent=4, separators=(",", ": "))
        )

def _download_UniProt_sequences(taxon, out_dir=out_dir):

    # Initialize
    faulty_profiles = {
        "MA0024.1": ["Q01094"],
        "MA0046.1": ["P20823"],
        "MA0052.1": ["Q02078"],
        "MA0058.1": ["P61244"],
        "MA0098.1": ["P14921"],
        "MA0110.1": ["P46667"],
        "MA0138.1": ["Q13127"],
        "MA0328.1": ["P0CY08"]
    }
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

    # Initialize
    profiles_json_file = os.path.join(out_dir, taxon + profiles_file_ext)

    # Skip if taxon uniprot JSON file already exists
    uniprot_json_file = os.path.join(out_dir, taxon + uniprot_file_ext)
    if not os.path.exists(uniprot_json_file):

        # Initialize
        uniaccs = {}

        # Load JSON file
        with open(profiles_json_file) as f:
            profiles = json.load(f)

        # For each profile...
        for profile in sorted(profiles):

            # Get profile detailed info
            url = os.path.join(jaspar_url, "api", "v1", "matrix", profile)
            response = client.get(url)
            json_obj = json.loads(codec.encode(response))

            # Fix faulty profiles
            if json_obj["matrix_id"] in faulty_profiles:
                json_obj["uniprot_ids"] = faulty_profiles[json_obj["matrix_id"]]

            # For each UniProt Accession...
            for uniacc in json_obj["uniprot_ids"]:

                # Initialize
                uniacc = uniacc.strip(" ")
                uniaccs.setdefault(uniacc, [[], None])

                # Add uniacc
                if profile not in uniaccs[uniacc][0]:
                    uniaccs[uniacc][0].append(profile)

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
            json.dumps(uniaccs, sort_keys=True, indent=4, separators=(",", ": "))
        )

def _get_Pfam_alignments(taxon, out_dir=out_dir):

    # Initialize
    global pfam_DBDs
    seq_file = ".seq.fasta"
    hmm_db = os.path.join(out_dir, "pfam-DBDs", "all_DBDs.hmm")
    pfam_DBD_file = os.path.join(out_dir, "pfam-DBDs.json")

    # Change dir
    os.chdir(out_dir)

    # Load JSON file
    with open(pfam_DBD_file) as f:
        pfam_DBDs = json.load(f)

    # For each taxon...
    for taxon in JTglobals.taxons:
        # Initialize
        uniprot_json_file = os.path.join(out_dir, taxon+uniprot_file_tail)
        # Skip if Pfam JSON file already exists
        pfam_json_file = os.path.join(out_dir, taxon+pfam_file_tail)
        if not os.path.exists(pfam_json_file):
            # Initialize
            pfams = {}
            # Load JSON file
            with open(uniprot_json_file) as f:
                uniaccs = json.load(f)
            # For each uniacc...
            for uniacc in uniaccs:
                # Initialize
                alignments = []
                pfams.setdefault(uniacc, [])
                # Make seq file
                seq = Seq(uniaccs[uniacc][1], IUPAC.protein)
                seq_record = SeqRecord(seq, id=uniacc, name=uniacc, description=uniacc)
                makeSeqFile(seq_record, seq_file)
                # For each DBD...
                for pfam_ac, start, end, evalue in hmmScan(seq_file, hmm_db, non_overlapping_domains=True):
                    # Initialize
                    hmm_file = os.path.join(out_dir, "pfam-DBDs", "%s.hmm" % pfam_ac)
                    # Make seq file
                    sub_seq = seq[start:end]
                    seq_record = SeqRecord(sub_seq, id=uniacc, name=uniacc, description=uniacc)
                    makeSeqFile(seq_record, seq_file)
                    # Add DBDs
                    alignment = hmmAlign(seq_file, hmm_file)
                    pfams[uniacc].append((pfam_ac, alignment, start+1, end, evalue))
            # Write
            JTglobals.write(pfam_json_file, json.dumps(
                pfams, sort_keys=True, indent=4, separators=(",", ": ")))

    # Change dir
    os.chdir(cwd)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()