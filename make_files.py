#!/usr/bin/env python2.7
import os, sys, re
import coreapi
import hashlib
import json
import optparse
import subprocess
import uniprot

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

    parser = optparse.OptionParser("./%prog -b <blast_dir> [-o <output_dir>]")

    parser.add_option("-b", action="store", type="string", dest="blast_dir", help="Full path to BLAST+ bin directory (i.e. where \"makeblastdb\" is located; e.g. $BLAST_PATH/bin)", metavar="<blast_dir>")
    parser.add_option("-o", default="./", action="store", type="string", dest="output_dir", help="Output directory (default = ./)", metavar="<output_dir>")

    (options, args) = parser.parse_args()

    if options.blast_dir is None:
        parser.error("missing arguments: type option \"-h\" for help")

    return options

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()

    # Create output dir #
    if not os.path.exists(os.path.abspath(options.output_dir)):
        os.makedirs(os.path.abspath(options.output_dir))

    # Initialize #
    cwd = os.getcwd()
    matrix_ids = set()
    codec = coreapi.codecs.CoreJSONCodec()
    taxons = ["fungi", "insects", "nematodes", "plants", "vertebrates"]

    # For each taxon... #
    for taxon in taxons:
        # Skip if taxon profiles JSON file already exists #
        profiles_json_file = os.path.join(os.path.abspath(options.output_dir), ".%s.profiles.json" % taxon)
        if not os.path.exists(profiles_json_file):
            try:
                # Initialize #
                profiles = {}
                client = coreapi.Client()        
                response = client.get("http://jaspar.genereg.net/api/v1/taxon/%s/" % taxon)
                json_obj = json.loads(codec.encode(response))
                # While there are more results... #
                while json_obj['next'] is not None:
                    # For each profile... #
                    for profile in json_obj['results']:
                        # If CORE collection profile... #
                        if profile['collection'] == "CORE":
                            # Add profile #
                            profiles.setdefault(profile['matrix_id'], profile['name'])
                    # Go to next #
                    response = client.get(json_obj['next'])
                    json_obj = json.loads(codec.encode(response))
                # Write #
                functions.write(profiles_json_file, json.dumps(profiles, sort_keys=True, indent=4, separators=(',', ': ')))
            except:
                raise ValueError("Could not fetch %s profiles from JASPAR" % taxon)
        # Skip if taxon uniprot JSON file already exists #
        uniprot_json_file = os.path.join(os.path.abspath(options.output_dir), ".%s.uniprot.json" % taxon)
        if not os.path.exists(uniprot_json_file):
            try:
                # Initialize #
                uniaccs = {}
                client = coreapi.Client()
                profiles = json.loads("\n".join([line for line in functions.parse_file(profiles_json_file)]))
                # For each profile... #
                for profile in sorted(profiles):
                    # Get profile detailed info #
                    response = client.get("http://jaspar.genereg.net/api/v1/matrix/%s/" % profile)
                    json_obj = json.loads(codec.encode(response))
                    # For each UniProt Accession... #
                    for uniacc in json_obj['uniprot_ids']:
                        # Add uniacc #
                        uniaccs.setdefault(uniacc, [[], None])
                        if profile not in uniaccs[uniacc][0]: uniaccs[uniacc][0].append(profile)
                # Get UniProt data
                uniprot_data = uniprot.batch_uniprot_metadata(sorted(uniaccs), 'cache')
                # For each UniProt Accession... #
                for uniacc in sorted(uniprot_data):
                    uniaccs[uniacc][1] = uniprot_data[uniacc]['sequence']
                # Write #
                functions.write(uniprot_json_file, json.dumps(uniaccs, sort_keys=True, indent=4, separators=(',', ': ')))
            except:
                raise ValueError("Could not fetch %s sequences from UniProt" % taxon)
        # Skip if taxon FASTA file already exists #
        fasta_file = os.path.join(os.path.abspath(options.output_dir), "%s.fa" % taxon)
        if not os.path.exists(fasta_file):
            # Initialize #
            uniaccs = json.loads("\n".join([line for line in functions.parse_file(uniprot_json_file)]))
            # For each UniProt Accession... #
            for uniacc in sorted(uniaccs):
                # Skip if no sequence #
                if uniaccs[uniacc][1] is not None:
                    # Write #
                    functions.write(fasta_file, ">%s\n%s" % (uniacc, uniaccs[uniacc][1]))
            # Format db #
            try:
                process = subprocess.check_output([os.path.join(os.path.abspath(options.blast_dir), "makeblastdb"), "-in", fasta_file, "-dbtype", "prot"], stderr=subprocess.STDOUT)
            except:
                raise ValueError("Could not format BLAST database: %s" % fasta_file)

    # Skip if FASTA file already exists #
    fasta_file = os.path.join(os.path.abspath(options.output_dir), "sequences.fa")
    if not os.path.exists(fasta_file):
        # Initialize #
        uniq_seqs = {}
        # For each taxon... #
        for taxon in taxons:
            # Skip if taxon FASTA file does not exist #
            taxon_fasta_file = os.path.join(os.path.abspath(options.output_dir), "%s.fa" % taxon)
            if os.path.exists(taxon_fasta_file):
                # For each header, sequence... #
                for header, sequence in functions.parse_fasta_file(taxon_fasta_file):
                    # Add sequence #
                    uniq_seqs.setdefault(header, sequence)
        # For each header... #
        for header in sorted(uniq_seqs):
            # Write #
            functions.write(fasta_file, ">%s\n%s" % (header, uniq_seqs[header]))
        # Format db #
        try:
            process = subprocess.check_output([os.path.join(os.path.abspath(options.blast_dir), "makeblastdb"), "-in", fasta_file, "-dbtype", "prot"], stderr=subprocess.STDOUT)
        except:
            raise ValueError("Could not format BLAST database: %s" % fasta_file)

    # Skip if Cis-BP JSON file already exists #
    cisbp_json_file = os.path.join(os.path.abspath(options.output_dir), ".cisbp.json")
    if not os.path.exists(cisbp_json_file):
        # Initialize
        proteins = {}
        prot_features = {}
        tfs = {}
        tf_families = {}        
        # Create Cis-BP dir #
        cisbp_dir = os.path.join(os.path.abspath(options.output_dir), "cisbp")
        if not os.path.exists(cisbp_dir):
            os.makedirs(cisbp_dir)
        # Change directory #
        os.chdir(cisbp_dir)
        # Skip if TFs file already exists #
        if not os.path.exists("cisbp_1.02.tfs.sql"):
            # Download SQL files #
            os.system("curl --silent -O http://cisbp.ccbr.utoronto.ca/data/1.02/DataFiles/SQLDumps/SQLArchive_cisbp_1.02.zip")
            # Unzip #
            os.system("unzip -qq SQLArchive_cisbp_1.02.zip")
            # Remove SQL files #
            os.remove("SQLArchive_cisbp_1.02.zip")
            # For each ZIP file... #
            for zip_file in frozenset(os.listdir(os.getcwd())):
                # Skip non-zip files #
                if not zip_file.endswith(".zip"): continue
                # Unzip #
                os.system("unzip -qq %s" % zip_file)
                os.remove(zip_file)
        # Return to original directory #
        os.chdir(cwd)
        # For each line... #
        for line in functions.parse_file(os.path.join(cisbp_dir, "cisbp_1.02.prot_features.sql")):
            m = re.search("\('.+', '(.+)', '.+', \d+, \d+, '(.+)'\)", line)
            if m:
                prot_features.setdefault(m.group(1), set())
                prot_features[m.group(1)].add(m.group(2))
        # For each line... #
        for line in functions.parse_file(os.path.join(cisbp_dir, "cisbp_1.02.tfs.sql")):
            m = re.search("\('(.+)', '(.+)', '.+', '.+', '.+', '.+', '.+'\)", line)
            if m:
                tfs.setdefault(m.group(1), m.group(2))
        # For each line... #
        for line in functions.parse_file(os.path.join(cisbp_dir, "cisbp_1.02.tf_families.sql")):
            m = re.search("\('(.+)', '.+', '.+', \d+, (.+)\)", line)
            if m:
                tf_families.setdefault(m.group(1), m.group(2))
        # For each line... #
        for line in functions.parse_file(os.path.join(cisbp_dir, "cisbp_1.02.proteins.sql")):
            m = re.search("\('(.+)', '(.+)', '.+', '.+', '([A-Z]+)\W*'\)", line)
            if m:
                if m.group(1) not in prot_features: continue
                # Digest to MD5 #
                h = hashlib.new('md5')
                h.update(m.group(3))
                md5 = h.hexdigest() + m.group(3)[:4] + m.group(3)[-4:]
                proteins.setdefault(md5, [tf_families[tfs[m.group(2)]], []])
                # For each domain... #
                for domain in prot_features[m.group(1)]:
                    if domain not in proteins[md5][1]:
                        proteins[md5][1].append(domain)
        # Write #
        functions.write(cisbp_json_file, json.dumps(proteins, sort_keys=True, indent=4, separators=(',', ': ')))

    # Skip if JSON files already exist #
    domains_json_file = os.path.join(os.path.abspath(options.output_dir), "domains.json")
    jaspar_json_file = os.path.join(os.path.abspath(options.output_dir), "jaspar.json")
    if not os.path.exists(domains_json_file) or not os.path.exists(jaspar_json_file):
        # Initialize #
        domains = {}
        jaspar = {}
        cisbp = json.loads("\n".join([line for line in functions.parse_file(cisbp_json_file)]))
        # Remove JSON files #
        if os.path.exists(domains_json_file): domains_json_file
        if os.path.exists(jaspar_json_file): jaspar_json_file
        # For each taxon... #
        for taxon in taxons:
            # Initialize #
            profiles_json_file = os.path.join(os.path.abspath(options.output_dir), ".%s.profiles.json" % taxon)
            uniprot_json_file = os.path.join(os.path.abspath(options.output_dir), ".%s.uniprot.json" % taxon)
            profiles = json.loads("\n".join([line for line in functions.parse_file(profiles_json_file)]))
            uniaccs = json.loads("\n".join([line for line in functions.parse_file(uniprot_json_file)]))
            # For each UniProt Accession... #
            for uniacc in sorted(uniaccs):
                # Skip if no sequence #
                if uniaccs[uniacc][1] is None: continue
                # Digest to MD5 #
                h = hashlib.new('md5')
                h.update(uniaccs[uniacc][1])
                md5 = h.hexdigest() + uniaccs[uniacc][1][:4] + uniaccs[uniacc][1][-4:]
                # If sequence in Cis-BP... #
                if md5 in cisbp:
                    # Add to domains #
                    domains.setdefault(uniacc, [cisbp[md5][1], cisbp[md5][0]])
                    # For each profile... #
                    for profile in uniaccs[uniacc][0]:
                        # Add to JASPAR #
                        jaspar.setdefault(uniacc, [])
                        jaspar[uniacc].append([profile, profiles[profile]])
        # Write #
        functions.write(domains_json_file, json.dumps(domains, sort_keys=True, indent=4, separators=(',', ': ')))
        functions.write(jaspar_json_file, json.dumps(jaspar, sort_keys=True, indent=4, separators=(',', ': ')))
