#!/usr/bin/env python

import os, re
import argparse
from Bio import SearchIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
#from bioservices import UniProt
import coreapi
#import hashlib
import json
# Download of Pfam/UniProt via RESTFUL API
from prody.database import pfam, uniprot
import shutil
import subprocess
import wget
from pathlib import Path
import zipfile # for unzip 
import pip # for wget 
import gzip
import shutil

# Import from JASPAR tools
from jaspartools import JTglobals
from Bio import SearchIO

#-------------#
# Globals     #
#-------------#

codec = coreapi.codecs.CoreJSONCodec()
cwd = os.getcwd()

# File extensions
pfam_DBD_file = "pfam-DBDs.json"
pfam_file_tail = ".pfam.json"
profiles_file_tail = ".profiles.json"
uniprot_file_tail = ".uniprot.json"

#-------------#
# Functions   #
#-------------#

def makeSeqFile(seq_record, file_name=".seq.fa"):

    # Remove seq file if exists...
    if os.path.exists(file_name):
        os.remove(file_name)

    # Write
    JTglobals.write(file_name, seq_record.format("fasta"))

def readDomainsTab(file_name):

    # Initialize
    doms = []
    # From PMID:22942020;
    # A hit has equal probability of being in the same clan as a different clan when the
    # E-value is 0.01 (log10 = −2). When the E-value is 10−5, the probability that a sequence
    # belongs to the same clan is >95%.
    cutoff_mod = 1e-5
    # From CIS-BP paper;
    # We scanned all protein sequences for putative DNA-binding domains (DBDs) using the 81
    # Pfam (Finn et al., 2010) models listed in (Weirauch and Hughes, 2011) and the HMMER tool
    # (Eddy, 2009), with the recommended detection thresholds of Per-sequence Eval < 0.01 and
    # Per-domain conditional Eval < 0.01.
    cutoff_dom = 0.01

    # For each result...
    for res in SearchIO.parse(file_name, "hmmscan3-domtab"):
        # For each model...
        for mod in res.iterhits():
            # Skip poor models
            if mod.evalue > cutoff_mod: continue
            # For each domain...
            for dom in mod.hsps:
                # Skip poor domains
                if dom.evalue_cond > cutoff_dom:
                    continue
                # Append domain
                doms.append((mod.id, dom.query_start, dom.query_end, dom.evalue_cond))

    return(doms)

def getNonOverlappingDomains(doms):

    # Initialize
    nov_doms = []

    # Sort domains by e-value
    for dom in sorted(doms, key=lambda x: x[-1]):
        # Initialize
        doms_ov = False
        # For each non-overlapping domain...
        for nov_dom in nov_doms:
            # domains 1 & 2 overlap?
            # ---------1111111---------
            # -------22222-------------  True
            # ----------22222----------  True
            # -------------22222-------  True
            # -----22222---------------  False
            # ---------------22222-----  False
            if dom[1] < nov_dom[2] and dom[2] > nov_dom[1]:
                # Domains overlap
                doms_ov = True
                break
        # If domain does not overlap previous domains
        if not doms_ov:
            nov_doms.append(dom)

    return(nov_doms)

def hmmScan(seq_file, hmm_file, non_overlapping_domains=False):

    # Initialize
    out_file = ".out.txt"
    cmd = "hmmscan --domtblout %s %s %s &> /dev/null" % (out_file, hmm_file, seq_file)

    # Scan
    process = subprocess.call(cmd, shell=True)
    # Read domains
    domains = readDomainsTab(out_file)
    # Filter overlapping domains
    if non_overlapping_domains:
        domains = getNonOverlappingDomains(domains)
    # Yield domains one by one
    for pfam_ac, start, end, evalue in sorted(domains, key=lambda x: x[1]):

        yield(pfam_ac, start, end, evalue)

def readPSIBLAST(psiblast_alignment): # called in hmmAlign

    # Initialize
    alignment = ""

    # For each chunk...
    for chunk in psiblast_alignment.split("\n"):
        m = re.search("\s+(\S+)$", chunk)
        if m:
            alignment += m.group(1)

    return(alignment)

def hmmAlign(seq_file, hmm_file):

    # Initialize
    cmd = "hmmalign --outformat PSIBLAST %s %s" % (hmm_file, seq_file)

    # Align
    process = subprocess.check_output([cmd], shell=True, universal_newlines=True)

    return(readPSIBLAST(process))

##### functions added by nicole #######################

## variables needed ###################################

seq_file_name = 'temp.seq'

hmmf_base = "jaspartools/files/Pfam/"

hmmAseq = "tmp.hmm.ali"

pfamAhmm = "jaspartools/files/Pfam-A.hmm" # to be changed when it is downloaded

hmms = "" # to be changed (temporary hmmalign output)

pfamReg = "jaspartools/files/Pfam-A.regions.uniprot.tsv"

pfam_file_tail = ".pfam.json"

pfam_aligned_file_tail = ".pfam.hmm.json"

cisbp_tfinfo = "jaspartools/files/TF_Information_all_motifs.txt" # update if update cisbp file path

dbd_pfam_file = "jaspartools/files/pfam_DBD_withName.json" # stores the pfam file with 

wrong_names = {"FLO_LFY": "SAM_LFY", "DUF260": "LOB"}

#######################################################

def name_corr(name):

    global wrong_names

    if name in wrong_names:
        real = wrong_names[name]
        return(real)
    else:
        return(name)

def installWget(): #called when make
    if hasattr(pip,'main'):
        pip.main(['install','wget'])
    else:
        pip._internal.main(['install','wget'])

def installBioS(): #called when make
    if hasattr(pip,'main'):
        pip.main(['install','bioservices'])
    else:
        pip._internal.main(['install','bioservices'])

def getPfamregions(): #get the pfam region file. ONLY CALL ONCE!!!
    global pfamReg
    config = Path(pfamReg)
    
    if config.is_file():
        return # if file already exist, don't download
    
    url = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.regions.uniprot.tsv.gz"
    
    try:
        zipped_pfamReg = wget.download(url)
    except:
        os.system("curl --silent -O %s" % url)
        zipped_pfamReg = "Pfam-A.regions.uniprot.tsv.gz"
    
    with gzip.open(zipped_pfamReg,'rb') as f_in:
        with open(pfamReg,'wb') as f_out:
            shutil.copyfileobj(f_in,f_out)

# commented out bc memory error
# def makePfamReg(): #called only by createPfamJson
#     with open(pfamReg,'r') as pf:
#         pfamR = pf.read()
#     return(pfamR)


def createPfamJson(unipJson,pfamOut): # Call when you have uniprot json names and want pfams
    # pass in the uniprot id json, will fetch the pfam id
    with open(unipJson,'r') as unip:
        unip_file = unip.read()
        obj = json.loads(unip_file)
    
    #pfamdic have uniprot ids as keys, and will be storing pfam ids as values
    pfamdic = {k:[] for k in obj} # create another dic with the same keys (but empty values)

    with open(pfamReg,'r') as pf: # now don't need to store, just need to read
        for line in pf:	
            line = line.replace("\n","")
            wdr = line.split("\t")
            #print(wdr)
            if wdr[0] in pfamdic and wdr[4] not in pfamdic[wdr[0]]:
                #print(wdr[0])
                #print(pfamdic[wdr[0]])
                pfamdic[wdr[0]].append(wdr[4])
                #print(pfamdic[wdr[0]])

#    with open (pfamOut,'w') as fp:
#        json.dump(pfamdic,fp)
    # Write
    JTglobals.write(pfamOut, json.dumps(
        pfamdic, sort_keys=True, indent=4, separators=(",", ": ")))

def getPfam(): # get the pfam and unzips ONLY CALL ONCE!!!!
    global pfamAhmm # assume know file name already 
    config = Path(pfamAhmm)
    
    if config.is_file():
        return # if file already exist, don't download

    # file download
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz'
    try:
        zipped_pfamAhmm = wget.download(url) #pfamAhmm contains the file name
    except:
        os.system("curl --silent -O %s" %url)
        zipped_pfamAhmm = "Pfam-A.hmm.gz"
    
    # file unzip
    with gzip.open(zipped_pfamAhmm,'rb') as f_in:
        with open(pfamAhmm,'wb') as f_out:
            shutil.copyfileobj(f_in,f_out)

def readHMMOutput(filename,seq,hmmfile): # pass complete file name; called in runhmm
    print("read HMM output") # test where it was run
    #this function takes the hmm file (in domain tbl) and returns the aligned seq
    seqList = list()
    hmmr = SearchIO.read(filename,"hmmscan3-domtab")
    for hsp in hmmr.hsps:
        st = hsp.query_start
        ed = hsp.query_end
        hitsq = seq[st:ed]
        aligned = hmmAlign(hitsq,hmmfile)
        seqList.append(aligned) # append sequence aligned already
    return(seqList) # return as a whole string
    
def runhmm(pfam,jseq): # called in runHMM_fromPfam
    global pfamAhmm
    # but the specific decimals need to be known
    try: # when the pfam id not found..
        hmm_name = subprocess.check_output(['grep',pfam,pfamAhmm],universal_newlines=True).split()[1]
    except: # this is more sustainable (can have shared and don't have to delete files everytime)
        print(pfam+' was not in the database')
        return([''])

    # this defines the file name (hmm file)
    global hmmf_base
    hmmf = hmmf_base + pfam + ".hmm"

    # if the file does not exist
    config = Path(hmmf)
    if config.is_file():
        pass
    else:
        # get pfam ID (and then fetch the hmm file)
    
        hmm = subprocess.run(['hmmfetch','-o',hmmf,pfamAhmm,hmm_name])

        # make hmmf ready for hmmscan
        subprocess.run(['hmmpress',hmmf])

    # hmmscan
    global hmms
    hmms = "temp.dom.out"
    global seq_file_name
    subprocess.run(['hmmscan','--cpu','10','--domtblout',hmms,hmmf,seq_file_name])

    # parse the result from hmmscan and pass to readHMMOutput
    return(readHMMOutput(hmms,jseq,hmmf))

def runHMM_fromPfam(uniprot_seqfile,pfam_jsonfile,pfam_aligned_file):
    global seq_file_name #declare global to use
    # read in the sequence
    with open(uniprot_seqfile,'r') as seq:
        seqd = seq.read()
        seqdic = json.loads(seqd)
    
    # read in the pfam
    with open(pfam_jsonfile,'r') as pfamJ:
        jsd = pfamJ.read()
        jdic = json.loads(jsd)
        for key,item in jdic.items():
            jseq = '>'+key+'\n'+seqdic[key][1]
            # write the sequence down
            fn1 = open(seq_file_name,"w+")
            fn1.write(jseq)
            fn1.close()
            count = 0 # for counting the pfams in item
            for i in item: #len(item), i is the pfam id
                align = runhmm(i,seqdic[key][1]) #pass in the Pfam string and the sequence string
                #print("working on "+seqdic[key][1])
                replacement = list()
                replacement.append(i)
                #print(replacement)
                #print(type(replacement))
                #print(align)
                #print(type(align))
                replacement = replacement + align
                #print(replacement)
                item[count] = replacement
                count += 1
                #print(psiblast).split()[1]
                # insert psiblast into the list, right after the 
    with open(pfam_aligned_file,'w+') as wf:
        json.dump(jdic,jwrite,sort_keys=True,indent=4,separators=(",",": "))
        jwrite.write("\n")

def hmmfetch(nid): # called in the get_DBD_TF(), to get the accession for Pfam ids
    global pfamAhmm
    try:
        out = subprocess.Popen(['hmmfetch',pfamAhmm,nid], stdout=subprocess.PIPE)
        pid = subprocess.check_output(['grep','ACC'], stdin=out.stdout, universal_newlines = True)
    except:
        #print(nid+" is not in the Pfam-A.hmm")
        return('')
    pid = pid.split()[1].replace('\n','')
    return(pid)

#pfamD = dict()
#speedName = dict()

def getCisbpTFinfo(): # download the Cisbp
    url = "http://cisbp.ccbr.utoronto.ca/data/2.00/DataFiles/Bulk_downloads/EntireDataset/TF_Information_all_motifs.txt.zip"
    zipped_info = wget.download(url)
    with zipfile.ZipFile(zipped_info,"r") as zip_ref:
        zip_ref.extractall("")

def get_DBD_TF(): # main function for TF DBD. Only call once

    global dbd_pfam_file # file name of the pfam dbd stored
    global cisbp_tfinfo # file name of the cisbp tsv
    
    # local variables 
    pfamD = dict()
    speedName = dict() # so don't repeat fetch

    with open(cisbp_tfinfo,'r') as cisbp:
        for line in cisbp:
            line = line.replace('\n','')
            wdr = line.split('\t')
            names = wdr[10]
            subN = names.split(',')
            for domain in subN:
                domain = name_corr(domain)
                if domain not in speedName:
                    pid = hmmfetch(domain)
                    if pid not in pfamD and pid != '':
                        pfamD[pid]=domain
                        speedName[domain] = ''
                    elif pid =='' and pid not in speedName:
                        speedName[domain] =''
                    else:
                        pass
                
                #print(pfamD)
                #exit(0)
            #print(wdr)
            #print(wdr[10])
            
#    with open(dbd_pfam_file,'w+') as pfamDoc: # write out DBD dict key as pfam ID, value as name
#        json.dump(pfamD , pfamDoc , sort_keys = True , indent = 4, separators = (",",": "))
    # Write
    JTglobals.write(dbd_pfam_file, json.dumps(
        pfamD, sort_keys=True, indent=4, separators=(",", ": ")))

##### functions added by nicole #######################
#
#   1. getPfam(), getPfamregions(), getCisbpTFinfo()
#
#   2. get pfam names from pfamRegions
#
#   3. runHMM_fromPfam() - in a loop with all the json 
#
#   4. TF DBD info from cisbp
#
#######################################################

def parse_args():
    """
    This function parses arguments provided via the command line using argparse.
    """

    parser = argparse.ArgumentParser(description="creates files for profile inference tool.")

    # Optional args
    files_dir = os.path.dirname(os.path.realpath(__file__))
    parser.add_argument("-o", metavar="OUTDIR", default=files_dir,
        help="Output directory (default = %s)" % files_dir)

    return parser.parse_args()

def make_files(out_dir=os.path.dirname(os.path.realpath(__file__))):
    
    # Initialize
    matrix_ids = set()


    # Globals
    global codec
    global cwd
    global profiles_file_tail
    global uniprot_file_tail

    # Create output dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Download JASPAR profiles
    _download_JASPAR_profiles(out_dir)

    # Download UniProt sequences
    _download_UniProt_sequences(out_dir)

    # Download Pfam DBDs
    _download_Pfam_DBDs(out_dir)

    # Get Pfam alignments
    _get_Pfam_alignments(out_dir)
    exit(0)

    # Skip if taxon FASTA file already exists
    fasta_file = os.path.join(out_dir, "%s.fa" % taxon)
    if not os.path.exists(fasta_file):
        # Load JSON file
        with open(uniprot_json_file) as f:
            uniaccs = json.load(f)
        # For each UniProt Accession...
        for uniacc in sorted(uniaccs):
            # Write
            JTglobals.write(fasta_file,
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

    # # Skip if Cis-BP JSON file already exists
    # cisbp_json_file = os.path.join(out_dir, "cisbp.json")
    # if not os.path.exists(cisbp_json_file):
    #     # Initialize
    #     proteins = {}
    #     prot_features = {}
    #     tfs = {}
    #     tf_families = {}
    #     # Create Cis-BP dir
    #     cisbp_dir = os.path.join(out_dir, "cisbp")
    #     if not os.path.exists(cisbp_dir):
    #         os.makedirs(cisbp_dir)
    #     # Change dir
    #     os.chdir(cisbp_dir)
    #     # Skip if TFs file already exists
    #     if not os.path.exists("cisbp_1.02.tfs.sql"):
    #         # Download SQL files
    #         os.system("curl --silent -O http://cisbp.ccbr.utoronto.ca/data/1.02/DataFiles/SQLDumps/SQLArchive_cisbp_1.02.zip")
    #         # Unzip
    #         os.system("unzip -qq SQLArchive_cisbp_1.02.zip")
    #         # Remove SQL files
    #         os.remove("SQLArchive_cisbp_1.02.zip")
    #         # For each ZIP file...
    #         for zip_file in frozenset(os.listdir(os.getcwd())):
    #             # Skip non-zip files
    #             if not zip_file.endswith(".zip"): continue
    #             # Unzip
    #             os.system("unzip -qq %s" % zip_file)
    #             os.remove(zip_file)
    #     # Return to original dir
    #     os.chdir(cwd)
    #     # Get protein features
    #     with open(os.path.join(cisbp_dir, "cisbp_1.02.prot_features.sql")) as f:
    #         # For each line...
    #         for line in f:
    #             m = re.search("\('.+', '(.+)', '.+', \d+, \d+, '(.+)'\)", line)
    #             if m:
    #                 prot_features.setdefault(m.group(1), set())
    #                 prot_features[m.group(1)].add(m.group(2))
    #     # Get TFs
    #     with open(os.path.join(cisbp_dir, "cisbp_1.02.tfs.sql")) as f:
    #         # For each line...
    #         for line in f:
    #             m = re.search("\('(.+)', '(.+)', '.+', '.+', '.+', '.+', '.+'\)", line)
    #             if m:
    #                 tfs.setdefault(m.group(1), m.group(2))
    #     # Get TF families
    #     with open(os.path.join(cisbp_dir, "cisbp_1.02.tf_families.sql")) as f:
    #         # For each line...
    #         for line in f:
    #             m = re.search("\('(.+)', '.+', '.+', \d+, (.+)\)", line)
    #             if m:
    #                 tf_families.setdefault(m.group(1), m.group(2))
    #     # Get proteins
    #     with open(os.path.join(cisbp_dir, "cisbp_1.02.proteins.sql")) as f:
    #         # For each line...
    #         for line in f:
    #             m = re.search("\('(.+)', '(.+)', '.+', '.+', '([A-Z]+)\W*'\)", line)
    #             if m:
    #                 if m.group(1) not in prot_features: continue
    #                 # Digest to MD5
    #                 h = hashlib.new("md5")
    #                 h.update(m.group(3).encode("utf-8"))
    #                 md5 = h.hexdigest() + m.group(3)[:4] + m.group(3)[-4:]
    #                 proteins.setdefault(md5, [tf_families[tfs[m.group(2)]], []])
    #                 # For each domain...
    #                 for domain in prot_features[m.group(1)]:
    #                     if domain not in proteins[md5][1]:
    #                         proteins[md5][1].append(domain)
    #     # Write
    #     JTglobals.write(cisbp_json_file, json.dumps(
    #         proteins, sort_keys=True, indent=4, separators=(",", ": ")))
    #     # Remove Cis-BP dir
    #     shutil.rmtree(cisbp_dir)

    # Skip if JSON files already exist
    #domains_json_file = os.path.join(out_dir, "domains.json")
    jaspar_json_file = os.path.join(out_dir, "jaspar.json")

    jaspar = {}

    if os.path.exists(jaspar_json_file): jaspar_json_file
        # For each taxon...
    for taxon in JTglobals.taxons:
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
            # #h = hashlib.new("md5")
            # h.update(uniaccs[uniacc][1].encode("utf-8"))
            # md5 = h.hexdigest() + uniaccs[uniacc][1][:4] + uniaccs[uniacc][1][-4:]
            # # If sequence in Cis-BP...
            # if md5 in cisbp:
            #     # Add to domains
            #     domains.setdefault(uniacc, [cisbp[md5][1], cisbp[md5][0]])
                # For each profile...
        for profile in uniaccs[uniacc][0]:
            # Add to JASPAR
            jaspar.setdefault(uniacc, [])
            jaspar[uniacc].append([profile, profiles[profile]])
    # Write
    #JTglobals.write(domains_json_file, json.dumps(
    #    domains, sort_keys=True, indent=4, separators=(",", ": ")))
    JTglobals.write(jaspar_json_file, json.dumps(
        jaspar, sort_keys=True, indent=4, separators=(",", ": ")))

    # For each taxon...
    for taxon in JTglobals.taxons:
        # Skip if taxon dir already exists
        taxon_dir = os.path.join(out_dir, taxon)
        if not os.path.exists(taxon_dir):
            # Create taxon dir
            os.makedirs(taxon_dir)
            # Change dir
            os.chdir(taxon_dir)
            # Download JASPAR profiles
            url = "http://jaspar.genereg.net/download/CORE/JASPAR2018_CORE_%s_redundant_pfms_jaspar.zip" % taxon
            os.system("curl --silent -O %s" % url)
            # Unzip #
            file_name = "JASPAR2018_CORE_%s_redundant_pfms_jaspar.zip" % taxon
            os.system("unzip -qq %s" % file_name)
            # Remove downloaded file
            os.remove("%s" % file_name)
            # Return to original dir
            os.chdir(cwd)

    # Nicole's added parts #####################################################

    # 1. get all the files needed
    getPfamregions() # used for getting pfam id for uniprot ids
    getPfam() # used for getting hmm files
    getCisbpTFinfo() # get the DBD TF info 

    global uniprot_file_tail
    global pfam_file_tail
    global pfam_aligned_file_tail
    global cisbp_tfinfo

    # get the pfam id files ready
    for taxon in JTglobals.taxons:
        # Initialize
        uniprot_json_file = os.path.join(
            out_dir, taxon+uniprot_file_tail)
        pfam_json_file = os.path.join(
            out_dir, taxon+pfam_file_tail)
        pfam_alignments = os.path.join(
            out_dir, taxon+pfam_aligned_file_tail)
        # create pfam jsons with uniprot id as key
        createPfamJson(uniprot_json_file, pfam_json_file) #unip Json, pfam out

        # run hmm stuff on the pfam files
        runHMM_fromPfam(uniprot_json_file, pfam_json_file, pfam_alignments) # uniprot sequence file, pfam json file
    
    # create json for DBD tf (using cisbp TF file info)
    get_DBD_TF() # automatically creates the file

def _download_JASPAR_profiles(out_dir=os.path.dirname(os.path.realpath(__file__))):

    # For each taxon...
    for taxon in JTglobals.taxons:
        # Initialize
        url = "http://jaspar.genereg.net/api/v1/taxon/%s/" % taxon
        # Skip if taxon profiles JSON file already exists
        profiles_json_file = os.path.join(out_dir, taxon+profiles_file_tail)
        if not os.path.exists(profiles_json_file):
            # Initialize
            profiles = {}
            client = coreapi.Client()        
            response = client.get(url)
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
            JTglobals.write(profiles_json_file, json.dumps(
                profiles, sort_keys=True, indent=4, separators=(",", ": ")))

def _download_UniProt_sequences(out_dir=os.path.dirname(os.path.realpath(__file__))):

    # Initialize
    failed_profiles = {
        "MA0024.1": ["Q01094"],
        "MA0046.1": ["P20823"],
        "MA0052.1": ["Q02078"],
        "MA0058.1": ["P61244"],
        "MA0098.1": ["P14921"],
        "MA0110.1": ["P46667"],
        "MA0138.1": ["Q13127"],
        "MA0328.1": ["P0CY08"]
    }
    failed_uniprot = {
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

    # For each taxon...
    for taxon in JTglobals.taxons:
        # Initialize
        profiles_json_file = os.path.join(out_dir, taxon+profiles_file_tail)
        # Skip if taxon uniprot JSON file already exists
        uniprot_json_file = os.path.join(out_dir, taxon+uniprot_file_tail)
        if not os.path.exists(uniprot_json_file):
            # Initialize
            uniaccs = {}
            client = coreapi.Client()
            # Load JSON file
            with open(profiles_json_file) as f:
                profiles = json.load(f)
            # For each profile...
            for profile in sorted(profiles):
                # Initialize
                url = "http://jaspar.genereg.net/api/v1/matrix/%s/" % profile
                # Get profile detailed info
                response = client.get(url)
                json_obj = json.loads(codec.encode(response))
                # Fix bugged cases
                if json_obj["matrix_id"] in failed_profiles:
                    json_obj["uniprot_ids"] = failed_profiles[json_obj["matrix_id"]]
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
                # Fix bugged cases
                if uniacc in failed_uniprot:
                    uniaccs[uniacc][1] = "".join(failed_uniprot[uniacc]) 
                    continue
                # Get UniProt sequence
                u = uniprot.queryUniprot(uniacc)
                uniaccs[uniacc][1] = "".join(u["sequence   0"].split("\n"))
            # Write
            JTglobals.write(uniprot_json_file, json.dumps(
                uniaccs, sort_keys=True, indent=4, separators=(",", ": ")))

def _download_Pfam_DBDs(out_dir=os.path.dirname(os.path.realpath(__file__))):

    # Initialize
    pfam_DBDs = {}
    pfam_ids = set()
    url = "http://cisbp.ccbr.utoronto.ca/data/2.00/DataFiles/Bulk_downloads/EntireDataset/"
    cisbp_file = "TF_Information_all_motifs.txt.zip"

    # Change dir
    os.chdir(out_dir)

    # Create Pfam dir
    pfam_dir = os.path.join(out_dir, "pfam-DBDs")
    if not os.path.exists(pfam_dir):
        os.makedirs(pfam_dir)

    # Skip if Pfam DBD file already exists
    pfam_DBD_file = os.path.join(out_dir, "pfam-DBDs.json")
    if not os.path.exists(pfam_DBD_file):
        # Skip if CIS-BP file already exists
        if not os.path.exists(cisbp_file):
            os.system("curl --silent -O %s%s" % (url, cisbp_file))
        # For each line...
        for line in JTglobals.parse_tsv_file(cisbp_file):
            pfam_ids.add(line[10])
        # Change dir
        os.chdir(pfam_dir)
        # For each Pfam ID...
        for pfam_id in sorted(pfam_ids):
            try:
                # Fetch MSA from Pfam
                msa_file = pfam.fetchPfamMSA(pfam_id, alignment="seed")
                # For each line...
                for line in JTglobals.parse_file(msa_file):
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
                process = subprocess.call(cmd, shell=True)
                # HMM press
                cmd = "hmmpress -f %s &> /dev/null" % hmm_file
                process = subprocess.call(cmd, shell=True)
                # Add Pfam
                pfam_DBDs.setdefault(pfam_ac, pfam_id_std)
                # Remove MSA file
                os.remove(msa_file)
            except:
                print("\nCould not fetch MSA for id: %s\n" % pfam_id)
        # Change dir
        os.chdir(out_dir)
        # Write
        JTglobals.write(pfam_DBD_file, json.dumps(
            pfam_DBDs, sort_keys=True, indent=4, separators=(",", ": ")))
        # Remove CIS-BP file
        os.remove(cisbp_file)

#    # Skip if HMM database of all DBDs already exists
#    hmm_db = os.path.join(pfam_dir, "all_DBDs.hmm")
#    if not os.path.exists(hmm_db):
#        # For each HMM file...
#        for hmm_file in os.listdir(pfam_dir):
#            # Skip if not HMM file
#            if not hmm_file.endswith(".hmm"): continue
#            # Initialize
#            file_name = os.path.join(pfam_dir, hmm_file)
#            # For each line...
#            for line in JTglobals.parse_file(file_name):
#                JTglobals.write(hmm_db, line)
#        # HMM press
#        cmd = "hmmpress -f %s &> /dev/null" % hmm_db
#        process = subprocess.call(cmd, shell=True)

    # Change dir
    os.chdir(cwd)

def _get_Pfam_alignments(out_dir=os.path.dirname(os.path.realpath(__file__))):

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

    # Parse arguments
    args = parse_args()

    # Make files
    make_files(args.o)

    # p = pfam.searchPfam("P49639")
    # msa = pfam.fetchPfamMSA("PF00046")