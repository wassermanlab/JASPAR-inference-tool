"""
JASPAR profile inference tool
"""

__author__ = "Oriol Fornes, Xi Zhang"
__email__ = "oriol.fornes@gmail.com, nzhang@cmmt.ubc.ca"
__organization__ = "The JASPAR Consortium"
__version__ = "2024.3.1"

from Bio import SeqIO
import gzip
import json
import numpy as np
import pandas
import sys
from zipfile import ZipFile

class Globals(object):
    """
    This class contains functions designed to work through the entire module.
    """

    #-------------#
    # Definitions #
    #-------------#

    version = 2024
    taxons = [
        # "cnidaria",
        # "diatoms",
        # "dictyostelium",
        "fungi",
        "insects",
        "nematodes",
        # "oomycota",
        "plants",
        # "trematodes",
        "urochordates",
        "vertebrates"
    ]

    #--------------#
    # Input/Output #
    #--------------#

    def _get_file_handle(self, file_name, mode="rt"):

        # Initialize
        raiseValueError = False

        # Open gzip file as handle
        if file_name.endswith(".gz"):
            try:
                handle = gzip.open(file_name, mode)
            except:
                raiseValueError = True
        # Open zip file as handle
        elif file_name.endswith(".zip"):
            try:
                zf = ZipFile(file_name, mode)
                for f in zf.infolist():

                    handle = zf.open(f, mode)
                    break
            except:
                raiseValueError = True
        # Open file as handle
        else:
            try:
                handle = open(file_name, mode)
            except:
                raiseValueError = True

        if raiseValueError:
            raise ValueError("Could not open file handle: %s" % file_name)

        return(handle)

    def parse_file(self, file_name):
        """
        Parses a file and yields lines one by one as a {str}.
        @input:
        file_name {str}
        @yield: {str}
        """

        # Get file handle
        handle = self._get_file_handle(file_name)

        # For each line...
        for line in handle:
            yield(line.strip("\n"))

        handle.close()

    def parse_csv_file(self, file_name, delimiter=","):
        """
        Parses a CSV file and yields lines one by one as a {list}.
        @input:
        file_name {str}
        delimiter {str} e.g. "\t"; default = ","
        @yield: {list}
        """

        # Get file handle
        handle = self._get_file_handle(file_name)

        # Read in chunks
        for chunk in pandas.read_csv(handle, header=None, encoding="utf8",
                                     sep=delimiter, chunksize=1024):
            for index, row in chunk.iterrows():
                yield(row.tolist())

        handle.close()

    def parse_tsv_file(self, file_name):
        """
        Parses a TSV file and yields lines one by one as a {list}.
        @input:
        file_name {str}
        @yield: {list}
        """

        # For each line...
        for line in self.parse_csv_file(file_name, delimiter="\t"):
            yield(line)

    def parse_fasta_file(self, file_name):
        """
        Parses a FASTA file and yields sequences one by one as a {SeqRecord}.
        @input:
        file_name {str}
        @yield: {SeqRecord}
        """

        # Get file handle
        handle = self._get_file_handle(file_name)

        # For each SeqRecord...
        for seq_record in SeqIO.parse(handle, "fasta"):
            yield(seq_record)

        handle.close()

    def write(self, file_name=None, content=None):
        """
        Writes content to a file or, if no file is provided, to STDOUT.
        Content will be appended at the end of the file.
        """

        if file_name:

            # Get file handle
            handle = self._get_file_handle(file_name, mode="at")

            # Write
            handle.write("%s\n" % content)
            handle.close()

        else:
            sys.stdout.write("%s\n" % content)

Jglobals = Globals()

Pfam2CisBP = {
    "AP2": "AP2",
    "KilA-N": "APSES",
    "ARID": "ARID/BRIGHT",
    "AT_hook": "AT hook",
    "HLH": "bHLH",
    "bZIP_1": "bZIP",
    "zf-C2H2": "C2H2 ZF",
    "DM": "DM",
    "zf-Dof": "Dof",
    "E2F_TDP": "E2F",
    "Ets": "Ets",
    "Forkhead": "Forkhead",
    "GATA": "GATA",
    "GCM": "GCM",
    "Homeobox": "Homeodomain",
    "Pou": "Homeodomain,POU",
    "HSF_DNA-bind": "HSF",
    "MADF_DNA_bdg": "MADF",
    "Myb_DNA-binding": "Myb/SANT",
    "NAM": "NAC/NAM",
    "zf-C4": "Nuclear receptor",
    "HTH_psq": "Pipsqueak",
    "RFX_DNA_binding": "RFX",
    "SAND": "SAND",
    "SBP": "SBP",
    "HMG_box": "Sox",
    "T-box": "T-box",
    "TCP": "TCP",
    "TCR": "TCR/CxC",
    "WRKY": "WRKY",
    "Zn_clus": "Zinc cluster",
    None: "NO_THRESHOLD"
}

CisBP2Pfam = {v: k for k, v in Pfam2CisBP.items()}

###########################################################
# Github: https://github.com/smlmbrt/SimilarityRegression #
# Script: ./similarityregression/PredictSimilarity.py     #
###########################################################

def ReadSRModel(filename):
    with open(filename) as SRModel:
        srmodel = json.load(SRModel)
    #Convert to NP arrarys
    if 'SR.Weights' in srmodel:
        srmodel['SR.Weights'] = np.asarray(srmodel['SR.Weights'])
        srmodel['SR.FeatureScales.mean'] = np.asarray(
            srmodel['SR.FeatureScales.mean']
        )
        #Convert 0's to NAs
        sd = np.asarray(srmodel['SR.FeatureScales.sd'])
        sd[sd == 0] = np.nan
        srmodel['SR.FeatureScales.sd'] = sd
    #Check for Amb/Dis threshold
    if np.isnan(srmodel['Threshold.Dis']):
        srmodel['Threshold.Dis'] = None
    return(srmodel)
    
def ScoreAlignmentResult(resultDict, scoreDict, applyidenticalRule = True):
    #Check if 100% identical (gets rid of proteins w/ truncations)
    if (applyidenticalRule == True) and (resultDict['PctID_L'] == 1):
        return(resultDict['PctID_L'], 'HSim')
    #Score The Sequence
    if scoreDict['Model.Class'] == 'SequenceIdentity':
        Score = resultDict[scoreDict['Model.Name']]
        threshold_hsim = scoreDict['Threshold.HSim']
        threshold_dis = scoreDict['Threshold.Dis']
        if Score >= threshold_hsim:
            Classification = 'HSim'
            return(Score, Classification)
        # else:
        #     Classification = 'Amb'
        # ##Check if Amb/Dis
        # if threshold_dis != None:
        #     if Score < threshold_dis:
        #         Classification = 'Dis'
    else:
        Score = resultDict[scoreDict['Baseline']['Name']]
        threshold_hsim = scoreDict['Baseline']['Threshold.HSim']
        # threshold_dis = scoreDict['Baseline']['Threshold.Dis']    
        if Score >= threshold_hsim:
            Classification = 'HSim'
            return(Score, Classification)
        # else:
        #     Classification = 'Amb'
        # ##Check if Amb/Dis
        # if threshold_dis != None:
        #     if Score < threshold_dis:
        #         Classification = 'Dis'
        SRweights = scoreDict['SR.Weights']
        #Get postional scores
        key = 'ByPos.'+scoreDict['SR.Features'].replace('_','.')
        ByPos = np.array(resultDict[key])
        # i.e. fix error when length of Pfam domain != from Cis-BP
        if len(ByPos) == len(scoreDict['SR.FeatureScales.mean']):
            #Normalize to features (f)
            f = (ByPos - scoreDict['SR.FeatureScales.mean']) / \
                scoreDict['SR.FeatureScales.sd']
            f[np.isnan(f)] = 0 #Cleanup NAs
            Score = scoreDict['SR.Intercept'] + np.dot(SRweights, f)
            if scoreDict['SR.LogisticTransform'] == True:
                logistic = lambda x: 1 / (1 + np.exp(-x))
                Score = logistic(Score)
            threshold_hsim = scoreDict['Threshold.HSim']
            threshold_dis = scoreDict['Threshold.Dis']
            if Score >= threshold_hsim:
                Classification = 'HSim'
                return(Score, Classification)
            # else:
            #     Classification = 'Amb'
            # ##Check if Amb/Dis
            # if threshold_dis != None:
            #     if Score < threshold_dis:
            #         Classification = 'Dis'

    return(np.nan, np.nan)