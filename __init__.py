"""
JASPAR profile inference tool
"""

__author__ = "Oriol Fornes, Xi Zhang"
__email__ = "oriol@cmmt.ubc.ca, nzhang@cmmt.ubc.ca"
__organization__ = "The JASPAR Consortium"
__version__ = "2.0.1"

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

    taxons = ["fungi", "insects", "nematodes", "plants", "vertebrates"]

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
    "Homeodomain": "Homeodomain",
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

class CisBP:

    def __init__(self, model_file):
        self.file = model_file
        self.__model = {}
        self.family = None
        self.__parse_model()

    def __parse_model(self):
        handle = Jglobals._get_file_handle(self.file)
        model = json.load(handle)
        self.family = CisBP2Pfam[model["Family_Name"]]
        if model["Model.Name"] == "PctID_L":
            self.__model.setdefault("pid", {})
            self.__model["pid"].setdefault("dis", model["Threshold.Dis"])
            self.__model["pid"].setdefault("hsim", model["Threshold.HSim"])
        else:
            self.__model.setdefault("pid", {})
            self.__model["pid"].setdefault("dis", model["Baseline"]["Threshold.Dis"])
            self.__model["pid"].setdefault("hsim", model["Baseline"]["Threshold.HSim"])
            self.__model.setdefault("sr", {})
            self.__model["sr"].setdefault("dis", model["Threshold.Dis"])
            self.__model["sr"].setdefault("hsim", model["Threshold.HSim"])
            self.__model["sr"].setdefault("mean", model["SR.FeatureScales.mean"])
            self.__model["sr"].setdefault("sd", model["SR.FeatureScales.sd"])
            self.__model["sr"].setdefault("intercept", model["SR.Intercept"])
            self.__model["sr"].setdefault("weights", model["SR.Weights"])
            if model["SR.Features"] == "AvgB62":
                self.__model["sr"].setdefault("similarity", "blosum62")
            else:
                self.__model["sr"].setdefault("similarity", "identity")

    def get_model_threshold(self, what_model, what_threshold):
        if what_model in self.__model:
            if what_threshold in self.__model[what_model]:
                threshold = self.__model[what_model][what_threshold]
                if np.isnan(threshold):
                    return(None)
                return(threshold)
        return(None)

    def get_model(self, what_model):
        if what_model in self.__model:
            return(self.__model[what_model])
        return(None)
