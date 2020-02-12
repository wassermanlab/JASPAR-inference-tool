"""
JASPAR profile inference tool
"""

__author__ = "Oriol Fornes, Xi Zhang"
__email__ = "oriol@cmmt.ubc.ca, nzhang@cmmt.ubc.ca"
__organization__ = "The JASPAR Consortium"
__version__ = "2.0.1"

from Bio import SeqIO
import gzip
import pandas
import sys
from zipfile import ZipFile

__all__ = ["get_files.py", "infer_profile.py"]

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