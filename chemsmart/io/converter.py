import logging
import os


from chemsmart.io.gaussian.output import Gaussian16Output, GaussianLogFolder
from chemsmart.utils.logger import create_logger

logger = logging.getLogger(__name__)
os.environ["OMP_NUM_THREADS"] = "1"


class FileConverter:
    """Class for converting files in different formats.
    Args:
        directory (str): Directory in which to convert files.
        type (str): Type of file to be converted, if directory is specified.
        filename (str): Input filename to be converted.
        output_filetype (str): Type of files to convert to, defaults to .xzy.
    """

    def __init__(self, directory, type, filename, output_filetype="xyz"):
        self.directory = directory
        self.type = type
        self.filename = filename
        self.output_filetype = output_filetype

    def convert_files(self):
        create_logger()
        if self.directory is not None:
            logger.info(f"Converting files in directory: {self.directory}")
            assert (
                self.type is not None
            ), "Type of file to be converted must be specified."
            self._convert_all_files(
                self.directory, self.type, self.output_filetype
            )
        else:
            logger.info(f"Converting file: {self.filename}")
            self._convert_single_file(self.filename, self.output_filetype)

    def _convert_all_files(self, directory, type, output_filetype):
        """Convert all files of specified type in the directory."""
        g16_folder = GaussianLogFolder(folder=directory)
        all_logpaths = g16_folder.all_logfiles

        for logpath in all_logpaths:
            outputfile = Gaussian16Output(filename=logpath)
            outputfile.write_xyz(output_filetype)
