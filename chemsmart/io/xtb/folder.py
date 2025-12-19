import logging
import os

from chemsmart.io.folder import BaseFolder
from chemsmart.utils.io import get_outfile_format

logger = logging.getLogger(__name__)


class XTBFolder(BaseFolder):
    """
    Folder containing all XTB-related output files for postprocessing.
    """

    def __init__(self, folder):
        """
        Initialize XTB output folder handler.

        Args:
            folder (str): Parent folder for all output files
        """
        self.folder = folder

    def _xtb_out(self):
        """Return the path to xtb.out (main output)."""
        outfiles = self.get_all_files_in_current_folder_by_suffix(".out")
        xtbout = [
            outfiles
            for outfile in outfiles
            if get_outfile_format(outfile) == "xtb"
        ]
        return xtbout[0] if xtbout else None

    def _xtb_err(self):
        """Return path to xtb.err (error output)."""
        return os.path.join(self.folder, "xtb.err")

    def _charge(self):
        """Return path to charge file."""
        return os.path.join(self.folder, "charge")

    def _xtbopt_log(self):
        """Return path to xtbopt.log (optimization trajectory)."""
        return os.path.join(self.folder, "xtbopt.log")

    def _gradient(self):
        """Return path to gradient file."""
        return os.path.join(self.folder, "gradient")
