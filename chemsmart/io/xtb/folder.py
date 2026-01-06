import logging
import os

from chemsmart.io.folder import BaseFolder
from chemsmart.utils.io import find_output_files_in_directory

logger = logging.getLogger(__name__)


class XTBFolder(BaseFolder):
    """
    Folder containing all XTB-related output files for postprocessing.
    """

    def _xtb_out(self):
        """Return the path to the main XTB output file."""
        xtbout = find_output_files_in_directory(self.folder, program="xtb")
        return xtbout[0] if xtbout else None

    def _xtb_err(self):
        """Return the path to XTB error output file."""
        errfiles = self.get_all_files_in_current_folder_by_suffix(".err")
        return errfiles[0] if errfiles else None

    def _charge(self):
        """Return the path to charge file."""
        charge_file = os.path.join(self.folder, "charges")
        return charge_file if os.path.exists(charge_file) else None

    def _wbo(self):
        """Return the path to wiberg bond order file."""
        wbo = os.path.join(self.folder, "wbo")
        return wbo if os.path.exists(wbo) else None

    def _xtbopt_log(self):
        """Return the path to optimization trajectory log file."""
        xtbopt_log = os.path.join(self.folder, "xtbopt.log")
        return xtbopt_log if os.path.exists(xtbopt_log) else None

    def _energy(self):
        """Return the path to energy file."""
        energy = os.path.join(self.folder, "energy")
        return energy if os.path.exists(energy) else None

    def _gradient(self):
        """Return the path to gradient file."""
        gradient = os.path.join(self.folder, "gradient")
        return gradient if os.path.exists(gradient) else None

    def _engrad(self):
        """Return the path to engrad file."""
        engrad = self.get_all_files_in_current_folder_by_suffix(".engrad")
        return engrad[0] if engrad else None

    def _hessian(self):
        """Return the path to hessian file."""
        hessian = os.path.join(self.folder, "hessian")
        return hessian if os.path.exists(hessian) else None

    def _g98_out(self):
        """Return the path to the GAUSSIAN-format vibrational frequencies file."""
        g98_out = os.path.join(self.folder, "g98.out")
        return g98_out if os.path.exists(g98_out) else None

    def _vibspectrum(self):
        """Return the path to vibrational spectrum file."""
        vibspectrum = os.path.join(self.folder, "vibspectrum")
        return vibspectrum if os.path.exists(vibspectrum) else None
