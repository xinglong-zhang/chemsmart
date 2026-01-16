import logging
import os

from chemsmart.io.folder import BaseFolder
from chemsmart.utils.io import find_output_files_in_directory

logger = logging.getLogger(__name__)


class XTBFolder(BaseFolder):
    """
    Folder containing all XTB-related output files for postprocessing.

    XTB calculations generate multiple output files in a single directory.
    This class provides methods to locate and access all these files.

    Typical XTB calculation folder structure:
        xtb_calculation/
        ├─ *.out            → main XTB output file (energy, frequencies, thermodynamics)
        ├─ xtbopt.log       → optimization trajectory log file (XYZ format)
        ├─ charges          → atomic partial charges
        ├─ energy           → total electronic energy
        ├─ *.engrad         → energy, gradient, and coordinates
        ├─ g98.out          → Gaussian-style vibrational analysis output
        ├─ gradient         → nuclear gradients (forces)
        ├─ hessian          → cartesian Hessian matrix
        ├─ vibspectrum      → vibrational frequencies and IR intensities
        ├─ wbo              → Wiberg bond orders
        ├─ xtbopt.*         → optimized geometry (XYZ, SDF, Turbomole coord, PDB, VASP POSCAR,
        │                     DFTB+ gen or Gaussian External format)
        ├─ xtbtopo.mol      → molecular topology and bonding information
        └─ ...              → other auxiliary files
    """

    def _xtb_out(self):
        """Return the path to the main XTB output file."""
        xtbout = find_output_files_in_directory(
            self.folder, program="xtb", recursive=False
        )
        return xtbout[0] if xtbout else None

    def _xtbopt_log(self):
        """Return the path to optimization trajectory log file."""
        xtbopt_log = os.path.join(self.folder, "xtbopt.log")
        return xtbopt_log if os.path.exists(xtbopt_log) else None

    def _charges(self):
        """Return the path to charge file."""
        charge_file = os.path.join(self.folder, "charges")
        return charge_file if os.path.exists(charge_file) else None

    def _energy(self):
        """Return the path to energy file."""
        energy = os.path.join(self.folder, "energy")
        return energy if os.path.exists(energy) else None

    def _engrad(self):
        """Return the path to engrad file."""
        engrad = self.get_all_files_in_current_folder_by_suffix(".engrad")
        return engrad[0] if engrad else None

    def _g98_out(self):
        """Return the path to the GAUSSIAN-format vibrational frequencies file."""
        g98_out = os.path.join(self.folder, "g98.out")
        return g98_out if os.path.exists(g98_out) else None

    def _gradient(self):
        """Return the path to gradient file."""
        gradient = os.path.join(self.folder, "gradient")
        return gradient if os.path.exists(gradient) else None

    def _hessian(self):
        """Return the path to hessian file."""
        hessian = os.path.join(self.folder, "hessian")
        return hessian if os.path.exists(hessian) else None

    def _vibspectrum(self):
        """Return the path to vibrational spectrum file."""
        vibspectrum = os.path.join(self.folder, "vibspectrum")
        return vibspectrum if os.path.exists(vibspectrum) else None

    def _wbo(self):
        """Return the path to wiberg bond order file."""
        wbo = os.path.join(self.folder, "wbo")
        return wbo if os.path.exists(wbo) else None

    def _xtbopt_geometry(self):
        """
        Return the path to optimized geometry file (xtbopt.*).

        XTB outputs the optimized geometry in the same format as the input.
        This method returns a parseable xtbopt.* file if available, otherwise
        falls back to an unsupported format with a warning.

        Search priority (parseable formats only):
            1. xtbopt.xyz     → XYZ format (supported)
            2. xtbopt.sdf     → SDF format (supported)

        Unsupported formats (will show warning):
            - xtbopt.coord   → Turbomole coord format
            - xtbopt.pdb     → PDB format
            - xtbopt.poscar  → VASP POSCAR format
            - xtbopt.gen     → DFTB+ gen format
            - xtbopt.EIn     → Gaussian External format

        Returns:
            str | None: Path to optimized geometry file, or None if
                        no xtbopt.* file exists
        """
        # Formats chemsmart can currently parse
        parseable_extensions = [".xyz", ".sdf"]
        # Formats not yet supported
        unsupported_extensions = [".coord", ".pdb", ".poscar", ".gen", ".EIn"]
        # Try parseable formats first
        for ext in parseable_extensions:
            filepath = os.path.join(self.folder, f"xtbopt{ext}")
            if os.path.exists(filepath):
                logger.debug(f"Found optimized geometry file: {filepath}")
                return filepath
        # Check if unsupported format exists
        for ext in unsupported_extensions:
            filepath = os.path.join(self.folder, f"xtbopt{ext}")
            if os.path.exists(filepath):
                logger.warning(
                    f"Found optimized geometry file {filepath}, but format {ext} "
                    "is not yet supported by chemsmart."
                )
                return filepath
        # No xtbopt.* file found
        return None

    def _xtbtopo_mol(self):
        """Return the path to molecular topology file."""
        xtbtopo_mol = os.path.join(self.folder, "xtbtopo.mol")
        return xtbtopo_mol if os.path.exists(xtbtopo_mol) else None
