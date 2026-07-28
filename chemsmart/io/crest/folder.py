import logging
import os

from chemsmart.io.folder import BaseFolder

logger = logging.getLogger(__name__)


class CRESTFolder(BaseFolder):
    """
    Folder containing all CREST-related output files for postprocessing.

    CREST calculations generate multiple output files in a single directory.
    This class provides methods to locate and access all these files.

    Typical CREST calculation folder structure:
        crest_calculation/
        ├─ *.out                → main CREST output log
        ├─ crest_conformers.xyz → conformer ensemble (sorted by energy)
        ├─ crest_best.xyz       → lowest-energy conformer
        ├─ crest_rotamers.xyz   → rotamer ensemble
        ├─ crest.energies       → energy list for all conformers
        ├─ coord                → TURBOMOLE format input geometry
        ├─ crestopt.log         → optimization trajectory log file (XYZ format)
        ├─ constraints.inp      → distance constraints (if constrained run)
        └─ ...                  → other auxiliary files
    """

    @property
    def is_crest_calculation_directory(self):
        """
        Check if this folder is a valid CREST calculation directory.

        A valid CREST directory must contain exactly one main output file
        and crest_conformers.xyz.

        Returns:
            bool: True if this is a valid CREST calculation directory.
        """
        if self._crest_out() is None:
            return False
        filepath = os.path.join(self.folder, "crest_conformers.xyz")
        if os.path.exists(filepath):
            return True
        logger.debug(
            f"Directory {self.folder} has CREST output file but lacks "
            f"crest_conformers.xyz."
        )
        return False

    def _crest_out(self):
        """
        Return the path to the main CREST output file.

        Returns:
            str | None: Path to the CREST output file, or None if not found.
        """
        crest_out = self.get_all_output_files_in_current_folder_by_program(
            program="crest"
        )
        if not crest_out:
            return None
        if len(crest_out) > 1:
            filenames = ", ".join(
                sorted(os.path.basename(path) for path in crest_out)
            )
            raise ValueError(
                f"Multiple CREST main output files found in '{self.folder}': "
                f"{filenames}. Each CREST calculation directory must contain "
                f"exactly one main output file."
            )
        return crest_out[0]

    def _conformers_xyz(self):
        """Return the path to the conformer ensemble file."""
        conformers_xyz = os.path.join(self.folder, "crest_conformers.xyz")
        return conformers_xyz if os.path.exists(conformers_xyz) else None

    def _best_xyz(self):
        """Return the path to the best (lowest energy) conformer file."""
        best_xyz = os.path.join(self.folder, "crest_best.xyz")
        return best_xyz if os.path.exists(best_xyz) else None

    def _rotamers_xyz(self):
        """Return the path to the rotamer ensemble file."""
        rotamers_xyz = os.path.join(self.folder, "crest_rotamers.xyz")
        return rotamers_xyz if os.path.exists(rotamers_xyz) else None

    def _energies(self):
        """Return the path to the energies file."""
        energies = os.path.join(self.folder, "crest.energies")
        return energies if os.path.exists(energies) else None

    def _constraints_inp(self):
        """Return the path to the constraints input file."""
        constraints_inp = os.path.join(self.folder, "constraints.inp")
        return constraints_inp if os.path.exists(constraints_inp) else None

    def _crestopt_log(self):
        """Return the path to optimization trajectory log file."""
        crestopt_log = os.path.join(self.folder, "crestopt.log")
        return crestopt_log if os.path.exists(crestopt_log) else None

    def _coord(self):
        """Return the path to TURBOMOLE coord file."""
        coord = os.path.join(self.folder, "coord")
        return coord if os.path.exists(coord) else None
