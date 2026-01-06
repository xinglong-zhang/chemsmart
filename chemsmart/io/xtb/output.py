import logging
from functools import cached_property

from chemsmart.io.xtb.file import XTBEngradFile, XTBMainOut, XTBOptLog
from chemsmart.io.xtb.folder import XTBFolder
from chemsmart.utils.io import create_molecule_list

logger = logging.getLogger(__name__)


class XTBOutput:
    """
    Comprehensive parser and coordinator for XTB calculation outputs.

    This class discovers and parses multiple XTB output files, then integrates
    data from different sources to construct complete Molecule objects with:
    - Coordinates and energies from xtbopt.log
    - Charge and multiplicity from main output file
    - Forces from .engrad file

    Args:
        folder: Path to folder containing XTB output files, or XTBFolder instance
    """

    def __init__(self, folder):
        self.folder = (
            folder if isinstance(folder, XTBFolder) else XTBFolder(folder)
        )

    @cached_property
    def _main_out(self):
        path = self.folder._xtb_out()
        return XTBMainOut(path) if path else None

    @cached_property
    def _engrad_file(self):
        path = self.folder._engrad()
        return XTBEngradFile(path) if path else None

    @cached_property
    def _optlog(self):
        path = self.folder._xtbopt_log()
        return XTBOptLog(path) if path else None

    @property
    def normal_termination(self):
        """Check if calculation terminated normally."""
        if self._main_out:
            return self._main_out.normal_termination
        return False

    @property
    def charge(self):
        """Get molecular charge from main output."""
        if self._main_out:
            return self._main_out.net_charge
        return None

    @property
    def multiplicity(self):
        """Get spin multiplicity from main output."""
        if self._main_out:
            return self._main_out.unpaired_electrons + 1
        return None

    @property
    def energies(self):
        """
        Get energies from all optimization steps.
        """
        if self._main_out:
            return self._main_out.energies
        return None

    @property
    def final_energy(self):
        """Get final converged energy."""
        if self._main_out:
            return self._main_out.total_energy
        return None

    @cached_property
    def all_structures(self):
        """
        Build Molecule objects from multiple data sources.

        This is where all data integration happens. Combines:
        - Coordinates and energies from xtbopt.log
        - Charge and multiplicity from main output
        - Forces from gradient file

        Returns:
            list[Molecule]: Molecule objects with all properties
        """
        if not self._optlog:
            logger.warning("No xtbopt.log file found")
            return []

        if not self._main_out:
            logger.warning("No main output file found")
            return []

        # 1. Get molecules from xtbopt.log (includes energies from comment lines)
        base_molecules = self._optlog.molecules
        n = len(base_molecules)

        if n == 0:
            logger.warning("No structures found in xtbopt.log")
            return []

        logger.debug(f"Found {n} structures in optimization trajectory")

        # Extract data from base molecules
        orientations = [mol.positions for mol in base_molecules]
        symbols = base_molecules[0].symbols
        energies = [
            mol.energy for mol in base_molecules
        ]  # Already extracted from comment lines!

        logger.debug(
            f"Extracted {len(energies)} energies from xtbopt.log comment lines"
        )

        # 2. PBC conditions (XTB doesn't use PBC for optimizations)
        orientations_pbc = [None] * n

        # 3. Get forces from .engrad file (usually only final structure)
        forces = [None] * n
        if self._engrad_file and self._engrad_file.force is not None:
            logger.debug(
                "Found gradient in .engrad file, assigning to last structure"
            )
            # The .engrad file typically contains only the final gradient
            # Assign it to the last structure
            forces[-1] = self._engrad_file.force
        else:
            logger.debug("No .engrad file or gradient data found")

        # 4. Build complete Molecule objects using create_molecule_list
        # This will override the basic properties and add charge/multiplicity/forces
        logger.info(f"Creating {n} Molecule objects with complete properties")
        molecules = create_molecule_list(
            orientations=orientations,
            orientations_pbc=orientations_pbc,
            energies=energies,  # From xtbopt.log comment lines!
            forces=forces,
            symbols=symbols,
            charge=self.charge,  # From main output
            multiplicity=self.multiplicity,  # From main output
            frozen_atoms=None,
            pbc_conditions=None,
        )

        return molecules

    @cached_property
    def optimized_structure(self):
        """
        Get the final optimized structure (only if calculation converged).

        Returns the last structure if calculation terminated normally,
        otherwise None.

        Returns:
            Molecule or None: Final optimized structure with all properties
        """
        if self.normal_termination:
            return self.all_structures[-1]
        return None

    @cached_property
    def last_structure(self):
        if self.all_structures:
            return self.all_structures[-1]
        return None
