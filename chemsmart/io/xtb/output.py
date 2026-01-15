import logging
from functools import cached_property

from chemsmart.io.file import SDFFile
from chemsmart.io.xtb.file import (
    XTBChargesFile,
    XTBEnergyFile,
    XTBEngradFile,
    XTBG98File,
    XTBGradientFile,
    XTBHessianFile,
    XTBMainOut,
    XTBVibSpectrumFile,
    XTBWibergBondOrderFile,
)
from chemsmart.io.xtb.folder import XTBFolder
from chemsmart.io.xyz.xyzfile import XYZFile
from chemsmart.utils.io import create_molecule_list
from chemsmart.utils.utils import string2index_1based

logger = logging.getLogger(__name__)


class XTBOutput:
    """
    Comprehensive parser and coordinator for xTB calculation outputs.

    This class discovers and parses multiple xTB output files, then integrates
    data from different sources to construct complete Molecule objects with:
    - Coordinates and energies from xtbopt.log
    - Charge and multiplicity from main output file
    - Forces from .engrad file

    Args:
        folder: Path to folder containing xTB output files, or XTBFolder instance
    """

    def __init__(self, folder):
        self.folder = (
            folder if isinstance(folder, XTBFolder) else XTBFolder(folder)
        )

    def __getattr__(self, name):
        """
        Delegate attribute access to internal file parsers if not found in XTBOutput.

        This allows XTBOutput to transparently expose all properties and methods
        from XTBMainOut, XTBEngradFile, XTBChargesFile, XTBWibergBondOrderFile,
        XTBHessianFile, XTBGradientFile, and any future parsers.

        The search order is:
        1. main_out (XTBMainOut) - main output file
        2. charges_file (XTBChargesFile) - atomic partial charges
        3. wbo_file (XTBWibergBondOrderFile) - Wiberg bond orders
        4. hessian_file (XTBHessianFile) - Hessian matrix
        5. engrad_file (XTBEngradFile) - energy and gradient file
        6. gradient_file (XTBGradientFile) - gradient file
        7. energy_file (XTBEnergyFile) - energy file
        8. g98_file (XTBG98File) - Gaussian 98 format vibrational analysis
        9. vibspectrum_file (XTBVibSpectrumFile) - vibrational spectrum

        Args:
            name: Attribute name being accessed

        Returns:
            The attribute value from the first parser that has it

        Raises:
            AttributeError: If the attribute doesn't exist in XTBOutput or any parser
        """
        # List of internal parsers to search, in order of priority
        file_parsers = [
            "main_out",
            "charges_file",
            "wbo_file",
            "hessian_file",
            "engrad_file",
            "gradient_file",
            "energy_file",
            "g98_file",
            "vibspectrum_file",
        ]
        # Try to get from each parser in order
        for file_parser in file_parsers:
            try:
                parser = object.__getattribute__(self, file_parser)
                if parser is not None:
                    try:
                        return getattr(parser, name)
                    except AttributeError:
                        # This parser doesn't have the attribute, try next one
                        continue
            except Exception as e:
                logger.debug(f"Failed to access parser {file_parser}: {e}")
                continue
        # If not found in any parser, raise AttributeError
        raise AttributeError(
            f"'{type(self).__name__}' object has no attribute '{name}'"
        )

    @cached_property
    def main_out(self):
        """Main xTB main output file parser"""
        path = self.folder._xtb_out()
        return XTBMainOut(path) if path else None

    @cached_property
    def charges_file(self):
        """Atomic partial charges file parser"""
        path = self.folder._charges()
        return XTBChargesFile(path) if path else None

    @cached_property
    def energy_file(self):
        """Energy file parser"""
        path = self.folder._energy()
        return XTBEnergyFile(path) if path else None

    @cached_property
    def engrad_file(self):
        """Energy and gradient file parser"""
        path = self.folder._engrad()
        return XTBEngradFile(path) if path else None

    @cached_property
    def g98_file(self):
        """GAUSSIAN-format vibrational frequencies file parser"""
        path = self.folder._g98_out()
        return XTBG98File(path) if path else None

    @cached_property
    def gradient_file(self):
        """Gradient file parser"""
        path = self.folder._gradient()
        return XTBGradientFile(path) if path else None

    @cached_property
    def hessian_file(self):
        """Hessian matrix file parser"""
        path = self.folder._hessian()
        return XTBHessianFile(path) if path else None

    @cached_property
    def vibspectrum_file(self):
        """Vibrational spectrum file parser"""
        path = self.folder._vibspectrum()
        return XTBVibSpectrumFile(path) if path else None

    @cached_property
    def wbo_file(self):
        """Wiberg bond order file parser"""
        path = self.folder._wbo()
        return XTBWibergBondOrderFile(path) if path else None

    @cached_property
    def xtbopt_geometry(self):
        """Read geometry from xtbopt.* file and return as Molecule object."""
        path = self.folder._xtbopt_geometry()
        if not path:
            return None
        try:
            if path.endswith(".xyz"):
                molecule = XYZFile(path).get_molecules(
                    index="-1", return_list=False
                )
            elif path.endswith(".sdf"):
                molecule = SDFFile(path).get_molecule()
            else:
                return None
            if molecule is None:
                return None

            # Enrich molecule with properties from main output
            if self.charge is not None:
                molecule.charge = self.charge
            if self.multiplicity is not None:
                molecule.multiplicity = self.multiplicity
            if self.final_energy is not None:
                molecule.energy = self.final_energy
            if self.final_forces is not None:
                molecule.forces = self.final_forces
            logger.info(f"Successfully read optimized geometry from {path}")
            return molecule

        except Exception as e:
            logger.error(f"Error reading optimized geometry from {path}: {e}")
            return None

    @property
    def charge(self):
        """Get molecular charge from main output."""
        if self.main_out:
            return self.main_out.net_charge
        if self.charges_file:
            return self.charges_file.total_charge
        return None

    @property
    def multiplicity(self):
        """Get spin multiplicity from main output."""
        if self.main_out and self.main_out.unpaired_electrons is not None:
            return self.main_out.unpaired_electrons + 1
        return None

    @property
    def final_energy(self):
        """Get final converged energy."""
        if self.main_out:
            return self.main_out.total_energy
        return None

    @property
    def final_forces(self):
        """Get final converged forces."""
        if self.engrad_file:
            return self.engrad_file.forces
        return None

    @cached_property
    def symbols(self):
        """Get atomic symbols."""
        if self.xtbopt_geometry:
            return list(self.xtbopt_geometry.symbols)
        if self.all_structures:
            return list(self.all_structures[0].symbols)
        return None

    @property
    def partial_charges(self):
        """Get atomic partial charges."""
        if (
            not self.charges_file
            or self.charges_file.partial_charges is None
            or self.symbols is None
        ):
            return None

        charges = self.charges_file.partial_charges
        symbols = self.symbols

        assert len(symbols) == len(charges), (
            f"Mismatch between number of symbols ({len(symbols)}) "
            f"and charges ({len(charges)})"
        )
        # Build dictionary with atom labels
        element_counts = {}
        labeled_charges = {}
        for symbol, charge in zip(symbols, charges):
            element_counts[symbol] = element_counts.get(symbol, 0) + 1
            label = f"{symbol}{element_counts[symbol]}"
            labeled_charges[label] = charge
        return labeled_charges

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
        # Get xtbopt.log path
        opt_log_path = self.folder._xtbopt_log()
        if not opt_log_path:
            logger.warning("No xtbopt.log file found")
            return []

        if not self.main_out:
            logger.warning("No main output file found")
            return []

        # 1. Get molecules from xtbopt.log (includes energies from comment lines)
        opt_log = XYZFile(opt_log_path)
        base_molecules = opt_log.get_molecules(index=":", return_list=True)
        n = len(base_molecules)

        if n == 0:
            logger.warning("No structures found in xtbopt.log")
            return []

        logger.debug(f"Found {n} structures in optimization trajectory")

        # Extract data from base molecules
        orientations = [mol.positions for mol in base_molecules]
        symbols = list(base_molecules[0].symbols)
        energies = [
            mol.energy for mol in base_molecules
        ]  # Already extracted from comment lines!

        logger.debug(
            f"Extracted {len(energies)} energies from xtbopt.log comment lines"
        )

        # 2. PBC conditions (xTB doesn't use PBC for optimizations)
        orientations_pbc = [None] * n

        # 3. Get forces from .engrad file (usually only final structure)
        has_final_forces = False
        final_forces = None
        if self.forces is not None:
            logger.debug(
                "Found gradient in .engrad file, assigning to last structure"
            )
            has_final_forces = True
            final_forces = self.forces
        else:
            logger.debug("No .engrad file or gradient data found")

        # 4. Build complete Molecule objects using create_molecule_list
        # This will override the basic properties and add charge/multiplicity/forces
        logger.info(f"Creating {n} Molecule objects with complete properties")
        molecules = create_molecule_list(
            orientations=orientations,
            orientations_pbc=orientations_pbc,
            energies=energies,  # From xtbopt.log comment lines
            forces=None,
            symbols=symbols,
            charge=self.charge,  # From main output
            multiplicity=self.multiplicity,  # From main output
            frozen_atoms=None,
            pbc_conditions=None,
        )
        # Assign forces to the final structure if available
        if has_final_forces and molecules:
            molecules[-1].forces = final_forces
        return molecules

    @cached_property
    def optimized_structure(self):
        """
        Get the final optimized structure (only if geometry optimization converged).

        This method attempts to read the optimized structure from multiple sources:
        1. First, try to read from optimized geometry file (xtbopt.xyz, xtbopt.sdf, etc.)

        2. If xtbopt.* file is not available, fall back to the last structure
           from the optimization trajectory (xtbopt.log)

        Returns:
            Molecule | None: Final optimized structure with all properties,
                            or None if optimization did not converge or no structure available
        """
        # Check if optimization converged and terminated normally
        if not (
            self.normal_termination
            and self.main_out
            and self.main_out.geometry_optimization_converged
        ):
            logger.debug("Optimization did not converge or terminate normally")
            return None

        # Strategy 1: Try to read from optimized geometry file (xtbopt.*)
        optimized_geometry = self.xtbopt_geometry
        if optimized_geometry is not None:
            logger.debug("Using optimized structure from xtbopt.* file")
            return optimized_geometry

        # Strategy 2: Fall back to last structure from trajectory
        if self.all_structures:
            logger.debug("Using last structure from optimization trajectory")
            return self.all_structures[-1]

        logger.warning("No optimized structure available")
        return None

    @cached_property
    def last_structure(self):
        if self.all_structures:
            return self.all_structures[-1]
        return None

    def get_molecule(self, index="-1"):
        """
        Get a specific molecule structure by index from the xTB folder.
        """
        index = string2index_1based(index)
        return self.all_structures[index]
