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
    Integrates and exposes all relevant data from an xTB calculation directory.

    This class discovers and parses all relevant xTB output files in a given folder,
    then integrates the data to construct complete Molecule objects with coordinates,
    energies, charges, multiplicity, and forces. It also provides convenient access
    to all underlying file parsers and their properties.

    Args:
        folder (str or XTBFolder): Path to xTB calculation folder or XTBFolder instance.
    """

    def __init__(self, folder):
        """
        Initialize XTBOutput with a folder path or XTBFolder instance.
        """
        self.folder = (
            folder if isinstance(folder, XTBFolder) else XTBFolder(folder)
        )

    def __getattr__(self, name):
        """
        Delegate attribute access to internal file parsers or cached properties.

        This allows XTBOutput to transparently expose all properties and methods
        from XTBMainOut, XTBEngradFile, XTBChargesFile, XTBWibergBondOrderFile,
        XTBHessianFile, XTBGradientFile, XTBVibSpectrumFile, and any future parsers.
        """
        # List of internal parsers to search, in order of priority
        file_parsers = [
            "main_out",
            "charges_file",
            "energy_file",
            "engrad_file",
            "g98_file",
            "gradient_file",
            "hessian_file",
            "vibspectrum_file",
            "wbo_file",
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
        raise AttributeError(
            f"'{type(self).__name__}' object has no attribute '{name}'"
        )

    @cached_property
    def main_out(self):
        """Main xTB output file parser (XTBMainOut)."""
        path = self.folder._xtb_out()
        return XTBMainOut(path) if path else None

    @cached_property
    def charges_file(self):
        """Atomic partial charges file parser (XTBChargesFile)."""
        path = self.folder._charges()
        return XTBChargesFile(path) if path else None

    @cached_property
    def energy_file(self):
        """Energy file parser (XTBEnergyFile)."""
        path = self.folder._energy()
        return XTBEnergyFile(path) if path else None

    @cached_property
    def engrad_file(self):
        """Energy and gradient file parser (XTBEngradFile)."""
        path = self.folder._engrad()
        return XTBEngradFile(path) if path else None

    @cached_property
    def g98_file(self):
        """GAUSSIAN-format vibrational frequencies file parser (XTBG98File)."""
        path = self.folder._g98_out()
        return XTBG98File(path) if path else None

    @cached_property
    def gradient_file(self):
        """Gradient file parser (XTBGradientFile)."""
        path = self.folder._gradient()
        return XTBGradientFile(path) if path else None

    @cached_property
    def hessian_file(self):
        """Hessian matrix file parser (XTBHessianFile)."""
        path = self.folder._hessian()
        return XTBHessianFile(path) if path else None

    @cached_property
    def vibspectrum_file(self):
        """Vibrational spectrum file parser (XTBVibSpectrumFile)."""
        path = self.folder._vibspectrum()
        return XTBVibSpectrumFile(path) if path else None

    @cached_property
    def wbo_file(self):
        """Wiberg bond order file parser (XTBWibergBondOrderFile)."""
        path = self.folder._wbo()
        return XTBWibergBondOrderFile(path) if path else None

    @cached_property
    def xtbopt_log_file(self):
        """Optimization trajectory log file parser (XYZFile)."""
        path = self.folder._xtbopt_log()
        return XYZFile(path) if path else None

    @cached_property
    def xtbopt_geometry_file(self):
        """Optimized geometry file parser (XYZFile or SDFFile)."""
        path = self.folder._xtbopt_geometry()
        if not path:
            return None
        if path.endswith(".xyz"):
            return XYZFile(path)
        elif path.endswith(".sdf"):
            return SDFFile(path)
        else:
            return None

    @cached_property
    def input_geometry_file(self):
        """Input geometry file parser (XYZFile or SDFFile)."""
        path = self.folder._input_geometry()
        if not path:
            return None
        if path.endswith(".xyz"):
            return XYZFile(path)
        elif path.endswith(".sdf"):
            return SDFFile(path)
        else:
            return None

    @property
    def normal_termination(self):
        """Check if the xTB calculation terminated normally."""
        if self.main_out:
            return self.main_out.normal_termination
        return False

    @property
    def geometry_optimization_converged(self):
        """Check if the geometry optimization converged."""
        if self.main_out:
            return self.main_out.geometry_optimization_converged
        return False

    @property
    def charge(self):
        """Get molecular charge from main output or charges file."""
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
    def mass(self):
        """Get molecular mass from main output."""
        if self.main_out:
            return self.main_out.molecular_mass
        return None

    @property
    def final_energy(self):
        """Get final converged energy in Hartree from main output, energy or .engrad file."""
        if self.main_out:
            return self.main_out.total_energy
        if self.energy_file:
            return self.energy_file.last_energy
        if self.engrad_file:
            return self.engrad_file.total_energy
        return None

    @property
    def final_forces(self):
        """Get final converged forces from .engrad file (if available)."""
        if self.engrad_file:
            return self.engrad_file.forces[-1]
        return None

    @cached_property
    def symbols(self):
        """Get atomic symbols from optimized geometry or first structure."""
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
    def xtbopt_geometry(self):
        """Read geometry from xtbopt.* file and return as Molecule object."""
        path = self.folder._xtbopt_geometry()
        molecule = self._read_geometry_file(path)
        if molecule:
            # Enrich molecule with properties from main output
            if self.charge is not None:
                molecule.charge = self.charge
            if self.multiplicity is not None:
                molecule.multiplicity = self.multiplicity
            if self.final_energy is not None:
                molecule.energy = self.final_energy
            if self.final_forces is not None:
                molecule.forces = self.final_forces
            logger.debug(f"Successfully read optimized geometry from {path}")
        return molecule

    @cached_property
    def input_geometry(self):
        """Read geometry from input geometry file and return as Molecule object."""
        path = self.folder._input_geometry()
        molecule = self._read_geometry_file(path)
        if molecule:
            if self.charge is not None:
                molecule.charge = self.charge
            if self.multiplicity is not None:
                molecule.multiplicity = self.multiplicity
            logger.debug(f"Successfully read input geometry from {path}")
        return molecule

    def _read_geometry_file(self, path):
        """Read geometry file (.xyz / .sdf) and return Molecule or None."""
        if not path:
            return None
        try:
            if path.endswith(".xyz"):
                return XYZFile(path).get_molecules(
                    index="-1", return_list=False
                )
            elif path.endswith(".sdf"):
                return SDFFile(path).get_molecule()
            else:
                return None
        except Exception as e:
            logger.error(f"Error reading geometry file {path}: {e}")
            return None

    @cached_property
    def all_structures(self):
        """
        Build Molecule objects from multiple data sources.

        This is where all data integration happens. Combines:
        - Coordinates and energies from available structure sources
        - Charge and multiplicity from main output
        - Forces from gradient file

        Returns:
            list[Molecule]: Molecule objects with all properties.
        """
        # 1. Choose orientations
        orientations, symbols, energies = self._choose_orientations()
        if not orientations:
            return []  # Nothing to build

        # 2. PBC conditions (xTB doesn't use PBC for optimizations)
        n = len(orientations)
        orientations_pbc = [None] * n

        # 3. Build complete Molecule objects using create_molecule_list
        logger.debug(f"Creating {n} Molecule objects with complete properties")
        molecules = create_molecule_list(
            orientations=orientations,
            orientations_pbc=orientations_pbc,
            energies=energies,
            forces=None,
            symbols=symbols,
            charge=self.charge,
            multiplicity=self.multiplicity,
            frozen_atoms=None,
            pbc_conditions=None,
        )
        # 4. Assign forces to the final structure if available
        if self.final_forces is not None and molecules:
            logger.debug(
                "Found gradient in .engrad file, assigning to last structure"
            )
            molecules[-1].forces = self.final_forces
        return molecules

    def _choose_orientations(self):
        """
        Choose the best available source of molecular orientations.

        Returns:
            tuple[list[np.ndarray], list[str], list[float] | None]
            orientations, symbols, energies
        """
        # Tier 1: xtbopt.log (trajectory)
        if self.xtbopt_log_file:
            logger.info(
                f"Reading molecule from file: {self.xtbopt_log_file.filepath}."
            )
            mols = self.xtbopt_log_file.get_molecules(":", return_list=True)
            if mols:
                return (
                    [mol.positions for mol in mols],
                    list(mols[0].symbols),
                    [mol.energy for mol in mols],
                )

        # Tier 2: xtbopt.xyz / xtbopt.sdf (optimized geometry)
        if self.xtbopt_geometry:
            logger.info(
                f"Reading molecule from file: {self.xtbopt_geometry_file.filepath}."
            )
            mol = self.xtbopt_geometry
            return [mol.positions], list(mol.symbols), [mol.energy]

        # Tier 3: g98 standard orientation (hess)
        if self.g98_file and self.g98_file.standard_orientation:
            logger.info(
                f"Reading molecule from file: {self.g98_file.filepath}."
            )
            orientation = list(self.g98_file.standard_orientation)
            symbols = list(self.g98_file.symbols)
            return [orientation], symbols, [self.final_energy]

        # Tier 4: input geometry (sp)
        if self.normal_termination and self.input_geometry:
            logger.info(
                f"Reading molecule from file: {self.input_geometry_file.filepath}."
            )
            mol = self.input_geometry
            return [mol.positions], list(mol.symbols), [self.final_energy]

        return [], None, None

    @property
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
        if (
            not self.normal_termination
            or not self.geometry_optimization_converged
        ):
            return None
        if self.xtbopt_geometry:
            return self.xtbopt_geometry
        if self.all_structures:
            return self.all_structures[-1]
        return None

    @cached_property
    def last_structure(self):
        if self.all_structures:
            return self.all_structures[-1]
        return None

    @property
    def molecule(self):
        """
        Alias for the last molecular structure.
        """
        return self.last_structure

    def get_molecule(self, index="-1"):
        """
        Get a specific molecule by index from the xTB folder.
        """
        index = string2index_1based(index)
        if not self.all_structures:
            raise ValueError(
                f"No molecule could be found in {self.folder.folderpath}."
            )
        return self.all_structures[index]
