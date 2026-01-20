import copy
import logging
import os
from functools import cached_property, lru_cache

import networkx as nx
import numpy as np
from ase import units
from ase.io import read as ase_read
from ase.symbols import Symbols
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Geometry import Point3D
from scipy.spatial.distance import cdist

from chemsmart.io.molecules import get_bond_cutoff
from chemsmart.utils.geometry import is_collinear
from chemsmart.utils.periodictable import PeriodicTable as pt
from chemsmart.utils.utils import file_cache, string2index_1based

p = pt()

logger = logging.getLogger(__name__)


class Molecule:
    """Class to represent a molcular structure.

    Parameters:

    symbols: a list of symbols.
    positions: a numpy array of atomic positions
        The shape of the array should be (n, 3) where n is the number
        of atoms in the molecule.
    charge: integer
        The charge of the molecule.
    multiplicity: integer
        The multiplicity of the molecule.
    frozen_atoms: list of integers, one for each atom, indicating which atoms are frozen.
        Follows Gaussian input file format where -1 denotes frozen atoms
        and 0 denotes relaxed atoms.
    pbc_conditions: list of integers
        The periodic boundary conditions for the molecule.
    translation_vectors: list of lists
        The translation vectors for the molecule.
    energy: float
        The energy of the molecule in eV.
    forces: numpy array
        The forces on the atoms in the molecule in eV/Å.
    velocities: numpy array
        The velocities of the atoms in the molecule.
    info: dict
        A dictionary containing additional information about the molecule.
    """

    def __init__(
        self,
        symbols=None,
        positions=None,
        charge=None,
        multiplicity=None,
        frozen_atoms=None,
        pbc_conditions=None,
        translation_vectors=None,
        energy=None,
        forces=None,
        velocities=None,
        vibrational_frequencies=None,
        vibrational_reduced_masses=None,
        vibrational_force_constants=None,
        vibrational_ir_intensities=None,
        vibrational_mode_symmetries=None,
        vibrational_modes=None,
        info=None,
    ):
        """
        Initialize molecular structure with atomic and quantum properties.
        """
        self.symbols = symbols
        self._positions = positions
        self.charge = charge
        self.multiplicity = multiplicity
        self.frozen_atoms = frozen_atoms
        self.pbc_conditions = pbc_conditions
        self.translation_vectors = translation_vectors
        self._energy = energy
        self.forces = forces
        self.velocities = velocities
        self.info = info
        self._num_atoms = len(self.symbols)

        # Define bond order classification multipliers (avoiding redundancy)
        # use the relationship between bond orders and bond lengths from J. Phys. Chem. 1959, 63, 8, 1346
        #                 1/Ls^3 : 1/Ld^3 : 1/Lt^3 = 2 : 3 : 4
        # From experimental data, the aromatic bond length in benzene is ~1.39 Å,
        # which is between single (C–C, ~1.54 Å) and double (C=C, ~1.34 Å) bonds.
        # Interpolate the aromatic bond length as L_ar ≈ L_s × (4/5)**1/3
        self.bond_length_multipliers = {
            "single": 1.0,
            "aromatic": (4 / 5) ** (1.0 / 3),  # Aromatic bond approximation
            "double": (2 / 3) ** (1.0 / 3),
            "triple": (1 / 2) ** (1.0 / 3),
        }

        # Validate essential molecular data
        if self.symbols is None or self.positions is None:
            raise ValueError(
                "Molecule must have symbols and positions defined."
            )

        # Ensure non-empty structure
        if len(self.symbols) == 0:
            raise ValueError("The number of symbols should not be empty!")

        # Validate symbols-positions consistency
        if len(self.symbols) != len(self.positions):
            logger.debug(f"Number of symbols: {len(self.symbols)}")
            logger.debug(f"Number of positions: {len(self.positions)}")
            logger.debug(f"Symbols: {self.symbols}")
            logger.debug(f"Positions: {self.positions}")
            raise ValueError(
                "The number of symbols and positions should be the same!"
            )

        self.vibrational_frequencies = (
            []
            if vibrational_frequencies is None
            else list(vibrational_frequencies)
        )
        self.vibrational_reduced_masses = (
            []
            if vibrational_reduced_masses is None
            else list(vibrational_reduced_masses)
        )
        self.vibrational_force_constants = (
            []
            if vibrational_force_constants is None
            else list(vibrational_force_constants)
        )
        self.vibrational_ir_intensities = (
            []
            if vibrational_ir_intensities is None
            else list(vibrational_ir_intensities)
        )
        self.vibrational_mode_symmetries = (
            []
            if vibrational_mode_symmetries is None
            else list(vibrational_mode_symmetries)
        )
        # modes is a list of (num_atoms x 3) arrays; keep as-is if provided
        self.vibrational_modes = (
            [] if vibrational_modes is None else vibrational_modes
        )

    def __len__(self):
        """
        Return the number of atoms in the molecule.
        """
        return len(self.chemical_symbols)

    def __getitem__(self, idx):
        """
        Get subset of molecule using 1-based indexing.

        Interprets the input idx as 1-based indices by index adjustment (i - 1).
        The method assumes that idx contains 1-based indices (e.g., [1, 2, 3]),
        so it subtracts 1 to convert them to Python's zero-based indexing.
        """
        symbols = [self.symbols[i - 1] for i in idx]
        positions = [self.positions[i - 1] for i in idx]
        return type(self)(symbols=symbols, positions=positions)

    @property
    def energy(self):
        """
        Total molecular energy in eV.
        """
        return self._energy

    @energy.setter
    def energy(self, value):
        """
        Set the molecular energy.
        """
        self._energy = value

    @property
    def positions(self):
        return self._positions

    @positions.setter
    def positions(self, value):
        arr = np.asarray(value, dtype=float)
        if arr.ndim != 2 or arr.shape[1] != 3:
            raise ValueError(f"positions must be (N, 3); got {arr.shape}")
        # store a copy to avoid accidental aliasing
        self._positions = arr.copy()

    @property
    def empirical_formula(self):
        """
        Empirical chemical formula using Hill notation.
        """
        return Symbols.fromsymbols(self.symbols).get_chemical_formula(
            mode="hill", empirical=True
        )

    @property
    def mass(self):
        """
        Total molecular mass using standard atomic masses.
        """
        return sum(p.to_atomic_mass(symbol) for symbol in self.symbols)

    @property
    def natural_abundance_weighted_mass(self):
        """
        Molecular mass weighted by natural isotope abundances.
        """
        return sum(
            p.to_weighted_atomic_mass_by_abundance(symbol)
            for symbol in self.symbols
        )

    @property
    def most_abundant_mass(self):
        """
        Molecular mass using most abundant isotopes.
        """
        return sum(
            p.to_most_abundant_atomic_mass(symbol) for symbol in self.symbols
        )

    @property
    def masses(self):
        """
        Numpy array of atomic masses of the molecule.
        """
        return np.array([p.to_atomic_mass(symbol) for symbol in self.symbols])

    @property
    def natural_abundance_weighted_masses(self):
        """
        Array of natural abundance weighted atomic masses.
        """
        return np.array(
            [
                p.to_weighted_atomic_mass_by_abundance(symbol)
                for symbol in self.symbols
            ]
        )

    @property
    def most_abundant_masses(self):
        """
        Array of most abundant isotope masses for each atom.
        """
        return np.array(
            [p.to_most_abundant_atomic_mass(symbol) for symbol in self.symbols]
        )

    @property
    def center_of_mass(self):
        """
        Compute the center of mass of the molecule.
        """
        return np.average(self.positions, axis=0, weights=self.masses)

    @property
    def chemical_formula(self):
        """
        Get the chemical formula in Hill notation.
        """
        return self.get_chemical_formula()

    @cached_property
    def chemical_symbols(self):
        """
        Return a list of chemical symbols strings
        """
        if self.symbols is not None:
            return list(self.symbols)

    @property
    def atomic_radii_list(self):
        """Return a list of atomic radii for each atom in the molecule."""
        return [p.covalent_radius(symbol) for symbol in self.symbols]

    @property
    def vdw_radii_list(self):
        """Return a list of van der Waals radii for each atom in the molecule."""
        return [p.vdw_radius(symbol) for symbol in self.symbols]

    @property
    def estimated_dispersion(self):
        """Estimated dispersion parameter for Voronoi tessellation.
        Returns:
        - dispersion (float): Estimated maximum distance for adjacent points.
        """
        n_points = self.distance_matrix.shape[0]
        max_distance = np.max(
            self.distance_matrix[np.triu_indices(n_points, k=1)]
        )

        max_radii_sum = 0.0
        radii = np.array(self.vdw_radii_list)
        if len(radii) != n_points:
            raise ValueError("Number of radii must match number of points.")
        for i in range(n_points):
            for j in range(i + 1, n_points):
                radii_sum = radii[i] + radii[j]
                max_radii_sum = max(max_radii_sum, radii_sum)

        # Use a factor of 1.5 to ensure sufficient dispersion for Voronoi tessellation
        dispersion = max(max_distance, max_radii_sum) * 1.5  # add 50% buffer
        return dispersion

    @property
    def voronoi_dirichlet_occupied_volume(self):
        """Calculate the occupied volume of the molecule using Voronoi-Dirichlet tessellation.

        Note: This method requires the pyvoro package, which can be installed with:
            pip install chemsmart[voronoi]
        Note: pyvoro requires Python < 3.12
        """
        from chemsmart.utils.geometry import (
            calculate_voronoi_dirichlet_occupied_volume,
        )

        return calculate_voronoi_dirichlet_occupied_volume(
            coords=self.positions,
            radii=self.vdw_radii_list,
            dispersion=self.estimated_dispersion,
        )

    @property
    def voronoi_dirichlet_polyhedra_occupied_volume(self):
        """Calculate the occupied volume of the molecule using Voronoi-Dirichlet Polyhedra (VDP)."""
        from chemsmart.utils.geometry import (
            calculate_molecular_volume_vdp,
        )

        return calculate_molecular_volume_vdp(
            coordinates=self.positions,
            vdw_radii=self.vdw_radii_list,
        )

    @property
    def crude_volume_by_atomic_radii(self):
        """Calculate the crude occupied volume of the molecule using atomic radii."""
        from chemsmart.utils.geometry import calculate_crude_occupied_volume

        return calculate_crude_occupied_volume(
            coords=self.positions, radii=self.atomic_radii_list
        )

    @property
    def crude_volume_by_vdw_radii(self):
        """Calculate the crude occupied volume of the molecule using van der Waals radii."""
        from chemsmart.utils.geometry import calculate_crude_occupied_volume

        return calculate_crude_occupied_volume(
            coords=self.positions, radii=self.vdw_radii_list
        )

    @property
    def vdw_volume(self):
        """Calculate the occupied volume of the molecule using van der Waals radii.

        Uses pairwise overlap correction. For more accurate results on complex
        molecules, consider using ``grid_vdw_volume`` instead.

        See Also
        --------
        grid_vdw_volume : Grid-based volume (more accurate for complex molecules)
        vdw_volume_from_rdkit : RDKit's grid-based implementation
        """
        from chemsmart.utils.geometry import calculate_vdw_volume

        return calculate_vdw_volume(
            coords=self.positions, radii=self.vdw_radii_list
        )

    @property
    def grid_vdw_volume(self):
        """Calculate the VDW volume using grid-based numerical integration.

        This method places the molecule in a 3D grid and counts grid points
        that fall inside any atomic VDW sphere. This approach correctly handles
        all orders of atomic overlaps and provides more accurate volume estimates
        for complex molecules compared to the pairwise method.

        This implementation is similar to RDKit's DoubleCubicLatticeVolume
        algorithm.

        Returns
        -------
        float
            Volume in cubic Ångstroms (Å³).

        See Also
        --------
        vdw_volume : Faster pairwise overlap method
        vdw_volume_from_rdkit : RDKit's native implementation
        """
        from chemsmart.utils.geometry import calculate_grid_vdw_volume

        return calculate_grid_vdw_volume(
            coords=self.positions, radii=self.vdw_radii_list
        )

    @property
    def vdw_volume_from_rdkit(self):
        """Calculate the van der Waals volume of the molecule using RDKit."""
        from rdkit.Chem.rdMolDescriptors import DoubleCubicLatticeVolume

        dclv = DoubleCubicLatticeVolume(self.to_rdkit())
        volume = dclv.GetVDWVolume()
        return volume

    @property
    def num_atoms(self):
        """
        Return the number of atoms in the molecule.
        """
        return self._num_atoms

    @num_atoms.setter
    def num_atoms(self, value):
        """
        Set the number of atoms in the molecule.
        """
        self._num_atoms = value

    @property
    def pbc(self):
        """
        Return the periodic boundary conditions.
        """
        if self.pbc_conditions is not None:
            return all(i == 0 for i in self.pbc_conditions)

    @property
    def is_chiral(self):
        """
        Check if molecule is chiral or not.
        """
        return Chem.FindMolChiralCenters(self.to_rdkit(), force=True) != []

    @property
    def is_aromatic(self):
        """
        Check if molecule is aromatic or not.
        """
        mol = self.to_rdkit()
        return any(atom.GetIsAromatic() for atom in mol.GetAtoms())

    @property
    def is_ring(self):
        """
        Check if molecule is a ring or not.
        """
        return Chem.GetSymmSSSR(self.to_rdkit()) != []

    @property
    def is_monoatomic(self):
        """
        Check if molecule is monoatomic or not.
        """
        return self.num_atoms == 1

    @property
    def is_diatomic(self):
        """
        Check if molecule is diatomic or not.
        """
        return self.num_atoms == 2

    @property
    def is_linear(self):
        """
        Check if molecule is linear or not.
        """
        if self.num_atoms <= 2:
            return True
        else:
            if self.num_atoms == 3:
                return is_collinear(self.positions)
            else:
                from sklearn.decomposition import PCA

                # Use PCA to check if all atoms lie on one principal axis
                pca = PCA(n_components=1)
                pca.fit(self.positions)
                reconstructed = pca.inverse_transform(
                    pca.transform(self.positions)
                )
                error = np.linalg.norm(
                    self.positions - reconstructed, axis=1
                ).max()
                return error < 1e-2

    @property
    def moments_of_inertia_tensor(self):
        """
        Calculate the moment of inertia tensor of the molecule.
        """
        moi_tensor, _, _ = self._get_moments_of_inertia
        return np.array(moi_tensor)

    @property
    def moments_of_inertia(self):
        """
        Obtain moments of inertia from molecular structure
        along principal axes as a list.
        """
        if self.is_monoatomic:
            return [0.0, 0.0, 0.0]
        else:
            _, eigenvalues, _ = self._get_moments_of_inertia
            return eigenvalues

    @property
    def moments_of_inertia_weighted_mass(self):
        """
        Get moments of inertia using natural abundance weighted masses.
        """
        if self.is_monoatomic:
            return [0.0, 0.0, 0.0]
        else:
            _, eigenvalues, _ = self._get_moments_of_inertia_weighted_mass
            return eigenvalues

    @property
    def moments_of_inertia_most_abundant_mass(self):
        """
        Get moments of inertia using most abundant isotope masses.
        """
        if self.is_monoatomic:
            return [0.0, 0.0, 0.0]
        else:
            _, eigenvalues, _ = self._get_moments_of_inertia_most_abundant_mass
            return eigenvalues

    @property
    def moments_of_inertia_principal_axes(self):
        """
        Obtain moments of inertia along principal axes from molecular structure.
        """
        _, _, eigenvectors = self._get_moments_of_inertia
        return eigenvectors

    @cached_property
    def _get_moments_of_inertia(self):
        """
        Calculate the moments of inertia of the molecule.
        Units of amu Å^2.
        """
        if self.num_atoms == 1:
            return np.zeros(3)
        else:
            from chemsmart.utils.geometry import calculate_moments_of_inertia

            return calculate_moments_of_inertia(self.masses, self.positions)

    @cached_property
    def _get_moments_of_inertia_weighted_mass(self):
        """
        Calculate the moments of inertia of the molecule. Use natural abundance weighted masses.
        Units of amu Å^2.
        """
        if self.num_atoms == 1:
            return np.zeros(3)
        else:
            from chemsmart.utils.geometry import calculate_moments_of_inertia

            return calculate_moments_of_inertia(
                self.natural_abundance_weighted_masses, self.positions
            )

    @cached_property
    def _get_moments_of_inertia_most_abundant_mass(self):
        """
        Calculate the moments of inertia of the molecule. Use most abundant masses.
        Units of amu Å^2.
        """
        if self.num_atoms == 1:
            return np.zeros(3)
        else:
            from chemsmart.utils.geometry import calculate_moments_of_inertia

            return calculate_moments_of_inertia(
                self.most_abundant_masses, self.positions
            )

    @cached_property
    def rotational_temperatures(self):
        """
        Obtain the rotational temperatures of the molecule in K.
        Θ_r,i = h^2 / (8 * pi^2 * I_i * k_B) for i = x, y, z
        """
        moi_in_SI_units = [
            float(i) * units._amu * (1 / units.m) ** 2
            for i in self.moments_of_inertia
        ]
        return [
            units._hplanck**2 / (8 * np.pi**2 * moi_in_SI_units[i] * units._k)
            for i in range(3)
        ]

    def get_chemical_formula(self, mode="hill", empirical=False):
        """
        Get chemical formula with flexible formatting options.
        """
        if self.symbols is not None:
            return Symbols.fromsymbols(self.symbols).get_chemical_formula(
                mode=mode, empirical=empirical
            )

    def get_distance(self, idx1, idx2):
        """
        Calculate the distance between two points.
        Use 1-based indexing for idx1 and idx2.
        """
        return np.linalg.norm(
            self.positions[idx1 - 1] - self.positions[idx2 - 1]
        )

    def get_angle(self, idx1, idx2, idx3):
        """
        Calculate the angle between three points.
        Use 1-based indexing for idx1, idx2, and idx3.
        """
        return self.get_angle_from_positions(
            self.positions[idx1 - 1],
            self.positions[idx2 - 1],
            self.positions[idx3 - 1],
        )

    def get_angle_from_positions(self, position1, position2, position3):
        """
        Calculate the angle between three points.
        """
        v1 = position1 - position2
        v2 = position3 - position2
        cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        return np.degrees(np.arccos(cos_theta))

    def get_dihedral(self, idx1, idx2, idx3, idx4):
        """
        Calculate the dihedral angle between four points, about bond formed by idx2 and idx3.
        Use 1-based indexing for idx1, idx2, idx3, and idx4.
        """
        return self.get_dihedral_from_positions(
            self.positions[idx1 - 1],
            self.positions[idx2 - 1],
            self.positions[idx3 - 1],
            self.positions[idx4 - 1],
        )

    def get_dihedral_from_positions(
        self, position1, position2, position3, position4
    ):
        """
        Calculate the dihedral angle between four points.
        """
        v1 = position1 - position2
        v2 = position3 - position2
        v3 = position4 - position3
        n1 = np.cross(v1, v2)
        n2 = np.cross(v2, v3)
        x = np.dot(n1, n2)
        y = np.dot(np.cross(n1, v2), n2)
        return np.degrees(np.arctan2(y, x))

    def copy(self):
        """
        Create a deep copy of the molecule.
        """
        return copy.deepcopy(self)

    @property
    def has_vibrations(self) -> bool:
        """True if any vibrational frequencies are present."""
        return bool(self.vibrational_frequencies)

    @property
    def num_vib_frequencies(self) -> int:
        """Number of vibrational frequencies available."""
        return len(self.vibrational_frequencies)

    @property
    def num_vib_modes(self) -> int:
        """Number of normal modes available."""
        return len(self.vibrational_modes)

    @classmethod
    def from_coordinate_block_text(cls, coordinate_block):
        """
        Create molecule from coordinate block text.
        """
        cb = CoordinateBlock(coordinate_block=coordinate_block)
        return cls(
            symbols=cb.symbols,
            positions=cb.positions,
            frozen_atoms=cb.constrained_atoms,
            translation_vectors=cb.translation_vectors,
        )

    @classmethod
    def from_symbols_and_positions_and_pbc_conditions(
        cls, list_of_symbols, positions, pbc_conditions=None
    ):
        """
        Create molecule from symbols, positions and periodic boundary conditions.
        """
        return cls(
            symbols=Symbols.fromsymbols(list_of_symbols),
            positions=positions,
            pbc_conditions=pbc_conditions,
        )

    @classmethod
    def from_filepath(cls, filepath, index="-1", return_list=False, **kwargs):
        """
        Create molecule from various file formats.
        """
        filepath = os.path.abspath(filepath)
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"{filepath} could not be found!")

        if os.path.getsize(filepath) == 0:
            return None

        molecule = cls._read_filepath(
            filepath, index=index, return_list=return_list, **kwargs
        )
        if return_list and not isinstance(molecule, list):
            return [molecule]
        else:
            return molecule

    @classmethod
    def _read_filepath(cls, filepath, index, return_list, **kwargs):
        """
        Internal method to read molecular data from various file formats.
        """
        basename = os.path.basename(filepath)

        if basename.endswith(".xyz"):
            logger.debug(f"Reading xyz file: {filepath}")
            return cls._read_xyz_file(
                filepath=filepath,
                index=index,
                return_list=return_list,
            )

        if basename.endswith(".sdf"):
            return cls._read_sdf_file(filepath)

        if basename.endswith((".com", ".gjf")):
            return cls._read_gaussian_inputfile(filepath)

        if basename.endswith(".log"):
            return cls._read_gaussian_logfile(filepath, index, **kwargs)

        if basename.endswith(".inp"):
            return cls._read_orca_inputfile(filepath, **kwargs)

        if basename.endswith(".out"):
            from chemsmart.utils.io import get_program_type_from_file

            program = get_program_type_from_file(filepath)
            if program == "orca":
                return cls._read_orca_outfile(filepath, index, **kwargs)
            if program == "gaussian":
                return cls._read_gaussian_logfile(filepath, index, **kwargs)
            raise ValueError(
                f"Unsupported .out file program type: {program}. "
                "Only Gaussian and ORCA are currently supported."
            )

        if basename.endswith(".gro"):
            return cls._read_gromacs_gro(filepath, index, **kwargs)

        if basename.endswith(".trr"):
            return cls._read_gromacs_trr(filepath, index, **kwargs)

        # if basename.endswith(".traj"):
        #     return cls._read_traj_file(filepath, index, **kwargs)

        if basename.endswith((".cdx", ".cdxml")):
            return cls._read_chemdraw_file(
                filepath=filepath,
                index=index,
                return_list=return_list,
            )

        if basename.endswith((".png", ".jpg", ".jpeg", ".tif", ".tiff")):
            return cls._read_image_file(filepath)

        return cls._read_other(filepath, index, **kwargs)

    @classmethod
    def _read_xyz_file(cls, filepath, index=":", return_list=False):
        """
        Read XYZ format molecular structure file.
        """
        from chemsmart.io.xyz.xyzfile import XYZFile

        xyz_file = XYZFile(filename=filepath)
        molecules = xyz_file.get_molecules(
            index=index, return_list=return_list
        )
        return molecules

    @staticmethod
    @file_cache()
    def _read_sdf_file(filepath):
        """
        Read SDF format molecular structure file.
        """
        from chemsmart.io.file import SDFFile

        sdf_file = SDFFile(filepath)
        return sdf_file.molecule

    @staticmethod
    @file_cache()
    def _read_gaussian_inputfile(filepath):
        """
        Read Gaussian input file (.com/.gjf) format.
        """
        from chemsmart.io.gaussian.input import Gaussian16Input

        try:
            g16_input = Gaussian16Input(filename=filepath)
            return g16_input.molecule
        except ValueError as e:
            # log the error or raise a more specific exception
            raise ValueError(
                f"Failed to read Gaussian input file {filepath}: {str(e)}"
            ) from e

    @staticmethod
    @file_cache()
    def _read_gaussian_logfile(filepath, index, **kwargs):
        """
        Returns a list of molecules.
        """
        from chemsmart.io.gaussian.output import Gaussian16Output

        g16_output = Gaussian16Output(filename=filepath, **kwargs)
        return g16_output.get_molecule(index=index)

    @staticmethod
    @file_cache()
    def _read_orca_inputfile(filepath):
        """
        Read ORCA input file (.inp) format.

        Args:
            filepath (str): Path to ORCA input file

        Returns:
            Molecule: Molecule object from ORCA input file
        """
        from chemsmart.io.orca.input import ORCAInput

        orca_input = ORCAInput(filename=filepath)
        return orca_input.molecule

    @staticmethod
    @file_cache()
    def _read_orca_outfile(filepath, index):
        """
        Read ORCA output file (.out) format.

        Args:
            filepath (str): Path to ORCA output file
            index (str or int): Index for multi-structure files

        Returns:
            Molecule: Molecule object from ORCA output file

        Note:
            TODO: Improve ORCAOutput object so that all structures
            can be obtained and returned via index
        """
        # TODO: to improve ORCAOutput object so that all the structures
        #  can be obtained and returned via index
        from chemsmart.io.orca.output import ORCAOutput

        orca_output = ORCAOutput(filename=filepath)
        return orca_output.get_molecule(index=index)

    @classmethod
    def _read_chemdraw_file(cls, filepath, index="-1", return_list=False):
        """
        Read ChemDraw file (.cdx or .cdxml) format.

        Args:
            filepath (str): Path to ChemDraw file (.cdx or .cdxml)
            index (str or int): Index for multi-structure files.
                Use "-1" for last molecule, ":" for all, or 1-based integer.
            return_list (bool): If True, return a list of molecules.

        Returns:
            Molecule or list[Molecule]: Molecule object(s) from ChemDraw file.

        Note:
            - .cdxml files are XML-based ChemDraw format.
            - .cdx files are binary ChemDraw format.
            - RDKit's MolsFromCDXMLFile supports both formats.
            - 3D coordinates are generated using RDKit's EmbedMolecule.
        """
        from chemsmart.io.file import CDXFile

        chemdraw_file = CDXFile(filename=filepath)
        return chemdraw_file.get_molecules(
            index=index, return_list=return_list
        )

    # @staticmethod
    # @file_cache()
    # def _read_gromacs_gro(filepath, index, **kwargs):
    #     # TODO: add handling of index
    #     from chemsmart.io.gromacs.outputs import GroGroOutput
    #
    #     gro_output = GroGroOutput(filename=filepath)
    #     return gro_output.get_gro()
    #
    # @staticmethod
    # @file_cache()
    # def _read_gromacs_trr(filepath, index, **kwargs):
    #     from chemsmart.io.gromacs.outputs import GroTrrOutput
    #
    #     trr_output = GroTrrOutput(filename=filepath)
    #     return trr_output.get_atoms(index=index)

    @classmethod
    def _read_image_file(cls, filepath):
        """
        Read molecular structure from an image file using DECIMER.

        Converts chemical structure images (PNG, JPG, JPEG, TIF, TIFF) to
        molecular structures by first extracting SMILES representation using
        DECIMER, then converting to a Molecule object via RDKit.
        
        Handles common chemical abbreviations (e.g., Ad for adamantyl, Ph for phenyl)
        by using OCR to detect text and expanding abbreviations to full structures.
        Abbreviation expansion requires pytesseract (optional dependency).

        Args:
            filepath (str): Path to the image file containing chemical structure

        Returns:
            Molecule: Molecule object created from the image

        Raises:
            ImportError: If DECIMER or opencv-python is not installed
            ValueError: If the image cannot be processed or converted
            
        Note:
            For better abbreviation support, install pytesseract:
            ``pip install chemsmart[image]`` or ``pip install pytesseract``
        """
        try:
            import cv2
            from DECIMER import predict_SMILES
        except ImportError as e:
            raise ImportError(
                "DECIMER and opencv-python are required to read image files. "
                "Install them with: `pip install decimer opencv-python`"
            ) from e

        # Try to import pytesseract for OCR (optional)
        try:
            import pytesseract
            has_ocr = True
        except ImportError:
            has_ocr = False
            logger.debug("pytesseract not available, abbreviation expansion disabled")

        # Common chemical abbreviations and their SMILES representations
        # These are neutral radicals/groups; connection point is implicit
        abbreviations = {
            'Ad': 'C1C2CC3CC1CC(C2)C3',  # Adamantyl (1-adamantyl)
            'Ph': 'c1ccccc1',  # Phenyl
            'Me': 'C',  # Methyl
            'Et': 'CC',  # Ethyl
            'nPr': 'CCC',  # n-Propyl
            'iPr': 'C(C)C',  # Isopropyl
            'Bu': 'CCCC',  # Butyl
            'iBu': 'CC(C)C',  # Isobutyl
            'sBu': 'C(C)CC',  # sec-Butyl
            'tBu': 'C(C)(C)C',  # tert-Butyl
            'Bn': 'Cc1ccccc1',  # Benzyl
            'Ac': 'C(=O)C',  # Acetyl
            'Bz': 'C(=O)c1ccccc1',  # Benzoyl
            'Ts': 'S(=O)(=O)c1ccc(C)cc1',  # Tosyl
            'Ms': 'S(=O)(=O)C',  # Mesyl
            'Tf': 'S(=O)(=O)C(F)(F)F',  # Triflyl
            'Cy': 'C1CCCCC1',  # Cyclohexyl
        }

        # Load the image
        img_orig = cv2.imread(filepath, cv2.IMREAD_GRAYSCALE)
        if img_orig is None:
            raise FileNotFoundError(filepath)

        # Detect abbreviations using OCR if available
        detected_abbrevs = {}
        detected_text = ""
        if has_ocr:
            try:
                # Use OCR to detect text in the image
                # PSM 6: Assume a single uniform block of text (good for chemical labels)
                detected_text = pytesseract.image_to_string(img_orig, config='--psm 6')
                logger.debug(f"OCR detected text: {detected_text}")
                
                # Check for known abbreviations in the detected text
                for abbrev, smiles in abbreviations.items():
                    # Case-sensitive match for chemical abbreviations
                    if abbrev in detected_text:
                        detected_abbrevs[abbrev] = smiles
                        logger.info(f"Detected abbreviation '{abbrev}' in image")
            except Exception as e:
                logger.debug(f"OCR failed: {e}")

        # Pre-process image for DECIMER (binarize + resize)
        img = cv2.resize(
            img_orig, None, fx=2.0, fy=2.0, interpolation=cv2.INTER_CUBIC
        )
        _, img = cv2.threshold(
            img, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU
        )

        # Write a temp preprocessed image (DECIMER takes a path)
        tmp = filepath + ".pre.png"
        cv2.imwrite(tmp, img)

        # Get SMILES from DECIMER
        smiles = predict_SMILES(tmp)  # DECIMER returns SMILES
        
        # If DECIMER fails or returns suspicious results when abbreviations are detected
        # (DECIMER typically fails with abbreviations), try to construct SMILES manually
        decimer_failed = smiles is None or not smiles.strip()
        should_use_abbrev = (
            detected_abbrevs 
            and (decimer_failed or len(smiles) < 3)  # Very short SMILES likely wrong
        )
        
        if should_use_abbrev:
            logger.warning(
                f"DECIMER returned incomplete/null SMILES: {smiles}. "
                f"Attempting to construct from detected abbreviations: {list(detected_abbrevs.keys())}"
            )
            
            # Store the original smiles (may be None)
            constructed_smiles = None
            
            # Try to construct molecule from abbreviations
            # For simple cases like "Ad-SH", we can combine parts
            if detected_text:
                # Simple heuristic: if we detected abbreviations, try common patterns
                text_upper = detected_text.upper()
                
                # Pattern: Ad-SH (adamantyl thiol)
                # Note: The dash might be a hyphen (-), en-dash (–), or em-dash (—)
                if 'Ad' in detected_abbrevs and 'SH' in text_upper:
                    constructed_smiles = detected_abbrevs['Ad'] + 'S'
                    logger.info(f"Constructed SMILES from Ad-SH pattern: {constructed_smiles}")
                # Pattern: Ph-X (phenyl with substituent)
                elif 'Ph' in detected_abbrevs and any(sep in detected_text for sep in ['-', '–', '—']):
                    # Try to get the substituent
                    parts = []
                    for sep in ['-', '–', '—']:
                        if sep in detected_text:
                            parts = detected_text.split(sep)
                            break
                    if len(parts) == 2:
                        substituent = parts[1].strip().upper()
                        # Map common substituents
                        if substituent == 'SH':
                            constructed_smiles = detected_abbrevs['Ph'] + 'S'
                        elif substituent == 'OH':
                            constructed_smiles = detected_abbrevs['Ph'] + 'O'
                        elif substituent in ['NH2', 'NH₂']:
                            constructed_smiles = detected_abbrevs['Ph'] + 'N'
                # If we have an abbreviation but couldn't construct, use it as-is
                elif len(detected_abbrevs) == 1:
                    constructed_smiles = list(detected_abbrevs.values())[0]
            
            # Use constructed SMILES if we managed to construct one
            if constructed_smiles:
                smiles = constructed_smiles
        
        if smiles is None:
            raise ValueError(
                f"Image {filepath} cannot be converted to SMILES."
            )

        return cls.from_smiles(smiles)

    @staticmethod
    @file_cache()
    def _read_other(filepath, index, **kwargs):
        """
        Reads a file using ASE and returns a Molecule object.
        """
        from .atoms import AtomsChargeMultiplicity

        # supplied index is 1-indexed, thus need to convert
        index = string2index_1based(index)

        ase_atoms = ase_read(filepath, index=index, **kwargs)
        logger.debug(f"Read ASE atoms: {ase_atoms} at index {index}")

        if isinstance(ase_atoms, list):
            logger.debug(f"Read {len(ase_atoms)} ASE atoms.")
            return [
                AtomsChargeMultiplicity.from_atoms(atoms).to_molecule()
                for atoms in ase_atoms
            ]
        return AtomsChargeMultiplicity.from_atoms(ase_atoms).to_molecule()

    @classmethod
    @lru_cache(maxsize=128)
    def from_pubchem(cls, identifier, return_list=False):
        """
        Create molecule object from PubChem database.

        Args:
            identifier (str): Compound identifier (name, CID, or SMILES string)
            return_list (bool): Whether to return list format. Default False

        Returns:
            Molecule or list or None: Molecule object from PubChem, None if not found

        Raises:
            requests.exceptions.RequestException: For network or HTTP-related issues
        """
        from chemsmart.io.molecules.pubchem import pubchem_search

        possible_attributes = (
            ["cid"]
            if identifier.isnumeric()
            else ["smiles", "name", "conformer"]
        )

        for attribute in possible_attributes:
            molecule = pubchem_search(**{attribute: identifier})
            if molecule is not None:
                logger.info(
                    f"Structure successfully created from pubchem with {attribute} = {identifier}"
                )
                if return_list:
                    return [molecule]
                return molecule

        logger.debug("Could not create structure from pubchem.")
        return None

    @classmethod
    def from_molecule(cls, molecule):
        """
        Create molecule from another molecule object.
        """
        return cls(**molecule.__dict__)

    @classmethod
    def from_ase_atoms(cls, atoms):
        """
        Creates a Molecule instance from an ASE Atoms object.
        """
        from .atoms import AtomsChargeMultiplicity

        return AtomsChargeMultiplicity.from_atoms(atoms).to_molecule()

    @classmethod
    def from_rdkit_mol(cls, rdMol: Chem.Mol) -> "Molecule":
        """
        Creates a Molecule instance from an RDKit Mol object, assuming a single conformer.
        """
        if rdMol is None:
            raise ValueError("Invalid RDKit molecule provided.")

        num_confs = rdMol.GetNumConformers()
        if num_confs == 0:
            raise ValueError(
                "RDKit molecule has no conformers. Ensure 3D coordinates are generated."
            )
        if num_confs > 1:
            raise ValueError(
                f"Expected a single conformer but found {num_confs}. Please provide a single conformer."
            )

        # Extract atomic symbols
        symbols = [atom.GetSymbol() for atom in rdMol.GetAtoms()]

        # Extract atomic positions from the first conformer (assuming single conformer)
        conf = rdMol.GetConformer(0)
        positions = np.array(
            [
                [
                    conf.GetAtomPosition(i).x,
                    conf.GetAtomPosition(i).y,
                    conf.GetAtomPosition(i).z,
                ]
                for i in range(rdMol.GetNumAtoms())
            ]
        )

        # Get molecular charge
        charge = Chem.GetFormalCharge(rdMol)

        # Estimate multiplicity (multiplicity = 2S + 1, where S is total spin)
        num_radical_electrons = sum(
            atom.GetNumRadicalElectrons() for atom in rdMol.GetAtoms()
        )
        multiplicity = (
            num_radical_electrons + 1 if num_radical_electrons > 0 else 1
        )  # Default to singlet

        # Store additional RDKit properties
        info = {
            "smiles": Chem.MolToSmiles(rdMol),
            "num_bonds": rdMol.GetNumBonds(),
        }

        rdkit_mol = cls(
            symbols=symbols,
            positions=positions,
            charge=charge,
            multiplicity=multiplicity,
            info=info,
        )

        rdkit_mol.num_atoms = rdMol.GetNumAtoms()

        return rdkit_mol

    @classmethod
    def from_smiles(cls, smiles):
        # Convert SMILES to RDKit molecule
        rdkit_mol = Chem.MolFromSmiles(smiles)
        if rdkit_mol is None:
            raise ValueError(
                f"Failed to create valid molecule from SMILES: {smiles}"
            )

        # Generate 3D coordinates
        from rdkit.Chem import AllChem

        rdkit_mol = Chem.AddHs(rdkit_mol)
        AllChem.EmbedMolecule(rdkit_mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(rdkit_mol)

        # Convert RDKit molecule to Molecule object using the classmethod
        molecule = cls.from_rdkit_mol(rdkit_mol)

        logger.info(
            f"Successfully created molecule from image with formula: "
            f"{molecule.chemical_formula}"
        )
        return molecule

    def write_coordinates(self, f, program=None):
        """
        Write molecular coordinates to file in specified format.

        Args:
            f (file): File object to write coordinates to
            program (str, optional): Format to use ('gaussian' or 'orca'). Default 'gaussian'

        Raises:
            ValueError: If program format is not supported

        Note:
            No empty end line at the end of the file
        """
        if program is None:
            program = "gaussian"  # use gaussian format by default
        if program.lower() == "gaussian":
            self._write_gaussian_coordinates(f)
            self._write_gaussian_pbc_coordinates(f)
        elif program.lower() == "orca":
            self._write_orca_coordinates(f)
            self._write_orca_pbc_coordinates(f)
        # elif program.lower() == "gromacs":
        # can implement other programs formats to write the coordinates for
        else:
            raise ValueError(
                f"Program {program} is not supported for writing coordinates."
            )

    def write(self, filename, format="xyz", mode="w", **kwargs):
        """
        Write molecule to file in specified format.

        Args:
            filename (str): Output file path
            format (str): File format ('xyz' or 'com'). Default 'xyz'
            mode (str): File write mode. Default 'w'
            **kwargs: Additional keyword arguments for format-specific writers

        Raises:
            ValueError: If format is not supported
        """
        if format.lower() == "xyz":
            self.write_xyz(filename, mode=mode, **kwargs)
        elif format.lower() == "com":
            self.write_com(filename, **kwargs)
        # elif format.lower() == "mol":
        #     self.write_mol(filename, **kwargs)
        else:
            raise ValueError(f"Format {format} is not supported for writing.")

    def write_xyz(self, filename, mode, **kwargs):
        """
        Write molecule to XYZ format file.

        Args:
            filename (str): Output XYZ file path
            mode (str): File write mode
            **kwargs: Additional keyword arguments (unused)
        """
        with open(filename, mode) as f:
            base_filename = os.path.basename(filename)
            if self.energy is not None:
                # energy found in file, e.g., .out, .log
                xyz_info = (
                    f"{base_filename}    Empirical formula: {self.chemical_formula}    "
                    f"Energy(Hartree): {self.energy:.6f}    "
                )
            else:
                # no energy found in file, e.g., .xyz or .com
                xyz_info = f"{base_filename}    Empirical formula: {self.chemical_formula}"

            logger.info(f"Writing outputfile to {filename}")
            f.write(f"{self.num_atoms}\n")
            f.write(f"{xyz_info}\n")
            self._write_orca_coordinates(f)

    def write_com(
        self,
        filename,
        charge=0,
        multiplicity=1,
        route="# opt freq m062x def2svp",
        **kwargs,
    ):
        """
        Write the molecule to a Gaussian input file.
        """
        with open(filename, "w") as f:
            basename = os.path.basename(filename).split(".")[0]
            f.write(f"%chk={basename}.chk\n")
            f.write("%mem=2GB\n")
            f.write("%nprocshared=32\n")
            f.write(f"{route}\n")  # example route
            f.write("\n")
            f.write(f"Generated from {filename}\n")
            f.write("\n")
            if self.charge is not None and self.multiplicity is not None:
                f.write(f"{self.charge} {self.multiplicity}\n")
            else:
                f.write(f"{charge} {multiplicity}\n")
            self.write_coordinates(f, program="gaussian")
            f.write("\n")

    def _write_gaussian_coordinates(self, f):
        """
        Write coordinates in Gaussian format.
        """
        assert self.symbols is not None, "Symbols to write should not be None!"
        assert (
            self.positions is not None
        ), "Positions to write should not be None!"
        if self.frozen_atoms is None or len(self.frozen_atoms) == 0:
            for i, (s, (x, y, z)) in enumerate(
                zip(self.chemical_symbols, self.positions)
            ):
                f.write(f"{s:5} {x:15.10f} {y:15.10f} {z:15.10f}\n")
        else:
            for i, (s, (x, y, z)) in enumerate(
                zip(self.chemical_symbols, self.positions)
            ):
                f.write(
                    f"{s:6} {self.frozen_atoms[i]:5} {x:15.10f} {y:15.10f} {z:15.10f}\n"
                )

    def _write_gaussian_pbc_coordinates(self, f):
        """
        Write the coordinates of the molecule with PBC conditions to a file.
        """
        if self.pbc_conditions is None or not any(self.pbc_conditions):
            # this happens when self.pbc_conditions = [False, False, False]
            # when the structure is read in from e.g., ASE database
            logger.debug("No PBC conditions to write.")
            return
        else:
            logger.debug(f"Writing PBC conditions: {self.pbc_conditions}")
            assert (
                self.translation_vectors is not None
            ), "Translation vectors should not be None when PBC conditions are given!"
            for i in range(len(self.translation_vectors)):
                f.write(
                    f"TV    {self.translation_vectors[i][0]:15.10f} "
                    f"{self.translation_vectors[i][1]:15.10f} "
                    f"{self.translation_vectors[i][2]:15.10f}\n"
                )

    def _write_orca_coordinates(self, f):
        """
        Write coordinates in ORCA format.
        """
        assert self.symbols is not None, "Symbols to write should not be None!"
        assert (
            self.positions is not None
        ), "Positions to write should not be None!"

        # if self.frozen_atoms is None:
        # commented above out since with frozen atom or not, the geometry is written the same way
        for i, (s, (x, y, z)) in enumerate(
            zip(self.chemical_symbols, self.positions)
        ):
            f.write(f"{s:5} {x:15.10f} {y:15.10f} {z:15.10f}\n")

    def _write_orca_pbc_coordinates(self, f):
        """
        Write PBC coordinates for ORCA.
        """
        # ORCA cannot do PBC calculations
        pass

    def __repr__(self):
        """
        Return string representation of molecule.
        """
        return f"{self.__class__.__name__}<{self.empirical_formula},energy: {self.energy}>"

    def __str__(self):
        """
        Return string representation of molecule.
        """
        return f"{self.__class__.__name__}<{self.empirical_formula},energy: {self.energy}>"

    @cached_property
    def distance_matrix(self):
        """
        Compute pairwise distance matrix.
        """
        return cdist(self.positions, self.positions)

    def determine_bond_order_one_bond(self, bond_length, bond_cutoff):
        """
        Determine bond order for single bond based on length and cutoff.
        """
        for bond_type, multiplier in [
            ("triple", 3),
            ("double", 2),
            ("aromatic", 1.5),
            ("single", 1),
        ]:
            if (
                bond_length
                < bond_cutoff * self.bond_length_multipliers[bond_type]
            ):
                return multiplier
        return 0  # No bond

    def determine_bond_order(self, bond_length, bond_cutoff):
        """
        Vectorized determination of bond order based on bond length and cutoff.
        """
        multipliers = np.array([3, 2, 1.5, 1])
        bond_types = ["triple", "double", "aromatic", "single"]

        # Compute bond order matrix
        bond_multiplier_matrix = np.array(
            [
                self.bond_length_multipliers[bond_type]
                for bond_type in bond_types
            ]
        )

        valid_bond = bond_length[..., np.newaxis] < (
            bond_cutoff[..., np.newaxis] * bond_multiplier_matrix
        )
        bond_order = np.where(valid_bond, multipliers, 0).max(axis=-1)

        return bond_order

    def bond_lengths(self):
        """
        Get all bond distances in the molecule.
        """
        # get all bond distances in the molecule
        return self.get_all_distances()

    def get_all_distances(self):
        """
        Calculate all pairwise atomic distances in the molecule.
        """
        bond_distances = []
        for i in range(self.num_atoms):
            for j in range(i + 1, self.num_atoms):
                bond_distances.append(
                    np.linalg.norm(self.positions[i] - self.positions[j])
                )
        return bond_distances

    def to_smiles(self):
        """
        Convert molecule to SMILES string.
        """
        # Create an RDKit molecule
        rdkit_mol = self.to_rdkit()

        # Convert RDKit molecule to SMILES
        return Chem.MolToSmiles(rdkit_mol)

    def to_rdkit(self, add_bonds=True, bond_cutoff_buffer=0.05, adjust_H=True):
        """Convert Molecule object to RDKit Mol with proper stereochemistry handling.
        Args:
            add_bonds (bool): Flag to add bonds to molecule or not.
            bond_cutoff_buffer (float): Additional buffer for bond cutoff distance.
            From testing, see test_resonance_handling, it seems that a value of 0.1Å
            works for ozone, acetone, benzene, and probably other molecules, too.
            adjust_Hs (bool): Adjust bond distances to H atoms.
        Returns:
            RDKit Mol: RDKit molecule object.
        """

        # Create molecule and add atoms
        rdkit_mol = Chem.RWMol()

        # Add atoms with correct element types
        for symbol in self.symbols:
            rdkit_mol.AddAtom(Chem.Atom(symbol))

        # Apply bond detection algorithm if requested
        if add_bonds:
            rdkit_mol = self._add_bonds_to_rdkit_mol(
                rdkit_mol, bond_cutoff_buffer, adjust_H
            )

        # Create conformer with 3D coordinates
        conformer = rdchem.Conformer(len(self.symbols))
        for i, pos in enumerate(self.positions):
            conformer.SetAtomPosition(i, Point3D(*pos))
        rdkit_mol.AddConformer(conformer)

        # Compute valences (Fix for getNumImplicitHs issue)
        rdkit_mol.UpdatePropertyCache(strict=False)

        # Partial sanitization for stereochemistry detection
        # Chem.SanitizeMol(rdkit_mol,
        #                  Chem.SANITIZE_ALL ^ Chem.SANITIZE_ADJUSTHS ^ Chem.SANITIZE_SETAROMATICITY)

        # I comment the following out since we do not want to modify the molecule
        # Validate the RDKit molecule
        # try:
        #     Chem.SanitizeMol(rdkit_mol)
        # except Chem.AtomValenceException as e:
        #     raise ValueError(f"Sanitization failed: {e}") from e

        # Detect stereochemistry from 3D coordinates
        Chem.AssignStereochemistryFrom3D(rdkit_mol, conformer.GetId())
        Chem.AssignAtomChiralTagsFromStructure(rdkit_mol, conformer.GetId())

        # Force update of stereo flags
        Chem.FindPotentialStereoBonds(rdkit_mol, cleanIt=True)

        return rdkit_mol.GetMol()

    def _add_bonds_to_rdkit_mol(
        self, rdkit_mol, bond_cutoff_buffer=0.05, adjust_H=True
    ):
        """
        Add bonds to the RDKit molecule.
        """
        for i in range(len(self.symbols)):
            for j in range(i + 1, len(self.symbols)):
                if adjust_H:
                    if self.symbols[i] == "H" and self.symbols[j] == "H":
                        # bond length of H-H is 0.74 Å
                        # covalent radius of H is 0.31 Å
                        cutoff_buffer = 0.2
                    elif self.symbols[i] == "H" or self.symbols[j] == "H":
                        # C-H bond distance of ~ 1.09 Å
                        # N-H bond distance of ~ 1.01 Å
                        # O-H bond distance of ~ 0.96 Å
                        # covalent radius of C is  0.76,
                        # covalent radius of N is 0.71,
                        # covalent radius of O is 0.66,
                        cutoff_buffer = 0.1
                    else:
                        cutoff_buffer = bond_cutoff_buffer
                else:
                    cutoff_buffer = 0.3  # default buffer
                cutoff = get_bond_cutoff(
                    self.symbols[i], self.symbols[j], cutoff_buffer
                )
                bond_order = self.determine_bond_order_one_bond(
                    bond_length=self.distance_matrix[i, j], bond_cutoff=cutoff
                )
                logger.debug(f"bond order: {bond_order}")
                if bond_order > 0:
                    bond_type = {
                        1: Chem.BondType.SINGLE,
                        1.5: Chem.BondType.AROMATIC,
                        2: Chem.BondType.DOUBLE,
                        3: Chem.BondType.TRIPLE,
                    }[bond_order]
                    rdkit_mol.AddBond(i, j, bond_type)
        return rdkit_mol

    def _add_bonds_to_rdkit_mol_vectorized(
        self, rdkit_mol, bond_cutoff_buffer=0.05, adjust_H=True
    ):
        """
        Add bonds to the RDKit molecule using a vectorized approach to compute bond orders.
        """
        num_atoms = len(self.symbols)

        # 1. Build an NxN cutoff matrix
        cutoff_matrix = np.zeros((num_atoms, num_atoms))
        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                # Decide on a pair-specific buffer
                if adjust_H:
                    if self.symbols[i] == "H" and self.symbols[j] == "H":
                        # bond length of H-H is ~0.74 Å
                        cutoff_buffer_ij = 0.2
                    elif self.symbols[i] == "H" or self.symbols[j] == "H":
                        # bond length to H is ~1.0–1.1 Å
                        cutoff_buffer_ij = 0.1
                    else:
                        cutoff_buffer_ij = bond_cutoff_buffer
                else:
                    cutoff_buffer_ij = 0.3  # some default if not adjusting H

                cutoff_ij = get_bond_cutoff(
                    self.symbols[i], self.symbols[j], cutoff_buffer_ij
                )
                cutoff_matrix[i, j] = cutoff_ij
                cutoff_matrix[j, i] = cutoff_ij

        # 2. Vectorized bond order calculation
        bond_orders = self.determine_bond_order(
            bond_length=self.distance_matrix,  # NxN distances
            bond_cutoff=cutoff_matrix,  # NxN cutoffs
        )

        # 3. Add edges (bonds) for non-zero bond orders
        bond_type_map = {
            1.0: Chem.BondType.SINGLE,
            1.5: Chem.BondType.AROMATIC,
            2.0: Chem.BondType.DOUBLE,
            3.0: Chem.BondType.TRIPLE,
        }

        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                bo = bond_orders[i, j]
                if bo > 0:
                    bond_type = bond_type_map.get(bo, Chem.BondType.SINGLE)
                    rdkit_mol.AddBond(i, j, bond_type)

        return rdkit_mol

    @cached_property
    def rdkit_fingerprints(self):
        """
        Return RDKit molecular fingerprints.
        """
        rdkit_mol = self.to_rdkit()
        return Chem.RDKFingerprint(rdkit_mol)

    @cached_property
    def bond_orders(self):
        """
        Return a list of bond orders from the molecular graph.
        Note that in conformers analysis, the bond orders should
        be the same for all conformers. In those cases, its best
        to use get_bond_orders_from_rdkit_mol(bond_cutoff_buffer=0.0)
        or get_bond_orders_from_graph(bond_cutoff_buffer=0.0) directly.
        """
        try:
            return self.get_bond_orders_from_graph()
        except Exception:
            return self.get_bond_orders_from_rdkit_mol()

    def get_bond_orders_from_rdkit_mol(self, **kwargs):
        """
        Return a list of bond orders from the RDKit molecule.
        """
        return [
            bond.GetBondTypeAsDouble()
            for bond in self.to_rdkit(**kwargs).GetBonds()
        ]

    def get_bond_orders_from_graph(self, **kwargs):
        """
        Return a list of bond orders from the molecular graph.
        """
        graph = self.to_graph(**kwargs)
        bond_orders = []
        for bond in graph.edges.values():
            bond_orders.append(bond["bond_order"])
        return bond_orders

    def to_graph(self, bond_cutoff_buffer=0.05, adjust_H=True) -> nx.Graph:
        """
        Convert a Molecule object to a connectivity graph with vectorized calculations.
        Bond cutoff value determines the maximum distance between two atoms
        to add a graph edge between them. Bond cutoff is obtained using Covalent
        Radii between the atoms via 𝑅_cutoff = 𝑅_𝐴 + 𝑅_𝐵 + tolerance_buffer.
        Args:
            bond_cutoff_buffer (float): Additional buffer for bond cutoff distance.
            adjust_H (bool): Whether to adjust hydrogen bond cutoffs.

        Returns:
            nx.Graph: A networkx graph object representing the molecule.
        """

        G = nx.Graph()
        positions = np.array(self.positions)
        symbols = np.array(self.chemical_symbols)
        num_atoms = len(symbols)

        # Add nodes
        for i, symbol in enumerate(symbols):
            G.add_node(i, element=symbol)

        # Compute distance matrix if not already available
        distance_matrix = np.linalg.norm(
            positions[:, np.newaxis] - positions, axis=2
        )

        # Compute bond cutoff matrix
        cutoff_matrix = np.zeros((num_atoms, num_atoms))

        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                element_i, element_j = str(symbols[i]), str(symbols[j])

                cutoff_buffer = bond_cutoff_buffer
                if adjust_H:
                    if element_i == "H" and element_j == "H":
                        cutoff_buffer = 0.12
                    elif element_i == "H" or element_j == "H":
                        cutoff_buffer = 0.05

                cutoff_matrix[i, j] = cutoff_matrix[j, i] = get_bond_cutoff(
                    element_i, element_j, cutoff_buffer
                )

        # Determine bond orders
        bond_orders = self.determine_bond_order(
            bond_length=distance_matrix, bond_cutoff=cutoff_matrix
        )

        # Add edges where bond order > 0
        i_indices, j_indices = np.where(bond_orders > 0)

        for i, j in zip(i_indices, j_indices):
            if i < j:  # Avoid duplicate edges
                G.add_edge(i, j, bond_order=bond_orders[i, j])

        return G

    def to_graph_non_vectorized(
        self, bond_cutoff_buffer=0.05, adjust_H=True
    ) -> nx.Graph:
        """Convert a Molecule object to a connectivity graph, non-vectorized.
        Bond cutoff value determines the maximum distance between two atoms
        to add a graph edge between them. Bond cutoff is obtained using Covalent
        Radii between the atoms via 𝑅_cutoff = 𝑅_𝐴 + 𝑅_𝐵 + tolerance_buffer.
        Args:
            bond_cutoff_buffer (float): Additional buffer for bond cutoff distance.
        Returns:
            nx.Graph: A networkx graph object representing the molecule.
        """
        G = nx.Graph()
        positions = self.positions

        # add nodes
        for i, symbol in enumerate(self.chemical_symbols):
            G.add_node(i, element=symbol)

        # Add edges (bonds) with bond order
        for i in range(len(positions)):
            for j in range(i + 1, len(positions)):
                element_i, element_j = (
                    self.chemical_symbols[i],
                    self.chemical_symbols[j],
                )

                cutoff_buffer = bond_cutoff_buffer

                if adjust_H:
                    if element_i == "H" and element_j == "H":
                        # bond length of H-H is 0.74 Å
                        # covalent radius of H is 0.31 Å
                        cutoff_buffer = 0.12
                    elif element_i == "H" or element_j == "H":
                        # C-H bond distance of ~ 1.09 Å
                        # N-H bond distance of ~ 1.01 Å
                        # O-H bond distance of ~ 0.96 Å
                        # covalent radius of C is  0.76,
                        # covalent radius of N is 0.71,
                        # covalent radius of O is 0.66,
                        cutoff_buffer = 0.05

                cutoff = get_bond_cutoff(
                    self.symbols[i], self.symbols[j], cutoff_buffer
                )
                bond_order = self.determine_bond_order_one_bond(
                    bond_length=self.distance_matrix[i, j], bond_cutoff=cutoff
                )
                if bond_order > 0:
                    G.add_edge(i, j, bond_order=bond_order)
        return G

    def to_ase(self):
        """
        Convert molecule object to ASE atoms object.
        """
        from .atoms import AtomsChargeMultiplicity

        return AtomsChargeMultiplicity(
            symbols=self.chemical_symbols,
            positions=self.positions,
            pbc=self.pbc,
            cell=self.translation_vectors,
            charge=self.charge,
            multiplicity=self.multiplicity,
            frozen_atoms=self.frozen_atoms,
            energy=self.energy,
            forces=self.forces,
            velocities=self.velocities,
            info=self.info,
        )

    def to_pymatgen(self):
        """
        Convert molecule object to pymatgen IStructure.
        """

        from pymatgen.io.ase import AseAtomsAdaptor

        return AseAtomsAdaptor.get_molecule(atoms=self.to_ase())

    def to_X_data(self, wbo=False):
        """
        Convert molecule object to X_data for ML models.
        """
        if self.positions is None:
            raise ValueError(
                "Positions are not available in the molecule object."
            )

        # Ensure energy is always included
        energy_array = np.array(
            [self.energy if self.energy is not None else 0.0]
        )  # Ensures shape (1,)
        positions_array = (
            np.array(self.positions).flatten().reshape(1, -1)
        )  # Ensures shape (1, num_atoms*3)

        # Concatenate energy and positions
        X = np.hstack(
            [energy_array.reshape(1, -1), positions_array]
        )  # Ensures (1, num_features)
        if wbo:
            # Add Wiberg bond orders if requested
            bond_orders = np.array(self.bond_orders).flatten()
            bond_orders = bond_orders.reshape(1, -1)
            X = np.hstack([X, bond_orders])

        return X

    def vibrationally_displaced(
        self,
        mode_idx,
        amp: float = 0.1,
        nframes: int | None = None,
        phase: float = np.pi / 2,
        normalize: bool = False,
        return_xyz: bool = False,
    ):
        """
        Create a geometry (or trajectory) displaced along a *mass-weighted* normal mode.

        Args:
            mode_idx (int): Mode index (1-based, negatives allowed)
            amp (float): Target maximum atomic displacement in Å after normalization.
            nframes (int | None): If provided, generate that many frames over one period
                [0, 2π); else return a single frame at `phase`.
            phase (float): Phase angle (radians) for the single-frame case.
                Defaults to pi/2 so that sin(phase) = 1.
            normalize (bool): If True, scale the (un-weighted) mode so its largest
                per-atom displacement is 1.0, making `amp` the max displacement.
            return_xyz (bool): If True and `nframes` is set, return a multi-frame XYZ string.

        Returns:
            Molecule | list[Molecule] | str
        """

        if isinstance(mode_idx, (int, np.integer)):
            idx = int(mode_idx)
            logger.debug(f"Displacing vibrational mode number #{idx}")
            if idx < 0:
                mw_mode = np.asarray(self.vibrational_modes[idx], dtype=float)
            else:
                # accept 1-based indices
                mw_mode = np.asarray(
                    self.vibrational_modes[idx - 1], dtype=float
                )
        else:
            raise ValueError("Vibrational mode number should be integer.")

        if mw_mode.shape != (self.num_atoms, 3):
            raise ValueError(
                f"vibrational mode must have shape ({self.num_atoms}, 3); "
                f"got {mw_mode.shape} instead!"
            )

        logger.debug(
            "Convert from mass-weighted mode to Cartesian mode: \n"
            "multiplying by M^{-1/2} * mode"
        )
        cart_mode = mw_mode / np.sqrt(self.masses)[:, None]  # masses in amu

        if normalize:
            logger.debug("normalize so max per-atom displacement = 1")
            per_atom = np.linalg.norm(cart_mode, axis=1)
            max_disp = float(per_atom.max())
            if max_disp == 0.0:
                raise ValueError("Provided vibrational mode has zero norm.")
            cart_mode = cart_mode / max_disp

        R0 = np.asarray(self.positions, float)
        logger.debug(f"Original positions: {R0}")
        symbols = list(self.symbols)

        def make_mol(theta: float):
            """Function to change the positions of the molecule by
            displacing along the vibrational mode."""
            Rt = R0 + (amp * np.sin(theta)) * cart_mode
            Rt_rounded = np.round(Rt, decimals=6)
            logger.debug(f"New positions: {Rt_rounded}")
            m = self.copy()
            m.positions = Rt_rounded
            return m

        if nframes is None:
            logger.debug(f"Return single frame/molecule with phase: {phase}")
            return make_mol(phase)

        thetas = np.linspace(0.0, 2.0 * np.pi, int(nframes), endpoint=False)
        frames = [make_mol(t) for t in thetas]
        logger.debug(
            f"Return {len(frames)} number of frames/molecules as a list."
        )

        if return_xyz:
            lines = []
            for t, m in zip(thetas, frames):
                lines.append(str(len(symbols)))
                lines.append(f"vib frame t={t:.6f} rad, amp={amp:.4f} Å")
                for s, (x, y, z) in zip(symbols, m.positions):
                    lines.append(f"{s:4s} {x: .6f} {y: .6f} {z: .6f}")
            return "\n".join(lines)

        return frames


class CoordinateBlock:
    """
    Class to create coordinate block object to abstract the geometry.
    """

    def __init__(self, coordinate_block):
        """
        Accepts a coordinate block either as text string or as a list of lines.
        If former, then convert to the latter before future usage.
        """
        coordinate_block_list = []
        if isinstance(coordinate_block, str):
            # Parse text format with newline separation
            for line in coordinate_block.split("\n"):
                coordinate_block_list.append(line.strip())
        elif isinstance(coordinate_block, list):
            # Use list format directly
            coordinate_block_list = coordinate_block
        else:
            raise TypeError(
                f"Coordinate block must be str or list, "
                f"got {type(coordinate_block)} instead!"
            )
        self.coordinate_block = coordinate_block_list

    @property
    def chemical_symbols(self):
        """
        Returns a list of chemical symbols for the molecule.
        """
        return self._get_symbols()

    @property
    def positions(self):
        """
        Returns a list of positions for the molecule.
        """
        return self._get_positions()

    @property
    def translation_vectors(self):
        """
        Return a list of translation vectors for systems with pbc.
        """
        return self._get_translation_vectors()

    @property
    def symbols(self) -> Symbols:
        """
        Returns a Symbols object.
        """
        return Symbols.fromsymbols(symbols=self.chemical_symbols)

    @property
    def molecule(self) -> Molecule:
        """
        Returns a molecule object.
        """
        return self.convert_coordinate_block_list_to_molecule()

    @property
    def constrained_atoms(self):
        """
        Returns a list of constraints in Gaussian format where 0 means unconstrained
        and -1 means constrained.
        """
        return self._get_constraints()

    def convert_coordinate_block_list_to_molecule(self):
        """
        Function to convert coordinate block supplied as text or as a list of lines into
        Molecule class.
        """
        return Molecule(
            symbols=self.symbols,
            positions=self.positions,
            frozen_atoms=self.constrained_atoms,
            pbc_conditions=self.pbc_conditions,
            translation_vectors=self.translation_vectors,
        )

    def _get_symbols(self):
        symbols = []
        for line in self.coordinate_block:
            line_elements = line.split()
            # assert len(line_elements) == 4, (
            # f'The geometry specification, `Symbol x y z` line should have 4 members \n'
            # f'but is {len(line_elements)} instead!')
            # not true for some cubes where the atomic number is repeated as a float:
            # 6    6.000000  -12.064399   -0.057172   -0.099010
            # also not true for Gaussian QM/MM calculations where "H" or "L" is
            # indicated at the end of the line

            if (
                len(line_elements) < 4 or len(line_elements) == 0
            ):  # skip lines that do not contain coordinates
                logger.debug(f"Line {line} has less than 4 line elements!")
                continue

            if (
                line_elements[0].upper() == "TV"
            ):  # cases where PBC system occurs in Gaussian
                logger.debug(f"Skipping line {line} with TV!")
                continue

            try:
                logger.debug(
                    f"Converting atomic number {line_elements[0]} to symbol."
                )
                atomic_number = int(
                    line_elements[0]
                )  # Could raise ValueError if not an integer
                chemical_symbol = p.to_symbol(
                    atomic_number=atomic_number
                )  # Could raise KeyError or similar
                logger.debug(
                    f"Successfully converted {line_elements[0]} to {chemical_symbol}."
                )
                symbols.append(chemical_symbol)
            except ValueError:
                # Handle case where line_elements[0] isn’t a valid integer
                logger.debug(
                    f"{line_elements[0]} is not a valid atomic number; treating as symbol."
                )
                try:
                    symbols.append(
                        p.to_element(element_str=str(line_elements[0]))
                    )
                except Exception as e:
                    logger.error(
                        f"Failed to convert {line_elements[0]} to element: {str(e)}"
                    )
            except Exception as e:
                # Catch any other unexpected errors
                logger.error(
                    f"Unexpected error processing {line_elements[0]}: {str(e)}"
                )
                try:
                    # Fallback attempt
                    symbols.append(
                        p.to_element(element_str=str(line_elements[0]))
                    )
                except Exception as fallback_e:
                    logger.error(
                        f"Fallback failed for {line_elements[0]}: {str(fallback_e)}"
                    )

        if len(symbols) == 0:
            raise ValueError(
                f"No symbols found in the coordinate block: {self.coordinate_block}!"
            )
        return symbols

    def _get_atomic_numbers_positions_and_constraints(self):
        atomic_numbers = []
        positions = []
        constraints = []
        for line in self.coordinate_block:
            if line.startswith(
                "TV"
            ):  # cases where PBC system occurs in Gaussian
                continue

            line_elements = line.split()
            if (
                len(line_elements) < 4 or len(line_elements) == 0
            ):  # skip lines that do not contain coordinates
                continue

            try:
                atomic_number = int(line_elements[0])
            except ValueError:
                atomic_number = p.to_atomic_number(
                    p.to_element(str(line_elements[0]))
                )
            atomic_numbers.append(atomic_number)

            second_value = float(line_elements[1])
            x_coordinate = 0.0
            y_coordinate = 0.0
            z_coordinate = 0.0
            if len(line_elements) > 4:
                if np.isclose(atomic_number, second_value, atol=10e-6):
                    # happens in cube file, where the second value is the same as
                    # the atomic number but in float format
                    x_coordinate = float(line_elements[2])
                    y_coordinate = float(line_elements[3])
                    z_coordinate = float(line_elements[4])
                elif np.isclose(second_value, -1, atol=10e-6) or np.isclose(
                    second_value, 0, atol=10e-6
                ):
                    # this is the case in frozen coordinates e.g.,
                    # C        -1      -0.5448210000   -1.1694570000    0.0001270000
                    # then ignore second value
                    constraints.append(int(second_value))
                    x_coordinate = float(line_elements[2])
                    y_coordinate = float(line_elements[3])
                    z_coordinate = float(line_elements[4])
            else:
                x_coordinate = float(line_elements[1])
                y_coordinate = float(line_elements[2])
                z_coordinate = float(line_elements[3])
            position = [x_coordinate, y_coordinate, z_coordinate]
            positions.append(position)
        if any(len(i) == 0 for i in [atomic_numbers, positions]):
            raise ValueError(
                f"No atomic numbers or positions found in the coordinate block: {self.coordinate_block}!"
            )
        return atomic_numbers, np.array(positions), constraints

    def _get_atomic_numbers(self):
        """
        Obtain a list of symbols as atomic numbers.
        """
        atomic_numbers, _, _ = (
            self._get_atomic_numbers_positions_and_constraints()
        )
        return atomic_numbers

    def _get_positions(self):
        """
        Obtain the coordinates of the molecule as numpy array.
        """
        _, positions, _ = self._get_atomic_numbers_positions_and_constraints()
        return positions

    def _get_constraints(self):
        _, _, constraints = (
            self._get_atomic_numbers_positions_and_constraints()
        )
        if len(constraints) == 0:
            return None
        if all(constraint == 0 for constraint in constraints):
            return None
        return constraints

    def _get_translation_vectors(self):
        tvs = []
        for line in self.coordinate_block:
            if line.startswith(
                "TV"
            ):  # cases where PBC system occurs in Gaussian
                line_elements = line.split()
                if len(line_elements) == 4:
                    x_coordinate = float(line_elements[1])
                    y_coordinate = float(line_elements[2])
                    z_coordinate = float(line_elements[3])
                else:
                    x_coordinate = float(line_elements[-3])
                    y_coordinate = float(line_elements[-2])
                    z_coordinate = float(line_elements[-1])
                tv = [x_coordinate, y_coordinate, z_coordinate]
                tvs.append(tv)
        if len(tvs) == 0:
            return None
        return tvs

    @property
    def pbc_conditions(self):
        """
        Obtain PBC conditions from given translation vectors.
        """
        if self.translation_vectors is not None:
            if len(self.translation_vectors) == 1:
                return [1, 0, 0]
            elif len(self.translation_vectors) == 2:
                return [1, 1, 0]
            elif len(self.translation_vectors) == 3:
                return [1, 1, 1]
        else:
            return None
