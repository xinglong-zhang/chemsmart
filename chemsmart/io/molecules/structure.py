import copy
import logging
import os
import re
from functools import cached_property, lru_cache

import ase
import networkx as nx
import numpy as np
from ase.symbols import Symbols
from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Geometry import Point3D
from scipy.spatial.distance import cdist

from chemsmart.io.molecules import get_bond_cutoff
from chemsmart.utils.mixins import FileMixin
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
    frozen_atoms: list of integers to freeze atoms in the molecule.
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
        info=None,
    ):
        self.symbols = symbols
        self.positions = positions
        self.charge = charge
        self.multiplicity = multiplicity
        self.frozen_atoms = frozen_atoms
        self.pbc_conditions = pbc_conditions
        self.translation_vectors = translation_vectors
        self.energy = energy
        self.forces = forces
        self.velocities = velocities
        self.info = info

        # Define bond order classification multipliers (avoiding redundancy)
        # use the relationship between bond orders and bond lengths from J. Phys. Chem. 1959, 63, 8, 1346
        #                 1/Ls^3 : 1/Ld^3 : 1/Lt^3 = 2 : 3 : 4
        # From experimental data, the aromatic bond length in benzene is ~1.39 Å,
        # which is between single (C–C, ~1.54 Å) and double (C=C, ~1.34 Å) bonds.
        # Interpolate the aromatic bond length as L_ar ≈ L_s × (4/5)**1/3
        self.bond_length_multipliers = {
            "single": 1.0,
            "aromatic": (4 / 5) ** (1.0 / 3),  # Aromatic bond approx.
            "double": (2 / 3) ** (1.0 / 3),
            "triple": (1 / 2) ** (1.0 / 3),
        }

        # Ensure symbols and positions are available
        if self.symbols is None or self.positions is None:
            raise ValueError(
                "Molecule must have symbols and positions defined."
            )

        # check that the number of symbols are not empty
        if len(self.symbols) == 0:
            raise ValueError("The number of symbols should not be empty!")

        # check that the number of symbols and positions are the same
        if len(self.symbols) != len(self.positions):
            raise ValueError(
                "The number of symbols and positions should be the same!"
            )

    def __len__(self):
        return len(self.chemical_symbols)

    def __getitem__(self, idx):
        """Interprets the input idx as 1-based indices by index Adjustment (i - 1).
        The method assumes that idx contains 1-based indices (e.g., [1, 2, 3]),
        so it subtracts 1 to convert them to Python's zero-based indexing.
        Retrieves the corresponding elements from self.symbols and self.positions using the provided indices.
        reates and returns a new instance of the same class, containing the selected symbols and positions.
        """
        symbols = [self.symbols[i - 1] for i in idx]
        positions = [self.positions[i - 1] for i in idx]
        return type(self)(symbols=symbols, positions=positions)

    @property
    def empirical_formula(self):
        return Symbols.fromsymbols(self.symbols).get_chemical_formula(
            mode="hill", empirical=True
        )

    @property
    def chemical_formula(self):
        return self.get_chemical_formula()

    @property
    def chemical_symbols(self):
        """Return a list of chemical symbols strings"""
        if self.symbols is not None:
            return list(self.symbols)

    @property
    def num_atoms(self):
        """Return the number of atoms in the molecule."""
        return len(self.chemical_symbols)

    @property
    def pbc(self):
        """Return the periodic boundary conditions."""
        return all(i == 0 for i in self.pbc_conditions)

    @property
    def is_chiral(self):
        """Check if molecule is chiral or not."""
        return Chem.FindMolChiralCenters(self.to_rdkit(), force=True) != []

    def get_chemical_formula(self, mode="hill", empirical=False):
        if self.symbols is not None:
            return Symbols.fromsymbols(self.symbols).get_chemical_formula(
                mode=mode, empirical=empirical
            )

    def copy(self):
        return copy.deepcopy(self)

    @classmethod
    def from_coordinate_block_text(cls, coordinate_block):
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
        return cls(
            symbols=Symbols.fromsymbols(list_of_symbols),
            positions=positions,
            pbc_conditions=pbc_conditions,
        )

    @classmethod
    def from_filepath(cls, filepath, index="-1", return_list=False, **kwargs):
        filepath = os.path.abspath(filepath)
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"{filepath} could not be found!")

        if os.path.getsize(filepath) == 0:
            return None

        molecule = cls._read_filepath(
            filepath, index=index, return_list=return_list, **kwargs
        )
        if return_list and isinstance(molecule, Molecule):
            return [molecule]
        else:
            return molecule

    @classmethod
    def _read_filepath(cls, filepath, index, return_list, **kwargs):
        basename = os.path.basename(filepath)

        if basename.endswith(".xyz"):
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
            return cls._read_gaussian_logfile(filepath, index)

        if basename.endswith(".inp"):
            return cls._read_orca_inputfile(filepath, **kwargs)

        if basename.endswith(".out"):
            return cls._read_orca_outfile(filepath, index, **kwargs)

        if basename.endswith(".gro"):
            return cls._read_gromacs_gro(filepath, index, **kwargs)

        if basename.endswith(".trr"):
            return cls._read_gromacs_trr(filepath, index, **kwargs)

        # if basename.endswith(".traj"):
        #     return cls._read_traj_file(filepath, index, **kwargs)

        return cls._read_other(filepath, index, **kwargs)

    @classmethod
    def _read_xyz_file(cls, filepath, index=":", return_list=False):
        xyz_file = XYZFile(filename=filepath)
        molecules = xyz_file.get_molecule(index=index, return_list=return_list)
        return molecules

    @staticmethod
    @file_cache()
    def _read_sdf_file(filepath):
        sdf_file = SDFFile(filepath)
        return sdf_file.molecule

    @staticmethod
    @file_cache()
    def _read_gaussian_inputfile(filepath):
        from chemsmart.io.gaussian.input import Gaussian16Input

        try:

            g16_input = Gaussian16Input(filename=filepath)
            return g16_input.molecule
        except ValueError:
            g16_input = Gaussian16Input(filename=filepath)

    @staticmethod
    @file_cache()
    def _read_gaussian_logfile(filepath, index):
        """Returns a list of molecules."""
        from chemsmart.io.gaussian.output import Gaussian16Output

        g16_output = Gaussian16Output(filename=filepath)
        return g16_output.get_molecule(index=index)

    @staticmethod
    @file_cache()
    def _read_orca_inputfile(filepath):
        from chemsmart.io.orca.input import ORCAInput

        orca_input = ORCAInput(filename=filepath)
        return orca_input.molecule

    @staticmethod
    @file_cache()
    def _read_orca_outfile(filepath, index):
        # TODO: to improve ORCAOutput object so that all the structures can be obtained and returned via index
        from chemsmart.io.orca.output import ORCAOutput

        orca_output = ORCAOutput(filename=filepath)
        return orca_output.molecule

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

    @staticmethod
    @file_cache()
    def _read_other(filepath, index, **kwargs):
        return ase.io.read(filepath, index=index, **kwargs)

    @classmethod
    @lru_cache(maxsize=128)
    def from_pubchem(cls, identifier, return_list=False):
        """Creates Molecule object from pubchem based on an identifier (CID, SMILES, or name).
        Args:
        identifier (str): The compound identifier (name, CID, or SMILES string).
        output_format (str): The desired format of the response. Default is "json".
                             Other options include "sdf" or "xml".
        Raises:
            requests.exceptions.RequestException: For network or HTTP-related issues.
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
        return cls(**molecule.__dict__)

    @classmethod
    def from_ase_atoms(cls, atoms):
        return cls(
            symbols=atoms.get_chemical_symbols(),
            positions=atoms.get_positions(),
            pbc_conditions=atoms.get_pbc(),
        )

    def write_coordinates(self, f, program=None):
        """Write the coordinates of the molecule to a file.
        No empty end line at the end of the file."""
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

    def _write_gaussian_coordinates(self, f):
        assert self.symbols is not None, "Symbols to write should not be None!"
        assert (
            self.positions is not None
        ), "Positions to write should not be None!"
        if self.frozen_atoms is None:
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
        """Write the coordinates of the molecule with PBC conditions to a file."""
        if self.pbc_conditions is not None:
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
        # ORCA cannot do PBC calculations
        pass

    def __repr__(self):
        return f"{self.__class__.__name__}<{self.empirical_formula}>"

    def __str__(self):
        return f"{self.__class__.__name__}<{self.empirical_formula}>"

    @cached_property
    def distance_matrix(self):
        """ "Compute pairwise distance matrix."""
        return cdist(self.positions, self.positions)

    def determine_bond_order(self, bond_length, bond_cutoff):
        """Determine the bond order based on bond length and cutoff."""
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

    def bond_lengths(self):
        # get all bond distances in the molecule
        return self.get_all_distances()

    def get_all_distances(self):
        bond_distances = []
        for i in range(self.num_atoms):
            for j in range(i + 1, self.num_atoms):
                bond_distances.append(
                    np.linalg.norm(self.positions[i] - self.positions[j])
                )
        return bond_distances

    def to_smiles(self):
        """Convert molecule to SMILES string."""
        # Create an RDKit molecule
        rdkit_mol = self.to_rdkit()

        # Convert RDKit molecule to SMILES
        return Chem.MolToSmiles(rdkit_mol)

    def to_rdkit_old(self):
        """Convert Molecule object to RDKit Mol."""

        # Ensure symbols and positions are available
        if self.symbols is None or self.positions is None:
            raise ValueError(
                "Molecule must have symbols and positions defined."
            )

        # Create an empty RDKit molecule
        rdkit_mol = Chem.RWMol()

        # Add atoms and store their indices
        atom_indices = []
        for symbol in self.symbols:
            atom = Chem.Atom(symbol)
            idx = rdkit_mol.AddAtom(atom)
            atom_indices.append(idx)

        # Create a conformer for the molecule and add 3D coordinates
        conformer = rdchem.Conformer(len(atom_indices))
        for idx, (x, y, z) in enumerate(self.positions):
            conformer.SetAtomPosition(idx, Point3D(x, y, z))

        # Add the conformer to the molecule
        rdkit_mol.AddConformer(conformer)

        # I comment the following out since we do not want to modify the molecule
        # Validate the RDKit molecule
        # try:
        #     Chem.SanitizeMol(rdkit_mol)
        # except Chem.AtomValenceException as e:
        #     raise ValueError(f"Sanitization failed: {e}") from e

        # # **Explicitly compute implicit valences** (Fix for getNumImplicitHs issue)
        # rdkit_mol.UpdatePropertyCache(strict=False)

        # # Detect chirality from 3D geometry
        # Chem.FindMolChiralCenters(rdkit_mol)
        #
        # # Assign bond stereochemistry (important for E/Z isomers)
        # Chem.AssignStereochemistry(rdkit_mol, force=True, cleanIt=True)

        Chem.AssignStereochemistryFrom3D(rdkit_mol)

        return rdkit_mol

    def to_rdkit(self, bond_cutoff_buffer=0.05, adjust_H=True):
        """Convert Molecule object to RDKit Mol with proper stereochemistry handling.
        Args:
            bond_cutoff_buffer (float): Additional buffer for bond cutoff distance.
            From testing, see test_resonance_handling, it seems that a value of 0.1Å
            works for ozone, acetone, benzene, and probably other molecules, too.
            adjust_Hs (bool): Adjust bond distances to H atoms.
        Returns:
            RDKit Mol: RDKit molecule object.
        """

        # Create molecule and add atoms
        rdkit_mol = Chem.RWMol()

        # Add atoms to the RDKit molecule
        for symbol in self.symbols:
            rdkit_mol.AddAtom(Chem.Atom(symbol))

        # add bonds
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
                print(f"bond cutoff buffer: {cutoff_buffer}")
                cutoff = get_bond_cutoff(
                    self.symbols[i], self.symbols[j], cutoff_buffer
                )
                print(
                    f"bond cutoff: {cutoff} between {self.symbols[i]} and {self.symbols[j]}"
                )
                bond_order = self.determine_bond_order(
                    bond_length=self.distance_matrix[i, j], bond_cutoff=cutoff
                )
                print(
                    f"bond order: {bond_order} between {self.symbols[i]} and {self.symbols[j]}"
                )
                if bond_order > 0:
                    bond_type = {
                        1: Chem.BondType.SINGLE,
                        1.5: Chem.BondType.AROMATIC,
                        2: Chem.BondType.DOUBLE,
                        3: Chem.BondType.TRIPLE,
                    }[bond_order]
                    rdkit_mol.AddBond(i, j, bond_type)
            # # reset bond cutoff buffer
            # bond_cutoff_buffer = bond_cutoff_buffer

        # Create a conformer and set 3D coordinates
        conformer = rdchem.Conformer(len(self.symbols))
        for i, pos in enumerate(self.positions):
            conformer.SetAtomPosition(i, Point3D(*pos))
        rdkit_mol.AddConformer(conformer)

        # Compute valences (Fix for getNumImplicitHs issue)
        rdkit_mol.UpdatePropertyCache(strict=False)

        # Partial sanitization for stereochemistry detection
        # Chem.SanitizeMol(rdkit_mol,
        #                  Chem.SANITIZE_ALL ^ Chem.SANITIZE_ADJUSTHS ^ Chem.SANITIZE_SETAROMATICITY)

        # Detect stereochemistry from 3D coordinates
        Chem.AssignStereochemistryFrom3D(rdkit_mol, conformer.GetId())
        Chem.AssignAtomChiralTagsFromStructure(rdkit_mol, conformer.GetId())

        # # **Explicitly set chiral tags**
        # for atom in rdkit_mol.GetAtoms():
        #     if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
        #         atom.SetProp("_CIPCode", atom.GetProp("_CIPCode"))  # Preserve RDKit's chirality perception

        # Force update of stereo flags
        Chem.FindPotentialStereoBonds(rdkit_mol, cleanIt=True)

        return rdkit_mol.GetMol()

    @cached_property
    def bond_orders(self):
        """Return a list of bond orders from the molecular graph."""

        # return [
        #     bond.GetBondTypeAsDouble() for bond in self.to_rdkit().GetBonds()
        # ]
        graph = self.to_graph()
        bond_orders = []
        for bond in graph.edges.values():
            bond_orders.append(bond["bond_order"])
        return bond_orders

    def to_graph(self, bond_cutoff_buffer=0.05, adjust_H=True) -> nx.Graph:
        """Convert a Molecule object to a connectivity graph.
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
                bond_order = self.determine_bond_order(
                    bond_length=self.distance_matrix[i, j], bond_cutoff=cutoff
                )
                if bond_order > 0:
                    G.add_edge(i, j, bond_order=bond_order)
        return G


class CoordinateBlock:
    """Class to create coordinate block object to abstract the geometry."""

    def __init__(self, coordinate_block):
        """Accepts a coordinate block either as text string or as a list of lines.
        If former, then convert to the latter before future usage."""
        coordinate_block_list = []
        if isinstance(coordinate_block, str):
            for line in coordinate_block.split("\n"):
                coordinate_block_list.append(line.strip())
        elif isinstance(coordinate_block, list):
            coordinate_block_list = coordinate_block
        else:
            raise TypeError(
                f"The given coordinate block should be str or list "
                f"but is {type(coordinate_block)} instead!"
            )
        self.coordinate_block = coordinate_block_list

    @property
    def chemical_symbols(self):
        """Returns a list of chemical symbols for the molecule."""
        return self._get_symbols()

    @property
    def positions(self):
        """Returns a list of positions for the molecule."""
        return self._get_positions()

    @property
    def translation_vectors(self):
        """Return a list of translation vectors for systems with pbc."""
        return self._get_translation_vectors()

    @property
    def symbols(self) -> Symbols:
        """Returns a Symbols object."""
        return Symbols.fromsymbols(symbols=self.chemical_symbols)

    @property
    def molecule(self) -> Molecule:
        """Returns a molecule object."""
        return self.convert_coordinate_block_list_to_molecule()

    @property
    def constrained_atoms(self):
        """Returns a list of contraints in Gaussian format where 0 means unconstrained
        and -1 means constrained."""
        return self._get_constraints()

    def convert_coordinate_block_list_to_molecule(self):
        """Function to convert coordinate block supplied as text or as a list of lines into
        Molecule class."""
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
                continue

            if (
                line_elements[0].upper() == "TV"
            ):  # cases where PBC system occurs in Gaussian
                continue

            try:
                atomic_number = int(line_elements[0])
                chemical_symbol = p.to_symbol(atomic_number=atomic_number)
                symbols.append(chemical_symbol)
            except ValueError:
                symbols.append(p.to_element(element_str=str(line_elements[0])))
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
        """Obtain a list of symbols as atomic numbers."""
        atomic_numbers, _, _ = (
            self._get_atomic_numbers_positions_and_constraints()
        )
        return atomic_numbers

    def _get_positions(self):
        """Obtain the coordinates of the molecule as numpy array."""
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
        """Obtain PBC conditions from given translation vectors."""
        if self.translation_vectors is not None:
            if len(self.translation_vectors) == 1:
                return [1, 0, 0]
            elif len(self.translation_vectors) == 2:
                return [1, 1, 0]
            elif len(self.translation_vectors) == 3:
                return [1, 1, 1]
        else:
            return None


class SDFFile(FileMixin):
    """SDF file object."""

    def __init__(self, filename):
        self.filename = filename

    @property
    def molecule(self):
        return self.get_molecule()

    def get_molecule(self):
        list_of_symbols = []
        cart_coords = []
        # sdf line pattern containing coordinates and element type
        from chemsmart.utils.repattern import sdf_pattern

        for line in self.contents:
            match = re.match(sdf_pattern, line)
            if match:
                x = float(match.group(1))
                y = float(match.group(2))
                z = float(match.group(3))
                atom_type = str(match.group(4))
                list_of_symbols.append(atom_type)
                cart_coords.append((x, y, z))

        cart_coords = np.array(cart_coords)

        if len(list_of_symbols) == 0 or len(cart_coords) == 0:
            raise ValueError("No coordinates found in the SDF file!")

        return Molecule.from_symbols_and_positions_and_pbc_conditions(
            list_of_symbols=list_of_symbols, positions=cart_coords
        )


class XYZFile(FileMixin):
    """xyz file object."""

    def __init__(self, filename):
        self.filename = filename

    @cached_property
    def num_atoms(self):
        return int(self.contents[0])

    def _get_molecules_and_comments(self, index=":", return_list=False):
        """Return a molecule object or a list of molecule objects from an xyz file.
        The xzy file can either contain a single molecule, as conventionally, or a list
        of molecules, such as those in crest_conformers.xyz file."""
        all_molecules = []
        comments = []
        i = 0
        while i < len(self.contents):
            # Read number of atoms
            num_atoms = int(self.contents[i].strip())
            if num_atoms == 0:
                raise ValueError("Number of atoms in the xyz file is zero!")
            i += 1
            # Read comment line
            comment = self.contents[i].strip()
            comments.append(comment)
            i += 1
            # Read the coordinate block
            coordinate_block = self.contents[i : i + num_atoms]
            i += num_atoms
            molecule = Molecule.from_coordinate_block_text(coordinate_block)

            # Store the molecule data
            all_molecules.append(molecule)

        molecules = all_molecules[string2index_1based(index)]
        comments = comments[string2index_1based(index)]
        if return_list and isinstance(molecules, Molecule):
            return [molecules], [comments]
        return molecules, comments

    def get_molecule(self, index=":", return_list=False):
        molecules, _ = self._get_molecules_and_comments(
            index=index, return_list=return_list
        )
        return molecules

    def get_comments(self, index=":", return_list=False):
        _, comments = self._get_molecules_and_comments(
            index=index, return_list=return_list
        )
        return comments
