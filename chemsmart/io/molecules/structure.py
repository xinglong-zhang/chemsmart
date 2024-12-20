import os
import re
import ase
import numpy as np
from ase.symbols import Symbols
from ase.io.formats import string2index
from functools import cached_property
from chemsmart.utils.utils import file_cache
from chemsmart.utils.utils import FileReadError
from chemsmart.utils.mixins import FileMixin
from chemsmart.utils.periodictable import PeriodicTable as pt

p = pt()


class Molecule:
    """Class to represent a molcular structure.

    Parameters:

    symbols: follows from :class:`ase.symbols.Symbols`
        str (formula) or list of str
        Can be a string formula, a list of symbols or a list of
        Atom objects.  Examples: 'H2O', 'COPt12', ['H', 'H', 'O'],
        [Atom('Ne', (x, y, z)), ...].
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

    @property
    def empirical_formula(self):
        return Symbols.fromsymbols(self.symbols).get_chemical_formula(
            mode="hill", empirical=True
        )

    @property
    def chemical_symbols(self):
        """Return a list of chemical symbols strings"""
        if self.symbols is not None:
            return list(self.symbols)

    @property
    def num_atoms(self):
        """Return the number of atoms in the molecule."""
        return len(self.chemical_symbols)

    def get_chemical_formula(self, mode="hill", empirical=False):
        if self.symbols is not None:
            return self.symbols.get_chemical_formula(
                mode=mode, empirical=empirical
            )

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

        try:
            molecule = cls._read_filepath(
                filepath, index=index, return_list=return_list, **kwargs
            )
            return molecule
        except Exception as e:
            raise FileReadError(
                f"Failed to create molecule from {filepath}."
            ) from e

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
        from chemsmart.io.gaussian.inputs import Gaussian16Input

        g16_input = Gaussian16Input(filename=filepath)
        return g16_input.molecule

    @staticmethod
    @file_cache()
    def _read_gaussian_logfile(filepath, index):
        from chemsmart.io.gaussian.output import Gaussian16Output

        g16_output = Gaussian16Output(filename=filepath)
        return g16_output.get_molecule(index=index)

    @staticmethod
    @file_cache()
    def _read_orca_inputfile(filepath):
        from chemsmart.io.orca.inputs import ORCAInput

        orca_input = ORCAInput(filename=filepath)
        return orca_input.molecule

    @staticmethod
    @file_cache()
    def _read_orca_outfile(filepath, index):
        # TODO: to improve ORCAOutput object so that all the structures can be obtained and returned via index
        from chemsmart.io.orca.outputs import ORCAOutput

        orca_output = ORCAOutput(filename=filepath)
        return orca_output.molecule

    # @staticmethod
    # @file_cache()
    # def _read_gromacs_gro(filepath, index, **kwargs):
    #     # TODO: add handling of index
    #     from pyatoms.io.gromacs.outputs import GroGroOutput
    #
    #     gro_output = GroGroOutput(filename=filepath)
    #     return gro_output.get_gro()
    #
    # @staticmethod
    # @file_cache()
    # def _read_gromacs_trr(filepath, index, **kwargs):
    #     from pyatoms.io.gromacs.outputs import GroTrrOutput
    #
    #     trr_output = GroTrrOutput(filename=filepath)
    #     return trr_output.get_atoms(index=index)

    @staticmethod
    @file_cache()
    def _read_other(filepath, index, **kwargs):
        return ase.io.read(filepath, index=index, **kwargs)

    def write(self, f):
        assert self.symbols is not None, "Symbols to write should not be None!"
        assert (
            self.positions is not None
        ), "Positions to write should not be None!"

        pass

    def __repr__(self):
        return f"{self.__class__.__name__}<{self.empirical_formula}>"

    def __str__(self):
        return f"{self.__class__.__name__}<{self.empirical_formula}>"

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

    def get_all_angles(self):
        # get all bond angles in the molecule
        bond_angles = []
        for i in range(self.num_atoms):
            for j in range(i + 1, self.num_atoms):
                for k in range(j + 1, self.num_atoms):
                    bond_angles.append(
                        self.get_angle(i, j, k)
                    )   # get the angle between atoms i, j, k


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

        molecules = all_molecules[string2index(index)]
        comments = comments[string2index(index)]
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
