import os
import numpy as np
import ase
from ase.symbols import Symbols
from ase.atoms import Atoms
from chemsmart.utils.utils import file_cache
from chemsmart.utils.utils import FileReadError
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
    """

    def __init__(
        self,
        symbols=None,
        positions=None,
        charge=None,
        multiplicity=None,
        constraint=None,
        energy=None,
        forces=None,
        velocities=None,
        info=None,
    ):
        self.symbols = symbols
        self.positions = positions
        self.charge = charge
        self.multiplicity = multiplicity
        self.constraint = constraint
        self.energy = energy
        self.forces = forces
        self.velocities = velocities
        self.info = info

    @property
    def empirical_formula(self):
        return self.symbols.get_chemical_formula(mode="hill", empirical=True)

    def get_chemical_formula(self, mode, empirical):
        if self.symbols is not None:
            return self.symbols.get_chemical_formula(mode="hill", empirical=False)

    @classmethod
    def from_coordinate_block_text(cls, coordinate_block):
        from chemsmart.utils.utils import CoordinateBlock

        c = CoordinateBlock(coordinate_block=coordinate_block)
        return cls(symbols=c.symbols, positions=c.positions)

    @classmethod
    def from_file(cls, **kwargs):
        return cls.from_filepath(**kwargs)

    @classmethod
    def from_filepath(cls, filepath, index="-1", **kwargs):
        filepath = os.path.abspath(filepath)
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"{filepath} could not be found!")

        if os.path.getsize(filepath) == 0:
            return None

        try:
            molecule = cls._read_file(filepath, index=index, **kwargs)
        except Exception as e:
            raise FileReadError(f"Failed to create molecule from {filepath}.") from e
        return molecule

    @classmethod
    def _read_file(cls, filepath, index, **kwargs):
        basename = os.path.basename(filepath)

        if basename.endswith(".db"):
            return cls._read_database(filepath, index, **kwargs)

        if basename.endswith(".log"):
            return cls._read_gaussian_logfile(filepath, index, **kwargs)

        if basename.endswith((".com", ".gjf")):
            return cls._read_gaussian_comfile(filepath, **kwargs)

        if basename.endswith(".inp"):
            return cls._read_orca_inputfile(filepath, **kwargs)

        if basename.endswith(".out"):
            return cls._read_orca_outfile(filepath, index, **kwargs)

        if basename.endswith(".gro"):
            return cls._read_gromacs_gro(filepath, index, **kwargs)

        if basename.endswith(".trr"):
            return cls._read_gromacs_trr(filepath, index, **kwargs)

        if basename.endswith(".traj"):
            return cls._read_traj_file(filepath, index, **kwargs)

        return cls._read_other(filepath, index, **kwargs)

    @staticmethod
    @file_cache()
    def _read_other(filepath, index, **kwargs):
        return ase.io.read(filepath, index=index, **kwargs)

    @staticmethod
    @file_cache()
    def _read_sdf_file(filepath):
        # read the file content into a string
        with open(filepath) as f:
            content = f.read()
        # convert to atoms object
        from pyatoms.utils.utils import sdf2atoms

        return sdf2atoms(content)

    @staticmethod
    @file_cache()
    def _read_gaussian_logfile(filepath, index, **kwargs):
        try:
            from pyatoms.io.gaussian.outputs import Gaussian16Output

            g16_output = Gaussian16Output(logfile=filepath)
            return g16_output.get_atoms(index=index, include_failed_logfile=True)
        except ValueError:
            from pyatoms.io.gaussian.outputs import Gaussian16OutputWithPBC

            g16_output = Gaussian16OutputWithPBC(logfile=filepath)
            return g16_output.get_atoms(index=index, include_failed_logfile=True)

    @staticmethod
    @file_cache()
    def _read_gromacs_gro(filepath, index, **kwargs):
        # TODO: add handling of index
        from pyatoms.io.gromacs.outputs import GroGroOutput

        gro_output = GroGroOutput(filename=filepath)
        return gro_output.get_gro()

    @staticmethod
    @file_cache()
    def _read_gromacs_trr(filepath, index, **kwargs):
        from pyatoms.io.gromacs.outputs import GroTrrOutput

        trr_output = GroTrrOutput(filename=filepath)
        return trr_output.get_atoms(index=index)

    @staticmethod
    @file_cache()
    def _read_gaussian_comfile(filepath, **kwargs):  # noqa: PLR0915
        from ase.constraints import FixAtoms

        from pyatoms.utils.periodictable import PERIODIC_TABLE
        from pyatoms.utils.utils import is_float

        symbols = []
        positions = []

        pbc_conditions = []
        translation_vectors = []

        frozen_coordinates_dict = {"frozen": []}
        frozen_coordinates_list = []
        coordinates_frozen_status = []

        with open(filepath) as f:
            for line in f.readlines():
                line_elements = line.strip().split()

                # read charge and multiplicity from Gaussian .com file
                if len(line_elements) == 2 and all(i.isdigit for i in line_elements):
                    line_elements[0]
                    line_elements[1]

                # read symbols and positions
                if (
                    len(line_elements) == 4
                    and line_elements[0] in PERIODIC_TABLE
                    and is_float(line_elements[1])
                    and is_float(line_elements[2])
                    and is_float(line_elements[3])
                ):
                    symbols.append(str(line_elements[0]))
                    each_coord = [
                        float(line_elements[1]),
                        float(line_elements[2]),
                        float(line_elements[3]),
                    ]
                    positions.append(each_coord)
                elif (
                    len(line_elements) == 5
                    and line_elements[0] in PERIODIC_TABLE
                    and is_float(line_elements[2])
                    and is_float(line_elements[3])
                    and is_float(line_elements[4])
                ):
                    symbols.append(str(line_elements[0]))
                    each_coord = [
                        float(line_elements[2]),
                        float(line_elements[3]),
                        float(line_elements[4]),
                    ]
                    positions.append(each_coord)
                    coordinates_frozen_status.append(int(line_elements[1]))

                # to be able to read PBC files
                if len(line_elements) == 4 and line_elements[0].upper() == "TV":
                    pbc_conditions.append(1)
                    tv = [
                        float(line_elements[1]),
                        float(line_elements[2]),
                        float(line_elements[3]),
                    ]
                    translation_vectors.append(tv)
                elif len(line_elements) == 5 and line_elements[0].upper() == "TV":
                    pbc_conditions.append(1)
                    tv = [
                        float(line_elements[2]),
                        float(line_elements[3]),
                        float(line_elements[4]),
                    ]
                    translation_vectors.append(tv)
                    coordinates_frozen_status.append(int(line_elements[1]))

        if translation_vectors and any(pbc_conditions):
            if len(translation_vectors) == 1:
                translation_vectors.append([0.0, 0.0, 0.0])
                translation_vectors.append([0.0, 0.0, 0.0])
                pbc_conditions = [1, 0, 0]
            elif len(translation_vectors) == 2:
                translation_vectors.append([0.0, 0.0, 0.0])
                pbc_conditions = [1, 1, 0]
            elif len(translation_vectors) == 3:
                pbc_conditions = [1, 1, 1]

            cells = np.array(translation_vectors)
        else:
            pbc_conditions = [0, 0, 0]
            cells = None

        # take care of frozen coordinates
        for i, status in enumerate(coordinates_frozen_status):
            if status == -1:
                frozen_coordinates_list.append(i)
                frozen_coordinates_dict["frozen"].append(i)

        c = FixAtoms(indices=frozen_coordinates_list)
        atoms = Atoms(
            symbols=symbols, positions=positions, pbc=pbc_conditions, cell=cells
        )
        atoms.set_constraint(c)

        return atoms

    @staticmethod
    @file_cache()
    def _read_orca_outfile(filepath, index, **kwargs):
        from pyatoms.io.orca.outputs import ORCAOutput

        orca_output = ORCAOutput(filename=filepath, **kwargs)
        return orca_output.get_atoms(index=index, include_failed_file=True)

    @staticmethod
    @file_cache()
    def _read_orca_inputfile(filepath, **kwargs):
        from pyatoms.io.orca.inputs import ORCAInput

        orca_input = ORCAInput(inpfile=filepath, **kwargs)
        return orca_input.atoms

    def write(self, f):
        assert self.symbols is not None, "Symbols to write should not be None!"
        assert self.positions is not None, "Positions to write should not be None!"

        pass


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
    def symbols(self) -> Symbols:
        """Returns a Symbols object."""
        return Symbols.fromsymbols(symbols=self.chemical_symbols)

    @property
    def molecule(self) -> Molecule:
        """Returns a molecule object."""
        return self.convert_coordinate_block_list_to_molecule()

    def convert_coordinate_block_list_to_molecule(self):
        """Function to convert coordinate block supplied as text or as a list of lines into
        Molecule class."""
        symbols = self._get_symbols(self.coordinate_block)
        positions = self._get_positions(self.coordinate_block)
        return Molecule(symbols=symbols, positions=positions)

    def _get_symbols(self):
        symbols = []
        for line in self.coordinate_block:
            line_elements = line.split()
            # assert len(line_elements) == 4, (
            # f'The geometry specification, `Symbol x y z` line should have 4 members \n'
            # f'but is {len(line_elements)} instead!')
            # not true for some cubes where the atomic number is repeated as a float:
            # 6    6.000000  -12.064399   -0.057172   -0.099010
            # also not true for Gaussian QM/MM calculations where "H" or "L" is indicated at the end of the line

            try:
                atomic_number = int(line_elements[0])
                chemical_symbol = p.to_symbol(atomic_number=atomic_number)
                symbols.append(chemical_symbol)
            except ValueError:
                symbols.append(p.to_element(element_str=str(line_elements[0])))
        return symbols

    def _get_positions(self):
        positions = []
        for line in self.coordinate_block:
            line_elements = line.split()

            try:
                atomic_number = int(line_elements[0])
            except ValueError:
                atomic_number = p.to_atomic_number(p.to_element(str(line_elements[0])))

            second_value = float(line_elements[1])
            if np.isclose(atomic_number, second_value, atol=10e-6):
                # happens in cube file, where the second value is the same as the atomic number but in float format
                x_coordinate = float(line_elements[2])
                y_coordinate = float(line_elements[3])
                z_coordinate = float(line_elements[4])
            elif second_value == -1 or second_value == 0:
                # this is the case in frozen coordinates e.g.,
                # C        -1      -0.5448210000   -1.1694570000    0.0001270000
                # then ignore second value
                x_coordinate = float(line_elements[2])
                y_coordinate = float(line_elements[3])
                z_coordinate = float(line_elements[4])
            else:
                x_coordinate = float(line_elements[1])
                y_coordinate = float(line_elements[2])
                z_coordinate = float(line_elements[3])
            position = [x_coordinate, y_coordinate, z_coordinate]
            positions.append(position)
        return np.array(positions)
