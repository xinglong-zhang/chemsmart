import numpy as np
from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.periodictable import PeriodicTable as pt

p = pt()


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
    def molecule(self):
        """Returns a molecule object."""
        return self.convert_coordinate_block_to_molecule()

    def convert_coordinate_block_list_to_molecule(self):
        """Function to convert coordinate block supplied as text or as a list of lines into
        Molecule class."""
        symbols = self._get_symbols(self.coordinate_block)
        positions = self._get_positions(self.coordinate_block)
        molecule = Molecule(symbols=symbols, positions=positions)
        return molecule

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
                atomic_number = p.to_atomic_number(
                    p.to_element(str(line_elements[0]))
                )

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
