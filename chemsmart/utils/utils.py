import numpy as np
from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.periodictable import PeriodicTable as p
periodic_table = p()
def convert_coordinate_block_to_molecule(coordinate_block):
    """ Function to convert coordinate block supplied as text or as a list of lines into
    Molecule class."""
    if isinstance(coordinate_block, str):
        molecule = _convert_coordinate_block_text_to_molecule(coordinate_block)
    elif isinstance(coordinate_block, list):
        molecule = _convert_coordinate_block_list_to_molecule(coordinate_block)
    else:
        raise TypeError(f'The given coordinate block should be str or list but is {type(coordinate_block)} instead!')
    return molecule

def _convert_coordinate_block_text_to_molecule(coordinate_block):
    lines = []
    for line in coordinate_block.split('\n'):
        lines.append(line.strip())
    moleclue = _convert_coordinate_block_list_to_molecule(lines)
    return moleclue

def _convert_coordinate_block_list_to_molecule(coordinate_block):
    symbols = _get_symbols(coordinate_block)
    positions = _get_positions(coordinate_block)
    molecule = Molecule(symbols=symbols, positions=positions)
    return molecule


def _get_symbols(coordinate_block_list):
    symbols = []
    for line in coordinate_block:
        line_elements = line.split()
        # assert len(line_elements) == 4, (f'The geometry specification, `Symbol x y z` line should have 4 members \n'
        #                                  f'but is {len(line_elements)} instead!')
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

def _get_positions(coordinate_block_list):
    positions = []
    for line in coordinate_block_list:
        line_elements = line.split()

        try:
            atomic_number = int(line_elements[0])
        except ValueError:
            atomic_number = p.to_atomic_number(p.to_element(str(line_elements[0])))

        second_value = int(line_elements[1])
        if atomic_number == second_value:
            x_coordinate = float(line_elements[2])
            y_coordinate = float(line_elements[3])
            z_coordinate = float(line_elements[4])
        elif second_value == -1 or second_value == 0:
            # what happens in frozen coordinates e.g.,
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








