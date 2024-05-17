from ase.symbols import Symbols, symbols2numbers

class Molecule:
    """ Class to represent a molcular structure.
    Follows from :class:`ase.symbols.Symbols` str (formula) or list of str
    Can be a string formula, a list of symbols or a list of Atom objects.
    Examples: 'H2O', 'COPt12', ['H', 'H', 'O'], [Atom('Ne', (x, y, z)), ...]. """
    def __init__(
            self, symbols=None, positions=None, charge=None, multiplicity=None, constraint=None, energy=None,
            forces=None, velocities=None, info=None
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

    def write(self, f):
        pass

