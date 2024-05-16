class Molecule:
    """ Class to represent a molcular structure."""
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

