import numpy as np
from ase import units

from chemsmart.io.molecules.structure import Molecule


class Thermochemistry:
    """Class for thermochemistry analysis.
    Requires filename from which thermochemistry data is extracted."""

    def __init__(self, filename):
        self.filename = filename
        self.molecule = Molecule.from_filepath(filename)

    def translational_entropy(self):
        """Obtain the translational entropy.
        Formula:
            S_t = R * ln(q_t + 1 + 3/2 ) for non-linear molecules.
            S_t = R * ln(q_t + 1 + 1) for linear molecules.
            S_t = R * ln(q_t + 1) for monoatomic molecules.
        """
        if self.molecule.is_monoatomic:
            return units.mol * units.kB * np.log(3 + 1)
        elif self.molecule.is_linear:
            return units.mol * units.kB * np.log(3 + 1 + 1)
        else:
            return units.mol * units.kB * np.log(3 + 1 + 3 / 2)

    def get_thermochemistry(self):
        pass


class GaussianThermochemistry(Thermochemistry):
    """Class for thermochemistry analysis of Gaussian output files."""

    def get_thermochemistry(self):
        pass


class OrcaThermochemistry(Thermochemistry):
    def get_thermochemistry(self):
        pass
