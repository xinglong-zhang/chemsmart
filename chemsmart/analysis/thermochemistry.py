import numpy as np
from chemsmart.utils.units import ase_units as units

from chemsmart.io.molecules.structure import Molecule


class Thermochemistry:
    """Class for thermochemistry analysis.
    Requires filename from which thermochemistry data is extracted.
    Args:
        filename: str. Filepath to the file from which thermochemistry data is extracted.
        temperature: float. Temperature of the system, in ºC.
        pressure: float. Pressure of the system, in atm.
    """

    def __init__(self, filename, temperature, pressure):
        self.filename = filename
        self.molecule = Molecule.from_filepath(filename)
        self.temperature = temperature + 273.15  # Convert to Kelvin
        self.pressure = pressure * units.atm_to_pa

    @property
    def translational_partition_function(self):
        """Obtain the translational partition function.
        Formula:
            q_t = (2 * pi * m * k_B * T / h^2)^(3/2) * (k_B * T / P)
        where:
            m = mass of the molecule
            k = Boltzmann constant
            T = temperature
            h = Planck constant
            P = pressure of the system
        """
        m = self.molecule.mass
        T = self.temperature
        P = self.pressure
        return ((2 * np.pi * m * units._k * T) / (units.hplanck) ** 2) ** (3 / 2) * (
            units._k * T / P
        )

    def translational_entropy(self):
        """Obtain the translational entropy.
        Formula:
            S_t = R * [ln(q_t) + 1 + d/2]
            where d = 3 for non-linear molecules;
                  d = 2 for linear molecules;
                  d = 1 for monoatomic molecules.
        """
        if self.molecule.is_monoatomic:
            d = 1
        elif self.molecule.is_linear:
            d = 2
        else:
            d = 3
        return units.kB * units._Nav * np.log(3 + 1 + d / 2)

    def translational_internal_energy(self):
        """Obtain the translational internal energy.
        Same for all types of molecules, whether linear, non-linear or monoatomic.
        Formula:
            U_t = 3/2 * R * T
        """
        return 3 / 2 * units.kB * units._Nav * self.temperature

    def get_thermochemistry(self):
        pass


class GaussianThermochemistry(Thermochemistry):
    """Class for thermochemistry analysis of Gaussian output files."""

    def get_thermochemistry(self):
        pass


class OrcaThermochemistry(Thermochemistry):
    def get_thermochemistry(self):
        pass
