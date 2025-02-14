import logging

import numpy as np
from ase import units

from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.constants import R, atm_to_pa
from chemsmart.utils.logger import create_logger

logger = logging.getLogger(__name__)

create_logger()


class Thermochemistry:
    """Class for thermochemistry analysis. Use SI units.
    Requires filename from which thermochemistry data is extracted.
    Args:
        filename: str. Filepath to the file from which thermochemistry data is extracted.
        temperature: float. Temperature of the system, in K.
        pressure: float. Pressure of the system, in atm.
    """

    def __init__(self, filename, temperature, pressure):
        self.filename = filename
        self.molecule = Molecule.from_filepath(filename)
        self.temperature = temperature  # in Kelvin
        self.pressure = pressure * atm_to_pa

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
        logger.debug(f"Mass: {m} g/mol, Temperature: {T} K, Pressure: {P} Pa")
        return (2 * np.pi * m * units._k * T / (units._hplanck**2)) ** (
            3 / 2
        ) * (units._k * T / P)

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
        return R * (np.log(self.translational_partition_function) + 1 + d / 2)

    def translational_internal_energy(self):
        """Obtain the translational internal energy.
        Same for all types of molecules, whether linear, non-linear or monoatomic.
        Formula:
            U_t = 3/2 * R * T
        """
        return 3 / 2 * R * self.temperature

    def get_thermochemistry(self):
        pass


class GaussianThermochemistry(Thermochemistry):
    """Class for thermochemistry analysis of Gaussian output files."""

    def get_thermochemistry(self):
        pass


class OrcaThermochemistry(Thermochemistry):
    def get_thermochemistry(self):
        pass
