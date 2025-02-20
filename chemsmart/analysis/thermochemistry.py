import logging
import math

import numpy as np
from ase import units

from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.orca.output import ORCAOutput
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

    def __init__(
        self,
        filename,
        temperature,
        pressure,
    ):
        self.filename = filename
        self.molecule = Molecule.from_filepath(filename)
        self.temperature = temperature
        self.pressure = pressure
        self.m = (
            self.molecule.mass
            * units._amu  # converts mass from g/mol to kg/molecule
            # units._amu is same as divide by Avogadro's number then by 1000 (g to kg)
        )  # convert the unit of mass of the molecule from amu to kg
        self.T = self.temperature  # temperature in K
        self.P = (
            self.pressure * atm_to_pa
        )  # convert the unit of pressure from atm to Pascal
        if self.moments_of_inertia is not None:
            self.I = self.moments_of_inertia * (
                units._amu * (units.Bohr / units.m) ** 2
            )  # convert the unit of moments of inertia from amu Å^2 to kg m^2
        if self.vibrational_frequencies is not None:
            self.v = [
                k * units._c * 1e2 for k in self.vibrational_frequencies
            ]  # convert the unit of vibrational frequencies from cm^-1 to Hz

            # Calculate the characteristic vibrational temperature, theta, for each vibrational mode
            self.theta = [units._hplanck * vk / units._k for vk in self.v]

    @property
    def file_object(self):
        """Open the file and return the file object."""
        if str(self.filename).endswith(".log"):
            # create a Gaussian16Output object if .log file
            return Gaussian16Output(self.filename)
        elif str(self.filename).endswith(".out"):
            # create an OrcaOutput object if .out file
            return ORCAOutput(self.filename)
        else:
            # can be added in future to parse other file formats
            raise ValueError(
                "Unsupported file format. Use .log or .out files."
            )

    @property
    def moments_of_inertia(self):
        """Obtain the moments of inertia of the molecule.
        Array of moments of inertia along the x, y, and z axes"""
        return self.file_object.moments_of_inertia

    @property
    def rotational_symmetry_number(self):
        """Obtain the rotational symmetry number."""
        return self.file_object.rotational_symmetry_number

    @property
    def vibrational_frequencies(self):
        """Obtain the vibrational frequencies of the molecule."""
        return self.file_object.vibrational_frequencies

    @property
    def translational_partition_function(self):
        """Obtain the translational partition function.
        Formula:
            q_t = (2 * pi * m * k_B * T / h^2)^(3/2) * (k_B * T / P)
        where:
            m = mass of the molecule (kg)
            k_B = Boltzmann constant (J K^-1)
            T = temperature (K)
            h = Planck constant (J s)
            P = pressure of the system (Pa)
        """
        return (
            2 * np.pi * self.m * units._k * self.T / units._hplanck**2
        ) ** (3 / 2) * (units._k * self.T / self.P)

    @property
    def translational_entropy(self):
        """Obtain the translational entropy.
        Formula:
            S_t = R * [ln(q_t) + 1 + d/2]
        where:
            R = gas constant (J mol^-1 K^-1)
            d = 3 for non-linear molecules;
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

    @property
    def translational_internal_energy(self):
        """Obtain the translational internal energy.
        Same for all types of molecules, whether linear, non-linear or monoatomic.
        Formula:
            U_t = 3/2 * R * T
        """
        return 3 / 2 * R * self.T

    @property
    def translational_heat_capacity(self):
        """Obtain the constant volume heat capacity.
        Formula:
            C_t = 3/2 * R
        """
        return 3 / 2 * R

    @property
    def electronic_partition_function(self):
        """Obtain the electronic partition function.
        Gaussian assumes first electronic excitation energy is much greater than k_B * T.
        Thus, first and higher excited states assumed inaccessible at any temperature.
        Further, energy of ground state is set to zero.
        Formula:
            q_e = ω_0
        where:
            ω_0 = degeneracy of the ground state,
            which is simply the electronic spin multiplicity of the molecule.
        """
        return self.file_object.multiplicity

    @property
    def electronic_entropy(self):
        """Obtain the electronic entropy.
        Formula:
            S_e = R * ln(q_e)
        """
        return R * np.log(self.electronic_partition_function)

    @property
    def electronic_internal_energy(self):
        """The internal thermal energy due to electronic motion is zero.
        Since there are no temperature dependent terms in electronic partition function
        """
        return 0

    @property
    def electronic_heat_capacity(self):
        """The electronic heat capacity is zero for all types of molecules.
        C_V = (\partial U_e / \partial T)_V = 0"""
        return 0

    def _calculate_rotational_partition_function_for_linear_molecule(self):
        """Calculate the rotational partition function of a linear molecule.
        Formula:
            q_r = 1 / σ_r * (T / Θ_r)
        where:
            σ_r = symmetry number for rotation
            Θ_r = h^2 / (8 * pi^2 * I * k_B)
            I = moment of inertia (kg m^2)
        """
        theta_r = units._hplanck**2 / (8 * np.pi**2 * self.I[0] * units._k)
        return (1 / self.rotational_symmetry_number) * (self.T / theta_r)

    def _calculate_rotational_partition_function_for_nonlinear_polyatomic_molecule(
        self,
    ):
        """Calculate the rotational partition function of a nonlinear polyatomic molecule.
        Formula:
            q_r = pi^(1/2) / σ_r * (T^(3/2) / (Θ_r,x * Θ_r,y * Θ_r,z)^(1/2))
        where:
            σ_r = symmetry number for rotation
            Θ_r,i = h^2 / (8 * pi^2 * I_i * k_B) for i = x, y, z
        """
        theta_ri = [
            units._hplanck**2 / (8 * np.pi**2 * i * units._k) for i in self.I
        ]
        return (
            np.pi ** (1 / 2)
            / self.rotational_symmetry_number
            * (self.T ** (3 / 2) / np.prod(theta_ri) ** (1 / 2))
        )

    @property
    def rotational_partition_function(self):
        """Obtain the rotational partition function.
        For a single atom, q_r = 1. Since q_r does not depend on temperature, contribution
        of rotation to internal thermal energy, heat capacity and entropy are all identically zero.
        """
        if self.molecule.is_monoatomic:
            return 1
        elif self.molecule.is_linear:
            return (
                self._calculate_rotational_partition_function_for_linear_molecule()
            )
        else:
            return (
                self._calculate_rotational_partition_function_for_nonlinear_polyatomic_molecule()
            )

    @property
    def rotational_entropy(self):
        """Obtain the rotational entropy.
        Formula:
            S_r = 0 for monoatomic molecules
                = R * (ln(q_r) + 1) for linear molecules
                = R * (ln(q_r) + 3/2) for nonlinear polyatomic molecules
        """
        if self.molecule.is_monoatomic:
            return 0
        elif self.molecule.is_linear:
            return R * (np.log(self.rotational_partition_function) + 1)
        else:
            return R * (np.log(self.rotational_partition_function) + 3 / 2)

    @property
    def rotational_internal_energy(self):
        """Obtain the rotational internal energy.
        Formula:
            U_r = 0 for monoatomic molecules
                = R * T for linear molecules
                = 3/2 * R * T for nonlinear polyatomic molecules
        """
        if self.molecule.is_monoatomic:
            return 0
        elif self.molecule.is_linear:
            return R * self.T
        else:
            return 3 / 2 * R * self.T

    @property
    def rotational_heat_capacity(self):
        """Obtain the rotational contribution to the heat capacity.
        Formula:
            C_r = 0 for monoatomic molecules
                = R for linear molecules
                = 3/2 * R for nonlinear polyatomic molecules
        """
        if self.molecule.is_monoatomic:
            return 0
        elif self.molecule.is_linear:
            return R
        else:
            return 3 / 2 * R

    @property
    def vibrational_partition_function_by_mode_BOT(self):
        """
        Obtain the partition function for each vibrational mode.
        The zero reference point is the bottom of the well (BOT).
        Formula:
            q_v,K = exp(-Θ_v,K / (2 * T)) / (1 - exp(-Θ_v,K / T))
        where:
            Θ_v,K = h * v_K / k_B
            v_K = vibrational frequency for mode K (Hz)
        """
        if self.molecule.is_monoatomic:
            return 1
        return [
            math.exp(-t / (2 * self.T)) / (1 - math.exp(-t / self.T))
            for t in self.theta
        ]

    @property
    def vibrational_partition_function_BOT(self):
        """Obtain the overall vibrational partition function with BOT.
        Formula:
            q_v = q_1 * q_2 * ... * q_vDOF
        where:
            vDOF = vibrational degrees of freedom
                 = 3 * N - 5 for linear molecules
                 = 3 * N - 6 for nonlinear polyatomic molecules
            N = number of atoms in molecule
        """
        if self.molecule.is_monoatomic:
            return 1
        return np.prod(self.vibrational_partition_function_by_mode_BOT)

    @property
    def vibrational_partition_function_by_mode_V0(self):
        """
        Obtain the partition function for each vibrational mode.
        The zero reference point is the first vibrational energy level (V=0).
        Formula:
            q_v,K = 1 / (1 - exp(-Θ_v,K / T))
        """
        if self.molecule.is_monoatomic:
            return 1
        return [1 / (1 - math.exp(-t / self.T)) for t in self.theta]

    @property
    def vibrational_partition_function_V0(self):
        """Obtain the overall vibrational partition function with V=0."""
        if self.molecule.is_monoatomic:
            return 1
        return np.prod(self.vibrational_partition_function_by_mode_V0)

    @property
    def vibrational_entropy(self):
        """Obtain the vibrational entropy.
        Formula:
            S_v = R * Σ((Θ_v,K / T) / (exp(Θ_v,K / T) - 1) - ln(1 - exp(-Θ_v,K / T)))
        """
        if self.molecule.is_monoatomic:
            return 0
        s = [
            (t / self.T) / (math.exp(t / self.T) - 1)
            - np.log(1 - math.exp(-t / self.T))
            for t in self.theta
        ]
        return R * sum(s)

    @property
    def vibrational_internal_energy(self):
        """Obtain the vibrational internal energy.
        Formula:
            U_v = R * Σ(Θ_v,K * (1/2 + 1 / (exp(Θ_v,K / T) - 1)))
        """
        if self.molecule.is_monoatomic:
            return 0
        u = [t * (1 / 2 + 1 / (math.exp(t / self.T) - 1)) for t in self.theta]
        return R * sum(u)

    @property
    def vibrational_heat_capacity(self):
        """Obtain the vibrational contribution to the heat capacity.
        Formula:
            C_v = R * Σ(exp(Θ_v,K / T) * ((Θ_v,K / T) / (exp(-Θ_v,K / T) - 1))^2)
        """
        if self.molecule.is_monoatomic:
            return 0
        c = [
            math.exp(t / self.T)
            * ((t / self.T) / (math.exp(-t / self.T) - 1)) ** 2
            for t in self.theta
        ]
        return R * sum(c)

    @property
    def total_partition_function(self):
        """Obtain the total partition function."""
        return (
            self.translational_partition_function
            * self.rotational_partition_function
            * self.electronic_partition_function
            * self.vibrational_partition_function_V0
        )

    @property
    def total_entropy(self):
        """Obtain the total entropy."""
        return (
            self.translational_entropy
            + self.rotational_entropy
            + self.electronic_entropy
            + self.vibrational_entropy
        )

    @property
    def total_internal_energy(self):
        """Obtain the total internal energy."""
        return (
            self.translational_internal_energy
            + self.rotational_internal_energy
            + self.electronic_internal_energy
            + self.vibrational_internal_energy
        )

    @property
    def total_heat_capacity(self):
        """Obtain the total heat capacity."""
        return (
            self.translational_heat_capacity
            + self.rotational_heat_capacity
            + self.electronic_heat_capacity
            + self.vibrational_heat_capacity
        )

    def get_thermochemistry(self):
        pass


class GaussianThermochemistry(Thermochemistry):
    """Class for thermochemistry analysis of Gaussian output files."""

    def get_thermochemistry(self):
        pass


class OrcaThermochemistry(Thermochemistry):
    def get_thermochemistry(self):
        pass
