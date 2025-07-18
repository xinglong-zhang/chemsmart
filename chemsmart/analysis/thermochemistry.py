import logging
import math
import os

import numpy as np
from ase import units

from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.orca.output import ORCAOutput
from chemsmart.utils.constants import (
    R,
    atm_to_pa,
    energy_conversion,
    hartree_to_joules,
)
from chemsmart.utils.references import (
    grimme_quasi_rrho_entropy_ref,
    head_gordon_damping_function_ref,
    head_gordon_quasi_rrho_enthalpy_ref,
    qrrho_header,
)

logger = logging.getLogger(__name__)


class Thermochemistry:
    """Class for thermochemistry analysis. Use SI units.
    Requires filename from which thermochemistry data is extracted.
    Args:
        filename: str. Filepath to the file from which thermochemistry data is extracted.
        temperature: float. Temperature of the system, in K.
        concentration: float. Concentration of the system, in mol/L.
        pressure: float. Pressure of the system, in atm.
        use_weighted_mass: bool. If True, use natural abundance weighted masses; otherwise, use most abundant masses.
        alpha: int. Interpolator exponent used in the quasi-RRHO approximation.
        s_freq_cutoff: float. The cutoff frequency of the damping function used in calculating entropy.
        h_freq_cutoff: float. The cutoff frequency of the damping function used in calculating enthalpy.
        energy_units: str. The energy units to use for output. Default is "hartree".
        outputfile: str. The output file to save the thermochemistry results.
        check_imaginary_frequencies: bool. If True, checks for imaginary frequencies in the vibrational analysis.
    """

    def __init__(
        self,
        filename,
        temperature=None,
        concentration=None,
        pressure=1.0,
        use_weighted_mass=False,
        alpha=4,
        s_freq_cutoff=None,
        h_freq_cutoff=None,
        energy_units="hartree",
        outputfile=None,
        check_imaginary_frequencies=True,
    ):
        self.filename = filename
        self.molecule = Molecule.from_filepath(filename)
        self.temperature = temperature
        self.pressure = pressure
        self.use_weighted_mass = use_weighted_mass
        self.m = (
            self.mass
            * units._amu  # converts mass from g/mol to kg/molecule
            # units._amu is same as divide by Avogadro's number then by 1000 (g to kg)
        )  # convert the unit of mass of the molecule from amu to kg
        self.T = self.temperature  # temperature in K
        self.P = (
            self.pressure * atm_to_pa
        )  # convert the unit of pressure from atm to Pascal
        self.concentration = concentration
        self.alpha = alpha
        self.s_freq_cutoff = (
            s_freq_cutoff * units._c * 1e2 if s_freq_cutoff else None
        )  # convert the unit of cutoff frequency from cm^-1 to Hz
        self.h_freq_cutoff = (
            h_freq_cutoff * units._c * 1e2 if h_freq_cutoff else None
        )  # convert the unit of cutoff frequency from cm^-1 to Hz
        self.c = (
            self.concentration * 1000 * units._Nav
            if self.concentration is not None
            else None
        )  # convert the unit of concentration from mol/L to Particle/m^3

        self.I = [
            i * (units._amu * (units.Ang / units.m) ** 2)
            for i in self.moments_of_inertia
        ]

        # convert the unit of moments of inertia from amu Ang^2 to kg m^2
        self.Bav = (
            (units._hplanck / self.average_rotational_constant)
            if self.average_rotational_constant
            else None
        )

        # convert the unit of vibrational frequencies from cm^-1 to Hz
        self.v = (
            [k * units._c * 1e2 for k in self.cleaned_frequencies]
            if self.cleaned_frequencies is not None
            else None
        )

        # Calculate the characteristic vibrational temperature, theta, for each vibrational mode
        self.theta = (
            [units._hplanck * vk / units._k for vk in self.v]
            if self.v is not None
            else None
        )

        self.energy_units = energy_units
        self.outputfile = outputfile
        self.check_imaginary_frequencies = check_imaginary_frequencies

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
    def job_type(self):
        return self.file_object.job_type

    @property
    def mass(self):
        """Obtain the molecular mass."""
        if self.use_weighted_mass:
            return self.molecule.natural_abundance_weighted_mass
        return self.molecule.most_abundant_mass

    @property
    def moments_of_inertia(self):
        """Obtain the moments of inertia of the molecule along principal axes.
        Direct calculation from molecular structure, since sometimes Gaussian
        output does not print it properly (prints as ***** if values too large)
        """
        if self.use_weighted_mass:
            return self.molecule.moments_of_inertia_weighted_mass
        return self.molecule.moments_of_inertia_most_abundant_mass

    @property
    def average_rotational_constant(self):
        if self.molecule.is_monoatomic:
            return None
        if self.molecule.is_linear:
            return units._hplanck / (8 * np.pi**2 * self.I[-1])

        assert (
            len(self.I) == 3
        ), "Number of moments of inertia should be 3 for nonlinear molecules."
        rotational_constants = [
            units._hplanck / (8 * np.pi**2 * i) for i in self.I
        ]
        return sum(rotational_constants) / len(rotational_constants)

    @property
    def rotational_symmetry_number(self):
        """Obtain the rotational symmetry number."""
        return self.file_object.rotational_symmetry_number

    @property
    def vibrational_frequencies(self):
        """Obtain the vibrational frequencies of the molecule."""
        if not self.file_object.freq:
            return None
        return self.file_object.vibrational_frequencies

    @property
    def real_frequencies(self):
        """Obtain the real vibrational frequencies of the molecule."""
        if self.vibrational_frequencies is None:
            return None
        return [k for k in self.vibrational_frequencies if k >= 0.0]

    @property
    def imaginary_frequencies(self):
        """Obtain the imaginary vibrational frequencies of the molecule."""
        if self.vibrational_frequencies is None:
            return None
        return [k for k in self.vibrational_frequencies if k < 0.0]

    @property
    def cleaned_frequencies(self):
        """Clean up the vibrational frequencies for thermochemical calculations.
        For transition states (job_type == "ts"), the first imaginary frequency is assumed to
        correspond to the reaction coordinate and is excluded from thermochemical calculation.
        For optimization, only geometries without imaginary frequencies are parsed for thermochemical calculations.
        """
        if self.vibrational_frequencies is None:
            return None
        if self.imaginary_frequencies:
            if self.job_type == "ts":
                if (
                    len(self.imaginary_frequencies) == 1
                    and self.vibrational_frequencies[0] < 0.0
                ):
                    return self.vibrational_frequencies[1:]
                else:
                    raise ValueError(
                        f"!! ERROR: Detected multiple imaginary frequencies in TS calculation for {self.filename}. "
                        f"Only one imaginary frequency is allowed for a valid TS. "
                        f"Please re-optimize the geometry to locate a true TS."
                    )
            else:
                raise ValueError(
                    f"!! ERROR: Detected imaginary frequencies in geometry optimization for {self.filename}. "
                    f"A valid optimized geometry should not contain imaginary frequencies. "
                    f"Please re-optimize the geometry to locate a true minimum."
                )
        return self.vibrational_frequencies

    @property
    def electronic_energy(self):
        """Obtain the total electronic energy in J mol^-1."""
        return self.file_object.energies[-1] * hartree_to_joules * units._Nav

    @property
    def multiplicity(self):
        """Obtain the multiplicity of the molecule."""
        return self.file_object.multiplicity

    @property
    def translational_partition_function(self):
        """Obtain the translational partition function.
        Formula in gas phase:
            q_t = (2 * pi * m * k_B * T / h^2)^(3/2) * (k_B * T / P)
        In solution, uses concentration instead of pressure. Formula:
            q_t = (2 * pi * m * k_B * T / h^2)^(3/2) * (1 / c)
        where:
            m = mass of the molecule (kg)
            k_B = Boltzmann constant (J K^-1)
            T = temperature (K)
            h = Planck constant (J s)
            P = pressure of the system (Pa)
            c = pressure of the system (m^-3)
        """
        if self.c is not None:
            return (
                2 * np.pi * self.m * units._k * self.T / units._hplanck**2
            ) ** (3 / 2) * (1 / self.c)
        else:
            return (
                2 * np.pi * self.m * units._k * self.T / units._hplanck**2
            ) ** (3 / 2) * (units._k * self.T / self.P)

    @property
    def translational_entropy(self):
        """Obtain the translational entropy in J mol^-1 K^-1.
        Formula:
            S_t = R * [ln(q_t) + 1 + 3/2]
        where:
            R = gas constant (J mol^-1 K^-1)
        """
        return R * (np.log(self.translational_partition_function) + 1 + 3 / 2)

    @property
    def translational_internal_energy(self):
        """Obtain the translational internal energy J mol^-1.
        Same for all types of molecules, whether linear, non-linear or monoatomic.
        Formula:
            E_t = 3/2 * R * T
        """
        return 3 / 2 * R * self.T

    @property
    def translational_heat_capacity(self):
        """Obtain the constant volume heat capacity in J mol^-1 K^-1.
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
        return self.multiplicity

    @property
    def electronic_entropy(self):
        """Obtain the electronic entropy in J mol^-1 K^-1.
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
        theta_r = units._hplanck**2 / (8 * np.pi**2 * self.I[-1] * units._k)
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
        """Obtain the rotational entropy in J mol^-1 K^-1.
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
        """Obtain the rotational internal energy J mol^-1.
        Formula:
            E_r = 0 for monoatomic molecules
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
        """Obtain the rotational contribution to the heat capacity in J mol^-1 K^-1.
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
    def vibrational_partition_function_by_mode_bot(self):
        """
        Obtain the partition function for each vibrational mode.
        The zero reference point is the bottom of the well (BOT).
        Formula:
            q_v,K = exp(-Θ_v,K / (2 * T)) / (1 - exp(-Θ_v,K / T))
        where:
            Θ_v,K = h * v_K / k_B
            v_K = vibrational frequency for mode K (Hz)
        """
        return (
            [
                math.exp(-t / (2 * self.T)) / (1 - math.exp(-t / self.T))
                for t in self.theta
            ]
            if self.theta is not None
            else None
        )

    @property
    def vibrational_partition_function_bot(self):
        """Obtain the overall vibrational partition function with BOT.
        Formula:
            q_v = q_1 * q_2 * ... * q_vDOF
        where:
            vDOF = vibrational degrees of freedom
                 = 3 * N - 5 for linear molecules
                 = 3 * N - 6 for nonlinear polyatomic molecules
            N = number of atoms in molecule
        """
        if self.vibrational_partition_function_by_mode_bot is None:
            return None
        return np.prod(self.vibrational_partition_function_by_mode_bot)

    @property
    def vibrational_partition_function_by_mode_v0(self):
        """
        Obtain the partition function for each vibrational mode.
        The zero reference point is the first vibrational energy level (V=0).
        Formula:
            q_v,K = 1 / (1 - exp(-Θ_v,K / T))
        """
        return (
            [1 / (1 - math.exp(-t / self.T)) for t in self.theta]
            if self.theta is not None
            else None
        )

    @property
    def vibrational_partition_function_v0(self):
        """Obtain the overall vibrational partition function with V=0."""
        if self.vibrational_partition_function_by_mode_v0 is None:
            return None
        return np.prod(self.vibrational_partition_function_by_mode_v0)

    @property
    def vibrational_entropy(self):
        """Obtain the vibrational entropy in J mol^-1 K^-1.
        Formula:
            S_v = R * Σ((Θ_v,K / T) / (exp(Θ_v,K / T) - 1) - ln(1 - exp(-Θ_v,K / T)))
        """
        if self.theta is None:
            return None
        s = [
            (t / self.T) / (math.exp(t / self.T) - 1)
            - np.log(1 - math.exp(-t / self.T))
            for t in self.theta
        ]
        return R * sum(s)

    @property
    def zero_point_energy(self):
        """Obtain the vibrational zero-point energy (ZPE) in J mol^-1.
        Formula:
            E_ZPE = R * Σ(1/2 * Θ_v,K)
        """
        if self.theta is None:
            return None
        u = [1 / 2 * t for t in self.theta]
        return R * sum(u)

    @property
    def vibrational_internal_energy(self):
        """Obtain the vibrational internal energy in J mol^-1.
        Formula:
            E_v = R * Σ(Θ_v,K * (1/2 + 1 / (exp(Θ_v,K / T) - 1)))
        """
        if self.theta is None:
            return None
        u = [t * (1 / 2 + 1 / (math.exp(t / self.T) - 1)) for t in self.theta]
        return R * sum(u)

    @property
    def vibrational_heat_capacity(self):
        """Obtain the vibrational contribution to the heat capacity in J mol^-1 K^-1.
        Formula:
            C_v = R * Σ(exp(-Θ_v,K / T) * ((Θ_v,K / T) / (exp(-Θ_v,K / T) - 1))^2)
        """
        if self.theta is None:
            return None
        c = [
            math.exp(-t / self.T)
            * ((t / self.T) / (math.exp(-t / self.T) - 1)) ** 2
            for t in self.theta
        ]
        return R * sum(c)

    @property
    def total_partition_function(self):
        """Obtain the total partition function.
        Formula:
            q_tot = q_t * q_r * q_v * q_e
        """
        if self.vibrational_partition_function_v0 is None:
            return None
        return (
            self.translational_partition_function
            * self.rotational_partition_function
            * self.electronic_partition_function
            * self.vibrational_partition_function_v0
        )

    @property
    def total_entropy(self):
        """Obtain the total entropy in J mol^-1 K^-1.
        Formula:
            S_tot = S_t + S_r + S_v + S_e
        """
        if self.vibrational_entropy is None:
            return None
        return (
            self.translational_entropy
            + self.rotational_entropy
            + self.electronic_entropy
            + self.vibrational_entropy
        )

    @property
    def total_internal_energy(self):
        """Obtain the total internal energy in J mol^-1.
        Formula:
            E_tot = E_t + E_r + E_v + E_e
        """
        if self.vibrational_internal_energy is None:
            return None
        return (
            self.translational_internal_energy
            + self.rotational_internal_energy
            + self.electronic_internal_energy
            + self.vibrational_internal_energy
        )

    @property
    def total_heat_capacity(self):
        """Obtain the total heat capacity in J mol^-1 K^-1.
        Formula:
            C_tot = C_t + C_r + C_v + C_e
        """
        if self.vibrational_heat_capacity is None:
            return None
        return (
            self.translational_heat_capacity
            + self.rotational_heat_capacity
            + self.electronic_heat_capacity
            + self.vibrational_heat_capacity
        )

    def _calculate_damping_function(self, freq_cutoff):
        """Calculate the damping function of Head-Gordon, which interpolates
        between the RRHO and the free rotor entropy.
        Formula:
            w(v_K) = 1 / (1 + (v_0 / v_K)^α)
        where:
            v_0 = cutoff frequency in Hz, default is 100 cm^-1 (already converted to Hz)
            α = dimensionless interpolator exponent, default value is 4
        """
        if freq_cutoff is None or self.v is None:
            return None
        damp = [1 / (1 + (freq_cutoff / vk) ** self.alpha) for vk in self.v]
        return damp

    @property
    def entropy_damping_function(self):
        return self._calculate_damping_function(self.s_freq_cutoff)

    @property
    def enthalpy_damping_function(self):
        return self._calculate_damping_function(self.h_freq_cutoff)

    @property
    def free_rotor_entropy(self):
        """Obtain the free rotor entropy in J mol^-1 K^-1, which is used to treat
        low frequency modes below cutoff.
        Formula:
            S_R,K = R * (1/2 + ln((8 * pi^3 * u'_K * k_B * T / h^2)^(1/2)))
        where:
            u'_K = u_K * B_av / (u_K + B_av)
            u_K = h / (8 * pi^2 * v_K)
            B_av = average molecular moment of inertia (kg m^2)
        """
        if self.v is None:
            return None
        bav = self.Bav
        mu = [units._hplanck / (8 * np.pi**2 * vk) for vk in self.v]
        mu_prime = [mu_k * bav / (mu_k + bav) for mu_k in mu]
        entropy = [
            R
            * (
                1 / 2
                + np.log(
                    (
                        8
                        * np.pi**3
                        * mu_prime_k
                        * units._k
                        * self.T
                        / units._hplanck**2
                    )
                    ** (1 / 2)
                )
            )
            for mu_prime_k in mu_prime
        ]
        return entropy

    @property
    def rrho_entropy(self):
        """Obtain the Harmonic Oscillator (within RRHO approximation)
        vibrational entropy in J mol^-1 K^-1.
        Formula:
            S^rrho_v,K = R * [(Θ_v,K / T) / (exp(Θ_v,K / T) - 1) - ln(1 - exp(-Θ_v,K / T))]
        """
        if self.theta is None:
            return None
        entropy = [
            R
            * (
                (t / self.T) / (math.exp(t / self.T) - 1)
                - np.log(1 - math.exp(-t / self.T))
            )
            for t in self.theta
        ]
        return entropy

    @property
    def qrrho_vibrational_entropy(self):
        """Obtain the vibrational entropy with quasi-RRHO approximation, in J mol^-1 K^-1.
        Formula:
            S^qrrho_v = Σ(w(v_K) * S^rrho_v,K + (1 - w(v_K)) * S_R,K)
        """
        if self.s_freq_cutoff is None or self.v is None:
            return None
        vib_entropy = []
        assert len(self.v) == len(self.entropy_damping_function), (
            f"The length of vibrational frequencies and damping function must be equal.\n"
            f"The damping function is {self.entropy_damping_function}.\n"
        )
        for j in range(0, len(self.v)):
            vib_entropy.append(
                self.entropy_damping_function[j] * self.rrho_entropy[j]
                + (1 - self.entropy_damping_function[j])
                * self.free_rotor_entropy[j]
            )
        return sum(vib_entropy)

    @property
    def rrho_internal_energy(self):
        """Obtain the Harmonic Oscillator (within RRHO approximation)
         vibrational internal energy in J mol^-1.
        Formula:
            E^rrho_v,K = R * Θ_v,K * (1/2 + 1 / (exp(Θ_v,K / T) - 1))
        """
        if self.theta is None:
            return None
        energy = [
            R * t * (1 / 2 + 1 / (math.exp(t / self.T) - 1))
            for t in self.theta
        ]
        return energy

    @property
    def qrrho_vibrational_internal_energy(self):
        """Obtain the vibrational internal energy with quasi-RRHO approximation, in J mol^-1.
        Formula:
            E^qrrho_v = Σ(w(v_K) * E^rrho_v,K + (1 - w(v_K)) * 1/2 * R * T)
        """
        if self.h_freq_cutoff is None or self.v is None:
            return None
        vib_energies = []
        assert len(self.v) == len(self.enthalpy_damping_function), (
            f"The length of vibrational frequencies and damping function must be equal.\n"
            f"The damping function is {self.enthalpy_damping_function}.\n"
        )
        for j in range(0, len(self.v)):
            vib_energies.append(
                self.enthalpy_damping_function[j]
                * self.rrho_internal_energy[j]
                + (1 - self.enthalpy_damping_function[j]) * 1 / 2 * R * self.T
            )
        return sum(vib_energies)

    @property
    def enthalpy(self):
        """Obtain the enthalpy in J mol^-1.
        Formula:
            H = E0 + E_tot + R * T
        where:
            E0 = the total electronic energy (J mol^-1)
        """
        if self.total_internal_energy is None:
            return None
        return self.electronic_energy + self.total_internal_energy + R * self.T

    @property
    def qrrho_total_entropy(self):
        """Obtain the quasi-RRHO total entropy in J mol^-1 K^-1.
        Formula:
            S^qrrho_tot = S_t + S_r + S^qrrho_v + S_e
        """
        if self.qrrho_vibrational_entropy is None:
            return None
        return (
            self.translational_entropy
            + self.rotational_entropy
            + self.electronic_entropy
            + self.qrrho_vibrational_entropy
        )

    @property
    def entropy_times_temperature(self):
        """Obtain the total entropy times temperature in J mol^-1.
        Formula:
            T * S_tot
        """
        if self.total_entropy is None:
            return None
        return self.T * self.total_entropy

    @property
    def qrrho_entropy_times_temperature(self):
        """Obtain the quasi-RRHO entropy times temperature in J mol^-1.
        Formula:
            T * S^qrrho_tot
        """
        if self.qrrho_total_entropy is None:
            return None
        return self.T * self.qrrho_total_entropy

    @property
    def gibbs_free_energy(self):
        """Obtain the Gibbs free energy in J mol^-1 .
        Formula:
            G = H - T * S_tot
        """
        if self.entropy_times_temperature is None or self.enthalpy is None:
            return None
        return self.enthalpy - self.entropy_times_temperature

    @property
    def qrrho_total_internal_energy(self):
        """Obtain the quasi-RRHO total internal energy in J mol^-1.
        Formula:
            E^qrrho_tot = E_t + E_r + E^qrrho_v + E_e
        """
        if self.qrrho_vibrational_internal_energy is None:
            return None
        return (
            self.translational_internal_energy
            + self.rotational_internal_energy
            + self.electronic_internal_energy
            + self.qrrho_vibrational_internal_energy
        )

    @property
    def qrrho_enthalpy(self):
        """Obtain the quasi-RRHO enthalpy in J mol^-1.
        Formula:
            H^qrrho = E0 + H^qrrho_corr
                    = E0 + E^qrrho_tot + R * T
        where:
            E0 = the total electronic energy (J mol^-1)
        """
        if self.qrrho_total_internal_energy is None:
            return None
        return (
            self.electronic_energy
            + self.qrrho_total_internal_energy
            + R * self.T
        )

    @property
    def qrrho_gibbs_free_energy(self):
        """Obtain the Gibbs free energy in J mol^-1, by quasi-RRHO corrections to both entropy and enthalpy.
        Formula:
            G^qrrho_q = H^qrrho - T * S^qrrho_tot
        """
        if (
            self.qrrho_enthalpy is None
            or self.qrrho_entropy_times_temperature is None
        ):
            return None
        return self.qrrho_enthalpy - self.qrrho_entropy_times_temperature

    @property
    def qrrho_gibbs_free_energy_qs(self):
        """Obtain the Gibbs free energy in J mol^-1, by a quasi-RRHO correction to entropy only.
        Formula:
            G^qrrho_qs = H - T * S^qrrho_tot
        """
        if (
            self.qrrho_entropy_times_temperature is None
            or self.enthalpy is None
        ):
            return None
        return self.enthalpy - self.qrrho_entropy_times_temperature

    @property
    def qrrho_gibbs_free_energy_qh(self):
        """Obtain the Gibbs free energy in J mol^-1, by a quasi-RRHO correction to enthalpy only.
        Formula:
            G^qrrho_qh = H^qrrho - T * S_tot
        """
        if (
            self.qrrho_enthalpy is None
            or self.entropy_times_temperature is None
        ):
            return None
        return self.qrrho_enthalpy - self.entropy_times_temperature

    def compute_thermochemistry(self):
        """Compute Boltzmann-averaged properties."""
        logger.debug(f"Computing thermochemistry for {self.filename}...")
        self._compute_thermochemistry()

    def _compute_thermochemistry(self):
        """Calculate thermochemical properties based on the parsed data."""
        # Check for imaginary frequencies if required
        if self.check_imaginary_frequencies:
            logger.debug("Checking imaginary frequencies.")
            self.check_frequencies()
        # convert energies to specified units
        logger.debug(f"Converting to energy units: {self.energy_units}")
        (
            electronic_energy,
            zero_point_energy,
            enthalpy,
            qrrho_enthalpy,
            entropy_times_temperature,
            qrrho_entropy_times_temperature,
            gibbs_free_energy,
            qrrho_gibbs_free_energy,
        ) = self.convert_energy_units()
        logger.debug(f"Finished converting energies to {self.energy_units}.")

        # Log the results to the output file or console
        structure = os.path.splitext(os.path.basename(self.filename))[0]
        self.log_results_to_file(
            structure,
            electronic_energy,
            zero_point_energy,
            enthalpy,
            qrrho_enthalpy,
            entropy_times_temperature,
            qrrho_entropy_times_temperature,
            gibbs_free_energy,
            qrrho_gibbs_free_energy,
        )

    def check_frequencies(self):
        """Check for imaginary frequencies and raise an error if found."""
        if self.imaginary_frequencies:
            if self.job_type == "ts":
                if len(self.imaginary_frequencies) == 1:
                    logger.info(
                        f"Correct Transition State detected: only 1 imaginary frequency\n"
                        f"Imaginary frequency excluded for thermochemistry calculation in {self.filename}."
                    )
                else:
                    raise ValueError(
                        f"Invalid number of imaginary frequencies for {self.filename}. "
                        f"Expected 0 for optimization or 1 for TS, but found "
                        f"{len(self.imaginary_frequencies)} for job: {self.job_type}!"
                    )
            else:
                raise ValueError(
                    f"Invalid geometry optimization for {self.filename}. "
                    f"A valid optimized geometry should not contain imaginary frequencies. "
                    f"Please re-optimize the geometry to locate a true minimum."
                )

    def convert_energy_units(self):
        """Convert all energies to the specified units."""
        electronic_energy = energy_conversion(
            "j/mol", self.energy_units, self.electronic_energy
        )
        zero_point_energy = energy_conversion(
            "j/mol", self.energy_units, self.zero_point_energy
        )
        enthalpy = energy_conversion("j/mol", self.energy_units, self.enthalpy)
        qrrho_enthalpy = (
            energy_conversion("j/mol", self.energy_units, self.qrrho_enthalpy)
            if self.qrrho_enthalpy
            else None
        )
        entropy_times_temperature = energy_conversion(
            "j/mol", self.energy_units, self.entropy_times_temperature
        )
        qrrho_entropy_times_temperature = (
            energy_conversion(
                "j/mol",
                self.energy_units,
                self.qrrho_entropy_times_temperature,
            )
            if self.qrrho_entropy_times_temperature
            else None
        )
        gibbs_free_energy = energy_conversion(
            "j/mol", self.energy_units, self.gibbs_free_energy
        )

        if self.s_freq_cutoff and self.h_freq_cutoff:
            qrrho_gibbs_free_energy = energy_conversion(
                "j/mol", self.energy_units, self.qrrho_gibbs_free_energy
            )
        elif self.s_freq_cutoff and not self.h_freq_cutoff:
            qrrho_gibbs_free_energy = energy_conversion(
                "j/mol", self.energy_units, self.qrrho_gibbs_free_energy_qs
            )
        elif not self.s_freq_cutoff and self.h_freq_cutoff:
            qrrho_gibbs_free_energy = energy_conversion(
                "j/mol", self.energy_units, self.qrrho_gibbs_free_energy_qh
            )
        else:
            qrrho_gibbs_free_energy = None

        return (
            electronic_energy,
            zero_point_energy,
            enthalpy,
            qrrho_enthalpy,
            entropy_times_temperature,
            qrrho_entropy_times_temperature,
            gibbs_free_energy,
            qrrho_gibbs_free_energy,
        )

    def __str__(self):
        """String representation of the thermochemistry results."""
        filename = getattr(self, 'filename', 'Unknown')
        temperature = getattr(self, 'temperature', None)
        concentration = getattr(self, 'concentration', None)
        pressure = getattr(self, 'pressure', None)
        use_weighted_mass = getattr(self, 'use_weighted_mass', False)
        energy_units = getattr(self, 'energy_units', 'Unknown')

        temperature_str = f"{temperature:.2f} K" if temperature is not None else "N/A"
        concentration_str = f"{concentration:.1f} mol/L" if concentration is not None else "N/A"
        pressure_str = f"{pressure:.1f} atm" if pressure is not None else "N/A"
        mass_weighted_str = (
            "Most Abundant Masses" if not use_weighted_mass else "Natural Abundance Weighted Masses"
        )

        return (
            f"Thermochemistry Results for {filename}:\n"
            f"Temperature: {temperature_str}\n"
            f"Concentration: {concentration_str}\n"
            f"Pressure: {pressure_str}\n"
            f"Mass Weighted: {mass_weighted_str}\n"
            f"Energy Unit: {energy_units}\n"
        )

    def log_results_to_file(
        self,
        structure,
        electronic_energy,
        zero_point_energy,
        enthalpy,
        qrrho_enthalpy,
        entropy_times_temperature,
        qrrho_entropy_times_temperature,
        gibbs_free_energy,
        qrrho_gibbs_free_energy,
    ):
        all_none = all(
            x is None
            for x in [
                zero_point_energy,
                enthalpy,
                qrrho_enthalpy,
                entropy_times_temperature,
                qrrho_entropy_times_temperature,
                gibbs_free_energy,
                qrrho_gibbs_free_energy,
            ]
        )
        no_freq = "{:39} {:13.6f}   {:<69}\n".format(
            structure,
            electronic_energy,
            "--- [NO FREQ INFO] Thermochemistry skipped. ---",
        )

        """Log the thermochemistry results to the output file."""
        if self.outputfile is None:
            # If no output file is specified, log to individual output file
            self.outputfile = os.path.splitext(self.filename)[0] + ".dat"

        # Write output to outputfile
        if not os.path.exists(self.outputfile):
            output_string = f"Temperature: {self.temperature:.2f} K\n"
            if self.concentration is not None:
                output_string += (
                    f"Concentration: {self.concentration:.1f} mol/L\n"
                )
            else:
                output_string += f"Pressure: {self.pressure:.1f} atm\n"

            if self.s_freq_cutoff:
                output_string += f"Entropy Frequency Cut-off: {(self.s_freq_cutoff/(units._c * 1e2)):.1f} cm^-1\n"

            if self.h_freq_cutoff:
                output_string += f"Enthalpy Frequency Cut-off: {(self.h_freq_cutoff/(units._c * 1e2)):.1f} cm^-1\n"
            if self.s_freq_cutoff or self.h_freq_cutoff:
                output_string += f"Damping Function Exponent: {self.alpha}\n"
            output_string += f"Mass Weighted: {'Most Abundant Masses' if not self.use_weighted_mass else 'Natural Abundance Weighted Masses'}\n"
            output_string += f"Energy Unit: {self.energy_units}\n\n"

            if self.h_freq_cutoff or self.s_freq_cutoff:
                output_string += qrrho_header
                output_string += head_gordon_damping_function_ref
            if self.s_freq_cutoff:
                output_string += grimme_quasi_rrho_entropy_ref
            if self.h_freq_cutoff:
                output_string += head_gordon_quasi_rrho_enthalpy_ref
            output_string += "\n"

            if self.h_freq_cutoff and self.s_freq_cutoff:
                output_string += "{:<39} {:>13} {:>10} {:>13} {:>13} {:>10} {:>10} {:>13} {:>13}\n".format(
                    "Structure",
                    "E",
                    "ZPE",
                    "H",
                    "qh-H",
                    "T.S",
                    "T.qh-S",
                    "G(T)",
                    "qh-G(T)",
                )

                output_string += "=" * 142 + "\n"
                output_string += (
                    "{:39} {:13.6f} {:10.6f} {:13.6f} {:13.6f} {:10.6f} {:10.6f} {:13.6f} {:13.6f}\n".format(
                        structure,
                        electronic_energy,
                        zero_point_energy,
                        enthalpy,
                        qrrho_enthalpy,
                        entropy_times_temperature,
                        qrrho_entropy_times_temperature,
                        gibbs_free_energy,
                        qrrho_gibbs_free_energy,
                    )
                    if not all_none
                    else no_freq
                )

            elif self.s_freq_cutoff and not self.h_freq_cutoff:
                output_string += "{:<39} {:>13} {:>10} {:>13} {:>10} {:>10} {:>13} {:>13}\n".format(
                    "Structure",
                    "E",
                    "ZPE",
                    "H",
                    "T.S",
                    "T.qh-S",
                    "G(T)",
                    "qh-G(T)",
                )
                output_string += "=" * 128 + "\n"
                output_string += (
                    "{:39} {:13.6f} {:10.6f} {:13.6f} {:10.6f} {:10.6f} {:13.6f} {:13.6f}\n".format(
                        structure,
                        electronic_energy,
                        zero_point_energy,
                        enthalpy,
                        entropy_times_temperature,
                        qrrho_entropy_times_temperature,
                        gibbs_free_energy,
                        qrrho_gibbs_free_energy,
                    )
                    if not all_none
                    else no_freq
                )
            elif self.h_freq_cutoff and not self.s_freq_cutoff:
                output_string += "{:<39} {:>13} {:>10} {:>13} {:>13} {:>10} {:>13} {:>13}\n".format(
                    "Structure",
                    "E",
                    "ZPE",
                    "H",
                    "qh-H",
                    "T.S",
                    "G(T)",
                    "qh-G(T)",
                )
                output_string += "=" * 131 + "\n"
                output_string += (
                    "{:39} {:13.6f} {:10.6f} {:13.6f} {:13.6f} {:10.6f} {:13.6f} {:13.6f}\n".format(
                        structure,
                        electronic_energy,
                        zero_point_energy,
                        enthalpy,
                        qrrho_enthalpy,
                        entropy_times_temperature,
                        gibbs_free_energy,
                        qrrho_gibbs_free_energy,
                    )
                    if not all_none
                    else no_freq
                )
            else:
                output_string += (
                    "{:<39} {:>13} {:>10} {:>13} {:>10} {:>13}\n".format(
                        "Structure", "E", "ZPE", "H", "T.S", "G(T)"
                    )
                )
                output_string += "=" * 103 + "\n"
                output_string += (
                    "{:39} {:13.6f} {:10.6f} {:13.6f} {:10.6f} {:13.6f}\n".format(
                        structure,
                        electronic_energy,
                        zero_point_energy,
                        enthalpy,
                        entropy_times_temperature,
                        gibbs_free_energy,
                    )
                    if not all_none
                    else no_freq
                )

            # If the output file does not exist, create it and write the header
            with open(self.outputfile, "w") as out:
                out.write(output_string)
            logger.info(f"Thermochemistry results saved to {self.outputfile}")
        else:
            # If the output file exists, append the results
            structure = os.path.splitext(os.path.basename(self.filename))[0]
            output_string = ""
            if self.h_freq_cutoff and self.s_freq_cutoff:
                output_string += (
                    "{:39} {:13.6f} {:10.6f} {:13.6f} {:13.6f} {:10.6f} {:10.6f} {:13.6f} {:13.6f}\n".format(
                        structure,
                        electronic_energy,
                        zero_point_energy,
                        enthalpy,
                        qrrho_enthalpy,
                        entropy_times_temperature,
                        qrrho_entropy_times_temperature,
                        gibbs_free_energy,
                        qrrho_gibbs_free_energy,
                    )
                    if not all_none
                    else no_freq
                )

            elif self.s_freq_cutoff and not self.h_freq_cutoff:
                output_string += (
                    "{:39} {:13.6f} {:10.6f} {:13.6f} {:10.6f} {:10.6f} {:13.6f} {:13.6f}\n".format(
                        structure,
                        electronic_energy,
                        zero_point_energy,
                        enthalpy,
                        entropy_times_temperature,
                        qrrho_entropy_times_temperature,
                        gibbs_free_energy,
                        qrrho_gibbs_free_energy,
                    )
                    if not all_none
                    else no_freq
                )
            elif self.h_freq_cutoff and not self.s_freq_cutoff:
                output_string += (
                    "{:39} {:13.6f} {:10.6f} {:13.6f} {:13.6f} {:10.6f} {:13.6f} {:13.6f}\n".format(
                        structure,
                        electronic_energy,
                        zero_point_energy,
                        enthalpy,
                        qrrho_enthalpy,
                        entropy_times_temperature,
                        gibbs_free_energy,
                        qrrho_gibbs_free_energy,
                    )
                    if not all_none
                    else no_freq
                )
            else:
                output_string += (
                    "{:39} {:13.6f} {:10.6f} {:13.6f} {:10.6f} {:13.6f}\n".format(
                        structure,
                        electronic_energy,
                        zero_point_energy,
                        enthalpy,
                        entropy_times_temperature,
                        gibbs_free_energy,
                    )
                    if not all_none
                    else no_freq
                )

            # Check if the header is already present
            with open(self.outputfile, "a") as out:
                out.write(output_string)
            logger.info(f"Thermochemistry results saved to {self.outputfile}")


class BoltzmannAverageThermochemistry(Thermochemistry):
    """Class to compute Boltzmann-averaged thermochemical properties from a list of files."""

    def __init__(self, files, energy_type="gibbs", **kwargs):
        super().__init__(
            filename=files[
                0
            ],  # No single file, we will take molecule from first filename
            **kwargs,
        )
        """
        Initialize with a list of Gaussian or ORCA output files.

        Parameters
        ----------
        files : list of str
            List of file paths (.log or .out) containing thermochemistry data for conformers.
        energy_type : str, optional
            Energy type to use for Boltzmann weighting ("electronic" or "gibbs"). Default is "gibbs".
        """
        if not files:
            raise ValueError("List of files cannot be empty.")
        if not all(
            isinstance(f, str) and f.endswith((".log", ".out")) for f in files
        ):
            raise ValueError("All files must be .log or .out files.")

        self.files = files
        self.energy_type = energy_type.lower()
        if self.energy_type not in ["electronic", "gibbs"]:
            raise ValueError("energy_type must be 'electronic' or 'gibbs'.")

        # Create Thermochemistry instances for each file
        self.thermochemistries = [
            Thermochemistry(filename=f, **kwargs) for f in files
        ]

    def compute_boltzmann_averages(self):
        """Compute Boltzmann-averaged properties."""
        # Check for imaginary frequencies if required
        if self.check_imaginary_frequencies:
            self.check_frequencies()

        self._calculate_boltzmann_averages()

        # convert energies to specified units
        (
            electronic_energy,
            zero_point_energy,
            enthalpy,
            qrrho_enthalpy,
            entropy_times_temperature,
            qrrho_entropy_times_temperature,
            gibbs_free_energy,
            qrrho_gibbs_free_energy,
        ) = self.convert_energy_units()
        # Log the results to the output file or console
        structure = (
            os.path.commonprefix(
                [os.path.splitext(os.path.basename(f))[0] for f in self.files]
            )
            + f"_boltzmann_avg_by_{self.energy_type}"
        )
        self.log_results_to_file(
            structure,
            electronic_energy,
            zero_point_energy,
            enthalpy,
            qrrho_enthalpy,
            entropy_times_temperature,
            qrrho_entropy_times_temperature,
            gibbs_free_energy,
            qrrho_gibbs_free_energy,
        )

    def _calculate_boltzmann_averages(self):
        """Compute Boltzmann-averaged thermochemical properties."""
        # Get temperature and units from settings
        temperature = self.temperature

        # Extract energies for Boltzmann weighting
        energies = []
        for thermo in self.thermochemistries:
            if self.energy_type == "electronic":
                energy = thermo.electronic_energy  # in J/mol
            else:  # gibbs
                if self.s_freq_cutoff and self.h_freq_cutoff:
                    energy = thermo.qrrho_gibbs_free_energy
                elif self.s_freq_cutoff and not self.h_freq_cutoff:
                    energy = thermo.qrrho_gibbs_free_energy_qs
                elif self.h_freq_cutoff and not self.s_freq_cutoff:
                    energy = thermo.qrrho_gibbs_free_energy_qh
                else:
                    energy = thermo.gibbs_free_energy
            if energy is None:
                raise ValueError(
                    f"Energy ({self.energy_type}) not available for file {thermo.filename}"
                )
            energies.append(energy)
        energies = np.array(energies)

        # Compute Boltzmann weights
        beta = 1.0 / (
            R * temperature
        )  #  beta = 1 / (R * temperature) in J^-1 mol
        energies_shifted = energies - np.min(
            energies
        )  # Shift to avoid overflow
        boltzmann_factors = np.exp(-beta * energies_shifted)
        partition_function = np.sum(boltzmann_factors)
        weights = boltzmann_factors / partition_function

        # Compute weighted averages for thermochemical properties
        self._electronic_energy = np.sum(
            [
                t.electronic_energy * w
                for t, w in zip(self.thermochemistries, weights)
            ]
        )

        self._zero_point_energy = np.sum(
            [
                t.zero_point_energy * w
                for t, w in zip(self.thermochemistries, weights)
            ]
        )

        self._qrrho_enthalpy = (
            np.sum(
                [
                    t.qrrho_enthalpy * w
                    for t, w in zip(self.thermochemistries, weights)
                ]
            )
            if self.h_freq_cutoff
            else None
        )

        self._enthalpy = np.sum(
            [t.enthalpy * w for t, w in zip(self.thermochemistries, weights)]
        )

        self._qrrho_entropy = (
            np.sum(
                [
                    t.qrrho_total_entropy * w
                    for t, w in zip(self.thermochemistries, weights)
                ]
            )
            if self.s_freq_cutoff
            else None
        )

        self._entropy = np.sum(
            [
                t.total_entropy * w
                for t, w in zip(self.thermochemistries, weights)
            ]
        )

        self._qrrho_gibbs_free_energy = (
            np.sum(
                [
                    (
                        t.qrrho_gibbs_free_energy
                        if (self.s_freq_cutoff and self.h_freq_cutoff)
                        else (
                            t.qrrho_gibbs_free_energy_qs
                            if (self.s_freq_cutoff and not self.h_freq_cutoff)
                            else (
                                t.qrrho_gibbs_free_energy_qh
                                if (
                                    self.h_freq_cutoff
                                    and not self.s_freq_cutoff
                                )
                                else t.gibbs_free_energy
                            )
                        )
                    )
                    * w
                    for t, w in zip(self.thermochemistries, weights)
                ]
            )
            if not (self.h_freq_cutoff is None and self.s_freq_cutoff is None)
            else None
        )

        self._gibbs_free_energy = np.sum(
            [
                t.gibbs_free_energy * w
                for t, w in zip(self.thermochemistries, weights)
            ]
        )

    @property
    def boltzmann_electronic_energy(self):
        """Boltzmann-averaged electronic energy."""
        return self._electronic_energy

    @property
    def boltzmann_zero_point_energy(self):
        """Boltzmann-averaged zero-point energy."""
        return self._zero_point_energy

    @property
    def boltzmann_qrrho_enthalpy(self):
        """Boltzmann-averaged enthalpy."""
        return self._qrrho_enthalpy

    @property
    def boltzmann_enthalpy(self):
        """Boltzmann-averaged enthalpy."""
        return self._enthalpy

    @property
    def boltzmann_entropy(self):
        """Boltzmann-averaged entropy."""
        return self._entropy

    @property
    def boltzmann_qrrho_entropy(self):
        """Boltzmann-averaged entropy."""
        return self._qrrho_entropy

    @property
    def boltzmann_entropy_times_temperature(self):
        """Boltzmann-averaged entropy."""
        return self._entropy * self.temperature

    @property
    def boltzmann_qrrho_entropy_times_temperature(self):
        """Boltzmann-averaged entropy times temperature."""
        return self._qrrho_entropy * self.temperature

    @property
    def boltzmann_gibbs_free_energy(self):
        """Boltzmann-averaged Gibbs free energy."""
        return self._gibbs_free_energy

    @property
    def boltzmann_qrrho_gibbs_free_energy(self):
        """Boltzmann-averaged Gibbs free energy."""
        return self._qrrho_gibbs_free_energy

    def convert_energy_units(self):
        """Convert all energies to the specified units."""
        boltzmann_electronic_energy = energy_conversion(
            "j/mol", self.energy_units, self.boltzmann_electronic_energy
        )
        boltzmann_zero_point_energy = energy_conversion(
            "j/mol", self.energy_units, self.boltzmann_zero_point_energy
        )
        boltzmann_enthalpy = energy_conversion(
            "j/mol", self.energy_units, self.boltzmann_enthalpy
        )
        boltzmann_qrrho_enthalpy = (
            energy_conversion(
                "j/mol", self.energy_units, self.boltzmann_qrrho_enthalpy
            )
            if self.h_freq_cutoff
            else None
        )
        boltzmann_entropy_times_temperature = energy_conversion(
            "j/mol",
            self.energy_units,
            self.boltzmann_entropy_times_temperature,
        )
        boltzmann_qrrho_entropy_times_temperature = (
            energy_conversion(
                "j/mol",
                self.energy_units,
                self.boltzmann_qrrho_entropy_times_temperature,
            )
            if self.s_freq_cutoff
            else None
        )
        boltzmann_gibbs_free_energy = energy_conversion(
            "j/mol", self.energy_units, self.boltzmann_gibbs_free_energy
        )
        boltzmann_qrrho_gibbs_free_energy = energy_conversion(
            "j/mol", self.energy_units, self.boltzmann_qrrho_gibbs_free_energy
        )

        return (
            boltzmann_electronic_energy,
            boltzmann_zero_point_energy,
            boltzmann_enthalpy,
            boltzmann_qrrho_enthalpy,
            boltzmann_entropy_times_temperature,
            boltzmann_qrrho_entropy_times_temperature,
            boltzmann_gibbs_free_energy,
            boltzmann_qrrho_gibbs_free_energy,
        )

    def __str__(self):
        """String representation of the Boltzmann-averaged thermochemistry."""
        energy_units = self.energy_units
        return (
            f"Boltzmann-Averaged Thermochemistry (using {self.energy_type} energy)\n"
            f"Temperature: {self.temperature:.2f} K\n"
            f"Electronic Energy: {self._electronic_energy:.6f} {energy_units}\n"
            f"Enthalpy: {self._enthalpy:.6f} {energy_units}\n"
            f"Entropy: {self._entropy:.6f} {energy_units}/K\n"
            f"Gibbs Free Energy: {self._gibbs_free_energy:.6f} {energy_units}"
        )
