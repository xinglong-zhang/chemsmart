import logging
import math

import numpy as np
from ase import units

from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.orca.output import ORCAOutput
from chemsmart.utils.constants import R, atm_to_pa, hartree_to_joules

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
    """

    def __init__(
        self,
        filename,
        temperature,
        concentration=1.0,
        pressure=1.0,
        use_weighted_mass=False,
        alpha=4,
        s_freq_cutoff=100.0,
        h_freq_cutoff=100.0,
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
        self.cutoff = max(min(s_freq_cutoff, h_freq_cutoff, 100.0), 1e-6)
        self.s_freq_cutoff = (
            s_freq_cutoff * units._c * 1e2
        )  # convert the unit of cutoff frequency from cm^-1 to Hz
        self.h_freq_cutoff = (
            h_freq_cutoff * units._c * 1e2
        )  # convert the unit of cutoff frequency from cm^-1 to Hz
        self.c = (
            self.concentration * 1000 * units._Nav
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
            if self.cleaned_frequencies
            else None
        )

        # Calculate the characteristic vibrational temperature, theta, for each vibrational mode
        self.theta = (
            [units._hplanck * vk / units._k for vk in self.v]
            if self.v
            else None
        )

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
            return self.molecule.mass
        return self.molecule.most_abundant_mass

    @property
    def moments_of_inertia(self):
        """Obtain the moments of inertia of the molecule along principal axes.
        Direct calculation from molecular structure, since sometimes Gaussian
        output does not print it properly (prints as ***** if values too large)
        """
        return self.molecule.moments_of_inertia

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
        return self.file_object.vibrational_frequencies

    @property
    def real_frequencies(self):
        """Obtain the real vibrational frequencies of the molecule."""
        return [k for k in self.vibrational_frequencies if k >= 0.0]

    @property
    def imaginary_frequencies(self):
        """Obtain the imaginary vibrational frequencies of the molecule."""
        return [k for k in self.vibrational_frequencies if k < 0.0]

    @property
    def cleaned_frequencies(self):
        """Replace residue imaginary frequencies with the cutoff value.
        For transition states (job_type == "ts"), the first imaginary frequency is assumed to
        correspond to the reaction coordinate and is excluded from thermochemical calculation.
        All other negative frequencies are replaced by the cutoff value.
        """
        if not self.vibrational_frequencies:
            return None
        if self.job_type == "ts" and self.vibrational_frequencies[0] < 0.0:
            return [
                self.cutoff if k < 0.0 else k
                for k in self.vibrational_frequencies[1:]
            ]
        return [
            self.cutoff if k < 0.0 else k for k in self.vibrational_frequencies
        ]

    @property
    def num_replaced_frequencies(self):
        if self.job_type == "ts" and self.vibrational_frequencies[0] < 0.0:
            return len(self.imaginary_frequencies) - 1
        return len(self.imaginary_frequencies)

    @property
    def energies(self):
        """Obtain the total electronic energy in J mol^-1."""
        return self.file_object.energies[-1] * hartree_to_joules * units._Nav

    @property
    def multiplicity(self):
        """Obtain the multiplicity of the molecule."""
        return self.file_object.multiplicity

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
        if self.molecule.is_monoatomic:
            return 1
        return (
            [
                math.exp(-t / (2 * self.T)) / (1 - math.exp(-t / self.T))
                for t in self.theta
            ]
            if self.theta
            else []
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
        if self.molecule.is_monoatomic:
            return 1
        return np.prod(self.vibrational_partition_function_by_mode_bot)

    @property
    def vibrational_partition_function_by_mode_v0(self):
        """
        Obtain the partition function for each vibrational mode.
        The zero reference point is the first vibrational energy level (V=0).
        Formula:
            q_v,K = 1 / (1 - exp(-Θ_v,K / T))
        """
        if self.molecule.is_monoatomic:
            return 1
        return (
            [1 / (1 - math.exp(-t / self.T)) for t in self.theta]
            if self.theta
            else []
        )

    @property
    def vibrational_partition_function_v0(self):
        """Obtain the overall vibrational partition function with V=0."""
        if self.molecule.is_monoatomic:
            return 1
        return np.prod(self.vibrational_partition_function_by_mode_v0)

    @property
    def vibrational_entropy(self):
        """Obtain the vibrational entropy in J mol^-1 K^-1.
        Formula:
            S_v = R * Σ((Θ_v,K / T) / (exp(Θ_v,K / T) - 1) - ln(1 - exp(-Θ_v,K / T)))
        """
        if self.molecule.is_monoatomic:
            return 0
        s = (
            [
                (t / self.T) / (math.exp(t / self.T) - 1)
                - np.log(1 - math.exp(-t / self.T))
                for t in self.theta
            ]
            if self.theta
            else []
        )
        return R * sum(s)

    @property
    def zero_point_energy(self):
        """Obtain the vibrational zero-point energy (ZPE) in J mol^-1.
        Formula:
            E_ZPE = R * Σ(1/2 * Θ_v,K)
        """
        if self.molecule.is_monoatomic:
            return 0
        u = [1 / 2 * t for t in self.theta] if self.theta else []
        return R * sum(u)

    @property
    def vibrational_internal_energy(self):
        """Obtain the vibrational internal energy in J mol^-1.
        Formula:
            E_v = R * Σ(Θ_v,K * (1/2 + 1 / (exp(Θ_v,K / T) - 1)))
        """
        if self.molecule.is_monoatomic:
            return 0
        u = (
            [t * (1 / 2 + 1 / (math.exp(t / self.T) - 1)) for t in self.theta]
            if self.theta
            else []
        )
        return R * sum(u)

    @property
    def vibrational_heat_capacity(self):
        """Obtain the vibrational contribution to the heat capacity in J mol^-1 K^-1.
        Formula:
            C_v = R * Σ(exp(-Θ_v,K / T) * ((Θ_v,K / T) / (exp(-Θ_v,K / T) - 1))^2)
        """
        if self.molecule.is_monoatomic:
            return 0
        c = (
            [
                math.exp(-t / self.T)
                * ((t / self.T) / (math.exp(-t / self.T) - 1)) ** 2
                for t in self.theta
            ]
            if self.theta
            else []
        )
        return R * sum(c)

    @property
    def total_partition_function(self):
        """Obtain the total partition function.
        Formula:
            q_tot = q_t * q_r * q_v * q_e
        """
        if self.molecule.is_monoatomic:
            return (
                self.translational_partition_function
                * self.electronic_partition_function
            )
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
        if self.molecule.is_monoatomic:
            return self.translational_entropy + self.electronic_entropy
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
        if self.molecule.is_monoatomic:
            return (
                self.translational_internal_energy
                + self.electronic_internal_energy
            )
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
        if self.molecule.is_monoatomic:
            return (
                self.translational_heat_capacity
                + self.electronic_heat_capacity
            )
        return (
            self.translational_heat_capacity
            + self.rotational_heat_capacity
            + self.electronic_heat_capacity
            + self.vibrational_heat_capacity
        )

    @property
    def zero_point_energy_hartree(self):
        """Obtain the ZPE in Hartree."""
        return self.zero_point_energy / (hartree_to_joules * units._Nav)

    def _calculate_damping_function(self, freq_cutoff):
        """Calculate the damping function of Head-Gordon, which interpolates
        between the RRHO and the free rotor entropy.
        Formula:
            w(v_K) = 1 / (1 + (v_0 / v_K)^α)
        where:
            v_0 = cutoff frequency in Hz, default is 100 cm^-1 (already converted to Hz)
            α = dimensionless interpolator exponent, default value is 4
        """
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
        entropy = (
            [
                R
                * (
                    (t / self.T) / (math.exp(t / self.T) - 1)
                    - np.log(1 - math.exp(-t / self.T))
                )
                for t in self.theta
            ]
            if self.theta
            else []
        )
        return entropy

    @property
    def qrrho_vibrational_entropy(self):
        """Obtain the vibrational entropy with quasi-RRHO approximation, in J mol^-1 K^-1.
        Formula:
            S^qrrho_v = Σ(w(v_K) * S^rrho_v,K + (1 - w(v_K)) * S_R,K)
        """
        vib_entropy = []
        if self.v:
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
        energy = (
            [
                R * t * (1 / 2 + 1 / (math.exp(t / self.T) - 1))
                for t in self.theta
            ]
            if self.theta
            else []
        )
        return energy

    @property
    def qrrho_vibrational_internal_energy(self):
        """Obtain the vibrational internal energy with quasi-RRHO approximation, in J mol^-1.
        Formula:
            E^qrrho_v = Σ(w(v_K) * E^rrho_v,K + (1 - w(v_K)) * 1/2 * R * T)
        """
        vib_energy = []
        if self.v:
            for j in range(0, len(self.v)):
                vib_energy.append(
                    self.enthalpy_damping_function[j]
                    * self.rrho_internal_energy[j]
                    + (1 - self.enthalpy_damping_function[j])
                    * 1
                    / 2
                    * R
                    * self.T
                )
            return sum(vib_energy)

    @property
    def enthalpy(self):
        """Obtain the enthalpy in J mol^-1.
        Formula:
            H = E0 + E_tot + R * T
        where:
            E0 = the total electronic energy (J mol^-1)
        """
        return self.energies + self.total_internal_energy + R * self.T

    @property
    def translational_partition_function_concentration(self):
        """Obtain the translational partition function.
        Uses concentration instead of pressure.
        Formula:
            q_t,c = (2 * pi * m * k_B * T / h^2)^(3/2) * (1 / c)
        where:
            c = pressure of the system (m^-3)
        """
        return (
            2 * np.pi * self.m * units._k * self.T / units._hplanck**2
        ) ** (3 / 2) * (1 / self.c)

    @property
    def translational_entropy_concentration(self):
        """Obtain the translational entropy in J mol^-1 K^-1.
        Uses concentration instead of pressure.
        Formula:
            S_t,c = R * [ln(q_t) + 1 + 3/2]
        """
        return R * (
            np.log(self.translational_partition_function_concentration)
            + 1
            + 3 / 2
        )

    @property
    def total_entropy_concentration(self):
        """Obtain the total entropy in J mol^-1 K^-1.
        Uses concentration instead of pressure.
        Formula:
            S_tot,c = S_t,c + S_r + S_v + S_e
        """
        if self.molecule.is_monoatomic:
            return (
                self.translational_entropy_concentration
                + self.electronic_entropy
            )
        return (
            self.translational_entropy_concentration
            + self.rotational_entropy
            + self.electronic_entropy
            + self.vibrational_entropy
        )

    @property
    def qrrho_total_entropy(self):
        """Obtain the quasi-RRHO total entropy in J mol^-1 K^-1.
        Formula:
            S^qrrho_tot = S_t,c + S_r + S^qrrho_v + S_e
        """
        if self.molecule.is_monoatomic:
            return (
                self.translational_entropy_concentration
                + self.electronic_entropy
            )
        return (
            self.translational_entropy_concentration
            + self.rotational_entropy
            + self.electronic_entropy
            + self.qrrho_vibrational_entropy
        )

    @property
    def entropy_times_temperature(self):
        """Obtain the total entropy times temperature in J mol^-1.
        Formula:
            T * S_tot,c
        """
        return self.T * self.total_entropy_concentration

    @property
    def qrrho_entropy_times_temperature(self):
        """Obtain the quasi-RRHO entropy times temperature in J mol^-1.
        Formula:
            T * S^qrrho_tot
        """
        return self.qrrho_total_entropy * self.T

    @property
    def gibbs_free_energy(self):
        """
        Formula:
            G = H - T * S_tot,c
        """
        return self.enthalpy - self.entropy_times_temperature

    @property
    def qrrho_total_internal_energy(self):
        """Obtain the quasi-RRHO total internal energy in J mol^-1.
        Formula:
            E^qrrho_tot = E_t + E_r + E^qrrho_v + E_e
        """
        if self.molecule.is_monoatomic:
            return (
                self.translational_internal_energy
                + self.electronic_internal_energy
            )
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
        return self.energies + self.qrrho_total_internal_energy + R * self.T

    @property
    def qrrho_gibbs_free_energy(self):
        """Obtain the Gibbs free energy in J mol^-1, by quasi-RRHO corrections to both entropy and enthalpy.
        Formula:
            G^qrrho_q = H^qrrho - T * S^qrrho_tot
        """
        return self.qrrho_enthalpy - self.qrrho_entropy_times_temperature

    @property
    def qrrho_gibbs_free_energy_qs(self):
        """Obtain the Gibbs free energy in J mol^-1, by a quasi-RRHO correction to entropy only.
        Formula:
            G^qrrho_qs = H - T * S^qrrho_tot
        """
        return self.enthalpy - self.qrrho_entropy_times_temperature

    @property
    def qrrho_gibbs_free_energy_qh(self):
        """Obtain the Gibbs free energy in J mol^-1, by a quasi-RRHO correction to enthalpy only.
        Formula:
            G^qrrho_qh = H^qrrho - T * S_tot,c
        """
        return self.qrrho_enthalpy - self.entropy_times_temperature
