import logging
import re
from functools import cached_property

import numpy as np

from chemsmart.utils.mixins import FileMixin, XTBFileMixin
from chemsmart.utils.repattern import normal_mode_pattern

logger = logging.getLogger(__name__)


class XTBMainOut(XTBFileMixin):
    """
    Parse xTB main output file (*.out) containing calculation results.

    This is the primary output file from xTB calculations, containing:
    - Setup information (charge, multiplicity, basis functions, etc.)
    - SCF convergence data
    - Molecular properties (dipole, quadrupole moments)
    - Thermodynamic properties (if frequency calculation)
    - Optimization trajectory (if geometry optimization)
    - Timing information

    Args:
        filename (str): Path to the xTB output file (e.g., xtb.out, water_ohess.out)
    """

    def __init__(self, filename):
        self.filename = filename

    @property
    def xtb_version(self):
        """xtb version used for the calculation."""
        for line in self.contents:
            if "xtb version" in line:
                return line.split()[3]
        return None

    @property
    def normal_termination(self):
        """Check if the xtb output file is terminated by checking each output file line from the last line onwards."""
        for line in reversed(self.contents):
            if "[ERROR]" in line:
                logger.info(f"File {self.filename} has error termination.")
                return False
            if "* finished run" in line:
                logger.info(f"File {self.filename} terminated normally.")
                return True
        return False

    def _get_route(self):
        """
        Get parsed route information.
        """
        for line in self.contents:
            if "program call" in line:
                route = line.split(":")[-1].strip()
                return route
        return None

    #           ...................................................
    #           :                      SETUP                      :
    #           :.................................................:
    def _get_setup_information(self, keyword):
        """Get xtb setup information."""
        # restrict to SETUP block
        for i, line in enumerate(self.contents):
            if "SETUP" in line:
                for j_line in self.contents[i + 2 :]:
                    if len(j_line) == 0:
                        break
                    if keyword in j_line:
                        return (
                            j_line.split(keyword)[-1]
                            .strip()
                            .split()[0]
                            .strip()
                        )
        return None

    @property
    def num_basis_functions(self):
        num_basis_functions = self._get_setup_information("# basis functions")
        return int(num_basis_functions) if num_basis_functions else None

    @property
    def num_atomic_orbital(self):
        num_atomic_orbital = self._get_setup_information("# atomic orbitals")
        return int(num_atomic_orbital) if num_atomic_orbital else None

    @property
    def num_shells(self):
        num_shells = self._get_setup_information("# shells")
        return int(num_shells) if num_shells else None

    @property
    def num_electrons(self):
        num_electrons = self._get_setup_information("# electrons")
        return int(num_electrons) if num_electrons else None

    @property
    def max_iter(self):
        max_iter = self._get_setup_information("max. iterations")
        return int(max_iter) if max_iter else None

    @property
    def hamiltonian(self):
        return self._get_setup_information("Hamiltonian")

    @property
    def restart(self):
        restart = self._get_setup_information("restarted?")
        return restart.lower() == "true"

    @property
    def solvent_on(self):
        # instead of checking for GBSA solvation,
        # we will search for the actual solvent model and solvent ID
        solvation = self._get_setup_information("GBSA solvation")
        return solvation.lower() == "true"

    # solvation model and solvent ID
    @property
    def solvent_model(self):
        for line in self.contents:
            if "Solvation model:" in line:
                return line.split()[-1]
        return None

    @property
    def solvent_id(self):
        for line in self.contents:
            if "Solvent" in line and len(line.split()) == 2:
                return line.split()[-1]
        return None

    @property
    def pc_potential(self):
        """Point Charge Potential, an external electrostatic potential
        applied to the system using classical point charges"""
        pc_pot = self._get_setup_information("PC potential")
        return pc_pot.lower() == "true"

    @property
    def electronic_temperature(self):
        electronic_temp = self._get_setup_information("electronic temp.")
        return float(electronic_temp) if electronic_temp else None

    @property
    def accuracy(self):
        accuracy = self._get_setup_information("accuracy")
        return float(accuracy) if accuracy else None

    @property
    def integral_cutoff(self):
        integral_cutoff = self._get_setup_information("integral cutoff")
        return float(integral_cutoff) if integral_cutoff else None

    @property
    def integral_neglect(self):
        integral_neglect = self._get_setup_information("integral neglect")
        return float(integral_neglect) if integral_neglect else None

    @property
    def scf_convergence(self):
        scf_convergence = self._get_setup_information("SCF convergence")
        return float(scf_convergence) if scf_convergence else None

    @property
    def wf_convergence(self):
        wf_convergence = self._get_setup_information("wf. convergence")
        return float(wf_convergence) if wf_convergence else None

    @property
    def broyden_damping(self):
        broyden_damping = self._get_setup_information("Broyden damping")
        return float(broyden_damping) if broyden_damping else None

    @property
    def net_charge(self):
        net_charge = self._get_setup_information("net charge")
        return int(net_charge) if net_charge else None

    @property
    def unpaired_electrons(self):
        unpaired_electrons = self._get_setup_information("unpaired electrons")
        return int(unpaired_electrons) if unpaired_electrons else None

    @property
    def optimization_level(self):
        """Optimization level."""
        opt_level = self._get_setup_information("optimization level")
        return opt_level

    @property
    def max_optcycles(self):
        max_optcycles = self._get_setup_information("max. optcycles")
        return int(max_optcycles) if max_optcycles else None

    @property
    def anc_microcycles(self):
        """ANC micro-cycles: Advanced Newton–Raphson Convergence (ANC) micro-iterations
        used in the Self-Consistent Field (SCF) algorithm to improve convergence
        """
        anc_microcycles = self._get_setup_information("ANC micro-cycles")
        return int(anc_microcycles) if anc_microcycles else None

    @property
    def degrees_of_freedom(self):
        dof = self._get_setup_information("degrees of freedom")
        return int(dof) if dof else None

    @property
    def rf_solver(self):
        """In xtb, the RF solver refers to the Relaxed Fock (RF) solver,
        an alternative approach for solving the Self-Consistent Field (SCF) equations.
        It is designed to improve SCF convergence, particularly for challenging
        systems where standard methods struggle."""
        return self._get_setup_information("RF solver")

    @property
    def write_all_intermediate_geometries(self):
        write_all = self._get_setup_information("write xtbopt.log")
        if write_all:
            return write_all.lower() == "true"

    @property
    def is_linear(self):
        """Determine if molecule is linear."""
        linear = self._get_setup_information("linear?")
        if linear:
            return linear.lower() == "true"

    @property
    def energy_convergence(self):
        energy_conv = self._get_setup_information("energy convergence")
        if energy_conv:
            return float(energy_conv)

    @property
    def gradient_convergence(self):
        """Gradient convergence threshold, in Eh/alpha."""
        gradient_conv = self._get_setup_information("grad. convergence")
        if gradient_conv:
            return float(gradient_conv)

    @property
    def max_rf_displacement(self):
        """Maximum displacement in the Relaxed Fock (RF) solver."""
        max_rf_disp = self._get_setup_information("maxmium RF displ.")
        if max_rf_disp:
            return float(max_rf_disp)

    @property
    def low_frequency_cutoff(self):
        """Low frequency cutoff in cm^-1."""
        low_freq_cutoff = self._get_setup_information("Hlow (freq-cutoff)")
        if low_freq_cutoff:
            return float(low_freq_cutoff)

    @property
    def max_frequency_cutoff(self):
        """Maximum frequency cutoff in cm^-1."""
        max_freq_cutoff = self._get_setup_information("Hmax (freq-cutoff)")
        if max_freq_cutoff:
            return float(max_freq_cutoff)

    @property
    def s6_in_model_hessian(self):
        """S6 parameter in the model Hessian. S6 is a scaling factor used for
        dispersion correction in GFN-xTB. It is part of D3 dispersion correction
        and related to the dispersion energy scaling.
        S6 is a global scaling factor that determines strength of dispersion correction
        applied in the model Hessian calculation, to account for long-range van der Waals
        interactions.
        The model Hessian is an approximate way to compute vibrational frequencies in xtb.
        Since dispersion affects intermolecular forces and bond stiffness, S6 influences
        force constants used in Hessian construction.
        Proper dispersion scaling ensures realistic vibrational spectra and thermodynamic properties.
        """
        s6 = self._get_setup_information("S6 in model hess.")
        if s6:
            return float(s6)

    @property
    def homo_energy(self):
        """Obtain HOMO energy of last optimized structure, in eV."""
        for line in reversed(self.contents):
            if "(HOMO)" in line:
                homo_energy = line.split()[-2]
                return float(homo_energy)
        return None

    @property
    def lumo_energy(self):
        """Obtain LUMO energy of last optimized structure, in eV."""
        for line in reversed(self.contents):
            if "(LUMO)" in line:
                lumo_energy = line.split()[-2]  # lumo energy in eV
                return float(lumo_energy)
        return None

    @property
    def c6_coefficient(self):
        """C6 dispersion coefficient in au·bohr⁶."""
        for line in reversed(self.contents):
            if "Mol. C6AA" in line:
                c6_coefficient = line.split()[-1]
                return float(c6_coefficient)
        return None

    @property
    def c8_coefficient(self):
        """C8 dispersion coefficient in au·bohr⁸."""
        for line in reversed(self.contents):
            if "Mol. C8AA" in line:
                c8_coefficient = line.split()[-1]
                return float(c8_coefficient)
        return None

    @property
    def alpha_coefficient(self):
        """Alpha coefficient α(0)."""
        for line in reversed(self.contents):
            if "Mol. α(0)" in line:
                alpha_coefficient = line.split()[-1]
                return float(alpha_coefficient)
        return None

    def get_all_summary_blocks(self):
        """Obtain all SUMMARY blocks from the output file."""
        summary_blocks = []
        for i, line in enumerate(self.contents):
            if "SUMMARY" in line:
                summary_block = []
                for j_line in self.contents[i + 2 :]:
                    if len(j_line) == 0:
                        break
                    if "::::::::" in j_line or "........" in j_line:
                        continue
                    else:
                        summary_block.append(j_line)
                summary_blocks.append(summary_block)
        if len(summary_blocks) == 0:
            return None
        return summary_blocks

    @property
    def vertical_ionization_potential(self):
        """Vertical Ionization Potential (VIP) in eV, using command line '--vip', '--vipea' or '--vomega'."""
        for line in reversed(self.contents):
            if "delta SCC IP (eV)" in line:
                vertical_ionization_potentials = line.split()[-1]
                return float(vertical_ionization_potentials)
        return None

    @property
    def vertical_electron_affinity(self):
        """Vertical electron Affinities (EA) in eV, using command line '--ea', '--vipea' or '--vomega'."""
        for line in reversed(self.contents):
            if "delta SCC EA (eV)" in line:
                vertical_electron_affinities = line.split()[-1]
                return float(vertical_electron_affinities)
        return None

    @property
    def global_electrophilicity_index(self):
        """Global Electrophilicity Indexes (GEI) in eV, using command line '--vomega'."""
        for line in reversed(self.contents):
            if "Global electrophilicity index (eV):" in line:
                global_electrophilicity_index = line.split()[-1]
                return float(global_electrophilicity_index)
        return None

    @property
    def fukui_index(self):
        """Return Fukui Index block, using command line '--vfukui'."""
        fukui_block = []
        for i, line in enumerate(self.contents):
            if "Fukui functions:" in line:
                for j_line in self.contents[i + 1 :]:
                    if "------" in j_line:
                        break
                    fukui_block.append(j_line)
                return fukui_block
        return None

    def _extract_summary_information(self, keyword):
        # make sure it only extracts from SUMMARY block
        last_summary_block = self.get_all_summary_blocks()[-1]
        if last_summary_block:
            for line in last_summary_block:
                if keyword in line:
                    return float(line.split()[-3])
        return None

    @property
    def total_energy_without_gsasa_hb(self):
        """Total free energy (electronic + thermal) of system without contributions
        from GSASA (Generalized Solvent Accessible Surface Area) and hydrogen
        bonding (hb) corrections."""
        return self._extract_summary_information("total w/o Gsasa/hb")

    @property
    def scc_energy(self):
        """Electronic energy from self-consistent charge (SCC) calculation."""
        return self._extract_summary_information("SCC energy")

    @property
    def isotropic_es(self):
        """Coulombic interactions assuming a spherical charge distribution
        (simplified electrostatics)."""
        return self._extract_summary_information("-> isotropic ES")

    @property
    def anisotropic_es(self):
        """Electrostatic interactions considering directional charge effects
        (higher-order multipoles)."""
        return self._extract_summary_information("-> anisotropic ES")

    @property
    def anisotropic_xc(self):
        """Direction-dependent exchange-correlation energy contribution,
        which accounts for electronic interactions beyond spherical approximation.
        """
        return self._extract_summary_information("-> anisotropic XC")

    @property
    def dispersion_energy(self):
        """Dispersion (van der Waals) correction, capturing weak
        interactions between non-bonded atoms."""
        return self._extract_summary_information("-> dispersion")

    @property
    def solvation_energy_gsolv(self):
        """Energy correction due to solvation effects (implicit solvent model).
        Gsolv in xtb."""
        return self._extract_summary_information("-> Gsolv")

    @property
    def electronic_solvation_energy_gelec(self):
        """Electronic solvation energy contribution, Gelec in xtb."""
        return self._extract_summary_information("-> Gelec")

    @property
    def surface_area_solvation_energy_gsasa(self):
        """SASA solvation energy contribution, Gsasa in xtb."""
        return self._extract_summary_information("-> Gsasa")

    @property
    def hydrogen_bonding_solvation_energy_ghb(self):
        """Hydrogen bonding solvation correction."""
        return self._extract_summary_information("-> Ghb")

    @property
    def empirical_shift_correction_gshift(self):
        """Empirical shift correction, Gshift in xtb."""
        return self._extract_summary_information("-> Gshift")

    @property
    def repulsion_energy(self):
        """Pauli repulsion energy between overlapping electron densities."""
        return self._extract_summary_information("repulsion energy")

    @property
    def additional_restraining_energy(self):
        """Additional restraining energy."""
        return self._extract_summary_information("add. restraining")

    @property
    def total_charge(self):
        return self._extract_summary_information("total charge")

    # Numerical Hessian
    @cached_property
    def numerical_hessian_block(self):
        """Information under Numerical Hessian block."""
        for i, line in enumerate(self.contents):
            if "Numerical Hessian" in line:
                numerical_hessian_block = []
                for j_line in self.contents[i + 2 :]:
                    if len(j_line) == 0:
                        break
                    numerical_hessian_block.append(j_line)
                return numerical_hessian_block
        return None

    @property
    def numfreq(self):
        """If numerical hessian is turned on."""
        for line in self.contents:
            if "Numerical Hessian" in line:
                return True
        return False

    @property
    def hessian_step_length(self):
        """Finite displacement step size used to compute the numerical Hessian.'
        A smaller step gives more accurate results but can be computationally expensive.
        Default is 0.005 Å, which balances accuracy and efficiency."""
        if self.numerical_hessian_block:
            for line in self.numerical_hessian_block:
                if "step length" in line:
                    return float(line.split()[-1])
        return None

    @property
    def scc_accuracy(self):
        """SCC accuracy, controls the self-consistent charge (SCC) cycle convergence
        for the Hessian calculation. Lower values (e.g., 0.10000) indicate stricter
        convergence, higher values (e.g., 0.50000) speed up calculation but may
        introduce slight inaccuracies."""
        if self.numerical_hessian_block:
            for line in self.numerical_hessian_block:
                if "SCC accuracy" in line:
                    return float(line.split()[-1])
        return None

    @property
    def hessian_scale_factor(self):
        """Hessian scaling factor, applied to adjust vibrational frequencies."""
        if self.numerical_hessian_block:
            for line in self.numerical_hessian_block:
                if "Hessian scale factor" in line:
                    return float(line.split()[-1])
        return None

    @property
    def rms_gradient(self):
        """Root mean square (RMS) gradient, a measure of the convergence of the
        numerical Hessian calculation. Lower values indicate better convergence.
        Generally, values below 0.001 Eh/a₀ suggest a well-converged structure.
        """
        if self.numerical_hessian_block:
            for line in self.numerical_hessian_block:
                if "RMS gradient" in line:
                    return float(line.split()[-1])
        return None

    @property
    def molecular_dipole_lines(self):
        dipole_lines = []
        for i, line in enumerate(self.contents):
            if line.startswith("molecular dipole:"):
                for j_line in self.contents[i + 2 : i + 4]:
                    # only get two lines
                    dipole_lines.append(j_line)
        if len(dipole_lines) == 0:
            return None
        return dipole_lines

    @property
    def qonly_molecular_dipole(self):
        """Charge only dipole, computed only from atomic partial charges
        (electrostatic contribution)."""
        if self.molecular_dipole_lines is not None:
            for line in self.molecular_dipole_lines:
                if line.startswith("q only:"):
                    return np.array([float(x) for x in line.split()[-3:]])
        return None

    @property
    def full_molecular_dipole(self):
        """Actual dipole moment including both charge distribution
        and electronic polarization contributions"""
        if self.molecular_dipole_lines is not None:
            for line in self.molecular_dipole_lines:
                if line.startswith("full:"):
                    return np.array([float(x) for x in line.split()[1:4]])
        return None

    @property
    def total_molecular_dipole_moment(self):
        """Total molecular dipole moment, in Debye."""
        if self.molecular_dipole_lines is not None:
            for line in self.molecular_dipole_lines:
                if line.startswith("full:"):
                    return float(line.split()[-1])
        return None

    @property
    def molecular_quadrupole_lines(self):
        quadrupole_lines = []
        for i, line in enumerate(self.contents):
            if line.startswith("molecular quadrupole (traceless):"):
                for j_line in self.contents[i + 1 : i + 9]:
                    quadrupole_lines.append(j_line)
        if len(quadrupole_lines) == 0:
            return None
        return quadrupole_lines

    @property
    def qonly_molecular_quadrupole(self):
        """Charge-only quadrupole moment, computed from atomic partial charges."""
        if self.molecular_quadrupole_lines is not None:
            for line in self.molecular_quadrupole_lines:
                if line.startswith("q only:"):
                    xx, xy, yy, xz, yz, zz = [
                        float(x) for x in line.split()[-6:]
                    ]
                    # convert to tensor
                    return np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])
        return None

    @property
    def q_dip_molecular_quadrupole(self):
        """Molecular quadrupole from both charge and dipole moment contributions."""
        if self.molecular_quadrupole_lines is not None:
            for line in self.molecular_quadrupole_lines:
                if line.startswith("q+dip:"):
                    xx, xy, yy, xz, yz, zz = [
                        float(x) for x in line.split()[-6:]
                    ]
                    # convert to tensor
                    return np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])
        return None

    @property
    def full_molecular_quadrupole(self):
        """Full molecular quadrupole moment, including charge, dipole, and quadrupole
        contributions."""
        if self.molecular_quadrupole_lines is not None:
            for line in self.molecular_quadrupole_lines:
                if line.startswith("full:"):
                    xx, xy, yy, xz, yz, zz = [
                        float(x) for x in line.split()[-6:]
                    ]
                    # convert to tensor
                    return np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])
        return None

    @property
    def dielectric_constant(self):
        if self.solvent_on:
            for line in self.contents:
                if "Dielectric constant" in line:
                    return float(line.split()[-1])
        return None

    @property
    def free_energy_shift(self):
        if self.solvent_on:
            for line in self.contents:
                if "Free energy shift" in line:
                    return float(line.split()[-4])  # free energy shift in Eh
        return None

    @property
    def temperature(self):
        if self.solvent_on:
            for line in self.contents:
                if "Temperature" in line:
                    return float(line.split()[-2])  # temperature in K
        return None

    @property
    def density(self):
        if self.solvent_on:
            for line in self.contents:
                if "Density" in line:
                    return float(line.split()[-2])  # density in kg/L
        return None

    @property
    def solvent_mass(self):
        if self.solvent_on:
            for line in self.contents:
                if "Solvent mass" in line:
                    return float(line.split()[-2])  # solvent mass in g/mol
        return None

    @property
    def h_bond_correction(self):
        if self.solvent_on:
            for line in self.contents:
                if "H-bond correction" in line:
                    return line.split()[-1] == "true"
        return None

    @property
    def ion_screening(self):
        if self.solvent_on:
            for line in self.contents:
                if "Ion screening" in line:
                    return line.split()[-1] == "true"
        return None

    @property
    def surface_tension(self):
        if self.solvent_on:
            for line in self.contents:
                if "Surface tension" in line:
                    return float(line.split()[-4])  # surface tension in Eh
        return None

    """
    GEOMETRY OPTIMIZATION
    """

    @cached_property
    def geometry_optimization_converged(self):
        for line in self.contents:
            if "GEOMETRY OPTIMIZATION CONVERGED" in line:
                return True
            elif "FAILED TO CONVERGE GEOMETRY OPTIMIZATION" in line:
                return False
        return False

    @property
    def optimized_structure_block(self):
        """Return optimized structure."""
        if self.geometry_optimization_converged:
            coordinates_blocks = []
            for i, line in enumerate(self.contents):
                if "final structure:" in line:
                    for j_line in self.contents[i + 2 :]:
                        if "Bond Distances (Angstroems)" in j_line:
                            break
                        coordinates_blocks.append(j_line)
            return coordinates_blocks
        return None

    @property
    def molecular_mass(self):
        if not self.geometry_optimization_converged:
            return None
        for line in self.contents:
            if "molecular mass/u" in line:
                molecular_mass = line.split(":")[1].strip()
                return float(molecular_mass)
        return None

    @property
    def center_of_mass(self):
        if not self.geometry_optimization_converged:
            return None
        for line in self.contents:
            if "center of mass at/Å" in line:
                return [float(x) for x in line.split(":")[1].strip().split()]
        return None

    @property
    def moments_of_inertia(self):
        if not self.geometry_optimization_converged:
            return None
        for line in self.contents:
            if "moments of inertia/u·Å²" in line:
                return [float(x) for x in line.split(":")[1].strip().split()]
        return None

    @property
    def rotational_constants(self):
        if not self.geometry_optimization_converged:
            return None
        for line in self.contents:
            if "rotational constants/cm⁻¹" in line:
                return [float(x) for x in line.split(":")[1].strip().split()]
        return None

    @property
    def optimizer_wall_time(self):
        if self.elapsed_walltime_by_jobs("ANC optimizer:"):
            return round(
                sum(self.elapsed_walltime_by_jobs("ANC optimizer:")), 6
            )
        return None

    @property
    def optimizer_cpu_time(self):
        if self.cpu_runtime_by_jobs("ANC optimizer:"):
            return round(sum(self.cpu_runtime_by_jobs("ANC optimizer:")), 6)
        return None

    # """
    # CALCULATION OF VIBRATIONAL FREQUENCIES
    # """

    @cached_property
    def all_vibrational_frequencies(self):
        """Read the vibrational frequencies from the xTB output file.
        The first six (for non-linear molecules) or five (for linear molecules)
        frequencies correspond to translations (3x) or rotations (3x/2x) of the molecule.
        """
        found_frequency_printout = False
        for i, line in enumerate(self.contents):
            if "Frequency Printout" in line:
                found_frequency_printout = True
                continue
            if found_frequency_printout and "vibrational frequencies" in line:
                frequencies = []
                for j_line in self.contents[i + 1 :]:
                    if "reduced masses (amu)" in j_line:
                        break
                    freq_line = j_line.split(":")[-1].strip().split()
                    for freq in freq_line:
                        frequencies.append(float(freq))
                return frequencies
        return None

    @property
    def vibrational_frequencies(self):
        """Return vibrational frequencies without translational and rotational modes."""
        if self.all_vibrational_frequencies is None:
            return []
        return [x for x in self.all_vibrational_frequencies if x != 0.0]

    @cached_property
    def reduced_masses(self):
        """Obtain list of reduced masses corresponding to the vibrational frequency."""
        for i, line in enumerate(self.contents):
            if line.startswith("reduced masses (amu)"):
                reduced_masses = []
                for j_line in self.contents[i + 1 :]:
                    if "IR intensities (km·mol⁻¹)" in j_line:
                        break
                    reduced_mass_line = j_line.split()[1::2]
                    for reduced_mass in reduced_mass_line:
                        reduced_masses.append(float(reduced_mass))
                return reduced_masses[-self.num_vib_frequencies :]
        return None

    @cached_property
    def ir_intensities(self):
        """Obtain list of IR intensities corresponding to the vibrational frequency."""
        for i, line in enumerate(self.contents):
            if line.startswith("IR intensities ("):
                ir_intensities = []
                for j_line in self.contents[i + 1 :]:
                    if "Raman intensities" in j_line:
                        break
                    ir_intensity_line = j_line.split()[1::2]
                    for ir_intensity in ir_intensity_line:
                        ir_intensities.append(float(ir_intensity))
                return ir_intensities[-self.num_vib_frequencies :]
        return None

    @cached_property
    def raman_intensities(self):
        """Obtain list of Raman intensities corresponding to the vibrational frequency."""
        for i, line in enumerate(self.contents):
            if line.startswith("Raman intensities ("):
                raman_intensities = []
                for j_line in self.contents[i + 1 :]:
                    if "output can be read by thermo" in j_line:
                        break
                    raman_intensity_line = j_line.split()[1::2]
                    for raman_intensity in raman_intensity_line:
                        try:
                            raman_intensities.append(float(raman_intensity))
                        except ValueError:
                            pass
                return raman_intensities[-self.num_vib_frequencies :]
        return None

    # Hessian SETUP information
    @property
    def num_vib_frequencies(self):
        num_freq = self._get_setup_information("# frequencies")
        return int(num_freq)

    @property
    def num_imaginary_frequencies(self):
        num_im_freq = self._get_setup_information("# imaginary freq.")
        if num_im_freq:
            return int(num_im_freq)

    @property
    def only_rot_calc(self):
        """compute only rotational contributions to Hessian,
        rather than the full vibrational analysis."""
        only_rot = self._get_setup_information("only rotational calc.")
        if only_rot:
            return only_rot.lower() == "true"

    @property
    def symmetry(self):
        """Molecular symmetry."""
        return self._get_setup_information("symmetry")

    @property
    def rotational_symmetry_number(self):
        rot_num = self._get_setup_information("rotational number")
        if rot_num:
            return int(rot_num)

    @property
    def scaling_factor(self):
        """Scaling factor used for vibrational frequencies."""
        scale_factor = self._get_setup_information("scaling factor")
        if scale_factor:
            return float(scale_factor)

    @property
    def rotor_cutoff(self):
        """Defines threshold below which low-frequency vibrational modes are treated
        as rotational modes (free internal rotations). Defaults to 50 cm^-1."""
        rotor_cut = self._get_setup_information("rotor cutoff")
        if rotor_cut:
            return float(rotor_cut)

    @property
    def imaginary_frequency_cutoff(self):
        """Imaginary frequency cutoff in cm^-1. Defaults to 20 cm^-1."""
        im_freq_cutoff = self._get_setup_information("imag. cutoff")
        if im_freq_cutoff:
            return float(im_freq_cutoff)

    # ^^ Hessian SETUP information

    @property
    def partition_function(self):
        partition_function = {}
        for line in self.contents:
            if "VIB" in line:
                partition_function["vibrational"] = float(line.split()[2])
            elif "ROT" in line:
                partition_function["rotational"] = float(line.split()[1])
            elif "INT" in line:
                partition_function["internal"] = float(line.split()[1])
            elif "TR" in line:
                partition_function["translational"] = float(line.split()[1])
        if partition_function:
            return partition_function
        return None

    def _get_thermodynamics_block(self):
        """Get thermodynamics block."""
        for i, line in enumerate(self.contents):
            if "THERMODYNAMIC" in line:
                thermodynamics_block = []
                for j_line in self.contents[i + 2 :]:
                    if len(j_line) == 0:
                        break
                    if "::::::::" in j_line or "........" in j_line:
                        continue
                    thermodynamics_block.append(j_line)
                return thermodynamics_block
        return None

    def _extract_thermodynamics_information(self, keyword):
        """Extract thermodynamic information from the output file."""
        thermodynamics_block = self._get_thermodynamics_block()
        if thermodynamics_block:
            for line in thermodynamics_block:
                if keyword in line:
                    return float(line.split()[-3])
        return None

    @property
    def zero_point_energy(self):
        """Zero point energy in Eh"""
        return self._extract_thermodynamics_information("zero point energy")

    @property
    def grrho_without_zpve(self):
        """Free energy in Eh within the rigid-rotor-harmonic-oscillator (RRHO) approximation,
        excluding zero-point vibrational energy (ZPVE)"""
        return self._extract_thermodynamics_information("G(RRHO) w/o ZPVE")

    @property
    def grrho_contribution(self):
        """Contribution of RRHO approximation to free energy in Eh"""
        return self._extract_thermodynamics_information("G(RRHO) contrib.")

    @cached_property
    def energies(self):
        """
        Return energies of the system from xTB output file.
        """
        return self._get_energies()

    def _get_energies(self):
        """
        Obtain a list of energies for each geometry optimization point.
        """
        energies = []
        for i, line in enumerate(self.contents):
            if "* total energy  :" in line:
                energy = float(line.split(":")[1].split()[0].strip())
                energies.append(energy)
        return energies or [self.total_energy]

    def _extract_final_information(self, keyword):
        for line in reversed(self.contents):
            if keyword in line:
                return float(
                    line.split(keyword)[-1].strip().split()[0].strip()
                )
        return None

    @property
    def total_energy(self):
        return self._extract_final_information(
            "TOTAL ENERGY"
        )  # total energy in Eh

    @property
    def gradient_norm(self):
        return self._extract_final_information(
            "GRADIENT NORM"
        )  # gradient norm in Eh/α

    @property
    def fmo_gap(self):
        """Extract HOMO-LUMO gap, in eV"""
        fmo_gap = self._extract_final_information("HOMO-LUMO GAP")
        assert np.isclose(
            self.lumo_energy - self.homo_energy, fmo_gap, atol=1e-3
        )
        return fmo_gap

    @property
    def total_enthalpy(self):
        """Total enthalpy in Eh"""
        return self._extract_final_information("TOTAL ENTHALPY")

    @property
    def total_free_energy(self):
        """Total free energy in Eh"""
        return self._extract_final_information("TOTAL FREE ENERGY")

    def sum_time_hours(self, line):
        n_days = float(line.split(" d,")[0].split()[-1])
        n_hours = float(line.split(" h,")[0].split()[-1])
        n_minutes = float(line.split(" min,")[0].split()[-1])
        n_seconds = float(line.split(" sec")[0].split()[-1])
        total_seconds = (
            n_days * 24 * 60 * 60
            + n_hours * 60 * 60
            + n_minutes * 60
            + n_seconds
        )
        total_hours = round(total_seconds / 3600, 6)
        return total_hours

    def elapsed_walltime_by_jobs(self, task_name):
        elapsed_walltime = []
        for i, line in enumerate(self.contents):
            if task_name in self.contents[i - 1] and "wall-time:" in line:
                total_hours = self.sum_time_hours(line)
                elapsed_walltime.append(total_hours)
        if elapsed_walltime:
            return elapsed_walltime
        return None

    def cpu_runtime_by_jobs(self, task_name):
        cpu_runtime = []
        for i, line in enumerate(self.contents):
            if task_name in self.contents[i - 2] and "cpu-time:" in line:
                total_hours = self.sum_time_hours(line)
                cpu_runtime.append(total_hours)
        if cpu_runtime:
            return cpu_runtime
        return None

    @property
    def total_wall_time(self):
        if self.elapsed_walltime_by_jobs("total:"):
            return round(sum(self.elapsed_walltime_by_jobs("total:")), 6)
        return None

    @property
    def total_cpu_time(self):
        if self.cpu_runtime_by_jobs("total:"):
            return round(sum(self.cpu_runtime_by_jobs("total:")), 6)
        return None

    @property
    def scf_wall_time(self):
        if self.elapsed_walltime_by_jobs("SCF:"):
            return round(sum(self.elapsed_walltime_by_jobs("SCF:")), 6)
        return None

    @property
    def scf_cpu_time(self):
        if self.cpu_runtime_by_jobs("SCF:"):
            return round(sum(self.cpu_runtime_by_jobs("SCF:")), 6)
        return None

    @property
    def hessian_wall_time(self):
        if self.elapsed_walltime_by_jobs("analytical hessian:"):
            return round(
                sum(self.elapsed_walltime_by_jobs("analytical hessian:")), 6
            )
        return None

    @property
    def hessian_cpu_time(self):
        if self.cpu_runtime_by_jobs("analytical hessian:"):
            return round(
                sum(self.cpu_runtime_by_jobs("analytical hessian:")), 6
            )
        return None


class XTBChargesFile(FileMixin):
    """
    Parse xTB charges file containing atomic partial charges.

    File format::
        -0.56472698
         0.28236349
         0.28236349

    Each line contains one floating point number representing the partial charge
    of an atom, in the order atoms appear in the structure.

    Args:
        filename (str): Path to the charges file
    """

    def __init__(self, filename):
        self.filename = filename

    @cached_property
    def partial_charges(self):
        """Get atomic partial charges from the file."""
        charges = []
        for line in self.contents:
            line = line.strip()
            if not line:
                continue
            try:
                charges.append(float(line))
            except ValueError:
                continue
        return charges if charges else None

    @property
    def total_charge(self):
        """Calculate total molecular charge from partial charges."""
        if self.partial_charges is not None:
            return sum(self.partial_charges)
        return None


class XTBEnergyFile(FileMixin):
    """
    Parse xTB energy file containing total energies from calculation.

    File format::

        $energy
             1    -5.07054444346    -5.07054444346    -5.07054444346
        $end

    The three energy values are identical and correspond to the final total
    energy in Hartree.

    Args:
        filename (str): Path to the energy file
    """

    def __init__(self, filename):
        self.filename = filename

    @cached_property
    def last_energy(self):
        """Final total energy in Hartree."""
        for line in self.contents:
            if line.strip().startswith("$") or not line.strip():
                continue
            return float(line.split()[1])
        return None


class XTBEngradFile(FileMixin):
    """
    Parse xTB energy and gradient file (.engrad).

    File format::

        #
        # Number of atoms
        #
                 3
        #
        # The current total energy in Eh
        #
             -5.070544443465
        #
        # The current gradient in Eh/bohr
        #
               0.000000000130
              -0.000000000000
               0.000057137032
              -0.000019065816
              ...
        #
        # The atomic numbers and current coordinates in Bohr
        #
           8    -0.0000025   -0.0000013   -0.7167751
           1     1.4592625   -0.0000070    0.3583826
           1    -1.4592600    0.0000083    0.3583928

    The file contains:
    - Number of atoms
    - Total energy in Hartree (Eh)
    - Gradient components (3N values) in Eh/bohr
    - Atomic numbers and Cartesian coordinates in Bohr

    Args:
        filename (str): Path to the .engrad file
    """

    def __init__(self, filename):
        self.filename = filename

    @cached_property
    def num_atoms(self):
        """Obtain the number of atoms in the system."""
        for i, line in enumerate(self.contents):
            if "Number of atoms" in line:
                # check following lines for the number
                for content in self.contents[i + 1 : i + 4]:
                    try:
                        return int(content.split()[0])
                    except (ValueError, IndexError):
                        pass
        return None

    @cached_property
    def total_energy(self):
        """Obtain the total energy in Hartree."""
        for i, line in enumerate(self.contents):
            if "current total energy" in line:
                # check following lines for the energy value
                for content in self.contents[i + 1 : i + 4]:
                    try:
                        energy_in_hartree = float(content.split()[0])
                        return energy_in_hartree
                    except (ValueError, IndexError):
                        pass
        return None

    @cached_property
    def forces(self):
        """Obtain forces on the atoms in Hartree/Bohr."""
        grad = self._get_gradient()
        if grad is not None and grad.shape == (self.num_atoms, 3):
            # Convert gradient to forces: F = -grad
            forces = -grad
            return [forces]
        return None

    def _get_gradient(self):
        """Get the gradient in Hartree/Bohr."""
        if self.num_atoms is None:
            return None
        for i, line in enumerate(self.contents):
            if "current gradient" in line:
                # Read 3N gradient components
                grad_data = []
                for content in self.contents[
                    i + 1 : i + 3 * self.num_atoms + 4
                ]:
                    try:
                        grad_value = float(
                            content.split()[0]
                        )  # in Hartree/Bohr
                        grad_data.append(grad_value)
                    except (ValueError, IndexError):
                        pass
                # Validate we have the correct number of gradient components
                if len(grad_data) == 3 * self.num_atoms:
                    gradient_array = np.array(grad_data).reshape(
                        (self.num_atoms, 3)
                    )
                    return gradient_array
        return None


class XTBG98File(FileMixin):
    """
    Parse xTB Gaussian 98 format vibrational analysis output (g98.out).

    This file contains vibrational frequencies and normal modes in Gaussian 98
    format, compatible with visualization programs like GaussView, Molden, etc.

    File format::

        Entering Gaussian System
        *********************************************
        Gaussian 98:
        frequency output generated by the xtb code
        *********************************************
                             Standard orientation:
        --------------------------------------------------------------------
         Center     Atomic     Atomic              Coordinates (Angstroms)
         Number     Number      Type              X           Y           Z
        --------------------------------------------------------------------
           1          8             0       -0.000001   -0.000001   -0.379301
           2          1             0        0.772209   -0.000004    0.189648
           3          1             0       -0.772207    0.000004    0.189653
        --------------------------------------------------------------------
        ...
         Harmonic frequencies (cm**-1), IR intensities (km*mol⁻¹),
         Raman scattering activities (A**4/amu), Raman depolarization ratios,
         reduced masses (AMU), force constants (mDyne/A) and normal coordinates:
                             1                      2                      3
                               a                      a                      a
        Frequencies --  1539.3017              3643.5149              3651.7072
        Red. masses --     2.1457                 1.5477                 2.1398
        Frc consts  --     0.0000                 0.0000                 0.0000
        IR Inten    --   133.2648                 6.7667                16.6406
        Raman Activ --     0.0000                 0.0000                 0.0000
        Depolar     --     0.0000                 0.0000                 0.0000
        Atom AN      X      Y      Z        X      Y      Z        X      Y      Z
          1   8     0.00   0.00   0.28    -0.00  -0.00  -0.19     0.27  -0.00  -0.00
          2   1     0.40  -0.00  -0.55     0.58  -0.00   0.38    -0.55   0.00  -0.40
          3   1    -0.40   0.00  -0.55    -0.58   0.00   0.38    -0.55   0.00   0.40

    Args:
        filename (str): Path to the g98.out file
    """

    def __init__(self, filename):
        self.filename = filename

    @cached_property
    def vibrational_frequencies(self):
        """Obtain list of vibrational frequencies in cm^-1."""
        frequencies = []
        for line in self.contents:
            if line.startswith("Frequencies --"):
                freq_string = line.split("--")[1].strip()
                for freq in freq_string.split():
                    frequencies.append(float(freq))
        return frequencies

    @cached_property
    def reduced_masses(self):
        """Obtain list of reduced masses corresponding to the vibrational frequency, in amu."""
        reduced_masses = []
        for line in self.contents:
            if line.startswith("Red. masses --"):
                reduced_masses_string = line.split("--")[1].strip()
                for mass in reduced_masses_string.split():
                    reduced_masses.append(float(mass))
        return reduced_masses

    @cached_property
    def force_constants(self):
        """Obtain list of force constants corresponding to the vibrational frequency, in mDyne/Å."""
        force_constants = []
        for line in self.contents:
            if line.startswith("Frc consts  --"):
                force_constants_string = line.split("--")[1].strip()
                for force in force_constants_string.split():
                    force_constants.append(float(force))
        return force_constants

    @cached_property
    def ir_intensities(self):
        """Obtain list of IR intensities corresponding to the vibrational frequency, in km/mol."""
        ir_intensities = []
        for line in self.contents:
            if line.startswith("IR Inten    --"):
                ir_intensities_string = line.split("--")[1].strip()
                for intensity in ir_intensities_string.split():
                    ir_intensities.append(float(intensity))
        return ir_intensities

    @cached_property
    def raman_activities(self):
        """Obtain list of Raman activities corresponding to the vibrational frequency, in A^4/amu."""
        raman_activities = []
        for line in self.contents:
            if line.startswith("Raman Activ --"):
                raman_string = line.split("--")[1].strip()
                for activity in raman_string.split():
                    raman_activities.append(float(activity))
        return raman_activities

    @cached_property
    def depolarization_ratios(self):
        """Obtain list of Raman depolarization ratios corresponding to the vibrational frequency."""
        depolar_ratios = []
        for line in self.contents:
            if line.startswith("Depolar     --"):
                depolar_string = line.split("--")[1].strip()
                for ratio in depolar_string.split():
                    depolar_ratios.append(float(ratio))
        return depolar_ratios

    @cached_property
    def vibrational_mode_symmetries(self):
        """Obtain list of vibrational mode symmetries corresponding to the vibrational frequency."""
        vibrational_mode_symmetries = []
        for i, line in enumerate(self.contents):
            if line.startswith("Frequencies --"):
                # go back one line to get the symmetries
                symmetries = self.contents[i - 1].split()
                for sym in symmetries:
                    vibrational_mode_symmetries.append(sym)
        return vibrational_mode_symmetries

    @cached_property
    def vibrational_modes(self):
        """
        Obtain list of vibrational normal modes corresponding
        to the vibrational frequency. Returns a list of normal modes,
        each of num_atoms x 3 (in dx, dy, and dz for each element)
        vibration.
        """
        list_of_vib_modes = []
        for i, line in enumerate(self.contents):
            if line.startswith("Frequencies --"):
                first_col_vib_modes = []
                second_col_vib_modes = []
                third_col_vib_modes = []
                for j_line in self.contents[i + 7 :]:
                    # if line match normal mode pattern
                    if re.match(normal_mode_pattern, j_line):
                        normal_mode = [float(val) for val in j_line.split()]
                        first_col_vib_mode = normal_mode[2:5]
                        second_col_vib_mode = normal_mode[5:8]
                        third_col_vib_mode = normal_mode[8:11]
                        first_col_vib_modes.append(first_col_vib_mode)
                        if second_col_vib_mode:
                            second_col_vib_modes.append(second_col_vib_mode)
                        if third_col_vib_mode:
                            third_col_vib_modes.append(third_col_vib_mode)
                    else:
                        break
                if first_col_vib_modes:
                    list_of_vib_modes.append(np.array(first_col_vib_modes))
                if second_col_vib_modes:
                    list_of_vib_modes.append(np.array(second_col_vib_modes))
                if third_col_vib_modes:
                    list_of_vib_modes.append(np.array(third_col_vib_modes))
        return list_of_vib_modes

    @cached_property
    def num_vib_modes(self):
        """Number of vibrational modes found."""
        return len(self.vibrational_modes)

    @cached_property
    def num_vib_frequencies(self):
        """Number of vibrational frequencies found."""
        return len(self.vibrational_frequencies)

    def _attach_vib_metadata(self, mol):
        """Attach vibrational data to a Molecule object as attributes."""
        vib = {
            "frequencies": self.vibrational_frequencies or [],
            "reduced_masses": self.reduced_masses or [],
            "force_constants": self.force_constants or [],
            "ir_intensities": self.ir_intensities or [],
            "mode_symmetries": self.vibrational_mode_symmetries or [],
            "modes": self.vibrational_modes or [],
        }

        setattr(mol, "vibrational_frequencies", vib["frequencies"])
        setattr(mol, "vibrational_reduced_masses", vib["reduced_masses"])
        setattr(mol, "vibrational_force_constants", vib["force_constants"])
        setattr(mol, "vibrational_ir_intensities", vib["ir_intensities"])
        setattr(mol, "vibrational_mode_symmetries", vib["mode_symmetries"])
        setattr(mol, "vibrational_modes", vib["modes"])

        return mol


class XTBGradientFile(FileMixin):
    """
    Parse XTB gradient file containing energy gradient information.

    File format:
        $grad
          cycle =      1    SCF energy =    -5.07054444346   |dE/dxyz| =  0.000075
           -0.00000250190431     -0.00000125099553     -0.71677514611431      O
            1.45926248160647     -0.00000701517601      0.35838257371508      H
           -1.45925997970150      0.00000826617187      0.35839276137184      H
           1.2982149851656E-10  -3.3293838441823E-18   5.7137032470948E-05
          -1.9065815509550E-05  -1.9010909790127E-17  -2.8568564065665E-05
           1.9065685688053E-05   2.2340293634309E-17  -2.8568468405277E-05
        $end

    Contains coordinates and gradient for each atom.

    Args:
        filename (str): Path to the gradient file
    """

    def __init__(self, filename):
        self.filename = filename

    # TODO: Add parsing methods for gradient.


class XTBHessianFile(FileMixin):
    """
    Parse xTB Hessian file containing the Cartesian Hessian matrix.

    File format::

        $hessian
                0.6095802844  -0.0000031918  -0.0000007369  -0.3047915657   0.0000012039
               -0.2245646512  -0.3047887187   0.0000019878   0.2245653880
               -0.0000031918   0.0000000000   0.0000006951   0.0000012763  -0.0000000000
                0.0000008283   0.0000019155  -0.0000000000  -0.0000015234
               ...

    The Hessian matrix (second derivatives of energy with respect to Cartesian
    coordinates) is written as a continuous stream of floating point numbers
    in Hartree/Bohr². The matrix is symmetric with size (3N × 3N) where N is
    the number of atoms.

    Args:
        filename (str): Path to the hessian file
    """

    def __init__(self, filename):
        self.filename = filename

    # TODO: Add parsing methods for hessian.


class XTBVibSpectrumFile(FileMixin):
    """
    Parse xTB vibrational spectrum file containing frequencies and IR intensities.

    File format::

        $vibrational spectrum
        #  mode     symmetry     wave number   IR intensity    selection rules
        #                         cm**(-1)      (km*mol⁻¹)        IR
             1                      -0.00         0.00000          -
             2                      -0.00         0.00000          -
             3                       0.00         0.00000          -
             4                       0.00         0.00000          -
             5                       0.00         0.00000          -
             6                       0.00         0.00000          -
             7        a           1539.30       133.26477         YES
             8        a           3643.51         6.76672         YES
             9        a           3651.71        16.64059         YES
        $end

    The first 6 modes (5 for linear molecules) are translations and rotations
    with zero frequency. Real vibrational modes have positive frequencies.

    Args:
        filename (str): Path to the vibspectrum file
    """

    def __init__(self, filename):
        self.filename = filename

    # TODO: Add parsing methods for vibrational spectrum.


class XTBWibergBondOrderFile(FileMixin):
    """
    Parse XTB Wiberg bond order (wbo) file.

    File format::
        1           2  0.92021379026732564
        1           3  0.92021379039282269

    Each line contains:
        atom1_index  atom2_index  bond_order_value

    Args:
        filename (str): Path to the wbo file
    """

    def __init__(self, filename):
        self.filename = filename

    # TODO: Add parsing methods for bond order.
