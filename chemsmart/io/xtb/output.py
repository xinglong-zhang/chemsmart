import logging
from functools import cached_property

import numpy as np

from chemsmart.utils.mixins import FileMixin

logger = logging.getLogger(__name__)


class XTBOutput(FileMixin):
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
        for line in self.contents:
            if "program call" in line:
                return line.split(":")[-1].strip()
        return None

    @property
    def route_string(self):
        return self._get_route()

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
        return int(num_basis_functions)

    @property
    def num_atomic_orbital(self):
        num_atomic_orbital = self._get_setup_information("# atomic orbitals")
        return int(num_atomic_orbital)

    @property
    def num_shells(self):
        num_shells = self._get_setup_information("# shells")
        return int(num_shells)

    @property
    def num_electrons(self):
        num_electrons = self._get_setup_information("# electrons")
        return int(num_electrons)

    @property
    def max_iter(self):
        max_iter = self._get_setup_information("max. iterations")
        return int(max_iter)

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
        return float(electronic_temp)

    @property
    def accuracy(self):
        accuracy = self._get_setup_information("accuracy")
        return float(accuracy)

    @property
    def integral_cutoff(self):
        integral_cutoff = self._get_setup_information("integral cutoff")
        return float(integral_cutoff)

    @property
    def integral_neglect(self):
        integral_neglect = self._get_setup_information("integral neglect")
        return float(integral_neglect)

    @property
    def scf_convergence(self):
        scf_convergence = self._get_setup_information("SCF convergence")
        return float(scf_convergence)

    @property
    def wf_convergence(self):
        wf_convergence = self._get_setup_information("wf. convergence")
        return float(wf_convergence)

    @property
    def broyden_damping(self):
        broyden_damping = self._get_setup_information("Broyden damping")
        return float(broyden_damping)

    @property
    def net_charge(self):
        net_charge = self._get_setup_information("net charge")
        return int(net_charge)

    @property
    def unpaired_electrons(self):
        unpaired_electrons = self._get_setup_information("unpaired electrons")
        return int(unpaired_electrons)

    @property
    def optimization_level(self):
        """Optimization level."""
        opt_level = self._get_setup_information("optimization level")
        return opt_level

    @property
    def max_optcycles(self):
        max_optcycles = self._get_setup_information("max. optcycles")
        if max_optcycles:
            return int(max_optcycles)

    @property
    def anc_microcycles(self):
        """ANC micro-cycles: Advanced Newton–Raphson Convergence (ANC) micro-iterations
        used in the Self-Consistent Field (SCF) algorithm to improve convergence
        """
        anc_microcycles = self._get_setup_information("ANC micro-cycles")
        if anc_microcycles:
            return int(anc_microcycles)

    @property
    def degrees_of_freedom(self):
        dof = self._get_setup_information("degrees of freedom")
        if dof:
            return int(dof)

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
        return summary_blocks

    def _extract_summary_information(self, keyword):
        # make sure it only extracts from SUMMARY block
        last_summary_block = self.get_all_summary_blocks()[-1]
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
        which accounts for electronic interactions beyond spherical approximation."""
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
        for line in self.numerical_hessian_block:
            if "Step length" in line:
                return float(line.split()[-1])

    @property
    def scc_accuracy(self):
        """SCC accuracy, controls the self-consistent charge (SCC) cycle convergence
        for the Hessian calculation. Lower values (e.g., 0.10000) indicate stricter
        convergence, higher values (e.g., 0.50000) speed up calculation but may
        introduce slight inaccuracies."""
        for line in self.numerical_hessian_block:
            if "SCC accuracy" in line:
                return float(line.split()[-1])

    @property
    def hessian_scale_factor(self):
        """Hessian scaling factor, applied to adjust vibrational frequencies."""
        for line in self.numerical_hessian_block:
            if "Hessian scale factor" in line:
                return float(line.split()[-1])

    @property
    def rms_gradient(self):
        """Root mean square (RMS) gradient, a measure of the convergence of the
        numerical Hessian calculation. Lower values indicate better convergence.
        Generally, values below 0.001 Eh/a₀ suggest a well-converged structure."""
        for line in self.numerical_hessian_block:
            if "RMS gradient" in line:
                return float(line.split()[-1])

    @property
    def molecular_dipole(self):
        dipole_lines = []
        for i, line in enumerate(self.contents):
            if line.startswith("molecular dipole:"):
                for j_line in self.contents[i + 2 : i + 4]:
                    dipole_lines.append(j_line.split(":")[1].strip().split())
        if len(dipole_lines) == 0:
            return None
        dipole_data = {
            "q_only": [float(x) for x in dipole_lines[0][0:3]],
            "full": [float(x) for x in dipole_lines[1][0:3]],
        }
        return dipole_data

    @property
    def total_dipole(self):
        for i, line in enumerate(self.contents):
            if line.startswith("molecular dipole:"):
                if "full:" in self.contents[i + 3]:
                    total_dipole = self.contents[i + 3].split()[
                        -1
                    ]  # total dipole in Debye
                    return float(total_dipole)
        return None

    @property
    def molecular_quadrupole(self):
        quadrupole_lines = []
        for i, line in enumerate(self.contents):
            if line.startswith("molecular quadrupole (traceless):"):
                for j_line in self.contents[i + 2 : i + 5]:
                    quadrupole_lines.append(
                        j_line.split(":")[1].strip().split()
                    )
        if len(quadrupole_lines) == 0:
            return None
        quadrupole_data = {
            "q_only": [float(x) for x in quadrupole_lines[0][0:6]],
            "q+dip": [float(x) for x in quadrupole_lines[1][0:6]],
            "full": [float(x) for x in quadrupole_lines[2][0:6]],
        }
        return quadrupole_data

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

    #"""
    #CALCULATION OF VIBRATIONAL FREQUENCIES
    #"""

    @cached_property
    def vibrational_frequencies(self):
        """Read the vibrational frequencies from the XTB output file.
        The first six (for non-linear molecules) or five (for linear molecules)
        frequencies correspond to translations (3x) or rotations (3x/2x) of the molecule.
        """
        found_frequency_printout = False
        for i, line in enumerate(self.contents):
            if "Frequency Printout" in line:
                found_frequency_printout = True
                continue
            if found_frequency_printout and line.startswith(
                "projected vibrational frequencies"
            ):
                frequencies = []
                for j_line in self.contents[i + 1 :]:
                    if "reduced masses (amu)" in j_line:
                        break
                    freq_line = j_line.split(":")[-1].strip().split()
                    for freq in freq_line:
                        frequencies.append(float(freq))
                return frequencies
        return None

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
                return reduced_masses
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
                return ir_intensities
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
                return raman_intensities
        return None

    # Hessian SETUP information
    @property
    def num_frequencies(self):
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

    @property
    def zero_point_energy(self):
        """Zero point energy in Eh"""
        return self._extract_summary_information("zero point energy")

    @property
    def grrho_without_zpve(self):
        """Free energy in Eh within the rigid-rotor-harmonic-oscillator (RRHO) approximation,
        excluding zero-point vibrational energy (ZPVE)"""
        return self._extract_summary_information("G(RRHO) w/o ZPVE")

    @property
    def grrho_contribution(self):
        """Contribution of RRHO approximation to free energy in Eh"""
        return self._extract_summary_information("G(RRHO) contrib.")

    @property
    def total_enthalpy(self):
        """Total enthalpy in Eh"""
        return self._extract_summary_information("TOTAL ENTHALPY")

    @property
    def total_free_energy(self):
        """Total free energy in Eh"""
        return self._extract_summary_information("TOTAL FREE ENERGY")

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
