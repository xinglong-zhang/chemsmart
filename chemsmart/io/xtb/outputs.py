import logging
import re
from functools import cached_property
from chemsmart.io.molecules.structure import CoordinateBlock
from chemsmart.utils.mixins import XTBFileMixin
import numpy as np

logger = logging.getLogger(__name__)

class XTBOutput(XTBFileMixin):
    def __init__(self, filename):
        self.filename = filename

    @property
    def contents(self):
        with open(self.filename) as f:
            return [line.strip() for line in f.readlines()]

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

    def _extract_setup_information(self, keyword):
        for line in self.contents:
            if keyword in line:
                return int(line.split()[-2])
        return None

    @property
    def num_basis_functions(self):
        return self._extract_setup_information("# basis functions")

    @property
    def num_atomic_orbital(self):
        return self._extract_setup_information("# atomic orbitals")

    @property
    def num_shells(self):
        return self._extract_setup_information("# shells")

    @property
    def num_electrons(self):
        return self._extract_setup_information("# electrons")

    @property
    def solvation(self):
        for line in self.contents:
            if "GBSA solvation" in line:
                line_elements = line.split()
                if "true" in line_elements:
                    return True
        return False

    @property
    def solvation_info(self):
        solvation_info = {}
        if self.solvation:
            lines = self.contents
            for i, line in enumerate(lines):
                if "* Solvation model:" in line:
                    solvation_info["solvation_model"] = line.split(":")[-1].strip()
                elif "Solvent" in line:
                    if "* Solvation model:" in lines[i - 1]:
                        solvation_info["solvent"] = line.split()[-1]
                elif "Dielectric constant" in line:
                    solvation_info["dielectric_constant"] = float(line.split()[-1])
                elif "Free energy shift" in line:
                    """Free energy shift in Eh"""
                    solvation_info["free_energy_shift"] = float(line.split()[-4])
                elif "Temperature" in line:
                    """Temperature in K"""
                    solvation_info["temperature"] = float(line.split()[-2])
                elif "Density" in line:
                    """Density in kg/L"""
                    solvation_info["density"] = float(line.split()[-2])
                elif "Solvent mass" in line:
                    """Solvent mass in g/mol"""
                    solvation_info["solvent_mass"] = float(line.split()[-2])
                elif "H-bond correction" in line:
                    solvation_info["H_bond_correction"] = line.split()[-1] == "true"
                elif "Ion screening" in line:
                    solvation_info["ion_screening"] = line.split()[-1] == "true"
                elif "Surface tension" in line:
                    """Surface tension in Eh"""
                    solvation_info["surface_tension"] = float(line.split()[-4])
            if solvation_info:
                return solvation_info
        return None

    @property
    def net_charge(self):
        return self._extract_setup_information("net charge")

    @property
    def unpaired_electrons(self):
        return self._extract_setup_information("unpaired electrons")

    @property
    def homo_energy(self):
        for line in self.contents:
            if "(HOMO)" in line:
                homo_energy = line.split()[-2]
                """homo energy in eV"""
                return float(homo_energy)
        return None

    @property
    def lumo_energy(self):
        for line in self.contents:
            if "(LUMO)" in line:
                lumo_energy = line.split()[-2]
                """lumo energy in eV"""
                return float(lumo_energy)
        return None

    def _extract_summary_information(self, keyword):
        for line in reversed(self.contents):
            if keyword in line:
                return float(line.split()[-3])
        return None

    @property
    def scc_energy(self):
        return self._extract_summary_information("SCC energy")

    @property
    def isotropic_es(self):
        return self._extract_summary_information("-> isotropic ES")

    @property
    def anisotropic_es(self):
        return self._extract_summary_information("-> anisotropic ES")

    @property
    def anisotropic_xc(self):
        return self._extract_summary_information("-> anisotropic XC")

    @property
    def dispersion(self):
        return self._extract_summary_information("-> dispersion")

    @property
    def gsolv(self):
        return self._extract_summary_information("-> Gsolv")

    @property
    def gelec(self):
        return self._extract_summary_information("-> Gelec")

    @property
    def gsasa(self):
        return self._extract_summary_information("-> Gsasa")

    @property
    def ghb(self):
        return self._extract_summary_information("-> Ghb")

    @property
    def gshift(self):
        return self._extract_summary_information("-> Gshift")

    @property
    def repulsion_energy(self):
        return self._extract_summary_information("repulsion energy")

    @property
    def total_charge(self):
        return self._extract_summary_information("total charge")

    @property
    def molecular_dipole(self):
        dipole_lines = []
        for i, line in enumerate(self.contents):
            if line.startswith("molecular dipole:"):
                for j_line in self.contents[i + 2: i + 4]:
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
                    total_dipole = self.contents[i + 3].split()[-1]
                    """Total dipole in Debye"""
                    return float(total_dipole)
        return None

    @property
    def molecular_quadrupole(self):
        quadrupole_lines = []
        for i, line in enumerate(self.contents):
            if line.startswith("molecular quadrupole (traceless):"):
                for j_line in self.contents[i + 2: i + 5]:
                    quadrupole_lines.append(j_line.split(":")[1].strip().split())
        if len(quadrupole_lines) == 0:
            return None
        quadrupole_data = {
            "q_only": [float(x) for x in quadrupole_lines[0][0:6]],
            "q+dip": [float(x) for x in quadrupole_lines[1][0:6]],
            "full": [float(x) for x in quadrupole_lines[2][0:6]],
        }
        return quadrupole_data

    @property
    def total_energy(self):
        """Total energy in Eh"""
        return self._extract_summary_information("TOTAL ENERGY")

    @property
    def gradient_norm(self):
        """Gradient norm in Eh/α"""
        return self._extract_summary_information("GRADIENT NORM")

    @property
    def fmo_gap(self):
        """Homo-Lumo gap in eV"""
        return self._extract_summary_information("HOMO-LUMO GAP")

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


    """
    GEOMETRY OPTIMIZATION
    """
    @cached_property
    def geometry_optimization_convergence(self):
        for line in self.contents:
            if "GEOMETRY OPTIMIZATION CONVERGED" in line:
                return True
            elif "FAILED TO CONVERGE GEOMETRY OPTIMIZATION" in line:
                return False
        return False

    @property
    def degrees_of_freedom(self):
        return self._extract_setup_information("degrees of freedom")

    @property
    def optimized_structure_block(self):
        """Return optimized structure."""
        if self.geometry_optimization_convergence:
            coordinates_blocks = []
            for i, line in enumerate(self.contents):
                if "final structure:" in line:
                    for j_line in self.contents[i + 2:]:
                        if "Bond Distances (Angstroems)" in j_line:
                            break
                        coordinates_blocks.append(j_line)
            return coordinates_blocks
        return None

    @property
    def molecular_mass(self):
        if not self.geometry_optimization_convergence:
            return None
        for line in self.contents:
            if "molecular mass/u" in line:
                molecular_mass = line.split(":")[1].strip()
                return float(molecular_mass)
        return None

    @property
    def center_of_mass(self):
        if not self.geometry_optimization_convergence:
            return None
        for line in self.contents:
            if "center of mass at/Å" in line:
                return [float(x) for x in line.split(":")[1].strip().split()]
        return None

    @property
    def moments_of_inertia(self):
        if not self.geometry_optimization_convergence:
            return None
        for line in self.contents:
            if "moments of inertia/u·Å²" in line:
                return [float(x) for x in line.split(":")[1].strip().split()]
        return None

    @property
    def rotational_constants(self):
        if not self.geometry_optimization_convergence:
            return None
        for line in self.contents:
            if "rotational constants/cm⁻¹" in line:
                return [float(x) for x in line.split(":")[1].strip().split()]
        return None

    @property
    def optimizer_wall_time(self):
        if self.elapsed_walltime_by_jobs("ANC optimizer:"):
            return round(sum(self.elapsed_walltime_by_jobs("ANC optimizer:")), 6)
        return None

    @property
    def optimizer_cpu_time(self):
        if self.cpu_runtime_by_jobs("ANC optimizer:"):
            return round(sum(self.cpu_runtime_by_jobs("ANC optimizer:")), 6)
        return None


    """
    CALCULATION OF VIBRATIONAL FREQUENCIES
    """
    @cached_property
    def vibrational_frequencies(self):
        """Read the vibrational frequencies from the XTB output file.
        The first six frequencies correspond to the rotations and translations of the molecule."""
        found_frequency_printout = False
        for i, line in enumerate(self.contents):
            if "Frequency Printout" in line:
                found_frequency_printout = True
                continue
            if found_frequency_printout and line.startswith("projected vibrational frequencies (cm⁻¹)"):
                frequencies = []
                for j_line in self.contents[i + 1:]:
                    if "reduced masses (amu)" in j_line:
                        break
                    freq_line = j_line.split(":")[1].strip().split()
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
                for j_line in self.contents[i + 1:]:
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
            if line.startswith("IR intensities (km·mol⁻¹)"):
                ir_intensities = []
                for j_line in self.contents[i + 1:]:
                    if "Raman intensities (Ä⁴*amu⁻¹)" in j_line:
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
            if line.startswith("Raman intensities (Ä⁴*amu⁻¹)"):
                raman_intensities = []
                for j_line in self.contents[i + 1:]:
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

    @property
    def num_frequencies(self):
        return self._extract_setup_information("# frequencies")

    @property
    def num_imaginary_frequencies(self):
        return self._extract_setup_information("# imaginary freq.")

    @property
    def symmetry(self):
        for line in self.contents:
            if line.startswith(":"):
                if "symmetry" in line:
                    return line.split()[-2]
        return None

    @property
    def rotational_number(self):
        return self._extract_setup_information("rotational number")

    @property
    def scaling_factor(self):
        for line in self.contents:
            if line.startswith(":"):
                if "scaling factor" in line:
                    return float(line.split()[-2])
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
            return round(sum(self.elapsed_walltime_by_jobs("analytical hessian:")), 6)
        return None

    @property
    def hessian_cpu_time(self):
        if self.cpu_runtime_by_jobs("analytical hessian:"):
            return round(sum(self.cpu_runtime_by_jobs("analytical hessian:")), 6)
        return None