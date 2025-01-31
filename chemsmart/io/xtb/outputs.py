import logging
import re
from functools import cached_property
from chemsmart.utils.mixins import XTBFileMixin

logger = logging.getLogger(__name__)

class XTBOutput(XTBFileMixin):
    def __init__(self, filename):
        self.filename = filename

    @cached_property
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
            if "finished run" in line:
                logger.info(f"File {self.filename} terminated normally.")
                return True
        return False

    @property
    def route_string(self):
        for line in self.contents:
            if "program call" in line:
                return line.split(":")[-1].strip()
        return None

    @property
    def num_basis_functions(self):
        for line in self.contents:
            if "# basis functions" in line:
                num_basis_functions = line.split()[-2]
                return int(num_basis_functions)
        return None

    @property
    def num_atomic_orbital(self):
        for line in self.contents:
            if "# atomic orbitals" in line:
                num_atomic_orbital = line.split()[-2]
                return int(num_atomic_orbital)
        return None

    @property
    def num_shells(self):
        for line in self.contents:
            if "# shells" in line:
                num_shells = line.split()[-2]
                return int(num_shells)
        return None

    @property
    def num_electrons(self):
        for line in self.contents:
            if "# electrons" in line:
                num_electrons = line.split()[-2]
                return int(num_electrons)
        return None

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
                    logger.info('Free energy shift in Eh')
                    solvation_info["free_energy_shift"] = float(line.split()[-4])
                elif "Temperature" in line:
                    logger.info('Temperature in K')
                    solvation_info["temperature"] = float(line.split()[-2])
                elif "Density" in line:
                    logger.info('Density in kg/L')
                    solvation_info["density"] = float(line.split()[-2])
                elif "Solvent mass" in line:
                    logger.info('Solvent mass in g/mol')
                    solvation_info["solvent_mass"] = float(line.split()[-2])
                elif "H-bond correction" in line:
                    solvation_info["H_bond_correction"] = line.split()[-1] == "true"
                elif "Ion screening" in line:
                    solvation_info["ion_screening"] = line.split()[-1] == "true"
                elif "Surface tension" in line:
                    logger.info('Surface tension in Eh')
                    solvation_info["surface_tension"] = float(line.split()[-4])
            if solvation_info:
                return solvation_info
        return None

    @property
    def net_charge(self):
        for line in self.contents:
            if "net charge" in line:
                net_charge = line.split()[-2]
                return int(net_charge)
        return None

    @property
    def unpaired_electrons(self):
        for line in self.contents:
            if "unpaired electrons" in line:
                unpaired_electrons = line.split()[-2]
                return int(unpaired_electrons)
        return None

    @property
    def total_energy(self):
        for line in self.contents:
            if "TOTAL ENERGY" in line:
                total_energy = line.split()[-3]
                logger.info('Total energy in Eh')
                return float(total_energy)
        return None

    @property
    def fmo_gap(self):
        for line in self.contents:
            if "HOMO-LUMO GAP" in line:
                fmo_gap = line.split()[-3]
                logger.info('homo-lumo gap in eV')
                return float(fmo_gap)
        return None

    @property
    def homo_energy(self):
        for line in self.contents:
            if "(HOMO)" in line:
                homo_energy = line.split()[-2]
                logger.info('homo energy in eV')
                return float(homo_energy)
        return None

    @property
    def lumo_energy(self):
        for line in self.contents:
            if "(LUMO)" in line:
                lumo_energy = line.split()[-2]
                logger.info('lumo energy in eV')
                return float(lumo_energy)
        return None

    @property
    def total_charge(self):
        for line in self.contents:
            if "total charge" in line:
                total_charge = line.split()[-3]
                return float(total_charge)
        return None

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
        if self.elapsed_walltime_by_jobs("total"):
            return round(sum(self.elapsed_walltime_by_jobs("total")), 6)
        return None

    @property
    def total_cpu_time(self):
        if self.cpu_runtime_by_jobs("total"):
            return round(sum(self.cpu_runtime_by_jobs("total")), 6)
        return None

    @property
    def scf_wall_time(self):
        if self.elapsed_walltime_by_jobs("SCF"):
            return round(sum(self.elapsed_walltime_by_jobs("SCF")), 6)
        return None

    @property
    def scf_cpu_time(self):
        if self.cpu_runtime_by_jobs("SCF"):
            return round(sum(self.cpu_runtime_by_jobs("SCF")), 6)
        return None

    @property
    def optimizer_wall_time(self):
        if self.elapsed_walltime_by_jobs("optimizer"):
            return round(sum(self.elapsed_walltime_by_jobs("optimizer")), 6)
        return None

    @property
    def optimizer_cpu_time(self):
        if self.cpu_runtime_by_jobs("optimizer"):
            return round(sum(self.cpu_runtime_by_jobs("optimizer")), 6)
        return None