import logging
import re
from functools import cached_property
from ase import units
from chemsmart.utils.mixins import FileMixin

class XTBOutput(FileMixin):
    def __init__(self, filename):
        self.filename = filename

    @property
    def normal_termination(self):
        """Check if the xtb output file is terminated by checking each output file line from the last line onwards."""
        for line in self.contents[::-1]:
            if "[ERROR]" in line:
                return False
            if "finished run" in line:
                return True
        return False

    @cached_property
    def optimized_output_lines(self):
        """Chunk of output file where the properties are calculated based on the final optimized structure.

        Returns empty for SP calculations.
        """
        optimized_output_lines = []
        for i, line in enumerate(self.contents):
            if "GEOMETRY OPTIMIZATION CONVERGED" in line:
                optimized_output_lines = list(self.contents[i:])
                break
        return optimized_output_lines

    @property
    def num_basis_functions(self):
        for line in self.contents:
            if "# basis functions" in line:
                line_elements = line.strip().split()
                num_basis_functions = line_elements[-2]
                return int(num_basis_functions)
        return None

    @property
    def num_atomic_orbital(self):
        for line in self.contents:
            if "# atomic orbitals" in line:
                line_elements = line.strip().split()
                num_atomic_orbital = line_elements[-2]
                return int(num_atomic_orbital)
        return None

    @property
    def num_shells(self):
        for line in self.contents:
            if "# shells" in line:
                line_elements = line.strip().split()
                num_shells = line_elements[-2]
                return int(num_shells)
        return None

    @property
    def num_electrons(self):
        for line in self.contents:
            if "# electrons" in line:
                line_elements = line.strip().split()
                num_electrons = line_elements[-2]
                return int(num_electrons)
        return None

    @property
    def get_total_energy(self):
        for line in self.contents:
            if "TOTAL ENERGY" in line:
                line_elements = line.strip().split()
                total_energy = line_elements[-3]
                return int(total_energy)
        return None

    @property
    def fmo_gap(self):
        for line in self.contents:
            if "HOMO-LUMO GAP" in line:
                line_elements = line.strip().split()
                fmo_gap = line_elements[-3]
                return int(fmo_gap)
        return None

    def sum_time_hours(self, line):
        n_days = float(line.split("d,")[0].split()[-1])
        n_hours = float(line.split("h,")[0].split()[-1])
        n_minutes = float(line.split("min,")[0].split()[-1])
        n_seconds = float(line.split("sec")[0].split()[-1])
        total_seconds = (
            n_days * 24 * 60 * 60
            + n_hours * 60 * 60
            + n_minutes * 60
            + n_seconds
        )
        total_hours = round(total_seconds / 3600, 4)
        return total_hours

    def elapsed_walltime_by_jobs(self, task_name):
        elapsed_walltime = []
        for i, line in enumerate(self.contents):
            if task_name in self.contents[i-1] and "wall-time:" in line:
                total_hours = self.sum_time_hours(line)
                elapsed_walltime.append(total_hours)
        return elapsed_walltime

    def cpu_runtime_by_jobs(self, task_name):
        cpu_runtime = []
        for i, line in enumerate(self.contents):
            if task_name in self.contents[i-1] and "cpu-time:" in line:
                total_hours = self.sum_time_hours(line)
                cpu_runtime.append(total_hours)
        return cpu_runtime


    @property
    def total_wall_time(self):
        return self.elapsed_walltime_by_jobs("total")

    @property
    def total_cpu_time(self):
        return self.cpu_runtime_by_jobs("total")

    @property
    def scf_wall_time(self):
        return self.elapsed_walltime_by_jobs("scf")

    @property
    def scf_cpu_time(self):
        return self.cpu_runtime_by_jobs("scf")

    @property
    def optimizer_wall_time(self):
        return self.elapsed_walltime_by_jobs("optimizer")

    @property
    def optimizer_cpu_time(self):
        return self.cpu_runtime_by_jobs("optimizer")