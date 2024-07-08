import re
import logging
from itertools import islice
from functools import cached_property
from chemsmart.utils.mixins import FileMixin
from chemsmart.utils.repattern import eV_pattern, nm_pattern, f_pattern, float_pattern

logger = logging.getLogger(__name__)

class Gaussian16Output(FileMixin):
    def __init__(self, filename):
        self.filename = filename

    @property
    def normal_termination(self):
        """Check for termination of gaussian file by checking the last line of the output file."""
        contents = self.contents
        if len(contents) == 0:
            return False

        last_line = contents[-1]
        if "Normal termination of Gaussian" in last_line:
            logger.info(f"File {self.filename} terminated normally.")
            return True

        logger.info(f"File {self.filename} has error termination.")
        return False

    @property
    def num_atoms(self):
        """Number of atoms in the molecule."""
        for line in self.contents:
            if line.startswith("NAtoms="):
                return int(line.split()[1])

    @property
    def charge(self):
        for line in self.contents:
            if "Charge" in line and "Multiplicity" in line:
                line_elem = line.split()
                return int(line_elem[2])

    @property
    def multiplicity(self):
        for line in self.contents:
            if "Charge" in line and "Multiplicity" in line:
                line_elem = line.split()
                return int(line_elem[-1])

    @property
    def spin(self):
        """Spin restricted vs spin unrestricted calculations."""
        for line in self.contents:
            if "SCF Done:" in line:
                line_elem = line.split("E(")
                theory = line_elem[1].split(")")[0].lower()
                # determine if the job is restricted or unrestricted
                if theory.startswith("r"):
                    spin = "restricted"
                elif theory.startswith("u"):
                    spin = "unrestricted"
                else:
                    spin = None
                return spin

    @property
    def num_basis_functions(self):
        for line in self.contents:
            if "basis functions," in line:
                line_elem = line.split(",")
                num_basis_functions = line_elem[0].strip().split()[0]
                return int(num_basis_functions)
        return None

    @property
    def num_primitive_gaussians(self):
        for line in self.contents:
            if "primitive gaussians," in line:
                line_elem = line.split(",")
                primitives = line_elem[1].strip().split()[0]
                return int(primitives)
        return None

    @property
    def num_cartesian_basis_functions(self):
        for line in self.contents:
            if (
                "cartesian basis functions" in line
                and "basis functions," in line
                and "primitive gaussians," in line
            ):
                line_elem = line.split(",")
                num_cartesian_basis_functions = line_elem[2].strip().split()[0]
                return int(num_cartesian_basis_functions)
        return None

    # Below gives computing time/resources used
    @property
    def cpu_runtime_by_jobs_core_hours(self):
        cpu_runtime = []
        for line in self.contents:
            if line.startswith("Job cpu time:"):
                n_days = float(line.split("days")[0].strip().split()[-1])
                n_hours = float(line.split("hours")[0].strip().split()[-1])
                n_minutes = float(line.split("minutes")[0].strip().split()[-1])
                n_seconds = float(line.split("seconds")[0].strip().split()[-1])
                total_seconds = (
                    n_days * 24 * 60 * 60
                    + n_hours * 60 * 60
                    + n_minutes * 60
                    + n_seconds
                )
                total_hours = round(total_seconds / 3600, 4)  # round to 1 nearest hour
                cpu_runtime.append(total_hours)
        return cpu_runtime

    @property
    def service_units_by_jobs(self):
        """SUs defined as the JOB CPU time in hours."""
        return self.cpu_runtime_by_jobs_core_hours

    @property
    def total_core_hours(self):
        return round(sum(self.cpu_runtime_by_jobs_core_hours), 4)

    @property
    def total_service_unit(self):
        return self.total_core_hours

    @property
    def elapsed_walltime_by_jobs(self):
        elapsed_walltime = []
        for line in self.contents:
            if line.startswith("Elapsed time:"):
                n_days = float(line.split("days")[0].strip().split()[-1])
                n_hours = float(line.split("hours")[0].strip().split()[-1])
                n_minutes = float(line.split("minutes")[0].strip().split()[-1])
                n_seconds = float(line.split("seconds")[0].strip().split()[-1])
                total_seconds = (
                    n_days * 24 * 60 * 60
                    + n_hours * 60 * 60
                    + n_minutes * 60
                    + n_seconds
                )
                total_hours = round(total_seconds / 3600, 4)
                elapsed_walltime.append(total_hours)
        return elapsed_walltime

    @property
    def total_elapsed_walltime(self):
        return round(sum(self.elapsed_walltime_by_jobs), 4)

    @cached_property
    def tddft_transitions(self):
        """
        Read a excitation energies after a TD-DFT calculation.

        Returns:
            A list: A list of tuple for each transition such as
                    [(energie (eV), lambda (nm), oscillatory strength), ... ]
        """
        tddft_transitions = []
        for line in self.contents:
            if line.startswith("Excited State"):
                eV_match = re.search(eV_pattern, line)
                nm_match = re.search(nm_pattern, line)
                f_match = re.search(f_pattern, line)
                if eV_match and nm_match and f_match:
                    # Extract and convert the matched values to float
                    excitation_energy_eV = float(eV_match.group(1))
                    absorption_wavelength = float(nm_match.group(1))
                    oscillatory_strength = float(f_match.group(1))

                    tddft_transitions.append(
                        (
                            excitation_energy_eV,
                            absorption_wavelength,
                            oscillatory_strength,
                        )
                    )
        return tddft_transitions

    @cached_property
    def excitation_energies_eV(self):
        """
        Read TDDFT transitions and return the transition energies in eV as a list.
        """
        excitation_energies_eV = []
        for i in self.tddft_transitions:
            excitation_energies_eV.append(i[0])
        return excitation_energies_eV

    @cached_property
    def absorptions_in_nm(self):
        """Read TDDFT transitions and return the absorbed wavelengths in nm as a list."""
        absorptions_in_nm = []
        for i in self.tddft_transitions:
            absorptions_in_nm.append(i[1])
        return absorptions_in_nm

    @cached_property
    def oscillatory_strengths(self):
        """Read TDDFT transitions and return the oscillatory strengths as a list."""
        oscillatory_strengths = []
        for i in self.tddft_transitions:
            oscillatory_strengths.append(i[2])
        return oscillatory_strengths

    @cached_property
    def transitions(self):
        """Read TDDFT transitions and return the MO transitions."""
        transitions, _ = self._read_transitions_and_contribution_coefficients()
        return transitions

    @cached_property
    def contribution_coefficients(self):
        """Read MO contribution coefficients."""
        _, cc = self._read_transitions_and_contribution_coefficients()
        return cc

    def _read_transitions_and_contribution_coefficients(self):
        transitions = []
        contribution_coefficients = []
        for i, line in enumerate(self.contents):
            if line.startswith("Excited State"):
                each_state_transitions = []
                each_state_contribution_coefficients = []
                # parse the lines that follow until an empty line is encountered
                j = 1
                while len(self.contents[i + j]) != 0:
                    line_element = self.contents[i + j].split()
                    if len(line_element) <= 4:
                        mo_transition = " ".join(
                            list(islice(line_element, len(line_element) - 1))
                        )
                        contribution_coefficient = float(line_element[-1])
                        each_state_transitions.append(mo_transition)
                        each_state_contribution_coefficients.append(
                            contribution_coefficient
                        )
                    j += 1
                transitions.append(each_state_transitions)
                contribution_coefficients.append(each_state_contribution_coefficients)

        return transitions, contribution_coefficients

    @cached_property
    def alpha_occ_eigenvalues(self):
        """Obtain all eigenenergies of the alpha occuplied orbitals."""
        alpha_occ_eigenvalues = []

        # Iterate through lines in reverse to find the last block of eigenvalues
        eigenvalue_blocks = []
        current_block = []
        found_first_block = False

        for line in reversed(self.contents):
            if line.startswith("Alpha  occ. eigenvalues"):
                # Add the line to the current block
                current_block.append(line)
                found_first_block = True
            elif found_first_block:
                # We've reached the end of the last block
                eigenvalue_blocks.append(current_block)
                current_block = []
                found_first_block = False

        if eigenvalue_blocks:
            # Extract the last block and process it
            last_block = eigenvalue_blocks[0]
            last_block.reverse()  # Reverse to original order

            # Flatten the last block and convert to list of floats
            last_block_values = []
            for line in last_block:
                # Find all floats in the line, including those without spaces
                values = re.findall(float_pattern, line)
                last_block_values.extend(map(float, values))
            return last_block_values

    @cached_property
    def alpha_virtual_eigenvalues(self):
        """Obtain all eigenenergies of the alpha unoccuplied orbitals."""

        # Iterate through lines in reverse to find the last block of eigenvalues
        eigenvalue_blocks = []
        current_block = []
        found_first_block = False

        for line in reversed(self.contents):
            if line.startswith("Alpha virt. eigenvalues"):
                # Add the line to the current block
                current_block.append(line)
                found_first_block = True
            elif found_first_block:
                # We've reached the end of the last block
                eigenvalue_blocks.append(current_block)
                current_block = []
                found_first_block = False

        if eigenvalue_blocks:
            # Extract the last block and process it
            last_block = eigenvalue_blocks[0]
            last_block.reverse()  # Reverse to original order

            # Flatten the last block and convert to list of floats
            last_block_values = []
            for line in last_block:
                # Find all floats in the line, including those without spaces
                values = re.findall(float_pattern, line)
                last_block_values.extend(map(float, values))

            # print(len(eigenvalue_blocks))  # number of eigenvalue blocks in the file

            return last_block_values

    @cached_property
    def beta_occ_eigenvalues(self):
        """Obtain all eigenenergies of the beta occuplied orbitals."""
        # Iterate through lines in reverse to find the last block of eigenvalues
        eigenvalue_blocks = []
        current_block = []
        found_first_block = False

        for line in reversed(self.contents):
            if line.startswith("Beta  occ. eigenvalues"):
                # Add the line to the current block
                current_block.append(line)
                found_first_block = True
            elif found_first_block:
                # We've reached the end of the last block
                eigenvalue_blocks.append(current_block)
                current_block = []
                found_first_block = False

        if eigenvalue_blocks:
            # Extract the last block and process it
            last_block = eigenvalue_blocks[0]
            last_block.reverse()  # Reverse to original order

            # Flatten the last block and convert to list of floats
            last_block_values = []
            for line in last_block:
                # Find all floats in the line, including those without spaces
                values = re.findall(float_pattern, line)
                last_block_values.extend(map(float, values))

            return last_block_values

    @cached_property
    def beta_virtual_eigenvalues(self):
        """Obtain all eigenenergies of the beta unoccuplied orbitals."""

        # Iterate through lines in reverse to find the last block of eigenvalues
        eigenvalue_blocks = []
        current_block = []
        found_first_block = False

        for line in reversed(self.contents):
            if line.startswith("Beta virt. eigenvalues"):
                # Add the line to the current block
                current_block.append(line)
                found_first_block = True
            elif found_first_block:
                # We've reached the end of the last block
                eigenvalue_blocks.append(current_block)
                current_block = []
                found_first_block = False

        if eigenvalue_blocks:
            # Extract the last block and process it
            last_block = eigenvalue_blocks[0]
            last_block.reverse()  # Reverse to original order

            # Flatten the last block and convert to list of floats
            last_block_values = []
            for line in last_block:
                # Find all floats in the line, including those without spaces
                values = re.findall(float_pattern, line)
                last_block_values.extend(map(float, values))

            return last_block_values

    @cached_property
    def homo_energy(self):
        if self.multiplicity == 1:
            assert (
                self.beta_occ_eigenvalues is None
                and self.beta_virtual_eigenvalues is None
            )
            return self.alpha_occ_eigenvalues[-1]

    @cached_property
    def num_unpaired_electrons(self):
        if self.multiplicity != 1:
            # the multiplicity is the number of unpaired electrons + 1
            assert (
                len(self.alpha_occ_eigenvalues) - len(self.beta_occ_eigenvalues) + 1
                == self.multiplicity
            )
            return len(self.alpha_occ_eigenvalues) - len(self.beta_occ_eigenvalues)

    @cached_property
    def somo_energy(self):
        if self.multiplicity != 1:
            # the multiplicity is the number of unpaired electrons + 1
            assert (
                len(self.alpha_occ_eigenvalues) - len(self.beta_occ_eigenvalues) + 1
                == self.multiplicity
            )
            return self.alpha_occ_eigenvalues[-1]

    @cached_property
    def lumo_energy(self):
        if self.multiplicity == 1:
            assert (
                self.beta_occ_eigenvalues is None
                and self.beta_virtual_eigenvalues is None
            )
            return self.alpha_virtual_eigenvalues[0]

    @cached_property
    def fmo_gap(self):
        if self.multiplicity == 1:
            return self.homo_energy - self.lumo_energy
        else:
            # to implement for radical systems
            pass
