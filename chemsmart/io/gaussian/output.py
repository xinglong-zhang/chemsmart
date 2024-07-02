import re
import logging
from itertools import islice
from functools import cached_property
from chemsmart.utils.mixins import FileMixin

logger = logging.getLogger(__name__)

# patterns for searching
eV_pattern = r"([\d\.]+) eV"
nm_pattern = r"([\d\.]+) nm"
f_pattern = r"f=([\d\.]+)"
float_pattern = r"[-]?\d*\.\d+|\d+"


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
