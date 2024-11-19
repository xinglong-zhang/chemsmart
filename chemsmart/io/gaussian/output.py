import re
import logging
from itertools import islice
from functools import cached_property
from chemsmart.utils.mixins import FileMixin

# patterns for searching
eV_pattern = r"([\d\.]+) eV"
nm_pattern = r"([\d\.]+) nm"
f_pattern = r"f=([\d\.]+)"

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
                total_hours = round(
                    total_seconds / 3600, 4
                )  # round to 1 nearest hour
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


class Gaussian16TDDFTOutput(Gaussian16Output):

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


class Gaussian16WBIOutput(Gaussian16Output):
    def __init__(self, filename):
        super().__init__(filename)

    @property
    def nbo_version(self):
        for line in self.contents:
            if "Gaussian NBO Version" in line:
                return line.split()[-1].split("*")[0]

    @cached_property
    def natural_atomic_orbitals(self):
        """Parse the NBO natural atomic orbitals."""
        nao = {}
        for i, line in enumerate(self.contents):
            if "NAO  Atom  No  lang   Type(AO)    Occupancy      Energy" in line:
                for j_line in self.contents[i + 2 :]:
                    if (
                        "WARNING" in j_line
                        or "Summary of Natural Population Analysis" in j_line
                    ):
                        break
                    if len(j_line) != 0:
                        columns = j_line.split()

                        # Extract values from each column
                        nao_number = int(columns[0])  # NAO Number (like 1, 2, etc.)
                        atom_type = columns[1]  # Atom type (e.g., 'Ni')
                        atom_number = columns[2]  # Atom number (e.g., '1')
                        lang = columns[3]  # Lang (e.g., 'S', 'px', 'py')
                        electron_type = columns[4].split("(")[
                            0
                        ]  # Electron type (e.g., 'Cor', 'Val', 'Ryd')
                        nao_type = columns[5].split(")")[
                            0
                        ]  # NAO Type (e.g., '1S', '2S', etc.)
                        nao_type += lang[1:]  # Append the lang to the NAO type
                        occupancy = float(columns[6])  # Occupancy
                        energy = float(columns[7])  # Energy

                        # Construct the atom key, e.g., "Ni1"
                        atom_key = f"{atom_type}{atom_number}"
                        # Construct the sub-key for each NAO entry, e.g., "NAO1", "NAO2", etc.
                        nao_key = f"NAO_{atom_type}{nao_number}"

                        # Initialize the atom dictionary if it doesn't exist
                        if atom_key not in nao:
                            nao[atom_key] = {}

                        # Populate the nested dictionary for each NAO entry
                        nao[atom_key][nao_key] = {
                            "nao_type": nao_type,
                            "electron_type": electron_type,
                            "occupancy": occupancy,
                            "energy": energy,
                        }
        return nao

    @cached_property
    def natural_population_analysis(self):
        """Parse the NBO natural population analysis."""
        npa = {}
        for i, line in enumerate(self.contents):
            if (
                "Atom  No    Charge         Core      Valence    Rydberg      Total"
                in line
            ):
                for j_line in self.contents[i + 2 :]:
                    if (
                        "======================================================================="
                        in j_line
                    ):
                        break
                    if len(j_line) != 0:
                        columns = j_line.split()

                        # Extract values from each column
                        atom_type = columns[0]
                        atom_number = columns[1]
                        charge = float(columns[2])
                        core = float(columns[3])
                        valence = float(columns[4])
                        rydberg = float(columns[5])
                        total = float(columns[6])

                        # Construct the atom key, e.g., "Ni1"
                        atom_key = f"{atom_type}{atom_number}"

                        # Initialize the atom dictionary if it doesn't exist
                        if atom_key not in npa:
                            npa[atom_key] = {}

                        # Populate the nested dictionary for each NAO entry
                        npa[atom_key] = {
                            "natural_charge": charge,
                            "core_electrons": core,
                            "valence_electrons": valence,
                            "rydberg_electrons": rydberg,
                            "total_electrons": total,
                        }
        return npa

    @cached_property
    def natural_charges(self):
        """Get natural charges corresponding to each atom as a dictionary."""
        natural_charges = {}
        for atom_key, atom_data in self.natural_population_analysis.items():
            natural_charges[atom_key] = atom_data["natural_charge"]
        return natural_charges

    @cached_property
    def total_electrons(self):
        """Get the total number of electrons corresponding to each atom as a dictionary."""
        total_electrons = {}
        for atom_key, atom_data in self.natural_population_analysis.items():
            total_electrons[atom_key] = atom_data["total_electrons"]
        return total_electrons

    @cached_property
    def electronic_configuration(self):
        """Get electronic configuration for each atom and store results in a dictionary."""
        electronic_configuration = {}
        for i, line in enumerate(self.contents):
            if "Natural Electron Configuration" in line:
                for j_line in self.contents[i + 2 :]:
                    if "Wiberg bond index matrix" in j_line:
                        break
                    if len(j_line) != 0:
                        columns = j_line.split()
                        atom_type = columns[0]
                        atom_number = columns[1]
                        configuration = "".join(columns[2:])
                        atom_key = f"{atom_type}{atom_number}"
                        electronic_configuration[atom_key] = configuration
        return electronic_configuration

    def get_num_naos(self, atom_key):
        """Get the number of NAOs for a given atom."""
        return len(self.natural_atomic_orbitals[atom_key])

    def get_total_electron_occ(self, atom_key):
        """Get the total electron occupancy for a given atom."""
        total_electron_occ = sum(
            entry["occupancy"]
            for entry in self.natural_atomic_orbitals[atom_key].values()
        )
        return total_electron_occ

    def get_electronic_configuration(self, atom_key):
        """Get the electronic configuration for a given atom."""
        return self.electronic_configuration[atom_key]
