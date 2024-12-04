import re
import logging
from itertools import islice
from functools import cached_property
import numpy as np
from ase import units
from chemsmart.utils.mixins import GaussianFileMixin
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.molecules.structure import CoordinateBlock
from chemsmart.utils.repattern import (
    eV_pattern,
    nm_pattern,
    f_pattern,
    float_pattern,
    normal_mode_pattern,
    frozen_coordinates_pattern,
    scf_energy_pattern,
    mp2_energy_pattern,
    oniom_energy_pattern,
)
from ase.io.formats import string2index

logger = logging.getLogger(__name__)


class Gaussian16Output(GaussianFileMixin):
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
    def num_steps(self):
        """Number of points scanned."""
        for line in self.contents:
            if line.startswith("Step number"):
                return int(line.split()[-1])
        return None

    @cached_property
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

    @cached_property
    def input_coordinates_block(self):
        """Obtain the coordinate block from the input that is printed in the outputfile."""
        coordinates_block_lines_list = []
        for i, line in enumerate(self.contents):
            if line.startswith("Symbolic Z-matrix:"):
                for j_line in self.contents[i + 2 :]:
                    if len(j_line) == 0:
                        break
                    coordinates_block_lines_list.append(j_line)
        cb = CoordinateBlock(coordinate_block=coordinates_block_lines_list)
        return cb

    @cached_property
    def symbols(self):
        return self.input_coordinates_block.chemical_symbols

    @cached_property
    def list_of_pbc_conditions(self):
        return self.input_coordinates_block.pbc_conditions

    # @cached_property
    # def all_structures(self):
    #     """Obtain all the structures from the output file.
    #     Use Standard orientations to get the structures; if not, use input orientations.
    #     Include their corresponding energy and forces if present.
    #     """
    #     all_structures = []
    #
    #     if self.normal_termination:
    #         # first clean up the number of structures and energies and forces by removing the structure without energy and forces
    #         if "opt" in self.route_string and "freq" in self.route_string:
    #             # both opt and freq jobs have been performed, then output structure is as follows:
    #             # 1. opt part
    #             # 		Input orientation
    #             # 		Standard orientation
    #             # 		Forces
    #             #
    #             # 		<after reaching optimized structure>
    #             # 		Input orientation
    #             # 		Standard orientation
    #             # 		<suspect this structure should be same as previously>
    #             #
    #             # 2. Freq part
    #             # 		Input orientation
    #             # 		Standard orientation
    #             # 		Forces
    #             assert np.allclose(
    #                 self.input_orientations[-1],
    #                 self.input_orientations[-2],
    #                 rtol=1e-5,
    #             ), "The last two input orientations should be the same."
    #             self.input_orientations.pop(-1)
    #             assert np.allclose(
    #                 self.standard_orientations[-1],
    #                 self.standard_orientations[-2],
    #                 rtol=1e-5,
    #             ), "The last two standard orientations should be the same."
    #             self.standard_orientations.pop(-1)
    #         # elif "opt" in self.route_string and not "freq" in self.route_string:
    #         #     # optimization job without freq calcs
    #
    #         if len(self.standard_orientations) == 0:
    #             # no standard orientation structures present, then use input orientation structures
    #
    #             for i, positions in enumerate(self.input_orientations):
    #                 all_structures.append(
    #                     Molecule(
    #                         symbols=self.symbols,
    #                         positions=positions,
    #                         charge=self.charge,
    #                         multiplicity=self.multiplicity,
    #                         frozen_atoms=self.frozen_atoms_masks,
    #                         pbc_conditions=self.list_of_pbc_conditions,
    #                         energy=self.energies_in_eV[i],
    #                         forces=self.forces_in_eV_per_A[i],
    #                     )
    #                 )
    #         else:
    #             # if standard orientation present, then use standard orientation structures
    #             for i, positions in enumerate(self.standard_orientations):
    #                 all_structures.append(
    #                     Molecule(
    #                         symbols=self.symbols,
    #                         positions=positions,
    #                         charge=self.charge,
    #                         multiplicity=self.multiplicity,
    #                         frozen_atoms=self.frozen_atoms_masks,
    #                         pbc_conditions=self.list_of_pbc_conditions,
    #                         energy=self.energies_in_eV[i],
    #                         forces=self.forces_in_eV_per_A[i],
    #                     )
    #                 )
    #         return all_structures
    #     else:
    #         # if the job is not terminated normally, then assemble structures according to least number of
    #         # energies and forces available
    #         num_structures_to_use = min(
    #             max(len(self.input_orientations), len(self.standard_orientations)),
    #             len(self.energies),
    #             len(self.forces)
    #         )
    #         if len(self.standard_orientations) == 0:
    #             # no standard orientation structures present, then use input orientation structures
    #
    #             for i, positions in enumerate(self.input_orientations[:num_structures_to_use]):
    #                 all_structures.append(
    #                     Molecule(
    #                         symbols=self.symbols,
    #                         positions=positions,
    #                         charge=self.charge,
    #                         multiplicity=self.multiplicity,
    #                         frozen_atoms=self.frozen_atoms_masks,
    #                         pbc_conditions=self.list_of_pbc_conditions,
    #                         energy=self.energies_in_eV[i],
    #                         forces=self.forces_in_eV_per_A[i],
    #                     )
    #                 )
    #         else:
    #             # if standard orientation present, then use standard orientation structures
    #             for i, positions in enumerate(self.standard_orientations[:num_structures_to_use]):
    #                 all_structures.append(
    #                     Molecule(
    #                         symbols=self.symbols,
    #                         positions=positions,
    #                         charge=self.charge,
    #                         multiplicity=self.multiplicity,
    #                         frozen_atoms=self.frozen_atoms_masks,
    #                         pbc_conditions=self.list_of_pbc_conditions,
    #                         energy=self.energies_in_eV[i],
    #                         forces=self.forces_in_eV_per_A[i],
    #                     )
    #                 )
    #         return all_structures

    @cached_property
    def all_structures(self):
        """
        Obtain all the structures from the output file.
        Use Standard orientations to get the structures; if not, use input orientations.
        Include their corresponding energy and forces if present.
        """

        def create_molecule_list(orientations, num_structures=None):
            """Helper function to create Molecule objects."""
            num_structures = num_structures or len(orientations)
            return [
                Molecule(
                    symbols=self.symbols,
                    positions=orientations[i],
                    charge=self.charge,
                    multiplicity=self.multiplicity,
                    frozen_atoms=self.frozen_atoms_masks,
                    pbc_conditions=self.list_of_pbc_conditions,
                    energy=self.energies_in_eV[i],
                    forces=self.forces_in_eV_per_A[i],
                )
                for i in range(num_structures)
            ]

        # If the job terminated normally
        if self.normal_termination:
            # Handle special case: both "opt" and "freq" present in route
            if "opt" in self.route_string and "freq" in self.route_string:
                assert np.allclose(
                    self.input_orientations[-1],
                    self.input_orientations[-2],
                    rtol=1e-5,
                ), "The last two input orientations should be the same."
                assert np.allclose(
                    self.standard_orientations[-1],
                    self.standard_orientations[-2],
                    rtol=1e-5,
                ), "The last two standard orientations should be the same."
                self.input_orientations.pop(-1)
                self.standard_orientations.pop(-1)

            # Use Standard orientations if available, otherwise Input orientations
            orientations = (
                self.standard_orientations
                if self.standard_orientations
                else self.input_orientations
            )
            return create_molecule_list(orientations)

        # If the job did not terminate normally
        num_structures_to_use = min(
            max(len(self.input_orientations), len(self.standard_orientations)),
            len(self.energies),
            len(self.forces),
        )
        orientations = (
            self.standard_orientations
            if self.standard_orientations
            else self.input_orientations
        )
        return create_molecule_list(
            orientations, num_structures=num_structures_to_use
        )

    @cached_property
    def optimized_structure(self):
        """Return optimized structure."""
        if self.normal_termination:
            return self.all_structures[-1]
        else:
            return None

    @cached_property
    def last_structure(self):
        """Return last structure, whether the output file has completed successfully or not."""
        return self.all_structures[-1]

    ######################### the following properties relate to intermediate geometry optimizations
    # for a constrained opt in e.g, scan/modred job

    @cached_property
    def intermediate_steps(self):
        """Return a list of intermediate steps."""
        initial_step = []
        final_step = []
        for line in self.contents:
            if line.startswith("Step number") and "on scan point" in line:
                line_elem = line.split()
                initial_step.append(int(line_elem[2]))
                final_step.append(int(line_elem[12]))
        steps_zip = zip(initial_step, final_step, strict=False)
        zipped_steps_list = list(steps_zip)
        if len(zipped_steps_list) != 0:
            return zipped_steps_list
        return None

    @cached_property
    def optimized_steps(self):
        """Return a list of optimized steps without intermediate steps."""
        steps = self.intermediate_steps
        if steps:
            optimized_steps = []
            for i in range(steps[-1][-1]):
                i_gaussian = i + 1  # gaussian uses 1-index
                each_steps = [step for step in steps if step[-1] == i_gaussian]
                optimized_steps.append(each_steps[-1])
            return optimized_steps
        return None

    @cached_property
    def optimized_steps_indices(self):
        if self.optimized_steps:
            return [
                self.intermediate_steps.index(i) for i in self.optimized_steps
            ]
        return None

    #########################

    def _get_route(self):
        lines = self.contents
        for i, line in enumerate(lines):
            if line.startswith("#"):
                if lines[i + 1].startswith("------"):
                    # route string in a single line
                    route = line.lower()
                elif not lines[i + 1].startswith("------") and lines[
                    i + 2
                ].startswith("------"):
                    # route string spans two lines
                    route = line.lower()
                    route += lines[i + 1].strip().lower()
                elif not lines[i + 1].startswith("------") and not lines[
                    i + 2
                ].startswith("------"):
                    # route string spans three lines
                    route = line.lower()
                    route += lines[i + 1].lower()
                    route += lines[i + 2].lower()
                else:
                    route = None
                return route
        return None

    @property
    def gen_genecp(self):
        """String specifying if gen or genecp is used in the calculation output file."""
        return self._get_gen_genecp()

    def _get_gen_genecp(self):
        if "gen" in self.basis:
            # return the string containing gen or genecp
            return self.basis
        return None

    @cached_property
    def num_basis_functions(self):
        for line in self.contents:
            if "basis functions," in line:
                line_elem = line.split(",")
                num_basis_functions = line_elem[0].strip().split()[0]
                return int(num_basis_functions)
        return None

    @cached_property
    def num_primitive_gaussians(self):
        for line in self.contents:
            if "primitive gaussians," in line:
                line_elem = line.split(",")
                primitives = line_elem[1].strip().split()[0]
                return int(primitives)
        return None

    @cached_property
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
    @cached_property
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

    @cached_property
    def service_units_by_jobs(self):
        """SUs defined as the JOB CPU time in hours."""
        return self.cpu_runtime_by_jobs_core_hours

    @cached_property
    def total_core_hours(self):
        return round(sum(self.cpu_runtime_by_jobs_core_hours), 4)

    @cached_property
    def total_service_unit(self):
        return self.total_core_hours

    @cached_property
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

    @cached_property
    def total_elapsed_walltime(self):
        return round(sum(self.elapsed_walltime_by_jobs), 4)

    #### FREQUENCY CALCULATIONS
    @cached_property
    def vibrational_frequencies(self):
        """Read the vibrational frequencies from the Gaussian output file."""
        frequencies = []
        for line in self.contents:
            if line.startswith("Frequencies --"):
                freq_string = line.split("--")[1].strip()
                for freq in freq_string.split():
                    frequencies.append(float(freq))
            else:
                continue
            if "Thermochemistry" in line:
                break
        return frequencies

    @cached_property
    def reduced_masses(self):
        """Obtain list of reduced masses corresponding to the vibrational frequency."""
        reduced_masses = []
        for line in self.contents:
            if line.startswith("Red. masses --"):
                reduced_masses_string = line.split("--")[1].strip()
                for mass in reduced_masses_string.split():
                    reduced_masses.append(float(mass))
            else:
                continue
            if "Thermochemistry" in line:
                break
        return reduced_masses

    @cached_property
    def force_constants(self):
        """Obtain list of force constants corresponding to the vibrational frequency."""
        force_constants = []
        for line in self.contents:
            if line.startswith("Frc consts  --"):
                force_constants_string = line.split("--")[1].strip()
                for force in force_constants_string.split():
                    force_constants.append(float(force))
            else:
                continue
            if "Thermochemistry" in line:
                break
        return force_constants

    @cached_property
    def ir_intensities(self):
        """Obtain list of IR intensities corresponding to the vibrational frequency."""
        IR_intensities = []
        for line in self.contents:
            if line.startswith("IR Inten    --"):
                IR_intensities_string = line.split("--")[1].strip()
                for intensity in IR_intensities_string.split():
                    IR_intensities.append(float(intensity))
            else:
                continue
            if "Thermochemistry" in line:
                break
        return IR_intensities

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
            else:
                continue
            if "Thermochemistry" in line:
                break
        return vibrational_mode_symmetries

    @cached_property
    def vibrational_modes(self):
        """Obtain list of vibrational normal modes corresponding to the vibrational frequency.
        Returns a list of normal modes, each of natoms x 3 (in dx, dy, and dz for each element) vibration.
        """
        list_of_vib_modes = []
        for i, line in enumerate(self.contents):
            if line.startswith("Frequencies --"):
                first_col_vib_modes = []
                second_col_vib_modes = []
                third_col_vib_modes = []
                for j_line in self.contents[i + 5 :]:
                    # if line match normal mode pattern
                    if re.match(normal_mode_pattern, j_line):
                        normal_mode = [float(val) for val in j_line.split()]
                        first_col_vib_mode = normal_mode[2:5]
                        second_col_vib_mode = normal_mode[5:8]
                        third_col_vib_mode = normal_mode[8:11]
                        first_col_vib_modes.append(first_col_vib_mode)
                        second_col_vib_modes.append(second_col_vib_mode)
                        third_col_vib_modes.append(third_col_vib_mode)
                    else:
                        break
                list_of_vib_modes.append(np.array(first_col_vib_modes))
                list_of_vib_modes.append(np.array(second_col_vib_modes))
                list_of_vib_modes.append(np.array(third_col_vib_modes))
            else:
                continue
            if "Thermochemistry" in line:
                break
        return list_of_vib_modes

    @cached_property
    def num_vib_modes(self):
        return len(self.vibrational_modes)

    @cached_property
    def num_vib_frequencies(self):
        return len(self.vibrational_frequencies)

    #### FREQUENCY CALCULATIONS
    @cached_property
    def has_frozen_coordinates(self):
        """Check if the output file has frozen coordinates."""
        has_frozen = []
        for i, line_i in enumerate(self.contents):
            if "Derivative Info." in line_i:
                for _j, line_j in enumerate(self.contents[i + 2 :]):
                    if "-----------------------------------" in line_j:
                        break
                    if line_j.split()[-2] == "Frozen":
                        has_frozen.append(True)
                    elif line_j.split()[-2] == "D2E/DX2":
                        has_frozen.append(False)
        return any(has_frozen)

    @cached_property
    def frozen_coordinate_indices(self):
        """Obtain list of frozen coordinate indices from the input format.
        Use 1-index to be the same as atom numbering."""
        frozen_coordinate_indices = []
        if self.has_frozen_coordinates:
            for i, line_i in enumerate(self.contents):
                if "Symbolic Z-matrix:" in line_i:
                    if len(line_i) == 0:
                        break
                    for j, line_j in enumerate(self.contents[i + 2 :]):
                        line_j_elem = line_j.split()
                        if (
                            re.match(frozen_coordinates_pattern, line_j)
                            and line_j_elem[1] == "-1"
                        ):
                            frozen_coordinate_indices.append(j + 1)
        return frozen_coordinate_indices

    @cached_property
    def free_coordinate_indices(self):
        """Obtain list of free coordinate indices from the input format by taking
        the complement of the frozen coordinates."""
        if self.has_frozen_coordinates:
            return [
                i
                for i in range(1, self.num_atoms + 1)
                if i not in self.frozen_coordinate_indices
            ]
        return None

    @cached_property
    def frozen_elements(self):
        frozen_atoms, _ = self._get_frozen_and_free_atoms()
        return frozen_atoms

    @cached_property
    def free_elements(self):
        _, free_atoms = self._get_frozen_and_free_atoms()
        return free_atoms

    def _get_frozen_and_free_atoms(self):
        """Obtain list of frozen and free atoms from the input format."""
        frozen_atoms = []
        free_atoms = []
        if self.has_frozen_coordinates:
            for i, line_i in enumerate(self.contents):
                if "Symbolic Z-matrix:" in line_i:
                    if len(line_i) == 0:
                        break
                    for j, line_j in enumerate(self.contents[i + 2 :]):
                        line_j_elem = line_j.split()
                        if (
                            re.match(frozen_coordinates_pattern, line_j)
                            and line_j_elem[1] == "-1"
                        ):
                            frozen_atoms.append(line_j_elem[0])
                        elif (
                            re.match(frozen_coordinates_pattern, line_j)
                            and line_j_elem[1] == "0"
                        ):
                            free_atoms.append(line_j_elem[0])
        return frozen_atoms, free_atoms

    @cached_property
    def frozen_atoms_masks(self):
        """Obtain list of frozen atoms masks from the input format.
        -1 is used for frozen atoms and 0 for free atoms."""
        frozen_atoms_masks = []
        if self.has_frozen_coordinates:
            for i, line_i in enumerate(self.contents):
                if "Symbolic Z-matrix:" in line_i:
                    if len(line_i) == 0:
                        break
                    for j, line_j in enumerate(self.contents[i + 2 :]):
                        line_j_elem = line_j.split()
                        if (
                            re.match(frozen_coordinates_pattern, line_j)
                            and line_j_elem[1] == "-1"
                        ):
                            frozen_atoms_masks.append(-1)
                        elif (
                            re.match(frozen_coordinates_pattern, line_j)
                            and line_j_elem[1] == "0"
                        ):
                            frozen_atoms_masks.append(0)
            return frozen_atoms_masks
        return None

    @cached_property
    def scf_energies(self):
        """Obtain SCF energies from the Gaussian output file. Default units of Hartree."""
        scf_energies = []
        for line in self.contents:
            match = re.match(scf_energy_pattern, line)
            if match:
                scf_energies.append(float(match[1]))
        return scf_energies

    @cached_property
    def mp2_energies(self):
        """Obtain MP2 energies from the Gaussian output file. Default units of Hartree."""
        mp2_energies = []
        for line in self.contents:
            match = re.search(mp2_energy_pattern, line)
            if match:
                mp2_energies.append(float(match[1].replace("D", "E")))
        return mp2_energies

    @cached_property
    def oniom_energies(self):
        """Obtain ONIOM energies from the Gaussian output file. Default units of Hartree."""
        oniom_energies = []
        for line in self.contents:
            match = re.match(oniom_energy_pattern, line)
            if match:
                oniom_energies.append(float(match[1]))
        return oniom_energies

    @cached_property
    def energies(self):
        """Return energies of the system."""
        if len(self.mp2_energies) == 0 and len(self.oniom_energies) == 0:
            return self.scf_energies
        elif len(self.mp2_energies) != 0:
            return self.mp2_energies
        elif len(self.oniom_energies) != 0:
            return self.oniom_energies

    @cached_property
    def energies_in_eV(self):
        """Convert energies from Hartree to eV."""
        return [energy * units.Hartree for energy in self.energies]

    @property
    def num_energies(self):
        return len(self.energies_in_eV)

    # check for convergence criterion not met (happens for some output files)
    @property
    def convergence_criterion_not_met(self):
        return any(
            ">>>>>>>>>> Convergence criterion not met." in line
            for line in self.contents
        )

    @cached_property
    def has_forces(self):
        """Check if the output file contains forces calculations."""
        for line in self.contents:
            if "Forces (Hartrees/Bohr)" in line:
                return True
        return False

    @cached_property
    def forces(self):
        """Obtain a list of cartesian forces.
        Each force is stored as a np array of shape (natoms, 3).
        Intrinsic units as used in Gaussian: Hartrees/Bohr."""
        list_of_all_forces = []
        for i, line in enumerate(self.contents):
            if "Forces (Hartrees/Bohr)" in line:
                forces = []
                for j_line in self.contents[i + 3 :]:
                    if "---------------------------" in j_line:
                        break
                    if j_line.startswith("-2"):
                        # remove those pbc forces
                        continue
                    forces.append([float(val) for val in j_line.split()[2:5]])
                list_of_all_forces.append(np.array(forces))
        return list_of_all_forces

    @cached_property
    def forces_in_eV_per_A(self):
        """Convert forces from Hartrees/Bohr to eV/Angstrom."""
        forces_in_eV_per_A = []
        for forces in self.forces:
            forces_in_eV_per_A.append(forces * units.Hartree / units.Bohr)
        return forces_in_eV_per_A

    @cached_property
    def num_forces(self):
        return len(self.forces_in_eV_per_A)

    @cached_property
    def input_orientations(self):
        """Obtain structures in Input Orientation from Gaussian output file."""
        input_orientations = []
        for i, line in enumerate(self.contents):
            if line.startswith("Input orientation:"):
                input_orientation = []
                for j_line in self.contents[i + 5 :]:
                    if "-----------------" in j_line:
                        break
                    if j_line.split()[1] == "-2":  # atomic number = -2 for TV
                        continue
                    input_orientation.append(
                        [float(val) for val in j_line.split()[3:6]]
                    )
                input_orientations.append(np.array(input_orientation))
        return input_orientations

    @cached_property
    def standard_orientations(self):
        """Obtain structures in Standard Orientation from Gaussian output file."""
        standard_orientations = []
        for i, line in enumerate(self.contents):
            if line.startswith("Standard orientation:"):
                standard_orientation = []
                for j_line in self.contents[i + 5 :]:
                    if "-----------------" in j_line:
                        break
                    if j_line.split()[1] == "-2":  # atomic number = -2 for TV
                        continue
                    standard_orientation.append(
                        [float(val) for val in j_line.split()[3:6]]
                    )
                standard_orientations.append(np.array(standard_orientation))
        return standard_orientations

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
                contribution_coefficients.append(
                    each_state_contribution_coefficients
                )

        return transitions, contribution_coefficients

    @cached_property
    def alpha_occ_eigenvalues(self):
        """Obtain all eigenenergies of the alpha occuplied orbitals and convert to eV."""
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

            alpha_occ_eigenvalues = [
                value * units.Hartree for value in last_block_values
            ]
        return alpha_occ_eigenvalues

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

            alpha_virtual_eigenvalues = [
                value * units.Hartree for value in last_block_values
            ]
            return alpha_virtual_eigenvalues

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

            beta_occ_eigenvalues = [
                value * units.Hartree for value in last_block_values
            ]
            return beta_occ_eigenvalues

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

            beta_virtual_eigenvalues = [
                value * units.Hartree for value in last_block_values
            ]
            return beta_virtual_eigenvalues

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
                len(self.alpha_occ_eigenvalues)
                - len(self.beta_occ_eigenvalues)
                + 1
                == self.multiplicity
            )
            return len(self.alpha_occ_eigenvalues) - len(
                self.beta_occ_eigenvalues
            )

    @cached_property
    def somo_energy(self):
        if self.multiplicity != 1:
            # the multiplicity is the number of unpaired electrons + 1
            assert (
                len(self.alpha_occ_eigenvalues)
                - len(self.beta_occ_eigenvalues)
                + 1
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
            return self.lumo_energy - self.homo_energy
        else:
            # to implement for radical systems
            pass

    def get_molecule(self, index="-1", include_failed_logfile=True):
        index = string2index(index)
        return self.all_structures[index]

    def to_dataset(self, **kwargs):
        """Convert Gaussian .log file to Dataset with all data points taken from the .log file.

        Returns:
            Dataset.
        """
        # TODO: to be implemented
        pass


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
            if (
                "NAO  Atom  No  lang   Type(AO)    Occupancy      Energy"
                in line
            ):
                for j_line in self.contents[i + 2 :]:
                    if (
                        "WARNING" in j_line
                        or "Summary of Natural Population Analysis" in j_line
                    ):
                        break
                    if len(j_line) != 0:
                        columns = j_line.split()

                        # Extract values from each column
                        nao_number = int(
                            columns[0]
                        )  # NAO Number (like 1, 2, etc.)
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


class Gaussian16OutputWithPBC(Gaussian16Output):
    """class for parsing and obtaining information from Gaussian output file with PBC."""

    def __init__(self, filename):
        super().__init__(filename=filename)

    def _parse(self, filename):
        pass

    @property
    def pbc(self):
        for line in self.contents:
            if "Periodicity:" in line:
                pbc_conditions = line.split("Periodicity:")[-1]
                pbc_conditions = pbc_conditions.split()
                assert (
                    len(pbc_conditions) == 3
                ), "Periodicity given for 3 dimensions."
                return np.array(
                    [
                        int(pbc_conditions[0]),
                        int(pbc_conditions[1]),
                        int(pbc_conditions[2]),
                    ]
                )
        return None

    @property
    def dim(self):
        d = 0
        for i in self.pbc:
            if i == 1:
                d += 1
        return d

    @property
    def input_translation_vectors(self):
        for i, line in enumerate(self.contents):
            if "Lengths of translation vectors:" in line:
                # get cells just once
                all_cells = []
                for tv_line in self.contents[i - 1 - self.num_atoms : i - 1]:
                    tv_line_elem = tv_line.split()
                    if "-----------------" in tv_line:
                        continue
                    if float(tv_line_elem[1]) == -2.0:
                        tv_vector = [
                            float(tv_line_elem[-3]),
                            float(tv_line_elem[-2]),
                            float(tv_line_elem[-1]),
                        ]
                        all_cells.append(tv_vector)
                return np.array(all_cells)
        return None

    @property
    def final_translation_vector(self):
        """Get final translation vectors from last step."""
        for i, line in enumerate(reversed(self.contents)):
            # read from backwards and get the last translation vector
            if "Lengths of translation vectors:" in line:
                start_idx = (
                    len(self.contents) - i
                )  # Line index in forward order
                # get cells just once
                all_cells = []
                for tv_line in self.contents[
                    start_idx - self.num_atoms - 4 : start_idx - 1
                ]:
                    if "-----------------" in tv_line:
                        continue
                    tv_line_elem = tv_line.split()
                    if float(tv_line_elem[1]) == -2.0:
                        tv_vector = [
                            float(tv_line_elem[-3]),
                            float(tv_line_elem[-2]),
                            float(tv_line_elem[-1]),
                        ]
                        all_cells.append(tv_vector)
                return np.array(all_cells)
        return None

    # def get_atoms(self, index="-1", include_failed_logfile=False):
    #     all_atoms = []
    #
    #     if not self.normal_termination and not include_failed_logfile:
    #         return None
    #
    #     all_energies = self.get_energies()
    #     all_structures = self.get_structures()
    #     all_forces = self.get_forces()
    #     chemical_symbols = self.get_symbols()
    #     assert len(all_energies) > 0, "No energies found!"
    #     assert len(all_structures) > 0, "No structures found!"
    #
    #     # order in the .log file: input orientation, then SCF Done, then Forces
    #     if len(all_forces) != 0:
    #         # forces are given
    #         min_num_structures = min(
    #             len(all_energies), len(all_structures), len(all_forces)
    #         )
    #         for i in range(min_num_structures):
    #             atoms = Atoms(
    #                 symbols=chemical_symbols,
    #                 positions=all_structures[i],
    #                 pbc=self.pbc,
    #                 cell=self.cells,
    #             )
    #             atoms.calc = SinglePointCalculator(
    #                 atoms=atoms, energy=all_energies[i], forces=all_forces[i]
    #             )
    #             all_atoms.append(atoms)
    #     else:
    #         # no forces information
    #         min_num_structures = min(len(all_energies), len(all_structures))
    #         for i in range(min_num_structures):
    #             atoms = Atoms(
    #                 symbols=chemical_symbols,
    #                 positions=all_structures[i],
    #                 pbc=self.pbc,
    #                 cell=self.cells,
    #             )
    #             atoms.calc = SinglePointCalculator(
    #                 atoms=atoms, energy=all_energies[i]
    #             )
    #             all_atoms.append(atoms)
    #
    #     index = string2index(index)
    #     return all_atoms[index]
