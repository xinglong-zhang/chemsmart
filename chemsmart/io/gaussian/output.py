import logging
import re
from functools import cached_property

import numpy as np
from ase import units

from chemsmart.io.molecules.structure import CoordinateBlock
from chemsmart.utils.io import clean_duplicate_structure, create_molecule_list
from chemsmart.utils.mixins import GaussianFileMixin
from chemsmart.utils.periodictable import PeriodicTable
from chemsmart.utils.repattern import (
    eV_pattern,
    f_pattern,
    float_pattern,
    frozen_coordinates_pattern,
    mp2_energy_pattern,
    nm_pattern,
    normal_mode_pattern,
    oniom_energy_pattern,
    oniom_gridpoint_pattern,
    scf_energy_pattern,
)
from chemsmart.utils.utils import (
    get_range_from_list,
    safe_min_lengths,
    string2index_1based,
)

p = PeriodicTable()
logger = logging.getLogger(__name__)


class Gaussian16Output(GaussianFileMixin):
    """Comprehensive parser for Gaussian 16 output files.

    This class provides extensive parsing capabilities for Gaussian output
    files, extracting molecular geometries, energies, vibrational frequencies,
    thermochemical data, and other computational results. It supports various
    calculation types including single-point energy calculations, geometry
    optimizations, frequency analyses, and potential energy surface scans.

    Args:
        filename (str): Path to the Gaussian output file to parse
        use_frozen (bool, optional): Whether to include frozen coordinates
            in molecular structures. Defaults to False. When False, frozen
            coordinates are excluded to avoid issues when using parsed
            molecules as input for subsequent calculations.
        include_intermediate (bool, optional): Whether to include intermediate
            optimization steps. Defaults to False. When False, only converged
            geometries are included (matching GaussView behavior). When True,
            all geometry steps are included, useful for detailed trajectory
            analysis or selecting specific points from scan calculations.
    """

    def __init__(self, filename, use_frozen=False, include_intermediate=False):
        """
        Initialize the Gaussian output parser.
        """
        self._energies = None
        self.filename = filename
        self.use_frozen = use_frozen
        self.include_intermediate = include_intermediate

    @property
    def normal_termination(self):
        """
        Check if the Gaussian calculation terminated normally.

        Examines the last line of the output file for the standard
        Gaussian termination message to determine if the calculation
        completed successfully.
        """
        contents = self.contents
        if len(contents) == 0:
            return False

        last_line = contents[-1]
        if "Normal termination of Gaussian" in last_line:
            logger.debug(f"File {self.filename} terminated normally.")
            return True

        logger.debug(f"File {self.filename} has error termination.")
        return False

    @property
    def heavy_elements(self):
        """TODO"""
        return None

    @property
    def heavy_elements_basis(self):
        """TODO"""
        return None

    @property
    def light_elements(self):
        """TODO"""
        return None

    @property
    def light_elements_basis(self):
        """TODO"""
        return None

    @property
    def custom_solvent(self):
        """TODO"""
        return None

    @cached_property
    def num_steps(self):
        """
        Number of points scanned.
        """
        for line in self.contents:
            if line.startswith("Step number"):
                return int(line.split()[-1])
        return None

    @cached_property
    def num_atoms(self):
        """
        Number of atoms in the molecular system.
        """
        for line in self.contents:
            if line.startswith("NAtoms="):
                return int(line.split()[1])

    @property
    def charge(self):
        """
        Charge of the molecule.
        """
        for line in self.contents:
            if "Charge" in line and "Multiplicity" in line:
                line_elem = line.split()
                return int(line_elem[2])

    @property
    def multiplicity(self):
        """
        Multiplicity of the molecule.
        """
        for line in self.contents:
            if "Charge" in line and "Multiplicity" in line:
                line_elem = line.split()
                # return int(line_elem[-1])
                # # in qmmm, not always last, e.g.,
                # # Charge = 1 Multiplicity = 2 for
                # low level calculation on real system.
                return int(line_elem[5])

    @property
    def spin(self):
        """
        Determine if calculation uses restricted or unrestricted spin.

        Analyzes the SCF method specification to determine whether
        the calculation uses restricted (R) or unrestricted (U) spin.
        """
        for line in self.contents:
            if "SCF Done:" in line:
                line_elem = line.split("E(")
                theory = line_elem[1].split(")")[0].lower()
                # Determine if the job is restricted or unrestricted
                if theory.startswith("r"):
                    spin = "restricted"
                elif theory.startswith("u"):
                    spin = "unrestricted"
                else:
                    spin = None
                return spin
        return None

    @cached_property
    def input_coordinates_block(self):
        """Obtain the coordinate block from the
        input that is printed in the outputfile."""
        coordinates_block_lines_list = []
        for i, line in enumerate(self.contents):
            if line.startswith("Symbolic Z-matrix:"):
                for j_line in self.contents[i + 2 :]:
                    if len(j_line) == 0:
                        break
                    if j_line.split()[0] not in p.PERIODIC_TABLE:
                        if j_line.startswith("Charge ="):
                            logger.debug(f"Skipping line: {j_line}")
                            # e.g., in QM/MM output files, the first
                            # element is not the coordinates information
                            # e.g., "Charge = 1 Multiplicity = 2 for
                            # low level calculation on real system."
                            continue
                        # elif j_line.startswith("TV"):
                        # we still add the line for PBC
                    coordinates_block_lines_list.append(j_line)
        cb = CoordinateBlock(coordinate_block=coordinates_block_lines_list)
        return cb

    @cached_property
    def all_structures(self):
        """
        Obtain all the structures from the output file.
        Use Standard orientations to get the
        structures; if not, use input orientations.
        Include their corresponding energy and forces if present.
        """
        return self._get_all_molecular_structures()

    def _get_all_molecular_structures(self):
        """
        Build and return the list of Molecule
        objects parsed from a calculation.

        Selection precedence for orientations:
          1) standard_orientations (+ PBC)
          2) input_orientations   (+ PBC)

        Special handling (as per design):
          - Non-link & normal termination:
          de-duplicate the terminal (e.g. freq) frame.
          - Non-link & abnormal termination: use safe_min_lengths.
          - Link & normal termination: drop the first frame, then behave like
            normal termination (incl. de-dup).
          - Link & abnormal termination:
              * if multiple frames: drop the first
              (carry-over), then use safe_min_lengths.
              * if single frame: return that only frame (no drop).

        Returns
        -------
        list
            List of Molecule objects (possibly empty).
        """
        # 1) Choose orientations (and their PBC)
        if self.standard_orientations:
            orientations = list(self.standard_orientations)
            orientations_pbc = list(self.standard_orientations_pbc or [])
        elif self.input_orientations:
            orientations = list(self.input_orientations)
            orientations_pbc = list(self.input_orientations_pbc or [])
        else:
            return []  # Nothing to build

        energies = list(self.energies) if self.energies else None
        forces = list(self.forces) if self.forces else None

        # Helper to drop the first item across all arrays (when present)
        def drop_first():
            nonlocal orientations, orientations_pbc, energies, forces
            if orientations:
                orientations = orientations[1:]
            if orientations_pbc:
                orientations_pbc = orientations_pbc[1:]
            if energies:
                energies = energies[1:]
            # Forces do not need to be dropped here because no
            # force computation occurs at the first link job.
            # This is intentional: the forces array is
            # already aligned with the relevant orientations.

        # Helper to keep only the last frame across all arrays
        def keep_last_only():
            nonlocal orientations, orientations_pbc, energies, forces
            orientations = orientations[-1:] if orientations else []
            orientations_pbc = (
                orientations_pbc[-1:] if orientations_pbc else []
            )
            if energies:
                energies = energies[-1:]
            if forces:
                forces = forces[-1:]

        # Right-trim auxiliaries to the number of
        # orientations (no data loss in orientations)
        def align_lengths_to_orientations():
            nonlocal orientations, orientations_pbc, energies, forces
            n = len(orientations)
            if orientations_pbc and len(orientations_pbc) > n:
                orientations_pbc = orientations_pbc[:n]
            if energies and len(energies) > n:
                energies = energies[:n]
            if forces and len(forces) > n:
                forces = forces[:n]

        # 2) Handle link jobs
        if self.is_link:
            if self.normal_termination:
                logger.debug(
                    "Link job with normal termination: dropping first frame."
                )
                if orientations:  # drop carried-over first frame if present
                    drop_first()

                # Single-point link jobs: keep only the last frame (after drop)
                if self.jobtype == "sp" and orientations:
                    keep_last_only()

                # Fall through to "normal termination" handling below
            else:
                logger.debug("Link job with error termination.")
                if len(orientations) > 1:
                    # Multiple frames: drop carried-over
                    # first, then treat as abnormal
                    drop_first()
                    # Fall through to abnormal
                    # handling below (safe_min_lengths)
                else:
                    # Single frame available: return it as-is (do NOT drop)
                    frozen_atoms = (
                        self.frozen_atoms_masks if self.use_frozen else None
                    )
                    return create_molecule_list(
                        orientations=orientations,
                        orientations_pbc=orientations_pbc,
                        energies=energies,
                        forces=forces,
                        symbols=self.symbols,
                        charge=self.charge,
                        multiplicity=self.multiplicity,
                        frozen_atoms=frozen_atoms,
                        pbc_conditions=self.list_of_pbc_conditions,
                    )

        # 3) De-dup only for normal-termination
        # paths (incl. link-normal after drop)
        if self.normal_termination:
            clean_duplicate_structure(orientations)
            # After dedup, ensure auxiliaries aren't longer than orientations
            align_lengths_to_orientations()

        # 4) Compute safe min length for logging and
        # for abnormal (non-link or link>1-after-drop)
        num_structures_to_use = safe_min_lengths(
            orientations, energies, forces
        )

        logger.debug(
            "Structures to use: %d | orientations=%d | energies=%d | forces=%d",
            num_structures_to_use,
            len(orientations),
            len(energies) if energies is not None else 0,
            len(forces) if forces is not None else 0,
        )

        frozen_atoms = self.frozen_atoms_masks if self.use_frozen else None

        # 5) Build Molecule list
        create_kwargs = dict(
            orientations=orientations,
            orientations_pbc=orientations_pbc,
            energies=energies,
            forces=forces,
            symbols=self.symbols,
            charge=self.charge,
            multiplicity=self.multiplicity,
            frozen_atoms=frozen_atoms,
            pbc_conditions=self.list_of_pbc_conditions,
        )

        if self.normal_termination:
            all_structures = create_molecule_list(**create_kwargs)
        else:
            # Abnormal (non-link, or link with >1
            # frame after drop): truncate safely
            all_structures = create_molecule_list(
                **create_kwargs, num_structures=num_structures_to_use
            )

        # 6) Keep only optimized steps if requested
        if (
            getattr(self, "optimized_steps_indices", None)
            and not self.include_intermediate
        ):
            logger.debug(
                "Ignoring intermediate optimization steps (constrained opt)."
            )
            all_structures = [
                all_structures[i] for i in self.optimized_steps_indices
            ]

        logger.debug(
            "Attaching vibrational data to the final structure if available..."
        )

        last_mol = all_structures[-1]
        # Attach vibrational data to the final structure if available
        if self.num_vib_frequencies:
            all_structures[-1] = self._attach_vib_metadata(last_mol)

        logger.debug("Total structures returned: %d", len(all_structures))
        return all_structures

    @cached_property
    def optimized_structure(self):
        """
        Return optimized structure.
        """
        if self.normal_termination:
            return self.all_structures[-1]
        else:
            return None

    @cached_property
    def last_structure(self):
        """
        Return the last molecular structure from the calculation.

        Returns the final structure regardless of whether the calculation
        completed successfully. Useful for analyzing partially converged
        optimizations or error cases.
        """
        return self.all_structures[-1]

    @property
    def molecule(self):
        """
        Alias for the last molecular structure.
        """
        return self.last_structure

    ###### the following properties relate to
    # intermediate geometry optimizations
    # for a constrained opt in e.g, scan/modred job

    @cached_property
    def intermediate_steps(self):
        """
        Return a list of intermediate steps.
        """
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
        """
        Return a list of optimized steps without intermediate steps.
        """
        steps = self.intermediate_steps
        if steps:
            optimized_steps = []
            for i in range(steps[-1][-1]):
                i_gaussian = i + 1  # gaussian uses 1-index
                each_steps = [step for step in steps if step[-1] == i_gaussian]
                optimized_steps.append(each_steps[-1])
            logger.debug(f"Optimized steps: {optimized_steps}")
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
    @property
    def modredundant_group(self):
        return self._get_modredundant_group()

    def _get_route(self):
        lines = self.contents
        for i, line in enumerate(lines):
            if line.startswith("#") and "stable=opt" in line:
                continue
            elif line.startswith("#"):
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

    def _get_modredundant_group(self):
        if "modred" in self.route_string:
            modredundant_group = []
            for i, line in enumerate(self.contents):
                if line.startswith(
                    "The following ModRedundant input section has been read:"
                ):
                    for j_line in self.contents[i + 1 :]:

                        if len(j_line) == 0:
                            break
                        modredundant_group.append(j_line)
            return modredundant_group
        return None

    @property
    def gen_genecp(self):
        """
        String specifying if gen or genecp is
        used in the calculation output file.
        """
        return self._get_gen_genecp()

    def _get_gen_genecp(self):
        if self.basis is None:
            # this happens for semi-empirical calculations
            return None
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
        cpu_runtimes = []
        for line in self.contents:
            if line.startswith("Job cpu time:"):
                cpu_runtime = []
                n_days = float(line.split("days")[0].strip().split()[-1])
                cpu_runtime.append(n_days * 24)
                n_hours = float(line.split("hours")[0].strip().split()[-1])
                cpu_runtime.append(n_hours)
                n_minutes = float(line.split("minutes")[0].strip().split()[-1])
                cpu_runtime.append(n_minutes / 60)
                n_seconds = float(line.split("seconds")[0].strip().split()[-1])
                cpu_runtime.append(n_seconds / 3600)
                cpu_runtimes.append(sum(cpu_runtime))
        return cpu_runtimes

    @cached_property
    def service_units_by_jobs(self):
        """
        SUs defined as the JOB CPU time in hours.
        """
        return self.cpu_runtime_by_jobs_core_hours

    @cached_property
    def total_core_hours(self):
        return round(sum(self.cpu_runtime_by_jobs_core_hours), 1)

    @cached_property
    def total_service_unit(self):
        return self.total_core_hours

    @cached_property
    def elapsed_walltime_by_jobs(self):
        elapsed_walltimes = []
        for line in self.contents:
            if line.startswith("Elapsed time:"):
                elapsed_walltime = []
                n_days = float(line.split("days")[0].strip().split()[-1])
                elapsed_walltime.append(n_days * 24)
                n_hours = float(line.split("hours")[0].strip().split()[-1])
                elapsed_walltime.append(n_hours)
                n_minutes = float(line.split("minutes")[0].strip().split()[-1])
                elapsed_walltime.append(n_minutes / 60)
                n_seconds = float(line.split("seconds")[0].strip().split()[-1])
                elapsed_walltime.append(n_seconds / 3600)
                elapsed_walltimes.append(sum(elapsed_walltime))
        return elapsed_walltimes

    @cached_property
    def total_elapsed_walltime(self):
        return round(sum(self.elapsed_walltime_by_jobs), 1)

    #### FREQUENCY CALCULATIONS
    @cached_property
    def vibrational_frequencies(self):
        """
        Read the vibrational frequencies from the Gaussian output file.
        """
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
        """
        Obtain list of reduced masses
        corresponding to the vibrational frequency.
        """
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
        """
        Obtain list of force constants
        corresponding to the vibrational frequency.
        """
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
        """
        Obtain list of IR intensities
        corresponding to the vibrational frequency.
        """
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
        """
        Obtain list of vibrational mode symmetries
        corresponding to the vibrational frequency.
        """
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
        """
        Number of vibrational frequencies found.
        """
        return len(self.vibrational_frequencies)

    @property
    def enthalpy(self):
        """
        Extract enthalpy from Gaussian output.

        This is the "Sum of electronic and thermal Enthalpies" value
        printed in the thermochemistry section of the output file.

        Returns:
            float: enthalpy in Hartree, or None if not found.
        """
        for line in self.contents:
            if "Sum of electronic and thermal Enthalpies=" in line:
                try:
                    return float(line.split()[-1])
                except (ValueError, IndexError):
                    pass
        return None

    @property
    def gibbs_free_energy(self):
        """
        Extract Gibbs free energy from Gaussian output.

        This is the "Sum of electronic and thermal Free Energies" value
        printed in the thermochemistry section of the output file.

        Returns:
            float: Gibbs free energy in Hartree, or None if not found.
        """
        for line in self.contents:
            if "Sum of electronic and thermal Free Energies=" in line:
                try:
                    return float(line.split()[-1])
                except (ValueError, IndexError):
                    pass
        return None

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

    #### FREQUENCY CALCULATIONS

    @cached_property
    def has_frozen_coordinates(self):
        """Check if the calculation includes frozen coordinates.

        Returns:
            bool: True if frozen coordinates are present
        """
        return self.frozen_coordinate_indices is not None

    @cached_property
    def frozen_coordinate_indices(self):
        """
        Obtain list of frozen coordinate indices from the input format.
        Use 1-index to be the same as atom numbering.
        """
        frozen_coordinate_indices = []
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
        if len(frozen_coordinate_indices) == 0:
            return None
        return frozen_coordinate_indices

    @cached_property
    def free_coordinate_indices(self):
        """
        Obtain list of free coordinate indices from the input format by taking
        the complement of the frozen coordinates.
        """
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
        """
        Obtain list of frozen and free atoms from the input format.
        """
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
        """
        Obtain list of frozen atoms masks (-1 = frozen, 0 = free)
        using precomputed frozen_coordinate_indices
        and free_coordinate_indices.
        """
        if not self.has_frozen_coordinates:
            return None

        # Initialize all as free (0), then mark frozen (-1)
        masks = [0] * self.num_atoms  # 0-based list
        for idx in self.frozen_coordinate_indices:
            masks[idx - 1] = -1  # Convert 1-based index to 0-based

        return masks

    @cached_property
    def scf_energies(self):
        """
        Obtain SCF energies from the Gaussian
        output file. Default units of Hartree.
        """
        scf_energies = []
        for line in self.contents:
            match = re.search(scf_energy_pattern, line)
            if match:
                scf_energies.append(float(match[1]))
        return scf_energies

    @cached_property
    def mp2_energies(self):
        """
        Obtain MP2 energies from the Gaussian
        output file. Default units of Hartree.
        """
        mp2_energies = []
        for line in self.contents:
            match = re.search(mp2_energy_pattern, line)
            if match:
                # Convert Gaussian's D notation to standard E notation
                mp2_energies.append(float(match[1].replace("D", "E")))
        return mp2_energies

    @cached_property
    def oniom_energies(self):
        """
        Obtain ONIOM energies from the Gaussian
        output file. Default units of Hartree.
        """
        oniom_energies = []
        for line in self.contents:
            match = re.match(oniom_energy_pattern, line)
            if match:
                oniom_energies.append(float(match[1]))
        return oniom_energies

    @cached_property
    def oniom_layer_energies(self):
        """Obtain ONIOM energies from the Gaussian
        output file. Default units of Hartree."""
        layer_energies = {}
        for line in self.contents:
            layer_match = re.match(oniom_gridpoint_pattern, line)
            if layer_match:
                formatted_layer = (
                    f"{line.split()[3]}  {line.split()[4]}, "
                    f"{line.split()[5]}  {line.split()[6]}"
                )
                layer_energies[formatted_layer] = float(line.split()[-1])
        return layer_energies

    @cached_property
    def energies(self):
        """
        Return energies of the system.
        """
        if len(self.mp2_energies) == 0 and len(self.oniom_energies) == 0:
            return self.scf_energies
        elif len(self.mp2_energies) != 0:
            return self.mp2_energies
        elif len(self.oniom_energies) != 0:
            return self.oniom_energies

    @cached_property
    def zero_point_energy(self):
        """
        Zero point energy in Hartree.
        """
        for line in self.contents:
            if "Zero-point correction=" in line:
                return float(line.split()[2])
        return None

    # check for convergence criterion not met (happens for some output files)
    @property
    def convergence_criterion_not_met(self):
        return any(
            ">>>>>>>>>> Convergence criterion not met." in line
            for line in self.contents
        )

    @cached_property
    def has_forces(self):
        """
        Check if the output file contains forces calculations.
        """
        for line in self.contents:
            if "Forces (Hartrees/Bohr)" in line:
                return True
        return False

    @cached_property
    def forces(self):
        list_of_all_forces, _ = self._get_forces_for_molecules_and_pbc()
        return list_of_all_forces

    def _get_forces_for_molecules_and_pbc(self):
        """
        Obtain a list of cartesian forces.
        Each force is stored as a np array of shape (num_atoms, 3).
        Intrinsic units as used in Gaussian: Hartrees/Bohr."""
        list_of_all_forces = []
        list_of_all_forces_pbc = []
        for i, line in enumerate(self.contents):
            if "Forces (Hartrees/Bohr)" in line:
                forces = []
                forces_pbc = []
                for j_line in self.contents[i + 3 :]:
                    if "---------------------------" in j_line:
                        break
                    if j_line.startswith("-2"):
                        # line indicates forces for pbc
                        forces_pbc.append(
                            [float(val) for val in j_line.split()[1:4]]
                        )
                    else:
                        forces.append(
                            [float(val) for val in j_line.split()[2:5]]
                        )
                list_of_all_forces.append(np.array(forces))
                list_of_all_forces_pbc.append(np.array(forces_pbc))
        if len(list_of_all_forces) == 0:
            return None, None
        return list_of_all_forces, list_of_all_forces_pbc

    @cached_property
    def pbc_forces(self):
        _, list_of_all_forces_pbc = self._get_forces_for_molecules_and_pbc()
        return list_of_all_forces_pbc

    @cached_property
    def num_forces(self):
        return len(self.forces)

    @cached_property
    def input_orientations(self):
        """Obtain structures in Input Orientation from Gaussian output file."""

        input_orientations, _ = self._get_input_orientations_and_pbc()
        return input_orientations

    @cached_property
    def input_orientations_pbc(self):
        """Obtain structures in Input Orientation
        with PBC from Gaussian output file."""
        _, input_orientations_pbc = self._get_input_orientations_and_pbc()
        return input_orientations_pbc

    def _get_input_orientations_and_pbc(self):
        input_orientations = []
        input_orientations_pbc = []
        for i, line in enumerate(self.contents):
            if line.startswith("Input orientation:"):
                input_orientation = []
                input_orientation_pbc = []
                for j_line in self.contents[i + 5 :]:
                    if "-----------------" in j_line:
                        break
                    if j_line.split()[1] == "-2":  # atomic number = -2 for TV
                        input_orientation_pbc.append(
                            [float(val) for val in j_line.split()[3:6]]
                        )
                    else:
                        input_orientation.append(
                            [float(val) for val in j_line.split()[3:6]]
                        )
                input_orientations.append(np.array(input_orientation))
                if len(input_orientation_pbc) != 0:
                    input_orientations_pbc.append(
                        np.array(input_orientation_pbc)
                    )
                else:
                    input_orientations_pbc.append(None)
        if len(input_orientations) == 0:
            return None, None
        return input_orientations, input_orientations_pbc

    @cached_property
    def standard_orientations(self):
        """Obtain structures in Standard
        Orientation from Gaussian output file."""
        standard_orientations, _ = self._get_standard_orientations_and_pbc()
        return standard_orientations

    @cached_property
    def standard_orientations_pbc(self):
        """Obtain structures in Standard Orientation
        with PBC from Gaussian output file."""
        _, standard_orientations_pbc = (
            self._get_standard_orientations_and_pbc()
        )
        return standard_orientations_pbc

    def _get_standard_orientations_and_pbc(self):
        standard_orientations = []
        standard_orientations_pbc = []
        for i, line in enumerate(self.contents):
            if line.startswith("Standard orientation:"):
                standard_orientation = []
                standard_orientation_pbc = []
                for j_line in self.contents[i + 5 :]:
                    if "-----------------" in j_line:
                        break
                    if j_line.split()[1] == "-2":  # atomic number = -2 for TV
                        standard_orientation_pbc.append(
                            [float(val) for val in j_line.split()[3:6]]
                        )
                    else:
                        standard_orientation.append(
                            [float(val) for val in j_line.split()[3:6]]
                        )
                standard_orientations.append(np.array(standard_orientation))
                if len(standard_orientation_pbc) != 0:
                    standard_orientations_pbc.append(
                        np.array(standard_orientation_pbc)
                    )
                else:
                    standard_orientations_pbc.append(None)
        if len(standard_orientations) == 0:
            return None, None
        return standard_orientations, standard_orientations_pbc

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
        Read TDDFT transitions and return the
        transition energies in eV as a list.
        """
        excitation_energies_eV = []
        for i in self.tddft_transitions:
            excitation_energies_eV.append(i[0])
        return excitation_energies_eV

    @cached_property
    def absorptions_in_nm(self):
        """
        Read TDDFT transitions and return the
        absorbed wavelengths in nm as a list.
        """
        absorptions_in_nm = []
        for i in self.tddft_transitions:
            absorptions_in_nm.append(i[1])
        return absorptions_in_nm

    @cached_property
    def oscillatory_strengths(self):
        """
        Read TDDFT transitions and return the oscillatory strengths as a list.
        """
        oscillatory_strengths = []
        for i in self.tddft_transitions:
            oscillatory_strengths.append(i[2])
        return oscillatory_strengths

    @cached_property
    def transitions(self):
        """
        Read TDDFT transitions and return the MO transitions.
        """
        transitions, _ = self._read_transitions_and_contribution_coefficients()
        return transitions

    @cached_property
    def contribution_coefficients(self):
        """
        Read MO contribution coefficients.
        """
        _, cc = self._read_transitions_and_contribution_coefficients()
        return cc

    def _read_transitions_and_contribution_coefficients(self):
        from chemsmart.utils.repattern import gaussian_tddft_transition_pattern

        transitions = []
        contribution_coefficients = []
        td_transition_pattern = re.compile(gaussian_tddft_transition_pattern)

        i = 0
        n = len(self.contents)

        while i < n:
            line = self.contents[i]
            if line.lstrip().startswith("Excited State"):
                each_state_transitions = []
                each_state_contribution_coefficients = []

                j = i + 1
                while j < n:
                    current = self.contents[j]

                    # stop on truly blank line
                    if not current.strip():
                        break

                    match = td_transition_pattern.match(current)
                    if match:
                        from_mo, arrow, to_mo, coeff = match.groups()
                        each_state_transitions.append(
                            f"{from_mo} {arrow} {to_mo}"
                        )
                        each_state_contribution_coefficients.append(
                            float(coeff)
                        )
                        j += 1
                        continue

                    # once transition lines have started, stop at first non-transition line
                    if each_state_transitions:
                        break

                    j += 1

                transitions.append(each_state_transitions)
                contribution_coefficients.append(
                    each_state_contribution_coefficients
                )

                i = j
            else:
                i += 1
        return transitions, contribution_coefficients

    @cached_property
    def contributions(self):
        """
        Return the contributions of each molecular orbital excitation
        to an excited state.

        The base value for each contribution coefficient is calculated as
        ``(coef**2) * 100``.

        The returned values then depend on ``self.spin``:

        - ``"restricted"``: the base percentages are doubled.
        - ``"unrestricted"``: the base percentages are unchanged.
        - any other value: a ``ValueError`` is raised.
        """
        contribution_percentage = [
            [(coef**2) * 100 for coef in cc]
            for cc in self.contribution_coefficients
        ]

        if self.spin == "restricted":
            logger.debug(
                "Closed-shell system: contribution percentage is doubled."
            )
            factor = 2
        elif self.spin == "unrestricted":
            logger.debug(
                "Unrestricted system: contribution percentage uses base values."
            )
            factor = 1
        else:
            raise ValueError(f"Unknown spin type: {self.spin!r}")

        return [
            [round(value * factor, 1) for value in percentage]
            for percentage in contribution_percentage
        ]

    @cached_property
    def alpha_occ_eigenvalues(self):
        """
        Obtain all eigenenergies of the alpha
        occupied orbitals and convert to eV.
        """
        alpha_occ_eigenvalues = []

        # Iterate through lines in reverse to
        # find the last block of eigenvalues
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
        """
        Obtain all eigenenergies of the alpha unoccupied orbitals.
        Units of eV, as for orbital energies.
        """

        # Iterate through lines in reverse to
        # find the last block of eigenvalues
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
        """
        Obtain all eigenenergies of the beta occupied orbitals.
        Units of eV, as for orbital energies.
        """
        # Iterate through lines in reverse to
        # find the last block of eigenvalues
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
        """
        Obtain all eigenenergies of the beta unoccupied orbitals.
        Units of eV, as for orbital energies.
        """

        # Iterate through lines in reverse to
        # find the last block of eigenvalues
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
    def mulliken_atomic_charges(self):
        mulliken_atomic_charges, _ = (
            self._get_mulliken_atomic_charges_and_spin_densities()
        )
        return mulliken_atomic_charges

    @cached_property
    def mulliken_spin_densities(self):
        _, mulliken_spin_densities = (
            self._get_mulliken_atomic_charges_and_spin_densities()
        )
        return mulliken_spin_densities

    def _get_mulliken_atomic_charges_and_spin_densities(self):
        """
        Obtain Mulliken charges from the output file.
        """
        all_mulliken_atomic_charges = []
        all_mulliken_spin_densities = []
        for i, line_i in enumerate(self.contents):
            mulliken_atomic_charges = {}
            mulliken_spin_densities = {}
            if line_i.startswith(
                ("Mulliken charges:", "Mulliken charges and spin densities:")
            ):
                for line_j in self.contents[i + 2 :]:
                    if "Sum of Mulliken charges" in line_j:
                        break
                    line_j_elements = line_j.split()
                    element = p.to_element(line_j_elements[1])
                    element_num = f"{element}{line_j_elements[0]}"
                    mulliken_atomic_charges[element_num] = float(
                        line_j_elements[2]
                    )
                    if len(line_j_elements) == 4:
                        mulliken_spin_densities[element_num] = float(
                            line_j_elements[3]
                        )
                all_mulliken_atomic_charges.append(mulliken_atomic_charges)
                if mulliken_spin_densities:
                    all_mulliken_spin_densities.append(mulliken_spin_densities)
        if all_mulliken_atomic_charges and all_mulliken_spin_densities:
            return (
                all_mulliken_atomic_charges[-1],
                all_mulliken_spin_densities[-1],
            )
        elif all_mulliken_atomic_charges:
            return all_mulliken_atomic_charges[-1], None
        else:
            return None, None

    @cached_property
    def mulliken_atomic_charges_heavy_atoms(self):
        mulliken_atomic_charges_heavy_atoms, _ = (
            self._get_mulliken_atomic_charges_and_spin_densities_heavy_atoms()
        )
        return mulliken_atomic_charges_heavy_atoms

    @cached_property
    def mulliken_spin_densities_heavy_atoms(self):
        _, mulliken_spin_densities_heavy_atoms = (
            self._get_mulliken_atomic_charges_and_spin_densities_heavy_atoms()
        )
        return mulliken_spin_densities_heavy_atoms

    def _get_mulliken_atomic_charges_and_spin_densities_heavy_atoms(self):
        """
        Obtain Mulliken charges with hydrogens summed into heavy atoms.
        """
        all_mulliken_atomic_charges_heavy_atoms = []
        all_mulliken_spin_densities_heavy_atoms = []
        for i, line_i in enumerate(self.contents):
            mulliken_atomic_charges_heavy_atoms = {}
            mulliken_spin_densities_heavy_atoms = {}
            if line_i.startswith(
                (
                    "Mulliken charges with hydrogens summed into heavy atoms:",
                    "Mulliken charges and spin densities with hydrogens summed into heavy atoms:",
                )
            ):
                for line_j in self.contents[i + 2 :]:
                    if "Electronic spatial extent" in line_j:
                        break
                    line_j_elements = line_j.split()
                    element = p.to_element(line_j_elements[1])
                    element_num = f"{element}{line_j_elements[0]}"
                    mulliken_atomic_charges_heavy_atoms[element_num] = float(
                        line_j_elements[2]
                    )
                    if len(line_j_elements) == 4:
                        mulliken_spin_densities_heavy_atoms[element_num] = (
                            float(line_j_elements[3])
                        )
                all_mulliken_atomic_charges_heavy_atoms.append(
                    mulliken_atomic_charges_heavy_atoms
                )
                if mulliken_spin_densities_heavy_atoms:
                    all_mulliken_spin_densities_heavy_atoms.append(
                        mulliken_spin_densities_heavy_atoms
                    )
        if (
            all_mulliken_atomic_charges_heavy_atoms
            and all_mulliken_spin_densities_heavy_atoms
        ):
            return (
                all_mulliken_atomic_charges_heavy_atoms[-1],
                all_mulliken_spin_densities_heavy_atoms[-1],
            )
        elif all_mulliken_atomic_charges_heavy_atoms:
            return all_mulliken_atomic_charges_heavy_atoms[-1], None
        else:
            # if spin densities present, charges must be present too
            return None, None

    @cached_property
    def hirshfeld_charges(self):
        hirshfeld_charges, _, _, _ = (
            self._get_hirshfeld_charges_spins_dipoles_cm5()
        )
        return hirshfeld_charges

    @cached_property
    def hirshfeld_spin_densities(self):
        _, spin_densities, _, _ = (
            self._get_hirshfeld_charges_spins_dipoles_cm5()
        )
        return spin_densities

    @cached_property
    def hirshfeld_dipoles(self):
        _, _, dipoles, _ = self._get_hirshfeld_charges_spins_dipoles_cm5()
        return dipoles

    @cached_property
    def hirshfeld_cm5_charges(self):
        _, _, _, cm5_charges = self._get_hirshfeld_charges_spins_dipoles_cm5()
        return cm5_charges

    @cached_property
    def hirshfeld_charges_heavy_atoms(self):
        hirshfeld_charges_heavy_atoms, _, _ = (
            self._get_hirshfeld_charges_spin_densities_cm5_charges_heavy_atoms()
        )
        return hirshfeld_charges_heavy_atoms

    @cached_property
    def hirshfeld_spin_densities_heavy_atoms(self):
        _, hirshfeld_spin_densities_heavy_atoms, _ = (
            self._get_hirshfeld_charges_spin_densities_cm5_charges_heavy_atoms()
        )
        return hirshfeld_spin_densities_heavy_atoms

    @cached_property
    def hirshfeld_cm5_charges_heavy_atoms(self):
        _, _, cm5_charges_heavy_atoms = (
            self._get_hirshfeld_charges_spin_densities_cm5_charges_heavy_atoms()
        )
        return cm5_charges_heavy_atoms

    def _get_hirshfeld_charges_spins_dipoles_cm5(self):
        """
        Obtain Hirshfeld charges, spin densities,
        dipoles, and CM5 charges from the output file.
        """
        all_hirshfeld_charges = []
        all_spin_densities = []
        all_dipoles = []
        all_cm5_charges = []
        for i, line_i in enumerate(self.contents):
            hirshfeld_charges = {}
            spin_densities = {}
            dipoles = {}
            cm5_charges = {}
            if (
                "Hirshfeld charges, spin densities, dipoles, and CM5 charges"
                in line_i
            ):
                for line_j in self.contents[i + 2 :]:
                    if line_j.startswith("Tot") or len(line_j) == 0:
                        break
                    line_j_elements = line_j.split()
                    element = p.to_element(line_j_elements[1])
                    element_num = f"{element}{line_j_elements[0]}"
                    hirshfeld_charges[element_num] = float(line_j_elements[2])
                    spin_densities[element_num] = float(line_j_elements[3])
                    dipoles[element_num] = np.array(
                        [
                            float(line_j_elements[4]),
                            float(line_j_elements[5]),
                            float(line_j_elements[6]),
                        ]
                    )
                    cm5_charges[element_num] = float(line_j_elements[7])
                all_hirshfeld_charges.append(hirshfeld_charges)
                all_spin_densities.append(spin_densities)
                all_dipoles.append(dipoles)
                all_cm5_charges.append(cm5_charges)
        return (
            all_hirshfeld_charges[-1],
            all_spin_densities[-1],
            all_dipoles[-1],
            all_cm5_charges[-1],
        )

    def _get_hirshfeld_charges_spin_densities_cm5_charges_heavy_atoms(self):
        """
        Obtain Hirshfeld charges, spin densities and
        CM5 with hydrogens summed into heavy atoms.
        """
        all_hirshfeld_charges_heavy_atoms = []
        all_hirshfeld_spin_densities_heavy_atoms = []
        all_cm5_charges_heavy_atoms = []
        for i, line_i in enumerate(self.contents):
            hirshfeld_charges_heavy_atoms = {}
            hirshfeld_spin_densities_heavy_atoms = {}
            cm5_charges_heavy_atoms = {}
            if line_i.startswith(
                (
                    "Hirshfeld charges with hydrogens summed into heavy atoms:",
                    "Hirshfeld charges and spin densities with hydrogens summed into heavy atoms:",
                )
            ):
                for line_j in self.contents[i + 2 :]:
                    if len(line_j) == 0 or line_j.startswith("Tot"):
                        break
                    line_j_elements = line_j.split()
                    element = p.to_element(line_j_elements[1])
                    element_num = f"{element}{line_j_elements[0]}"
                    hirshfeld_charges_heavy_atoms[element_num] = float(
                        line_j_elements[2]
                    )
                    if len(line_j_elements) == 4:
                        cm5_charges_heavy_atoms[element_num] = float(
                            line_j_elements[3]
                        )
                    elif len(line_j_elements) == 5:
                        hirshfeld_spin_densities_heavy_atoms[element_num] = (
                            float(line_j_elements[3])
                        )
                        cm5_charges_heavy_atoms[element_num] = float(
                            line_j_elements[4]
                        )
                all_hirshfeld_charges_heavy_atoms.append(
                    hirshfeld_charges_heavy_atoms
                )
                all_cm5_charges_heavy_atoms.append(cm5_charges_heavy_atoms)
                if hirshfeld_spin_densities_heavy_atoms:
                    all_hirshfeld_spin_densities_heavy_atoms.append(
                        hirshfeld_spin_densities_heavy_atoms
                    )

        if (
            all_hirshfeld_charges_heavy_atoms
            and all_hirshfeld_spin_densities_heavy_atoms
        ):
            return (
                all_hirshfeld_charges_heavy_atoms[-1],
                all_hirshfeld_spin_densities_heavy_atoms[-1],
                all_cm5_charges_heavy_atoms,
            )
        elif (
            all_hirshfeld_charges_heavy_atoms
            and not all_hirshfeld_spin_densities_heavy_atoms
        ):
            return (
                all_hirshfeld_charges_heavy_atoms[-1],
                None,
                all_cm5_charges_heavy_atoms[-1],
            )
        else:
            return None, None, None

    def get_molecule(self, index="-1"):
        index = string2index_1based(index)
        return self.all_structures[index]

    @property
    def mass(self):
        for line in self.contents:
            if "Molecular mass:" and "amu." in line:
                return float(line.split()[2])

    @cached_property
    def moments_of_inertia(self):
        """
        Obtain moments of inertia from the
        output file which are in atomic units
        (amu * Bohr^2) and convert to SI units (kg * m^2).
        """
        moments_of_inertia, _ = (
            self._get_moments_of_inertia_and_principal_axes()
        )
        return moments_of_inertia

    @cached_property
    def moments_of_inertia_principal_axes(self):
        _, principal_axes = self._get_moments_of_inertia_and_principal_axes()
        return principal_axes

    def _get_moments_of_inertia_and_principal_axes(self):
        """
        Obtain moments of inertia along principal axes from the output file
        (amu * Bohr^2 in Gaussian) and convert to units of (amu * Å^2).
        """
        for i, line in enumerate(self.contents):
            if "Principal axes and moments of inertia" in line:
                moments_of_inertia = []
                moments_of_inertia_principal_axes = []
                for j_line in self.contents[i + 2 :]:
                    if j_line.startswith("This molecule"):
                        break
                    if j_line.startswith("Eigenvalue"):
                        for eigenval in j_line.split("Eigenvalues --")[
                            -1
                        ].split():
                            try:
                                moments_of_inertia.append(
                                    float(eigenval) * units.Bohr**2
                                )
                            except ValueError:
                                logger.warning(
                                    f"Could not convert '{j_line}' due to "
                                    f"Gaussian incorrect printing."
                                )
                                moments_of_inertia.append(
                                    np.array([np.inf] * 3)
                                )
                    else:
                        if len(j_line.split()) == 4:
                            moments_of_inertia_principal_axes.append(
                                np.array(j_line.split()[1:4], dtype=float)
                            )
                moments_of_inertia_principal_axes = np.array(
                    moments_of_inertia_principal_axes
                ).transpose()
                return np.array(moments_of_inertia), np.array(
                    moments_of_inertia_principal_axes
                )

    @cached_property
    def rotational_symmetry_number(self):
        """
        Obtain the rotational symmetry number from the output file.
        """
        for line in self.contents:
            if "Rotational symmetry number" in line:
                return int(line.split()[-1].split(".")[0])

    @cached_property
    def rotational_temperatures(self):
        """
        Rotational temperatures in Kelvin, as a list.
        """
        rot_temps = []
        for line in reversed(self.contents):
            # take from the end of outputfile
            if "Rotational temperature" in line and "(Kelvin)" in line:
                for rot_temp in line.split("(Kelvin)")[-1].split():
                    # linear molecules may have only one rot temp,
                    # non-linear has three
                    rot_temps.append(float(rot_temp))
                return rot_temps

    @cached_property
    def rotational_constants_in_Hz(self):
        """
        Rotational constants in Hz, as a list.
        """
        rot_consts = []
        for line in reversed(self.contents):
            # take from the end of outputfile
            if "Rotational constant" in line and "(GHZ):" in line:
                for rot_const in line.split("(GHZ):")[-1].split():
                    rot_consts.append(float(rot_const) * 1e9)
                return rot_consts

    def to_dataset(self, **kwargs):
        """
        Convert Gaussian .log file to Dataset with
        all data points taken from the .log file.

        Returns:
            Dataset.
        """
        # TODO: to be implemented
        pass

    @cached_property
    def oniom_partition(self):
        """Obtain the atomic indices of each layer in the ONIOM calculation.
        Returns:
            indices of each layer as a dictionary"""
        high_level = []
        medium_level = []
        low_level = []
        for i, line in enumerate(self.contents):
            if "Symbolic Z-matrix:" in line:
                # First ONIOM format: coordinates
                # start 4 lines after the header
                if "Charge" not in self.contents[i + 4]:
                    atom_index = 1
                    for j_line in self.contents[i + 4 :]:
                        if len(j_line) == 0:
                            break
                        if len(j_line) > 4:
                            tokens = j_line.split()
                            if tokens[1] == "-1" or tokens[1] == "0":
                                layer = str(tokens[5])
                            else:
                                layer = str(tokens[4])
                            if layer == "H":
                                high_level.append(atom_index)
                            elif layer == "M":
                                medium_level.append(atom_index)
                            elif layer == "L":
                                low_level.append(atom_index)
                            atom_index += 1
                else:
                    # Alternative ONIOM format: coordinates
                    # start 7 lines after the header
                    atom_index = 1
                    for j_line in self.contents[i + 7 :]:
                        if len(j_line) == 0:
                            break
                        if len(j_line) > 4:
                            tokens = j_line.split()
                            if tokens[1] == "-1" or tokens[1] == "0":
                                layer = str(tokens[5])
                            else:
                                layer = str(tokens[4])
                            if layer == "H":
                                high_level.append(atom_index)
                            elif layer == "M":
                                medium_level.append(atom_index)
                            elif layer == "L":
                                low_level.append(atom_index)
                            atom_index += 1
        partition = {}
        for level_name, level_list in [
            ("high level atoms", high_level),
            ("medium level atoms", medium_level),
            ("low level atoms", low_level),
        ]:
            if len(level_list) != 0:
                partition[level_name] = get_range_from_list(level_list)
        return partition

    @cached_property
    def oniom_cutting_bonds(self):
        """Obtain the cutting bonds in the ONIOM calculation.
        Returns:
            cutting bonds as a dictionary"""
        cutting_bonds = {}
        for i, line in enumerate(self.contents):
            if "Cut between" in line:
                atom1 = int(self.contents[i].split()[5])
                atom2 = int(self.contents[i].split()[8])
                factor1 = float(self.contents[i].split()[10])
                factor2 = float(self.contents[i].split()[11])
                cutting_bonds[(atom1, atom2)] = (factor1, factor2)
        return cutting_bonds

    @cached_property
    def oniom_get_charge_and_multiplicity(self):
        """Obtain the charge and multiplicity
        of the system in the ONIOM calculation.
        Returns:
            charge and multiplicity as a dictionary"""
        charge_multiplicity = {}
        for line in self.contents:
            if "Charge" in line and "low   level calculation on real" in line:
                charge_multiplicity["low-level, real system"] = (
                    int(line.split()[2]),
                    int(line.split()[5]),
                )
            if "Charge" in line and "med   level calculation on mid" in line:
                charge_multiplicity["medium-level, mid system"] = (
                    int(line.split()[2]),
                    int(line.split()[5]),
                )
            if "Charge" in line and "low   level calculation on mid" in line:
                charge_multiplicity["low-level, mid system"] = (
                    int(line.split()[2]),
                    int(line.split()[5]),
                )
            if "Charge" in line and "high  level calculation on model" in line:
                charge_multiplicity["high-level, model system"] = (
                    int(line.split()[2]),
                    int(line.split()[5]),
                )
            if "Charge" in line and "med   level calculation on model" in line:
                charge_multiplicity["medium-level, model system"] = (
                    int(line.split()[2]),
                    int(line.split()[5]),
                )
            if "Charge" in line and "low   level calculation on model" in line:
                charge_multiplicity["low-level, model system"] = (
                    int(line.split()[2]),
                    int(line.split()[5]),
                )
        return charge_multiplicity


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
        """
        Parse the NBO natural atomic orbitals.
        """
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
                        # Construct the sub-key for each NAO
                        # entry, e.g., "NAO1", "NAO2", etc.
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
        """
        Parse the NBO natural population analysis.
        """
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
        """
        Get natural charges corresponding to each atom as a dictionary.
        """
        natural_charges = {}
        for atom_key, atom_data in self.natural_population_analysis.items():
            natural_charges[atom_key] = atom_data["natural_charge"]
        return natural_charges

    @cached_property
    def total_electrons(self):
        """
        Get the total number of electrons
        corresponding to each atom as a dictionary.
        """
        total_electrons = {}
        for atom_key, atom_data in self.natural_population_analysis.items():
            total_electrons[atom_key] = atom_data["total_electrons"]
        return total_electrons

    @cached_property
    def electronic_configuration(self):
        """
        Get electronic configuration for each
        atom and store results in a dictionary.
        """
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
        """
        Get the number of NAOs for a given atom.
        """
        return len(self.natural_atomic_orbitals[atom_key])

    def get_total_electron_occ(self, atom_key):
        """
        Get the total electron occupancy for a given atom.
        """
        total_electron_occ = sum(
            entry["occupancy"]
            for entry in self.natural_atomic_orbitals[atom_key].values()
        )
        return total_electron_occ

    def get_electronic_configuration(self, atom_key):
        """
        Get the electronic configuration for a given atom.
        """
        return self.electronic_configuration[atom_key]


class Gaussian16OutputWithPBC(Gaussian16Output):
    """
    class for parsing and obtaining information
    from Gaussian output file with PBC.
    """

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
        """
        Get final translation vectors from last step.
        """
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


class Gaussian16pKaOutput(Gaussian16Output):
    """
    Extended Gaussian16Output for pKa calculations with thermochemistry support.

    This class provides methods to extract electronic energy and quasi-harmonic
    Gibbs free energy from Gaussian optimization output files, which are essential
    for pKa calculations using thermodynamic cycles.

    The thermochemistry calculations use Grimme's quasi-RRHO method for entropy
    and Head-Gordon's quasi-RRHO method for enthalpy corrections, matching the
    behavior of:
        chemsmart run thermochemistry -f <file> -T <temp> -c <conc> -csg <cutoff> -ch <cutoff>

    Attributes:
        filename (str): Path to the Gaussian output file.
        temperature (float): Temperature in Kelvin for thermochemistry. Default 298.15 K.
        concentration (float): Concentration in mol/L. Default 1.0 mol/L.
        pressure (float): Pressure in atm. Default 1.0 atm.
        cutoff_entropy_grimme (float): Cutoff frequency for entropy (cm^-1). Default 100.0.
        cutoff_enthalpy (float): Cutoff frequency for enthalpy (cm^-1). Default 100.0.
        energy_units (str): Energy units for output. Default 'hartree'.

    Example:
        output = Gaussian16pKaOutput(
            "acetic_acid_opt.log",
            temperature=333.15,
            concentration=1.0,
            cutoff_entropy_grimme=100,
            cutoff_enthalpy=100
        )
        E = output.electronic_energy_in_units  # E in hartree
        G = output.qh_gibbs_free_energy  # qh-G(T) in hartree
    """

    def __init__(
        self,
        filename,
        temperature=298.15,
        concentration=1.0,
        pressure=1.0,
        cutoff_entropy_grimme=100.0,
        cutoff_enthalpy=100.0,
        energy_units="hartree",
    ):
        """
        Initialize Gaussian16pKaOutput with thermochemistry settings.

        Args:
            filename (str): Path to Gaussian output file.
            temperature (float): Temperature in Kelvin. Default 298.15 K.
            concentration (float): Concentration in mol/L. Default 1.0 mol/L.
            pressure (float): Pressure in atm. Default 1.0 atm.
            cutoff_entropy_grimme (float): Cutoff frequency for entropy
                in cm^-1 using Grimme's quasi-RRHO method. Default 100.0.
            cutoff_enthalpy (float): Cutoff frequency for enthalpy
                in cm^-1 using Head-Gordon's method. Default 100.0.
            energy_units (str): Energy units for output values.
                Options: 'hartree', 'eV', 'kcal/mol', 'kJ/mol'. Default 'hartree'.
        """
        super().__init__(filename=filename)
        self.temperature = temperature
        self.concentration = concentration
        self.pressure = pressure
        self.cutoff_entropy_grimme = cutoff_entropy_grimme
        self.cutoff_enthalpy = cutoff_enthalpy
        self.energy_units = energy_units.lower()
        self._thermochemistry = None

    @property
    def thermochemistry(self):
        """
        Get or create the Thermochemistry analysis object.

        Returns:
            Thermochemistry: Configured thermochemistry analysis object.

        Raises:
            ValueError: If the output file did not terminate normally.
        """
        if self._thermochemistry is None:
            from chemsmart.analysis.thermochemistry import Thermochemistry

            self._thermochemistry = Thermochemistry(
                filename=self.filename,
                temperature=self.temperature,
                concentration=self.concentration,
                pressure=self.pressure,
                use_weighted_mass=False,
                alpha=4,
                s_freq_cutoff=self.cutoff_entropy_grimme,
                entropy_method="grimme",
                h_freq_cutoff=self.cutoff_enthalpy,
                energy_units=self.energy_units,
                check_imaginary_frequencies=True,
            )
        return self._thermochemistry

    @property
    def electronic_energy_in_units(self):
        """
        Get the electronic energy (E) in specified units.

        This is the raw SCF energy from the Gaussian calculation,
        converted to the specified energy units.

        Returns:
            float: Electronic energy in specified units (default: hartree).
        """
        from chemsmart.utils.constants import energy_conversion

        # Get electronic energy in J/mol from thermochemistry
        electronic_energy_j_mol = self.thermochemistry.electronic_energy
        # Convert to specified units
        return energy_conversion(
            "j/mol", self.energy_units, electronic_energy_j_mol
        )

    @property
    def qh_gibbs_free_energy(self):
        """
        Get the quasi-harmonic Gibbs free energy qh-G(T) in specified units.

        This uses Grimme's quasi-RRHO method for entropy and Head-Gordon's
        quasi-RRHO method for enthalpy corrections, which is equivalent to
        running:
            chemsmart run thermochemistry -f <file> -T <temp> -c <conc> -csg <cutoff> -ch <cutoff>

        The qh-G(T) value corresponds to the quasi-RRHO corrected Gibbs free energy
        that accounts for low-frequency vibrations using interpolation to free rotor
        entropy and enthalpy.

        Returns:
            float: Quasi-harmonic Gibbs free energy in specified units (default: hartree).

        Raises:
            ValueError: If the file doesn't contain frequency data.
        """
        from chemsmart.utils.constants import energy_conversion

        # Get qh-G in J/mol from thermochemistry
        qh_gibbs_j_mol = self.thermochemistry.qrrho_gibbs_free_energy
        if qh_gibbs_j_mol is None:
            raise ValueError(
                f"Cannot compute qh-Gibbs free energy for {self.filename}. "
                "The file may not contain frequency calculation data."
            )
        # Convert to specified units
        return energy_conversion("j/mol", self.energy_units, qh_gibbs_j_mol)

    @property
    def zero_point_energy_in_units(self):
        """
        Get the zero-point energy (ZPE) in specified units.

        Returns:
            float: Zero-point energy in specified units (default: hartree).
        """
        from chemsmart.utils.constants import energy_conversion

        zpe_j_mol = self.thermochemistry.zero_point_energy
        if zpe_j_mol is None:
            raise ValueError(
                f"Cannot compute zero-point energy for {self.filename}. "
                "The file may not contain frequency calculation data."
            )
        return energy_conversion("j/mol", self.energy_units, zpe_j_mol)

    @property
    def enthalpy_in_units(self):
        """
        Get the enthalpy (H) in specified units.

        Returns:
            float: Enthalpy in specified units (default: hartree).
        """
        from chemsmart.utils.constants import energy_conversion

        enthalpy_j_mol = self.thermochemistry.enthalpy
        if enthalpy_j_mol is None:
            raise ValueError(
                f"Cannot compute enthalpy for {self.filename}. "
                "The file may not contain frequency calculation data."
            )
        return energy_conversion("j/mol", self.energy_units, enthalpy_j_mol)

    @property
    def qh_enthalpy_in_units(self):
        """
        Get the quasi-harmonic enthalpy qh-H(T) in specified units.

        Returns:
            float: Quasi-harmonic enthalpy in specified units (default: hartree).
        """
        from chemsmart.utils.constants import energy_conversion

        qh_enthalpy_j_mol = self.thermochemistry.qrrho_enthalpy
        if qh_enthalpy_j_mol is None:
            raise ValueError(
                f"Cannot compute qh-enthalpy for {self.filename}. "
                "The file may not contain frequency calculation data."
            )
        return energy_conversion("j/mol", self.energy_units, qh_enthalpy_j_mol)

    @property
    def gibbs_free_energy_in_units(self):
        """
        Get the standard Gibbs free energy G(T) in specified units.

        This is the uncorrected Gibbs free energy without quasi-RRHO corrections.

        Returns:
            float: Gibbs free energy in specified units (default: hartree).
        """
        from chemsmart.utils.constants import energy_conversion

        gibbs_j_mol = self.thermochemistry.gibbs_free_energy
        if gibbs_j_mol is None:
            raise ValueError(
                f"Cannot compute Gibbs free energy for {self.filename}. "
                "The file may not contain frequency calculation data."
            )
        return energy_conversion("j/mol", self.energy_units, gibbs_j_mol)

    @property
    def thermochemical_properties(self):
        """
        Compute and return all thermochemical properties.

        Returns:
            dict: Dictionary containing thermochemical properties:
                - electronic_energy: Electronic energy in specified units
                - zero_point_energy: Zero-point energy in specified units
                - enthalpy: Enthalpy in specified units
                - qh_enthalpy: Quasi-harmonic enthalpy in specified units
                - gibbs_free_energy: Gibbs free energy in specified units
                - qh_gibbs_free_energy: Quasi-harmonic Gibbs free energy in specified units
        """
        return {
            "electronic_energy": self.electronic_energy_in_units,
            "zero_point_energy": self.zero_point_energy_in_units,
            "enthalpy": self.enthalpy_in_units,
            "qh_enthalpy": self.qh_enthalpy_in_units,
            "gibbs_free_energy": self.gibbs_free_energy_in_units,
            "qh_gibbs_free_energy": self.qh_gibbs_free_energy,
        }

    def compute_thermochemistry(self):
        """
        Compute all thermochemistry properties.

        Returns:
            dict: Dictionary containing all thermochemistry values:
                - structure: Base filename
                - electronic_energy: E in specified units
                - zero_point_energy: ZPE in specified units
                - enthalpy: H in specified units
                - qh_enthalpy: qh-H(T) in specified units
                - entropy_times_temperature: T*S in specified units
                - qh_entropy_times_temperature: T*qh-S in specified units
                - gibbs_free_energy: G(T) in specified units
                - qh_gibbs_free_energy: qh-G(T) in specified units
        """
        import os

        from chemsmart.utils.constants import energy_conversion

        thermo = self.thermochemistry
        structure = os.path.splitext(os.path.basename(self.filename))[0]

        return {
            "structure": structure,
            "electronic_energy": self.electronic_energy_in_units,
            "zero_point_energy": self.zero_point_energy_in_units,
            "enthalpy": self.enthalpy_in_units,
            "qh_enthalpy": self.qh_enthalpy_in_units,
            "entropy_times_temperature": (
                energy_conversion(
                    "j/mol",
                    self.energy_units,
                    thermo.entropy_times_temperature,
                )
                if thermo.entropy_times_temperature
                else None
            ),
            "qh_entropy_times_temperature": (
                energy_conversion(
                    "j/mol",
                    self.energy_units,
                    thermo.qrrho_entropy_times_temperature,
                )
                if thermo.qrrho_entropy_times_temperature
                else None
            ),
            "gibbs_free_energy": self.gibbs_free_energy_in_units,
            "qh_gibbs_free_energy": self.qh_gibbs_free_energy,
        }

    @classmethod
    def from_settings(cls, filename, settings):
        """
        Create a Gaussian16pKaOutput from GaussianpKaJobSettings.

        Args:
            filename (str): Path to Gaussian output file.
            settings (GaussianpKaJobSettings): pKa job settings containing
                thermochemistry parameters.

        Returns:
            Gaussian16pKaOutput: Configured output object.
        """
        return cls(
            filename=filename,
            temperature=settings.temperature,
            concentration=settings.concentration,
            pressure=settings.pressure,
            cutoff_entropy_grimme=settings.cutoff_entropy_grimme,
            cutoff_enthalpy=settings.cutoff_enthalpy,
            energy_units=settings.energy_units,
        )

    # =========================================================================
    # Multi-file pKa thermochemistry support
    # =========================================================================

    @classmethod
    def compute_pka_thermochemistry(
        cls,
        ha_file=None,
        a_file=None,
        href_file=None,
        ref_file=None,
        temperature=298.15,
        concentration=1.0,
        pressure=1.0,
        cutoff_entropy_grimme=100.0,
        cutoff_enthalpy=100.0,
        energy_units="hartree",
    ):
        """
        Compute thermochemistry for pKa species (HA, A-, HB, B-).

        This helper method extracts thermochemical properties from gas-phase
        frequency calculation files for all species involved in pKa calculations.

        Args:
            ha_file (str, optional): Path to HA (protonated acid) output file.
            a_file (str, optional): Path to A- (conjugate base) output file.
            href_file (str, optional): Path to HB (reference acid) output file.
            ref_file (str, optional): Path to B- (reference conjugate base) output file.
            temperature (float): Temperature in Kelvin. Default 298.15 K.
            concentration (float): Concentration in mol/L. Default 1.0 mol/L.
            pressure (float): Pressure in atm. Default 1.0 atm.
            cutoff_entropy_grimme (float): Cutoff for entropy (cm^-1). Default 100.0.
            cutoff_enthalpy (float): Cutoff for enthalpy (cm^-1). Default 100.0.
            energy_units (str): Energy units for output. Default 'hartree'.

        Returns:
            dict: Dictionary containing thermochemistry for each species and settings:
                - 'settings': Dict of calculation settings
                - 'HA': Dict with E, qh_G, ZPE, H, qh_H, G, name for protonated acid
                - 'A': Dict with same keys for conjugate base
                - 'HRef': Dict for reference acid (if provided)
                - 'Ref': Dict for reference conjugate base (if provided)
        """
        results = {
            "settings": {
                "temperature": temperature,
                "concentration": concentration,
                "pressure": pressure,
                "cutoff_entropy_grimme": cutoff_entropy_grimme,
                "cutoff_enthalpy": cutoff_enthalpy,
                "energy_units": energy_units,
            }
        }

        common_kwargs = {
            "temperature": temperature,
            "concentration": concentration,
            "pressure": pressure,
            "cutoff_entropy_grimme": cutoff_entropy_grimme,
            "cutoff_enthalpy": cutoff_enthalpy,
            "energy_units": energy_units,
        }

        def get_species_thermo(file_path, name):
            """Extract thermochemistry for a species."""
            if file_path is None:
                return None
            output = cls(filename=file_path, **common_kwargs)
            return {
                "name": name,
                "E": output.electronic_energy_in_units,
                "qh_G": output.qh_gibbs_free_energy,
                "ZPE": output.zero_point_energy_in_units,
                "H": output.enthalpy_in_units,
                "qh_H": output.qh_enthalpy_in_units,
                "G": output.gibbs_free_energy_in_units,
            }

        if ha_file is not None:
            results["HA"] = get_species_thermo(ha_file, "HA")
        if a_file is not None:
            results["A"] = get_species_thermo(a_file, "A-")
        if href_file is not None:
            results["HRef"] = get_species_thermo(href_file, "HRef")
        if ref_file is not None:
            results["Ref"] = get_species_thermo(ref_file, "Ref-")

        return results

    @classmethod
    def from_pka_settings(
        cls,
        settings,
        ha_file=None,
        a_file=None,
        hb_file=None,
        b_file=None,
    ):
        """
        Create Gaussian16pKaOutput objects from GaussianpKaJobSettings.

        Factory method that creates output objects for pKa species using
        settings from a GaussianpKaJobSettings object.

        Args:
            settings: GaussianpKaJobSettings object containing thermochemistry
                parameters (temperature, concentration, cutoffs, etc.).
            ha_file (str, optional): Path to HA (protonated acid) output file.
            a_file (str, optional): Path to A- (conjugate base) output file.
            hb_file (str, optional): Path to HB (reference acid) output file.
            b_file (str, optional): Path to B- (reference conjugate base) output file.

        Returns:
            dict: Dictionary with Gaussian16pKaOutput objects for each species.
        """
        return cls.for_pka_species(
            ha_file=ha_file,
            a_file=a_file,
            hb_file=hb_file,
            b_file=b_file,
            temperature=settings.temperature,
            concentration=settings.concentration,
            pressure=settings.pressure,
            cutoff_entropy_grimme=settings.cutoff_entropy_grimme,
            cutoff_enthalpy=settings.cutoff_enthalpy,
            energy_units=settings.energy_units,
        )

    @classmethod
    def for_pka_species(
        cls,
        ha_file=None,
        a_file=None,
        hb_file=None,
        b_file=None,
        temperature=298.15,
        concentration=1.0,
        pressure=1.0,
        cutoff_entropy_grimme=100.0,
        cutoff_enthalpy=100.0,
        energy_units="hartree",
    ):
        """
        Create Gaussian16pKaOutput objects for multiple pKa species.

        Factory method that creates output objects for all species involved
        in pKa calculations: HA, A-, and optionally HB, B- for proton exchange.

        Args:
            ha_file (str, optional): Path to HA (protonated acid) output file.
            a_file (str, optional): Path to A- (conjugate base) output file.
            hb_file (str, optional): Path to HB (reference acid) output file.
            b_file (str, optional): Path to B- (reference conjugate base) output file.
            temperature (float): Temperature in Kelvin. Default 298.15 K.
            concentration (float): Concentration in mol/L. Default 1.0 mol/L.
            pressure (float): Pressure in atm. Default 1.0 atm.
            cutoff_entropy_grimme (float): Cutoff frequency for entropy (cm^-1).
            cutoff_enthalpy (float): Cutoff frequency for enthalpy (cm^-1).
            energy_units (str): Energy units for output.

        Returns:
            dict: Dictionary with Gaussian16pKaOutput objects:
                - 'HA': Output for protonated acid (if ha_file provided)
                - 'A': Output for conjugate base (if a_file provided)
                - 'HB': Output for reference acid (if hb_file provided)
                - 'B': Output for reference conjugate base (if b_file provided)

        Example:
            outputs = Gaussian16pKaOutput.for_pka_species(
                ha_file="acid_opt.log",
                a_file="base_opt.log",
                temperature=298.15
            )
            print(f"E(HA) = {outputs['HA'].electronic_energy_in_units}")
            print(f"qh-G(A-) = {outputs['A'].qh_gibbs_free_energy}")
        """
        outputs = {}
        common_kwargs = {
            "temperature": temperature,
            "concentration": concentration,
            "pressure": pressure,
            "cutoff_entropy_grimme": cutoff_entropy_grimme,
            "cutoff_enthalpy": cutoff_enthalpy,
            "energy_units": energy_units,
        }

        if ha_file is not None:
            outputs["HA"] = cls(filename=ha_file, **common_kwargs)
        if a_file is not None:
            outputs["A"] = cls(filename=a_file, **common_kwargs)
        if hb_file is not None:
            outputs["HB"] = cls(filename=hb_file, **common_kwargs)
        if b_file is not None:
            outputs["B"] = cls(filename=b_file, **common_kwargs)

        return outputs

    @classmethod
    def compute_pka(
        cls,
        ha_gas_file,
        a_gas_file,
        href_gas_file,
        ref_gas_file,
        ha_solv_file,
        a_solv_file,
        href_solv_file,
        ref_solv_file,
        pka_reference,
        temperature=298.15,
        concentration=1.0,
        pressure=1.0,
        cutoff_entropy_grimme=100.0,
        cutoff_enthalpy=100.0,
    ):
        """
        Compute pKa using the Dual-level Proton Exchange scheme.

        This is the default and recommended method for pKa calculations.
        It implements a dual-level approach using the proton exchange
        (isodesmic) thermodynamic cycle:

            HA + B⁻ → A⁻ + HB

        where:
            - HA: Target acid (protonated form, e.g., 5PQ_Me_ts1)
            - A⁻: Conjugate base of target acid
            - HB: Reference acid (protonated form, e.g., collidine-H, pKa=6.75)
            - B⁻: Conjugate base of reference acid

        The dual-level approach separates:
        1. **Thermal corrections (G_corr)**: Extracted from gas-phase frequency
           calculations using quasi-harmonic Gibbs free energy:
               G_corr = qh-G(T) - E_gas  [Hartree]

        2. **Solvent energies (E_solv)**: High-level single-point electronic
           energies calculated in implicit solvent (e.g., SMD). [Hartree]

        3. **Total free energy in solution** for each species:
               G_soln = E_solv + G_corr  [Hartree]

        4. **Reaction free energy** using proton exchange:
               ΔG_soln = [G(A⁻)_soln + G(HB)_soln] - [G(HA)_soln + G(B⁻)_soln]
               Converted to kcal/mol for pKa calculation.

        5. **pKa calculation**:
               pKa = pKa_ref + ΔG_soln / (RT × ln10)

        Note: All internal energies are stored in Hartree (au).
        Only ΔG_soln is converted to kcal/mol for the pKa formula.
        Conversion factor: 1 Hartree = 627.5094740631 kcal/mol

        Args:
            ha_gas_file (str): Path to HA gas-phase optimization+freq output file.
            a_gas_file (str): Path to A⁻ gas-phase optimization+freq output file.
            href_gas_file (str): Path to HB gas-phase optimization+freq output file.
            ref_gas_file (str): Path to B⁻ gas-phase optimization+freq output file.
            ha_solv_file (str): Path to HA solvent single-point output file.
            a_solv_file (str): Path to A⁻ solvent single-point output file.
            href_solv_file (str): Path to HB solvent single-point output file.
            ref_solv_file (str): Path to B⁻ solvent single-point output file.
            pka_reference (float): Experimental pKa of the reference acid HB.
                Default reference: collidine (pKa = 6.75).
            temperature (float): Temperature in Kelvin. Default 298.15 K.
            concentration (float): Concentration in mol/L. Default 1.0 mol/L.
            pressure (float): Pressure in atm. Default 1.0 atm.
            cutoff_entropy_grimme (float): Cutoff frequency for entropy
                in cm^-1 using Grimme's quasi-RRHO method. Default 100.0.
            cutoff_enthalpy (float): Cutoff frequency for enthalpy
                in cm^-1 using Head-Gordon's method. Default 100.0.

        Returns:
            dict: Dictionary containing pKa calculation results.
                All energies in Hartree (au) unless noted otherwise:
                - 'pKa': Computed pKa value of target acid HA
                - 'pKa_reference': Experimental pKa of reference acid HB
                - 'delta_G_soln_kcal_mol': ΔG_soln in kcal/mol (for pKa calc)
                - 'delta_G_soln_au': ΔG_soln in Hartree (au)
                - 'temperature': Temperature in Kelvin
                - 'G_soln_HA_au': Solution free energy of HA (Hartree)
                - 'G_soln_A_au': Solution free energy of A⁻ (Hartree)
                - 'G_soln_HB_au': Solution free energy of HB (Hartree)
                - 'G_soln_B_au': Solution free energy of B⁻ (Hartree)
                - 'E_solv_HA_au': Solvent SP energy of HA (Hartree)
                - 'E_solv_A_au': Solvent SP energy of A⁻ (Hartree)
                - 'E_solv_HB_au': Solvent SP energy of HB (Hartree)
                - 'E_solv_B_au': Solvent SP energy of B⁻ (Hartree)
                - 'G_corr_HA_au': Thermal correction for HA (Hartree)
                - 'G_corr_A_au': Thermal correction for A⁻ (Hartree)
                - 'G_corr_HB_au': Thermal correction for HB (Hartree)
                - 'G_corr_B_au': Thermal correction for B⁻ (Hartree)
                - 'E_gas_HA_au': Gas-phase electronic energy of HA (Hartree)
                - 'E_gas_A_au': Gas-phase electronic energy of A⁻ (Hartree)
                - 'E_gas_HB_au': Gas-phase electronic energy of HB (Hartree)
                - 'E_gas_B_au': Gas-phase electronic energy of B⁻ (Hartree)

        Example:
            # Calculate pKa of 5PQ_Me_ts1 using collidine as reference (pKa=6.75)
            result = Gaussian16pKaOutput.compute_pka(
                ha_gas_file="5PQ_Me_ts1_no_pd_opt.log",
                a_gas_file="5PQ_Me_ts1_b_no_pd_opt.log",
                href_gas_file="collidine-H_opt.log",
                ref_gas_file="collidine_opt.log",
                ha_solv_file="5PQ_Me_ts1_no_pd_opt_sp_smd.log",
                a_solv_file="5PQ_Me_ts1_b_no_pd_opt_sp_smd.log",
                href_solv_file="collidine-H_opt_sp_smd.log",
                ref_solv_file="collidine_opt_sp_smd.log",
                pka_reference=6.75,
                temperature=298.15
            )
            print(f"pKa(HA) = {result['pKa']:.2f}")
        """
        from chemsmart.utils.constants import HARTREE_TO_KCAL_MOL

        # Helper function to get thermochemistry from gas-phase file
        # All energies returned in Hartree (au)
        def get_gas_phase_data(gas_file):
            """Extract E_gas and qh-G(T) from gas-phase frequency calculation.

            Returns:
                tuple: (E_gas_au, qh_G_au, G_corr_au) all in Hartree
            """
            output = cls(
                filename=gas_file,
                temperature=temperature,
                concentration=concentration,
                pressure=pressure,
                cutoff_entropy_grimme=cutoff_entropy_grimme,
                cutoff_enthalpy=cutoff_enthalpy,
            )
            E_gas_au = output.electronic_energy_in_units  # Hartree
            qh_G_au = output.qh_gibbs_free_energy  # Hartree
            G_corr_au = qh_G_au - E_gas_au  # Thermal correction in Hartree
            return E_gas_au, qh_G_au, G_corr_au

        # Helper function to get E_solv from solvent SP file
        # Returns energy in Hartree (au)
        def get_solvent_energy(solv_file):
            """Extract electronic energy from solvent single-point calculation.

            Returns:
                float: E_solv in Hartree (au)
            """
            output = Gaussian16Output(filename=solv_file)
            # Get the last SCF energy (solvent SP) in Hartree
            E_solv_au = output.energies[-1] if output.energies else None
            if E_solv_au is None:
                raise ValueError(
                    f"Could not extract SCF energy from solvent file: {solv_file}"
                )
            return E_solv_au

        # Step 1: Get gas-phase data (E_gas, qh-G, G_corr) for all species
        # All values in Hartree (au)
        E_gas_HA_au, qh_G_HA_au, G_corr_HA_au = get_gas_phase_data(ha_gas_file)
        E_gas_A_au, qh_G_A_au, G_corr_A_au = get_gas_phase_data(a_gas_file)
        E_gas_HRef_au, qh_G_HRef_au, G_corr_HRef_au = get_gas_phase_data(
            href_gas_file
        )
        E_gas_Ref_au, qh_G_Ref_au, G_corr_Ref_au = get_gas_phase_data(
            ref_gas_file
        )

        # Step 2: Get solvent SP energies (E_solv) for all species
        # All values in Hartree (au)
        E_solv_HA_au = get_solvent_energy(ha_solv_file)
        E_solv_A_au = get_solvent_energy(a_solv_file)
        E_solv_HRef_au = get_solvent_energy(href_solv_file)
        E_solv_Ref_au = get_solvent_energy(ref_solv_file)

        # Step 3: Calculate G_soln = E_solv + G_corr
        # Solution free energy in Hartree (au)
        G_soln_HA_au = E_solv_HA_au + G_corr_HA_au
        G_soln_A_au = E_solv_A_au + G_corr_A_au
        G_soln_HRef_au = E_solv_HRef_au + G_corr_HRef_au
        G_soln_Ref_au = E_solv_Ref_au + G_corr_Ref_au

        # Step 4: Calculate ΔG_soln in Hartree (au)
        # ΔG_soln = [G(A⁻)_soln + G(HB)_soln] - [G(HA)_soln + G(B⁻)_soln]
        delta_G_soln_au = (G_soln_A_au + G_soln_HRef_au) - (
            G_soln_HA_au + G_soln_Ref_au
        )

        # Step 5: Convert ΔG_soln to kcal/mol for pKa calculation
        delta_G_soln_kcal_mol = delta_G_soln_au * HARTREE_TO_KCAL_MOL

        # Step 6: Calculate pKa
        # pKa = pKa_ref + ΔG_soln / (RT × ln10)
        # R = 1.987204 cal/(mol·K) = 0.001987204 kcal/(mol·K)
        R_kcal = 0.001987204  # kcal/(mol·K)
        ln10 = 2.302585093

        pka = pka_reference + delta_G_soln_kcal_mol / (
            R_kcal * temperature * ln10
        )

        return {
            "pKa": pka,
            "pKa_reference": pka_reference,
            # ΔG_soln in both units (kcal/mol needed for pKa formula)
            "delta_G_soln_kcal_mol": delta_G_soln_kcal_mol,
            "delta_G_soln_au": delta_G_soln_au,
            "temperature": temperature,
            # Solution free energies in Hartree (au)
            "G_soln_HA_au": G_soln_HA_au,
            "G_soln_A_au": G_soln_A_au,
            "G_soln_HRef_au": G_soln_HRef_au,
            "G_soln_Ref_au": G_soln_Ref_au,
            # Solvent SP energies in Hartree (au)
            "E_solv_HA_au": E_solv_HA_au,
            "E_solv_A_au": E_solv_A_au,
            "E_solv_HRef_au": E_solv_HRef_au,
            "E_solv_Ref_au": E_solv_Ref_au,
            # Thermal corrections in Hartree (au)
            "G_corr_HA_au": G_corr_HA_au,
            "G_corr_A_au": G_corr_A_au,
            "G_corr_HRef_au": G_corr_HRef_au,
            "G_corr_Ref_au": G_corr_Ref_au,
            # Gas-phase electronic energies in Hartree (au)
            "E_gas_HA_au": E_gas_HA_au,
            "E_gas_A_au": E_gas_A_au,
            "E_gas_HRef_au": E_gas_HRef_au,
            "E_gas_Ref_au": E_gas_Ref_au,
        }

    @classmethod
    def print_pka_summary(
        cls,
        ha_gas_file,
        a_gas_file,
        href_gas_file,
        ref_gas_file,
        ha_solv_file,
        a_solv_file,
        href_solv_file,
        ref_solv_file,
        pka_reference,
        temperature: float = 298.15,
        concentration: float = 1.0,
        pressure: float = 1.0,
        cutoff_entropy_grimme: float = 100.0,
        cutoff_enthalpy: float = 100.0,
    ):
        """
        Print a formatted summary of Dual-level Proton Exchange pKa calculation.

        All energies are displayed in Hartree (au) except ΔG_soln which is
        shown in kcal/mol as required for the pKa formula.

        Args:
            ha_gas_file (str): Path to HA gas-phase optimization+freq output file.
            a_gas_file (str): Path to A⁻ gas-phase optimization+freq output file.
            href_gas_file (str): Path to HB gas-phase optimization+freq output file.
            ref_gas_file (str): Path to B⁻ gas-phase optimization+freq output file.
            ha_solv_file (str): Path to HA solvent single-point output file.
            a_solv_file (str): Path to A⁻ solvent single-point output file.
            href_solv_file (str): Path to HB solvent single-point output file.
            ref_solv_file (str): Path to B⁻ solvent single-point output file.
            pka_reference (float): Experimental pKa of the reference acid HB.
            temperature (float): Temperature in Kelvin. Default 298.15 K.
            concentration (float): Concentration in mol/L. Default 1.0 mol/L.
            pressure (float): Pressure in atm. Default 1.0 atm.
            cutoff_entropy_grimme (float): Cutoff for entropy (cm⁻¹). Default 100.0.
            cutoff_enthalpy (float): Cutoff for enthalpy (cm⁻¹). Default 100.0.

        Example:
            Gaussian16pKaOutput.print_pka_summary(
                ha_gas_file="5PQ_Me_ts1_no_pd_opt.log",
                a_gas_file="5PQ_Me_ts1_b_no_pd_opt.log",
                href_gas_file="collidine-H_opt.log",
                ref_gas_file="collidine_opt.log",
                ha_solv_file="5PQ_Me_ts1_no_pd_opt_sp_smd.log",
                a_solv_file="5PQ_Me_ts1_b_no_pd_opt_sp_smd.log",
                href_solv_file="collidine-H_opt_sp_smd.log",
                ref_solv_file="collidine_opt_sp_smd.log",
                pka_reference=6.75,
                temperature=298.15
            )
        """
        result = cls.compute_pka(
            ha_gas_file=ha_gas_file,
            a_gas_file=a_gas_file,
            href_gas_file=href_gas_file,
            ref_gas_file=ref_gas_file,
            ha_solv_file=ha_solv_file,
            a_solv_file=a_solv_file,
            href_solv_file=href_solv_file,
            ref_solv_file=ref_solv_file,
            pka_reference=pka_reference,
            temperature=temperature,
            concentration=concentration,
            pressure=pressure,
            cutoff_entropy_grimme=cutoff_entropy_grimme,
            cutoff_enthalpy=cutoff_enthalpy,
        )

        print("=" * 78)
        print("pKa Calculation - Dual-level Proton Exchange Scheme")
        print("=" * 78)
        print("Reaction: HA + Ref⁻ → A⁻ + HRef")
        print(f"Temperature: {temperature} K")
        print()
        print("Method:")
        print("  G_corr = qh-G(T) - E_gas  (from gas-phase freq calculation)")
        print("  G_soln = E_solv + G_corr  (solution free energy)")
        print(
            "  ΔG_soln = [G(A⁻)_soln + G(HRef)_soln] - [G(HA)_soln + G(Ref⁻)_soln]"
        )
        print("  pKa = pKa_ref + ΔG_soln / (RT × ln10)")
        print("-" * 78)
        print()
        print("Gas-Phase Electronic Energies (E_gas, au):")
        print(f"  HA:  {result['E_gas_HA_au']:.10f}")
        print(f"  A⁻:  {result['E_gas_A_au']:.10f}")
        print(f"  HRef:  {result['E_gas_HRef_au']:.10f}")
        print(f"  Ref⁻:  {result['E_gas_Ref_au']:.10f}")
        print()
        print("Thermal Corrections (G_corr = qh-G - E_gas, au):")
        print(f"  HA:  {result['G_corr_HA_au']:.10f}")
        print(f"  A⁻:  {result['G_corr_A_au']:.10f}")
        print(f"  HRef:  {result['G_corr_HRef_au']:.10f}")
        print(f"  Ref⁻:  {result['G_corr_Ref_au']:.10f}")
        print()
        print("Solvent Single-Point Energies (E_solv, au):")
        print(f"  HA:  {result['E_solv_HA_au']:.10f}")
        print(f"  A⁻:  {result['E_solv_A_au']:.10f}")
        print(f"  HRef:  {result['E_solv_HRef_au']:.10f}")
        print(f"  Ref⁻:  {result['E_solv_Ref_au']:.10f}")
        print()
        print("Solution Free Energies (G_soln = E_solv + G_corr, au):")
        print(f"  HA:  {result['G_soln_HA_au']:.10f}")
        print(f"  A⁻:  {result['G_soln_A_au']:.10f}")
        print(f"  HRef:  {result['G_soln_HRef_au']:.10f}")
        print(f"  Ref⁻:  {result['G_soln_Ref_au']:.10f}")
        print("-" * 78)
        print()
        print("pKa Calculation:")
        print(f"  ΔG_soln = {result['delta_G_soln_au']:.10f} au")
        print(f"         = {result['delta_G_soln_kcal_mol']:.4f} kcal/mol")
        print(f"  pKa(HRef)_ref = {pka_reference:.2f}")
        print()
        print(f"  *** Computed pKa(HA) = {result['pKa']:.2f} ***")
        print("=" * 78)

    # Alias for backward compatibility
    print_pka_dual_level_summary = print_pka_summary
    compute_pka_dual_level = compute_pka
