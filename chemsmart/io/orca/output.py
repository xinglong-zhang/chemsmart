import logging
import math
import os
import re
from functools import cached_property

import numpy as np
from ase import units

from chemsmart.io.molecules.structure import CoordinateBlock, Molecule
from chemsmart.utils.io import (
    clean_duplicate_structure,
    create_molecule_list,
    increment_numbers,
)
from chemsmart.utils.mixins import ORCAFileMixin
from chemsmart.utils.periodictable import PeriodicTable
from chemsmart.utils.repattern import (
    orca_constrained_coordinates_pattern,
    orca_frozen_atoms_output_pattern,
    orca_input_coordinate_in_output,
    orca_nproc_used_line_pattern,
    standard_coord_pattern,
)
from chemsmart.utils.utils import (
    get_range_from_list,
    is_float,
    string2index_1based,
)

p = PeriodicTable()

logger = logging.getLogger(__name__)


class ORCAOutput(ORCAFileMixin):
    """ORCA output file with .out extension.
    ORCA does NOT support any PBC calculations as of version 6.0.1."""

    def __init__(self, filename):
        self.filename = filename

    @property
    def normal_termination(self):
        """Check if ORCA job has completed successfully.

        Checks each of the output file line from the last line onwards.
        """

        def _line_contains_success_indicators(line):
            success_indicators = [
                "****ORCA TERMINATED NORMALLY****",
                "ORCA TERMINATED NORMALLY",
            ]
            return any(indicator in line for indicator in success_indicators)

        return any(
            _line_contains_success_indicators(line)
            for line in self.contents[::-1]
        )

    @cached_property
    def has_forces(self):
        """Check if the output file contains forces."""
        for line in self.contents:
            if "CARTESIAN GRADIENT" in line:
                return True
        return False

    @cached_property
    def forces(self):
        return self._get_forces_for_molecules()

    def _get_forces_for_molecules(self):
        """Obtain a list of cartesian forces.
        Each force is stored as a np array of shape (num_atoms, 3).
        Intrinsic units as used in ORCA: Hartrees/Bohr."""
        list_of_all_forces = []
        for i, line in enumerate(self.contents):
            if "CARTESIAN GRADIENT" in line:
                forces = []
                for line_j in self.contents[i + 3 :]:
                    if len(line_j) == 0:
                        break
                    line_elements = line_j.split()
                    if len(line_elements) == 6:
                        forces.append(
                            [
                                float(line_elements[-3]),
                                float(line_elements[-2]),
                                float(line_elements[-1]),
                            ]
                        )
                list_of_all_forces.append(np.array(forces))
        if len(list_of_all_forces) == 0:
            return None
        return list_of_all_forces

    @cached_property
    def energies(self):
        """Return energies of the system from ORCA output file."""
        return self._get_energies()

    def _get_energies(self):
        """Obtain a list of energies for each geometry optimization point."""
        energies = []
        for i, line in enumerate(self.contents):
            if "FINAL SINGLE POINT ENERGY" in line:
                energy = float(line.split()[-1])
                energies.append(energy)
        return energies

    @cached_property
    def input_coordinates_block(self):
        """Obtain the coordinate block from the input that is printed in the outputfile."""
        return self._get_first_structure_coordinates_block_in_output()

    def _get_input_structure_coordinates_block_in_output(self):
        """In ORCA output file, the input structure is rewritten and for single points,
        is same as the output structure.

        An example of the relevant part of the output describing the structure is:
        | 20> * xyz 0 1
        | 21>   O   -0.00000000323406      0.00000000000000      0.08734060152197
        | 22>   H   -0.75520523910536      0.00000000000000     -0.50967029975151
        | 23>   H   0.75520524233942      0.00000000000000     -0.50967030177046
        | 24> *.
        # this will not work if the input file is supplied separately as .xyz file
        """
        coordinates_block_lines_list = []
        pattern = re.compile(orca_input_coordinate_in_output)
        for i, line in enumerate(self.contents):
            if "INPUT FILE" in line:
                for line in self.contents[i:]:
                    match = pattern.match(line)
                    if match:
                        coord_line = "  ".join(
                            line.split()[-4:]
                        )  # get the last 4 elements
                        coordinates_block_lines_list.append(coord_line)
                    if len(line) == 0:
                        break
        cb = CoordinateBlock(coordinate_block=coordinates_block_lines_list)
        return cb

    @cached_property
    def get_frozen_atoms(self):
        return self._get_constraints()["frozen_atoms"]

    @cached_property
    def get_constrained_bond_angles(self):
        return self._get_constraints()["constrained_bond_angles"]

    @cached_property
    def get_constrained_bond_lengths(self):
        return self._get_constraints()["constrained_bond_lengths"]

    @cached_property
    def get_constrained_dihedrals(self):
        return self._get_constraints()["constrained_dihedrals"]

    def _get_constraints(self):
        """Extract frozen atoms and constrained internal coordinates from ORCA output."""
        constraints = {
            "frozen_atoms": [],
            "constrained_bond_lengths": [],
            "constrained_bond_angles": [],
            "constrained_dihedrals": [],
        }

        # Regex patterns
        constrained_bond_length_pattern = (
            r"\|\s*(\d+)>.*\{\s*B\s+(\d+)\s+(\d+)\s+C\s*\}"
        )
        constrained_bond_angles_pattern = (
            r"\|\s*(\d+)>.*\{\s*A\s+(\d+)\s+(\d+)\s+(\d+)\s+C\s*\}"
        )
        constrained_dihedrals_pattern = (
            r"\|\s*(\d+)>.*\{\s*D\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+C\s*\}"
        )
        orca_frozen_atoms_output_pattern = (
            r"Will constrain atom \d+ coordinate \d"
        )

        for line in self.contents:
            frozen_match = re.match(orca_frozen_atoms_output_pattern, line)
            if frozen_match:
                atom_index = int(line.split()[3])
                if atom_index not in constraints["frozen_atoms"]:
                    constraints["frozen_atoms"].append(atom_index)
            elif re.match(constrained_bond_length_pattern, line):
                constraints["constrained_bond_lengths"].append(
                    f"{line.split()[4]}-{line.split()[5]}"
                )
            elif re.match(constrained_bond_angles_pattern, line):
                constraints["constrained_bond_angles"].append(
                    f"{line.split()[4]}-{line.split()[5]}-{line.split()[6]}"
                )
            elif re.match(constrained_dihedrals_pattern, line):
                constraints["constrained_dihedrals"].append(
                    f"{line.split()[4]}-{line.split()[5]}-{line.split()[6]}-{line.split()[7]}"
                )
        return constraints

    def _get_first_structure_coordinates_block_in_output(self):
        """Obtain the first structure coordinates block in the output file."""
        coordinates_block_lines_list = []
        pattern = re.compile(standard_coord_pattern)
        found_header = False

        for line in self.contents:
            # Start collecting after finding the header
            if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
                found_header = True
                continue  # Skip the header line itself

            # If we've found the header, process lines
            if found_header:
                match = pattern.match(line)
                if match:
                    # Extract the last 4 elements (symbol, x, y, z) and join with double spaces
                    coord_line = "  ".join(line.split()[-4:])
                    coordinates_block_lines_list.append(coord_line)
                elif (
                    coordinates_block_lines_list
                ):  # Stop if we hit a non-matching line after collecting coords
                    break

        # Return CoordinateBlock instance
        cb = CoordinateBlock(coordinate_block=coordinates_block_lines_list)
        return cb

    @property
    def frozen_atoms(self):
        """Get frozen atoms from the ORCA output file."""
        frozen_atoms = []
        for line in self.contents:
            if re.match(orca_frozen_atoms_output_pattern, line):
                atom_index = int(line.split()[3])
                atom_index += 1  # Convert to 1-based index
                if atom_index not in frozen_atoms:
                    frozen_atoms.append(atom_index)
        return frozen_atoms

    @property
    def constrained_bond_lengths(self):
        constrained_bond_lengths, _, _ = self._get_constraints
        return constrained_bond_lengths

    @property
    def constrained_bond_angles(self):
        _, constrained_bond_angles, _ = self._get_constraints
        return constrained_bond_angles

    @property
    def constrained_dihedral_angles(self):
        _, _, constrained_dihedral_angles = self._get_constraints
        return constrained_dihedral_angles

    @cached_property
    def _get_constraints(self):
        """Extract constrained internal coordinates from ORCA output.
        Reads from Redundant Internal Coordinates block where, if DOF is constrained,
        there is a "C" at the end of the line.
        Returns:
             a dict similar to optimized parameters, specifying the constraints and
             the associated values, e.g.,
                optimized_geometry == {
                    "B(H1,O0)": 0.9627,
                    "B(H2,O0)": 0.9627,
                    "A(H1,O0,H2)": 103.35,
                }
        """
        constrained_bond_lengths = {}
        constrained_bond_angles = {}
        constrained_dihedral_angles = {}

        for i, line in enumerate(self.contents):
            if "Redundant Internal Coordinates" in line:
                for j, line_j in enumerate(self.contents[i + 5 :]):
                    if "------------------------------------" in line_j:
                        continue
                    if len(line_j) == 0:
                        break
                    if re.match(orca_constrained_coordinates_pattern, line_j):
                        line_elements = line_j.split()
                        if line_elements[1].lower().startswith("b"):  # bond
                            parameter = f"{line_elements[1]}{line_elements[2]}{line_elements[3]}"
                            parameter = increment_numbers(parameter, 1)
                            constrained_bond_lengths[parameter] = float(
                                line_elements[-3]
                            )
                        elif line_elements[1].lower().startswith("a"):  # angle
                            parameter = f"{line_elements[1]}{line_elements[2]}{line_elements[3]}{line_elements[4]}"
                            parameter = increment_numbers(parameter, 1)
                            constrained_bond_angles[parameter] = float(
                                line_elements[-3]
                            )
                        elif (
                            line_elements[1].lower().startswith("d")
                        ):  # dihedral
                            parameter = f"{line_elements[1]}{line_elements[2]}{line_elements[3]}{line_elements[4]}{line_elements[5]}"
                            parameter = increment_numbers(parameter, 1)
                            constrained_dihedral_angles[parameter] = float(
                                line_elements[-3]
                            )
                        else:
                            raise ValueError(
                                f"Unknown parameter type in line: {line_j}"
                            )
                return (
                    constrained_bond_lengths,
                    constrained_bond_angles,
                    constrained_dihedral_angles,
                )

    @cached_property
    def optimized_output_lines(self):
        """Chunk of outputfile where the properties are calculated based on the final optimized structure.

        FOR SP CALCULATION, THIS WILL BE EMPTY!
        """
        optimized_output_lines = []
        for i, line_i in enumerate(self.contents):
            if "THE OPTIMIZATION HAS CONVERGED" in line_i:
                optimized_output_lines = list(self.contents[i:])
        return optimized_output_lines

    @property
    def route_string(self):
        """Route string for ORCA file, convert to lower case."""
        for line in self.contents:
            if line.startswith("|  1> !"):
                return line.lower().split("1> ")[-1]
        return None

    @property
    def natoms(self):
        for line in self.contents:
            if "Number of atoms" in line:
                return int(line.split()[-1])
        return None

    @property
    def num_basis_functions(self):
        for line in self.contents:
            if "Number of basis functions" in line:
                return int(line.split()[-1])
        return None

    @property
    def num_shells(self):
        for line in self.contents:
            if "Number of shells" in line:
                return int(line.split()[-1])
        return None

    @property
    def max_ang_mom(self):
        """Max. angular momentum."""
        for line in self.contents:
            if "Maximum angular momentum" in line:
                return int(line.split()[-1])
        return None

    @property
    def contraction_scheme(self):
        for line in self.contents:
            if "Contraction scheme used" in line:
                return line.split("...")[-1].strip()
        return None

    @property
    def coulomb_range_seperation(self):
        for line in self.contents:
            if "Coulomb Range Separation" in line:
                return line.split("...")[-1].strip()
        return None

    @property
    def exchange_range_seperation(self):
        for line in self.contents:
            if "Exchange Range Separation" in line:
                return line.split("...")[-1].strip()
        return None

    @property
    def finite_nucleus_model(self):
        for line in self.contents:
            if "Finite Nucleus Model" in line:
                return line.split("...")[-1].strip()
        return None

    @property
    def aux_j_fitting_basis(self):
        """Auxiliary Coulomb fitting basis."""
        for line in self.contents:
            if "Auxiliary Coulomb fitting basis" in line:
                return line.split("...")[-1].strip()
        return None

    @property
    def aux_j_num_basis_functions(self):
        """# of basis functions in Aux-J."""
        for line in self.contents:
            if "# of basis functions in Aux-J" in line:
                return int(line.split("...")[-1].strip())
        return None

    @property
    def aux_j_num_shells(self):
        """# of shells in Aux-J."""
        for line in self.contents:
            if "# of shells in Aux-J" in line:
                return int(line.split("...")[-1].strip())
        return None

    @property
    def aux_j_max_ang_mom(self):
        """# of angular momentum in Aux-J."""
        for line in self.contents:
            if "Maximum angular momentum in Aux-J" in line:
                return int(line.split("...")[-1].strip())
        return None

    @property
    def aux_jk_fitting_basis(self):
        """Auxiliary J/K fitting basis."""
        for line in self.contents:
            if "Auxiliary J/K fitting basis" in line:
                return line.split("...")[-1].strip()
        return None

    @property
    def aux_k_fitting_basis(self):
        """Auxiliary Correlation fitting basis."""
        for line in self.contents:
            if "Auxiliary Correlation fitting basis" in line:
                return line.split("...")[-1].strip()
        return None

    @property
    def aux_external_fitting_basis(self):
        """Auxiliary 'external' fitting basis."""
        for line in self.contents:
            if "Auxiliary 'external' fitting basis" in line:
                return line.split("...")[-1].strip()
        return None

    @property
    def integral_threshold(self):
        """Integral threshold."""
        for line in self.contents:
            if "Integral threshold" in line:
                return float(line.split()[-1])
        return None

    @property
    def primitive_cutoff(self):
        """Primitive cut-off."""
        for line in self.contents:
            if "Primitive cut-off" in line:
                return float(line.split()[-1])
        return None

    @property
    def primitive_pair_threshold(self):
        """Primitive pair pre-selection threshold."""
        for line in self.contents:
            if "Primitive pair pre-selection threshold" in line:
                return float(line.split()[-1])
        return None

    @property
    def ri_approx(self):
        """RI-approximation to the Coulomb term."""
        for line in self.contents:
            if "RI-approximation to the Coulomb term is turned" in line:
                return line.split()[-1] == "on"
        return None

    @property
    def rij_cosx(self):
        """RIJ-COSX(HFX calculated with COS-X))."""
        for line in self.contents:
            if "RIJ-COSX (HFX calculated with COS-X)" in line:
                return line.split()[-1] == "on"
        return None

    @property
    def charge(self):
        pattern = re.compile(r"Total Charge\s+Charge")
        for line in self.contents:
            if pattern.search(line):
                return int(line.split()[-1])
        return None

    @property
    def multiplicity(self):
        pattern = re.compile(r"Multiplicity\s+Mult")
        for line in self.contents:
            if pattern.search(line):
                return int(line.split()[-1])
        return None

    @property
    def num_electrons(self):
        pattern = re.compile(r"Number of Electrons\s+NEL")
        for line in self.contents:
            if pattern.search(line):
                return int(line.split()[-1])
        return None

    @property
    def basis_dim(self):
        pattern = re.compile(r"Basis Dimension\s+Dim")
        for line in self.contents:
            if pattern.search(line):
                return int(line.split()[-1])
        return None

    @property
    def diis_acceleration(self):
        """Convergence Acceleration using DIIS."""
        pattern = re.compile(r"DIIS\s+CNVDIIS")
        for line in self.contents:
            if pattern.search(line):
                return line.split()[-1] == "on"
        return None

    @property
    def scf_maxiter(self):
        pattern = re.compile(r"Maximum # iterations\s+MaxIter")
        for line in self.contents:
            if pattern.search(line):
                return int(line.split()[-1])
        return None

    @property
    def scf_convergence(self):
        pattern = re.compile(r"\|.*>.*convergence", re.IGNORECASE)
        for line in self.contents:
            if pattern.search(line):
                # return string directly after 'convergence' -- convergence criteria
                # ['sloppy', 'loose', 'medium', 'strong', 'tight', 'verytight', 'extreme']
                return line.lower().split("convergence")[-1].strip().split()[0]
        return None

    @property
    def dipole(self):
        pattern = re.compile(r"\|.*>.*dipole", re.IGNORECASE)
        for line in self.contents:
            if pattern.search(line):
                # return string directly after 'dipole'
                # ['True', 'False']
                return line.lower().split("dipole")[-1].strip().split()[0]
        return None

    @property
    def quadrupole(self):
        pattern = re.compile(r"\|.*>.*quadrupole", re.IGNORECASE)
        # same as
        # `if '|' in line and '>' in line and 'quadrupole' in line:`
        for line in self.contents:
            if pattern.search(line):
                # return string directly after 'quadrupole'
                # ['True', 'False']
                return line.lower().split("quadrupole")[-1].strip().split()[0]
        return None

    @property
    def converged(self):
        for line in self.contents:
            if "THE OPTIMIZATION HAS CONVERGED" in line:
                return True
        return None

    @cached_property
    def all_structures(self):
        """Obtain all structures in ORCA output file, including intermediate points if present.
        Include corresponding energies and forces where available."""

        # Extract all raw structure data
        orientations = self._get_all_orientations()
        if not orientations:
            return []  # No structures found

        # Clean duplicate structures (e.g., last structure might repeat in some cases)
        clean_duplicate_structure(orientations)

        # Handle PBC (default to None if not present in ORCA output)
        orientations_pbc = self._get_pbc_conditions() or [None] * len(
            orientations
        )

        # Prepare energies and forces (adjust for length mismatches)
        energies = (
            self.energies_in_eV
            if self.energies_in_eV is not None
            else [None] * len(orientations)
        )
        forces = (
            self.forces_in_eV_per_angstrom
            if self.forces_in_eV_per_angstrom is not None
            else [None] * len(orientations)
        )

        # Handle job termination
        if self.normal_termination:
            num_structures = len(orientations)
        else:
            # For abnormal termination, exclude the last structure (likely incomplete)
            num_structures = max(0, len(orientations) - 1)
            if num_structures == 0:
                return []

        # Truncate energies/forces to match number of structures
        num_structures_to_use = min(num_structures, len(energies), len(forces))
        orientations = orientations[:num_structures_to_use]
        orientations_pbc = orientations_pbc[:num_structures_to_use]
        energies = energies[:num_structures_to_use]
        forces = forces[:num_structures_to_use]

        # Create molecule list
        all_structures = create_molecule_list(
            orientations=orientations,
            orientations_pbc=orientations_pbc,
            energies=energies,
            forces=forces,
            symbols=self.symbols,
            charge=self.charge,
            multiplicity=self.multiplicity,
            frozen_atoms=self.frozen_atoms,
            pbc_conditions=(
                self.list_of_pbc_conditions
                if hasattr(self, "list_of_pbc_conditions")
                else None
            ),
            num_structures=num_structures_to_use,
        )

        # Filter optimized steps if requested (e.g., for geometry optimization)
        if (
            hasattr(self, "optimized_steps_indices")
            and self.optimized_steps_indices
            and not hasattr(self, "include_intermediate")
        ):
            all_structures = [
                all_structures[i] for i in self.optimized_steps_indices
            ]

        logger.info(
            f"Total number of structures located: {len(all_structures)}"
        )
        return all_structures

    def _get_all_orientations(self):
        """Extract all Cartesian coordinate blocks from the ORCA output."""
        orientations = []
        for i, line in enumerate(self.contents):
            if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
                coordinate_lines = []
                for line_j in self.contents[i + 2 :]:  # Skip header lines
                    if (
                        not line_j.strip() or "----" in line_j
                    ):  # End of coordinate block
                        break
                    pattern = re.compile(standard_coord_pattern)
                    if pattern.match(line_j):
                        coordinate_lines.append(
                            [float(val) for val in line_j.split()[1:]]
                        )
                if coordinate_lines:
                    orientations.append(np.array(coordinate_lines))
        return orientations

    def _get_pbc_conditions(self):
        """Extract periodic boundary conditions if present (rare in ORCA outputs)."""
        # ORCA rarely includes PBC in standard outputs; this is a placeholder
        # If your ORCA output includes lattice vectors, implement parsing here
        return None  # Default: no PBC unless explicitly parsed

    def _get_all_structures(self):
        """Extract all Cartesian coordinate blocks from the ORCA output.
        This does not however include energy and forces."""
        structures = []
        for i, line in enumerate(self.contents):
            if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
                coordinate_lines = []
                for line_j in self.contents[i:]:
                    pattern = re.compile(standard_coord_pattern)
                    if len(line_j) == 0:
                        break
                    if pattern.match(line_j):
                        coordinate_lines.append(line_j)
                cb = CoordinateBlock(coordinate_block=coordinate_lines)
                structures.append(cb.molecule)
        return structures

    ################################################
    #######  GET OPTIMIZED PARAMETERS ##############
    ################################################
    def get_optimized_parameters(self):
        """Obtain a list of optimized geometry parameters in ORCA output.

        --- Optimized Parameters --- (Angstrom and degrees)
        Example in the output file:
                          --- Optimized Parameters ---
                           (Angstroem and degrees).

            Definition                    OldVal   dE/dq     Step     FinalVal
        ----------------------------------------------------------------------------
         1. B(H   1,O   0)                0.9627 -0.000014  0.0000    0.9627
         2. B(H   2,O   0)                0.9627 -0.000014  0.0000    0.9627
         3. A(H   1,O   0,H   2)          103.34 -0.000009    0.00    103.35
        ----------------------------------------------------------------------------
        #TODO: need to convert to 1-indexing
        """
        optimized_geometry = {}
        for i, line_i in enumerate(self.optimized_output_lines):
            if "--- Optimized Parameters ---" in line_i:
                for line_j in self.optimized_output_lines[i + 5 :]:
                    parameter = None
                    if "---------------" in line_j:
                        break
                    line_elements = line_j.split()
                    if line_elements[1].lower().startswith("b"):  # bond
                        parameter = f"{line_elements[1]}{line_elements[2]}{line_elements[3]}"
                    if line_elements[1].lower().startswith("a"):  # angle
                        parameter = f"{line_elements[1]}{line_elements[2]}{line_elements[3]}{line_elements[4]}"
                    if line_elements[1].lower().startswith("d"):  # dihedral
                        parameter = (
                            f"{line_elements[1]}{line_elements[2]}{line_elements[3]}{line_elements[4]}"
                            f"{line_elements[5]}"
                        )
                    if parameter is not None:
                        # Convert to 1-indexing
                        parameter = increment_numbers(parameter, 1)
                    optimized_geometry[parameter] = float(line_elements[-1])
        ## the above result will return a dictionary storing the optimized parameters:
        ## optimized_geometry = { b(h1,o0) : 0.9627,  b(h2,o0) : 0.9627,  a(h1,o0, h2) : 103.35 }
        return optimized_geometry

    @cached_property
    def optimized_structure(self):
        if self.normal_termination:
            return self.all_structures[-1]
        else:
            return self._get_optimized_final_structure()

    def _get_optimized_final_structure(self):
        """Obtain the final optimized structure from ORCA geometry optimization job.

        An example of the output for this portion will look like:

                 *** FINAL ENERGY EVALUATION AT THE STATIONARY POINT ***
                 ***               (AFTER    5 CYCLES)               ***
                 *******************************************************
        ---------------------------------
        CARTESIAN COORDINATES (ANGSTROEM)
        ---------------------------------
          O     -0.000000    0.000000    0.087341
          H     -0.755205    0.000000   -0.509670
          H      0.755205    0.000000   -0.509670

        ----------------------------
        CARTESIAN COORDINATES (A.U.)
        ----------------------------
        NO LB      ZA    FRAG     MASS         X           Y           Z
        0 O     8.0000    0    15.999   -0.000000    0.000000    0.165050
        """
        pattern = re.compile(standard_coord_pattern)
        coordinate_lines = []
        if len(self.optimized_output_lines) != 0:
            for i, line_i in enumerate(self.optimized_output_lines):
                if "FINAL ENERGY EVALUATION AT THE STATIONARY POINT" in line_i:
                    # only start getting the structures that appear after this line (stationary point)
                    for line_j in self.optimized_output_lines[i + 4 :]:
                        if len(line_j) == 0:
                            # stop when an empty line appears, else it will continue parsing non-geometry
                            break
                        # start reading 4 lines after
                        if pattern.match(line_j):
                            coordinate_lines.append(line_j)
        cb = CoordinateBlock(coordinate_block=coordinate_lines)
        return cb.molecule

    @property
    def final_structure(self):
        if self.optimized_output_lines is not None:
            return self.optimized_structure
        try:
            return (
                self.last_structure
            )  # for job that does not terminate normally
        except (ValueError, IndexError):
            return (
                self._get_molecule_from_sp_output_file()
            )  # no structure can be created from output thus use input structure

    @cached_property
    def last_structure(self):
        """Return last structure, whether the output file has completed successfully or not."""
        return self.all_structures[-1]

    @property
    def molecule(self):
        return self.final_structure

    def get_molecule(self, index="-1"):
        index = string2index_1based(index)
        return self.all_structures[index]

    def _get_molecule_from_sp_output_file(self):
        molecule = None

        # if sp output file contains line read from .xyz
        for line in self.contents:
            if "coordinates will be read from file:" in line:
                xyz_file = line.strip().split("file: ")[
                    -1
                ]  # the lines here have all been converted to lower case
                xyz_filepath = os.path.join(self.folder, xyz_file)
                assert os.path.exists(
                    xyz_filepath
                ), f".xyz file read from {xyz_filepath} does not exist!"
                if os.path.exists(xyz_filepath):
                    molecule = Molecule.from_filepath(filepath=xyz_filepath)
                    break
            else:
                # If molecule is not found, get it from the input lines in the output file
                molecule = self._get_input_structure_in_output()
        return molecule

    def _get_input_structure_in_output(self):
        """In ORCA output file, the input structure is rewritten and for single points,
        is same as the output structure.

        An example of the relevant part of the output describing the structure is:
        | 20> * xyz 0 1
        | 21>   O   -0.00000000323406      0.00000000000000      0.08734060152197
        | 22>   H   -0.75520523910536      0.00000000000000     -0.50967029975151
        | 23>   H   0.75520524233942      0.00000000000000     -0.50967030177046
        | 24> *.
        """
        final_symbols = []
        final_positions = []
        pattern = re.compile(orca_input_coordinate_in_output)
        for line in self.contents:
            match = pattern.match(line)
            if match:
                line_num, atom_type, x, y, z = match.groups()
                final_symbols.append(str(atom_type))
                each_coord = [float(x), float(y), float(z)]
                final_positions.append(each_coord)

        final_positions = np.array(final_positions)
        if len(final_symbols) == 0:
            raise ValueError("No structure found!")

        return Molecule.from_symbols_and_positions_and_pbc_conditions(
            list_of_symbols=final_symbols, positions=final_positions
        )

    @property
    def empirical_formula(self):
        return self.molecule.get_chemical_formula(empirical=True)

    @property
    def optimized_geometry(self):
        return self.molecule.positions

    def _get_optimized_scf_energy(self):
        """Get the final SCF energy in Hartree."""
        for i, line_i in enumerate(self.optimized_output_lines):
            if "TOTAL SCF ENERGY" in line_i:
                for line_j in self.optimized_output_lines[i:]:
                    if "Total Energy       :" in line_j:
                        line_j_elements = line_j.split()
                        energy_in_hartree = float(line_j_elements[-4])
                        energy_in_eV = float(line_j_elements[-2])
                        assert math.isclose(
                            energy_in_hartree * units.Hartree,
                            energy_in_eV,
                            rel_tol=1e-4,
                        )
                        return energy_in_hartree

    @property
    def _get_sp_scf_energy(self):
        """Get the final SCF energy in Hartree."""
        if self.optimized_output_lines is None:
            for line in self.contents:
                if "Total Energy       :" in line:
                    line_elements = line.split()
                    energy_in_hartree = float(line_elements[-4])
                    energy_in_eV = float(line_elements[-2])
                    assert math.isclose(
                        energy_in_hartree * units.Hartree,
                        energy_in_eV,
                        rel_tol=1e-4,
                    )
                    return energy_in_hartree

    @property
    def final_scf_energy(self):
        if self.optimized_output_lines is not None:
            return self._get_optimized_scf_energy()
        return self._get_sp_scf_energy()

    @property
    def final_energy(self):
        if self.final_scf_energy is not None:
            return self.final_scf_energy
        return self.single_point_energy

    @property
    def single_point_energy(self):
        for line in self.contents:
            if "FINAL SINGLE POINT ENERGY" in line:
                sp_energy_in_hartree = float(line.split()[-1])
                # convert hartree to eV
                return sp_energy_in_hartree

    @property
    def single_point_energy_eV(self):
        """Get the single point energy in eV."""
        return self.single_point_energy * units.Hartree

    @property
    def final_nuclear_repulsion(self):
        """Get the final nuclear repulsion energy in Hartree."""
        final_nuclear_repulsion_hartree, _ = (
            self._get_final_nuclear_repulsion()
        )
        return final_nuclear_repulsion_hartree

    @property
    def final_nuclear_repulsion_eV(self):
        """Get the final nuclear repulsion energy in eV."""
        _, final_nuclear_repulsion_eV = self._get_final_nuclear_repulsion()
        return final_nuclear_repulsion_eV

    def _get_final_nuclear_repulsion(self):
        final_nuclear_repulsion_hartree = []
        final_nuclear_repulsion_eV = []
        for line in self.contents:
            if "Nuclear Repulsion  :" in line:
                line_elements = line.split()
                energy_in_hartree = float(line_elements[-4])
                energy_in_eV = float(line_elements[-2])
                assert np.isclose(
                    energy_in_hartree * units.Hartree,
                    energy_in_eV,
                    rtol=1e-4,
                )
                final_nuclear_repulsion_hartree.append(energy_in_hartree)
                final_nuclear_repulsion_eV.append(energy_in_eV)
        return (
            final_nuclear_repulsion_hartree[-1],
            final_nuclear_repulsion_eV[-1],
        )

    @property
    def final_electronic_energy(self):
        """Get the final electronic energy in Hartree."""
        final_electronic_energy_hartree, _ = (
            self._get_final_electronic_energy()
        )
        return final_electronic_energy_hartree

    @property
    def final_electronic_energy_eV(self):
        """Get the final electronic energy in eV."""
        _, final_electronic_energy_eV = self._get_final_electronic_energy()
        return final_electronic_energy_eV * units.Hartree

    def _get_final_electronic_energy(self):
        final_electronic_energy_hartree = []
        final_electronic_energy_eV = []
        for line in self.contents:
            if "Electronic Energy  :" in line:
                line_elements = line.split()
                energy_in_hartree = float(line_elements[-4])
                energy_in_eV = float(line_elements[-2])
                assert np.isclose(
                    energy_in_hartree * units.Hartree,
                    energy_in_eV,
                    rtol=1e-4,
                )
                final_electronic_energy_hartree.append(energy_in_hartree)
                final_electronic_energy_eV.append(energy_in_eV)
        return (
            final_electronic_energy_hartree[-1],
            final_electronic_energy_eV[-1],
        )

    @property
    def one_electron_energy(self):
        """Get the one-electron energy in Hartree."""
        one_electron_energy_hartree, _ = self._get_one_electron_energy()
        return one_electron_energy_hartree

    @property
    def one_electron_energy_eV(self):
        """Get the one-electron energy in eV."""
        _, one_electron_energy_eV = self._get_one_electron_energy()
        return one_electron_energy_eV * units.Hartree

    def _get_one_electron_energy(self):
        one_electron_energy_hartree = []
        one_electron_energy_eV = []
        for line in self.contents:
            if "One Electron Energy:" in line:
                line_elements = line.split()
                energy_in_hartree = float(line_elements[-4])
                energy_in_eV = float(line_elements[-2])
                assert np.isclose(
                    energy_in_hartree * units.Hartree,
                    energy_in_eV,
                    rtol=1e-4,
                )
                one_electron_energy_hartree.append(energy_in_hartree)
                one_electron_energy_eV.append(energy_in_eV)
        return one_electron_energy_hartree[-1], one_electron_energy_eV[-1]

    @property
    def two_electron_energy(self):
        """Get the two-electron energy in Hartree."""
        two_electron_energy_hartree, _ = self._get_two_electron_energy()
        return two_electron_energy_hartree

    @property
    def two_electron_energy_eV(self):
        """Get the two-electron energy in eV."""
        _, two_electron_energy_eV = self._get_two_electron_energy()
        return two_electron_energy_eV * units.Hartree

    def _get_two_electron_energy(self):
        two_electron_energy_hartree = []
        two_electron_energy_eV = []
        for line in self.contents:
            if "Two Electron Energy:" in line:
                line_elements = line.split()
                energy_in_hartree = float(line_elements[-4])
                energy_in_eV = float(line_elements[-2])
                assert np.isclose(
                    energy_in_hartree * units.Hartree,
                    energy_in_eV,
                    rtol=1e-4,
                )
                two_electron_energy_hartree.append(energy_in_hartree)
                two_electron_energy_eV.append(energy_in_eV)
        return two_electron_energy_hartree[-1], two_electron_energy_eV[-1]

    @property
    def max_cosx_asymmetry_energy(self):
        """Get the max COSX asymmetry energy in Hartree."""
        max_cosx_asymmetry_energy_hartree = (
            self._get_max_cosx_asymmetry_energy()
        )
        if max_cosx_asymmetry_energy_hartree is not None:
            return max_cosx_asymmetry_energy_hartree[-1]

    @property
    def max_cosx_asymmetry_energy_eV(self):
        """Get the max COSX asymmetry energy in eV."""
        max_cosx_asymmetry_energy_eV = self._get_max_cosx_asymmetry_energy_eV()
        if max_cosx_asymmetry_energy_eV is not None:
            return max_cosx_asymmetry_energy_eV[-1]

    def _get_max_cosx_asymmetry_energy(self):
        max_cosx_asymmetry_energy_hartree = []
        for line in self.contents:
            if "Max COSX asymmetry :" in line:
                energy_in_hartree = float(line.split()[-4])
                max_cosx_asymmetry_energy_hartree.append(energy_in_hartree)
        if len(max_cosx_asymmetry_energy_hartree) != 0:
            return max_cosx_asymmetry_energy_hartree

    def _get_max_cosx_asymmetry_energy_eV(self):
        max_cosx_asymmetry_energy_hartree = (
            self._get_max_cosx_asymmetry_energy()
        )
        if len(max_cosx_asymmetry_energy_hartree) != 0:
            max_cosx_asymmetry_energy_eV = [
                value * units.Hartree
                for value in max_cosx_asymmetry_energy_hartree
            ]
            return max_cosx_asymmetry_energy_eV

    @property
    def potential_energy(self):
        """Get the potential energy in Hartree."""
        potential_energy_hartree = self._get_potential_energy_hartree()
        if potential_energy_hartree is not None:
            return potential_energy_hartree[-1]

    @property
    def potential_energy_eV(self):
        """Get the potential energy in eV."""
        potential_energy_eV = self._get_potential_energy_eV()
        if potential_energy_eV is not None:
            return potential_energy_eV[-1]

    def _get_potential_energy_hartree(self):
        potential_energy_hartree = []
        for line in self.contents:
            if "Potential Energy   :" in line:
                energy_in_hartree = float(line.split()[-4])
                potential_energy_hartree.append(energy_in_hartree)
        if len(potential_energy_hartree) != 0:
            return potential_energy_hartree

    def _get_potential_energy_eV(self):
        potential_energy_hartree = self._get_potential_energy_hartree()
        if len(potential_energy_hartree) != 0:
            potential_energy_eV = [
                value * units.Hartree for value in potential_energy_hartree
            ]
            return potential_energy_eV

    @property
    def kinetic_energy(self):
        """Get the kinetic energy in Hartree."""
        kinetic_energy_in_hartree = self._get_kinetic_energy_hartree()
        if kinetic_energy_in_hartree is not None:
            return kinetic_energy_in_hartree[-1]

    @property
    def kinetic_energy_eV(self):
        """Get the kinetic energy in eV."""
        kinetic_energy_eV = self._get_kinetic_energy_eV()
        if kinetic_energy_eV is not None:
            return kinetic_energy_eV[-1]

    def _get_kinetic_energy_hartree(self):
        kinetic_energy_hartree = []
        for line in self.contents:
            if "Kinetic Energy     :" in line:
                energy_in_hartree = float(line.split()[-4])
                kinetic_energy_hartree.append(energy_in_hartree)
        if len(kinetic_energy_hartree) != 0:
            return kinetic_energy_hartree

    def _get_kinetic_energy_eV(self):
        kinetic_energy_eV = [
            value * units.Hartree
            for value in self._get_kinetic_energy_hartree()
        ]
        if len(kinetic_energy_eV) != 0:
            return kinetic_energy_eV

    @property
    def virial_ratio(self):
        virial_ratios = []
        for line in self.contents:
            if "Virial Ratio       :" in line:
                virial_ratios.append(float(line.split()[-1]))
        if len(virial_ratios) != 0:
            return virial_ratios[-1]

    @property
    def xc_energy(self):
        """Get the XC energy in Hartree."""
        xc_energy_hartree = self._get_xc_energy_hartree()
        if xc_energy_hartree is not None:
            return xc_energy_hartree[-1]

    @property
    def xc_energy_eV(self):
        """Get the XC energy in eV."""
        xc_energy_eV = self._get_xc_energy_eV()
        if xc_energy_eV is not None:
            return xc_energy_eV[-1]

    def _get_xc_energy_hartree(self):
        xc_energy_hartree = []
        for line in self.contents:
            if "E(XC)              :" in line:
                xc_energy_hartree.append(float(line.split()[-2]))

        return xc_energy_hartree

    def _get_xc_energy_eV(self):
        xc_energy_hartree = self._get_xc_energy_hartree()
        if xc_energy_hartree is not None:
            xc_energy_eV = [
                value * units.Hartree for value in xc_energy_hartree
            ]
            return xc_energy_eV

    @property
    def dfet_embed_energy(self):
        """Get the DFET-embed energy in Hartree."""
        dfet_embed_energy_hartree = self._get_dfet_embed_energy()
        if dfet_embed_energy_hartree is not None:
            print(dfet_embed_energy_hartree)
            return dfet_embed_energy_hartree[-1]

    @property
    def dfet_embed_energy_eV(self):
        """Get the DFET-embed energy in eV."""
        dfet_embed_energy_eV = self._get_dfet_embed_energy_eV()
        if dfet_embed_energy_eV is not None:
            return dfet_embed_energy_eV[-1]

    def _get_dfet_embed_energy(self):
        """Get the DFET-embed energy in Hartree."""
        dfet_embed_energy_hartree = []
        for line in self.contents:
            if "DFET-embed. en.    :" in line:
                dfet_embed_energy_hartree.append(float(line.split()[-2]))
        if len(dfet_embed_energy_hartree) != 0:
            return dfet_embed_energy_hartree

    def _get_dfet_embed_energy_eV(self):
        dfet_embed_energy_hartree = self._get_dfet_embed_energy()
        if len(dfet_embed_energy_hartree) != 0:
            dfet_embed_energy_eV = [
                value * units.Hartree for value in dfet_embed_energy_hartree
            ]
            return dfet_embed_energy_eV

    @property
    def orbital_occupancy(self):
        _, orbital_occupancy = self._get_orbital_energies_and_occupancy()
        return orbital_occupancy

    @property
    def orbital_energies(self):
        orbital_energies, _ = self._get_orbital_energies_and_occupancy()
        return orbital_energies

    def _get_orbital_energies_and_occupancy(self):
        """Get the orbital energies and occupancy from the ORCA output file.
        Orbital energies are in eV."""
        orbital_occupancy = []
        orbital_energies = []
        for line in self._get_last_orbital_energies_section()[2:]:
            # ignore the lines '----------------' and one empty line that follows
            line_elements = line.split()
            if len(line_elements) == 0:
                break
            if len(line_elements) != 4 or line.startswith("NO"):
                continue
            else:
                occ = int(float(line_elements[1]))
                orbital_occupancy.append(occ)
                energy_in_hartree = float(line_elements[2])
                orbital_energies.append(energy_in_hartree * units.Hartree)
        return orbital_energies, orbital_occupancy

    def _get_last_orbital_energies_section(self):
        """Get the last section of orbital energies"""
        reversed_lines = []
        for line in reversed(self.contents):
            if "ORBITAL ENERGIES" not in line:
                reversed_lines.append(line)
            else:
                break
        return reversed_lines[::-1]

    @property
    def homo_energy(self):
        # get all filled orbitals
        orbitals = list(zip(self.orbital_energies, self.orbital_occupancy))
        occupied_energies = [
            energy for energy, occupancy in orbitals if occupancy == 2
        ]
        # return the highest occupied MO energy
        return max(occupied_energies)

    @property
    def lumo_energy(self):
        # get all empty orbitals
        orbitals = list(zip(self.orbital_energies, self.orbital_occupancy))
        unoccupied_energies = [
            energy for energy, occupancy in orbitals if occupancy == 0
        ]
        # return the lowest unoccupied MO energy
        return min(unoccupied_energies)

    @cached_property
    def fmo_gap(self):
        if self.multiplicity == 1:
            return self.lumo_energy - self.homo_energy
        else:
            # to implement for radical systems
            pass

    @property
    def mulliken_atomic_charges(self):
        all_mulliken_atomic_charges = []
        for i, line_i in enumerate(self.contents):
            mulliken_atomic_charges = {}
            if "MULLIKEN ATOMIC CHARGES" in line_i:
                for line_j in self.contents[i + 2 :]:
                    if "Sum of atomic charges" in line_j:
                        break
                    line_j_elements = line_j.split()
                    element = p.to_element(line_j_elements[1])
                    element_num = f"{element}{line_j_elements[0]}"
                    mulliken_atomic_charges[element_num] = float(
                        line_j_elements[-1]
                    )
                all_mulliken_atomic_charges.append(mulliken_atomic_charges)
        return all_mulliken_atomic_charges[-1]

    @property
    def loewdin_atomic_charges(self):
        all_loewdin_atomic_charges = []
        for i, line_i in enumerate(self.contents):
            loewdin_atomic_charges = {}
            if "LOEWDIN ATOMIC CHARGES" in line_i:
                for line_j in self.contents[i + 2 :]:
                    if "LOEWDIN REDUCED ORBITAL CHARGES" in line_j:
                        break
                    line_j_elements = line_j.split()
                    if len(line_j_elements) == 4:
                        element = p.to_element(line_j_elements[1])
                        element_num = f"{element}{line_j_elements[0]}"
                        loewdin_atomic_charges[element_num] = float(
                            line_j_elements[-1]
                        )
                all_loewdin_atomic_charges.append(loewdin_atomic_charges)
        return all_loewdin_atomic_charges[-1]

    # ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
    # * MAYER POPULATION ANALYSIS *
    # ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
    @property
    def mayer_mulliken_gross_atomic_population(self):
        all_mayer_mulliken_gross_atomic_population = []
        for i, line_i in enumerate(self.contents):
            mayer_mulliken_gross_atomic_population = {}
            if "MAYER POPULATION ANALYSIS" in line_i:
                for line_j in self.contents[i:]:
                    if "TIMINGS" in line_j:
                        break
                    line_j_elements = line_j.split()
                    if len(line_j_elements) == 8:
                        element = p.to_element(line_j_elements[1])
                        element_num = f"{element}{line_j_elements[0]}"
                        mayer_mulliken_gross_atomic_population[element_num] = (
                            float(line_j_elements[2])
                        )
                all_mayer_mulliken_gross_atomic_population.append(
                    mayer_mulliken_gross_atomic_population
                )
        return all_mayer_mulliken_gross_atomic_population[-1]

    @property
    def mayer_total_nuclear_charge(self):
        all_mayer_total_nuclear_charge = []
        for i, line_i in enumerate(self.contents):
            mayer_total_nuclear_charge = {}
            if "MAYER POPULATION ANALYSIS" in line_i:
                for line_j in self.contents[i:]:
                    if "TIMINGS" in line_j:
                        break
                    line_j_elements = line_j.split()
                    if len(line_j_elements) == 8:
                        element = p.to_element(line_j_elements[1])
                        element_num = f"{element}{line_j_elements[0]}"
                        mayer_total_nuclear_charge[element_num] = float(
                            line_j_elements[3]
                        )
                all_mayer_total_nuclear_charge.append(
                    mayer_total_nuclear_charge
                )
        return all_mayer_total_nuclear_charge[-1]

    @property
    def mayer_mulliken_gross_atomic_charge(self):
        all_mayer_mulliken_gross_atomic_charge = []
        for i, line_i in enumerate(self.contents):
            mayer_mulliken_gross_atomic_charge = {}
            if "MAYER POPULATION ANALYSIS" in line_i:
                for line_j in self.contents[i:]:
                    if "TIMINGS" in line_j:
                        break
                    line_j_elements = line_j.split()
                    if len(line_j_elements) == 8:
                        element = p.to_element(line_j_elements[1])
                        element_num = f"{element}{line_j_elements[0]}"
                        mayer_mulliken_gross_atomic_charge[element_num] = (
                            float(line_j_elements[4])
                        )
                all_mayer_mulliken_gross_atomic_charge.append(
                    mayer_mulliken_gross_atomic_charge
                )
        return all_mayer_mulliken_gross_atomic_charge[-1]

    @property
    def mayer_total_valence(self):
        all_mayer_total_valence = []
        for i, line_i in enumerate(self.contents):
            mayer_total_valence = {}
            if "MAYER POPULATION ANALYSIS" in line_i:
                for line_j in self.contents[i:]:
                    if "TIMINGS" in line_j:
                        break
                    line_j_elements = line_j.split()
                    if len(line_j_elements) == 8:
                        element = p.to_element(line_j_elements[1])
                        element_num = f"{element}{line_j_elements[0]}"
                        mayer_total_valence[element_num] = float(
                            line_j_elements[5]
                        )
                all_mayer_total_valence.append(mayer_total_valence)
        return all_mayer_total_valence[-1]

    @property
    def mayer_bonded_valence(self):
        all_mayer_bonded_valence = []
        for i, line_i in enumerate(self.contents):
            mayer_bonded_valence = {}
            if "MAYER POPULATION ANALYSIS" in line_i:
                for line_j in self.contents[i:]:
                    if "TIMINGS" in line_j:
                        break
                    line_j_elements = line_j.split()
                    if len(line_j_elements) == 8:
                        element = p.to_element(line_j_elements[1])
                        element_num = f"{element}{line_j_elements[0]}"
                        mayer_bonded_valence[element_num] = float(
                            line_j_elements[6]
                        )
                all_mayer_bonded_valence.append(mayer_bonded_valence)
        return all_mayer_bonded_valence[-1]

    @property
    def mayer_free_valence(self):
        all_mayer_free_valence = []
        for i, line_i in enumerate(self.contents):
            mayer_free_valence = {}
            if "MAYER POPULATION ANALYSIS" in line_i:
                for line_j in self.contents[i + 1 :]:
                    if "TIMINGS" in line_j:
                        break
                    line_j_elements = line_j.split()
                    if len(line_j_elements) == 8:
                        element = p.to_element(line_j_elements[1])
                        element_num = f"{element}{line_j_elements[0]}"
                        mayer_free_valence[element_num] = float(
                            line_j_elements[-1]
                        )
                all_mayer_free_valence.append(mayer_free_valence)
        return all_mayer_free_valence[-1]

    @property
    def mayer_bond_orders_larger_than_zero_point_one(self):
        """Find the bonds with bond orders.

        Given input example:
         Mayer bond orders larger than 0.100000
        B(  0-O ,  1-H ) :   0.9959 B(  0-O ,  2-H ) :   0.9959.
        """

        from chemsmart.utils.repattern import mayer_bond_order_segment_pattern

        all_mayer_bond_orders_larger_than_zero_point_one = []

        for i, line_i in enumerate(self.contents):
            mayer_bond_orders_larger_than_zero_point_one = {}
            if "Mayer bond orders larger than 0.100000" in line_i:
                pattern = re.compile(mayer_bond_order_segment_pattern)

                for line_j in self.contents[i + 1 :]:
                    if len(line_j) == 0:
                        break
                    matches = pattern.findall(line_j)
                    for match in matches:
                        formatted_key = (
                            f"B({match[1]}{match[0]},{match[3]}{match[2]})"
                        )
                        formatted_value = float(match[4])
                        mayer_bond_orders_larger_than_zero_point_one[
                            formatted_key
                        ] = formatted_value
                all_mayer_bond_orders_larger_than_zero_point_one.append(
                    mayer_bond_orders_larger_than_zero_point_one
                )
        return all_mayer_bond_orders_larger_than_zero_point_one[-1]

    # ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
    # * HIRSHFELD ANALYSIS *
    # ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
    @property
    def total_integrated_alpha_density(self):
        all_hirshfeld_alpha_density = []
        for i, line_i in enumerate(self.contents):
            if "HIRSHFELD ANALYSIS" in line_i:
                for line_j in self.contents[i + 3 :]:
                    if len(line_j) == 0:
                        break
                    if "Total integrated alpha density" in line_j:
                        all_hirshfeld_alpha_density.append(
                            float(line_j.split()[-1])
                        )
        return all_hirshfeld_alpha_density[-1]

    @property
    def total_integrated_beta_density(self):
        all_hirshfeld_beta_density = []
        for i, line_i in enumerate(self.contents):
            if "HIRSHFELD ANALYSIS" in line_i:
                for line_j in self.contents[i + 3 :]:
                    if len(line_j) == 0:
                        break
                    if "Total integrated beta density" in line_j:
                        all_hirshfeld_beta_density.append(
                            float(line_j.split()[-1])
                        )
        return all_hirshfeld_beta_density[-1]

    def _get_hirshfeld_charges_and_spins(self):
        all_hirshfeld_charges = []
        all_hirshfeld_spins = []
        for i, line_i in enumerate(self.contents):
            hirshfeld_charges = {}
            hirshfeld_spins = {}
            if "HIRSHFELD ANALYSIS" in line_i:
                for line_j in self.contents[i + 6 :]:
                    if len(line_j) == 0:
                        break
                    line_j_elements = line_j.split()
                    if len(line_j_elements) == 4:
                        dict_key = (
                            f"{line_j_elements[1]}{int(line_j_elements[0])+1}"
                        )
                        charge_value = float(line_j_elements[2])
                        hirshfeld_charges[dict_key] = charge_value

                        spin_value = float(line_j_elements[3])
                        hirshfeld_spins[dict_key] = spin_value

                all_hirshfeld_charges.append(hirshfeld_charges)
                all_hirshfeld_spins.append(hirshfeld_spins)
        return all_hirshfeld_charges[-1], all_hirshfeld_spins[-1]

    @property
    def hirshfeld_charges(self):
        hirshfeld_charges, _ = self._get_hirshfeld_charges_and_spins()
        return hirshfeld_charges

    @property
    def hirshfeld_spin_densities(self):
        _, hirshfeld_spins = self._get_hirshfeld_charges_and_spins()
        return hirshfeld_spins

    # ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
    # *     ORCA property calculations      *
    # ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
    @property
    def dipole_moment_electric_contribution(self):
        all_dipole_moment_electric_contribution = []
        for i, line_i in enumerate(self.contents):
            dipole_moment = np.zeros((3, 1))
            if line_i == "DIPOLE MOMENT":
                for line_j in self.contents[i + 3 :]:
                    if len(line_j) == 0:
                        break
                    if "Electronic contribution:" in line_j:
                        line_j_elements = line_j.split()
                        dipole_moment = np.array(
                            [
                                float(line_j_elements[-3]),
                                float(line_j_elements[-2]),
                                float(line_j_elements[-1]),
                            ]
                        )
                all_dipole_moment_electric_contribution.append(dipole_moment)
        return all_dipole_moment_electric_contribution[-1]

    @property
    def dipole_moment_nuclear_contribution(self):
        all_dipole_moment_nuclear_contribution = []
        for i, line_i in enumerate(self.contents):
            dipole_moment = np.zeros((3, 1))
            if line_i == "DIPOLE MOMENT":
                for line_j in self.contents[i + 3 :]:
                    if len(line_j) == 0:
                        break
                    if "Nuclear contribution   :" in line_j:
                        line_j_elements = line_j.split()
                        dipole_moment = np.array(
                            [
                                float(line_j_elements[-3]),
                                float(line_j_elements[-2]),
                                float(line_j_elements[-1]),
                            ]
                        )
                all_dipole_moment_nuclear_contribution.append(dipole_moment)
        return all_dipole_moment_nuclear_contribution[-1]

    @property
    def total_dipole_moment(self):
        all_dipole_moment = []
        for i, line_i in enumerate(self.contents):
            dipole_moment = np.zeros((3, 1))
            if line_i == "DIPOLE MOMENT":
                for line_j in self.contents[i + 3 :]:
                    if len(line_j) == 0:
                        break
                    if "Total Dipole Moment    :" in line_j:
                        line_j_elements = line_j.split()
                        dipole_moment = np.array(
                            [
                                float(line_j_elements[-3]),
                                float(line_j_elements[-2]),
                                float(line_j_elements[-1]),
                            ]
                        )
                all_dipole_moment.append(dipole_moment)
        return all_dipole_moment[-1]

    @property
    def dipole_moment_in_au(self):
        all_dipole_moment = []
        for i, line_i in enumerate(self.contents):
            dipole_moment = 0.0
            if line_i == "DIPOLE MOMENT":
                for line_j in self.contents[i + 3 :]:
                    if len(line_j) == 0:
                        break
                    if "Magnitude (a.u.)       :" in line_j:
                        line_j_elements = line_j.split()
                        dipole_moment = float(line_j_elements[-1])
                all_dipole_moment.append(dipole_moment)
        return all_dipole_moment[-1]

    @property
    def dipole_moment_in_debye(self):
        all_dipole_moment = []
        for i, line_i in enumerate(self.contents):
            dipole_moment = 0.0
            if line_i == "DIPOLE MOMENT":
                for line_j in self.contents[i + 3 :]:
                    if len(line_j) == 0:
                        break
                    if "Magnitude (Debye)      :" in line_j:
                        line_j_elements = line_j.split()
                        dipole_moment = float(line_j_elements[-1])
                all_dipole_moment.append(dipole_moment)
        return all_dipole_moment[-1]

    @property
    def dipole_moment_along_axis_in_au(self):
        all_dipole_moment = []
        for i, line_i in enumerate(self.contents):
            dipole_moment = np.zeros((3, 1))
            if line_i == "Dipole components along the rotational axes:":
                for line_j in self.contents[i:]:
                    if len(line_j) == 0:
                        break
                    if "x,y,z [a.u.] :" in line_j:
                        line_j_elements = line_j.split()
                        dipole_moment = np.array(
                            [
                                float(line_j_elements[-3]),
                                float(line_j_elements[-2]),
                                float(line_j_elements[-1]),
                            ]
                        )
                all_dipole_moment.append(dipole_moment)
        return all_dipole_moment[-1]

    @property
    def dipole_moment_along_axis_in_debye(self):
        all_dipole_moment = []
        for i, line_i in enumerate(self.contents):
            dipole_moment = np.zeros((3, 1))
            if line_i == "Dipole components along the rotational axes:":
                for line_j in self.contents[i:]:
                    if len(line_j) == 0:
                        break
                    if "x,y,z [Debye]:" in line_j:
                        line_j_elements = line_j.split()
                        dipole_moment = np.array(
                            [
                                float(line_j_elements[-3]),
                                float(line_j_elements[-2]),
                                float(line_j_elements[-1]),
                            ]
                        )
                all_dipole_moment.append(dipole_moment)
        return all_dipole_moment[-1]

    @property
    def rotational_symmetry_number(self):
        """Obtain the rotational symmetry number from the output file."""
        for i, line_i in enumerate(self.contents):
            if line_i == "ENTHALPY":
                for line_j in self.contents[i:]:
                    if (
                        "Point Group:" in line_j
                        and "Symmetry Number:" in line_j
                    ):
                        line_j_elements = line_j.split()
                        rotational_symmetry_number = int(
                            line_j_elements[-1].strip()
                        )
                        return rotational_symmetry_number
        return None

    @property
    def rotational_constants_in_wavenumbers(self):
        """Rotational constants in wavenumbers."""
        all_rotational_constants_in_wavenumbers = []
        for i, line_i in enumerate(self.contents):
            rotational_constants_in_wavenumbers = []
            if line_i == "Rotational spectrum":
                for line_j in self.contents[i + 3 :]:
                    if len(line_j) == 0:
                        break
                    if "Rotational constants in cm-1:" in line_j:
                        line_j_elements = line_j.split()
                        for line_j_element in line_j_elements:
                            if is_float(line_j_element):
                                rotational_constants_in_wavenumbers.append(
                                    float(line_j_element)
                                )  # noqa: PERF401
                        all_rotational_constants_in_wavenumbers.append(
                            rotational_constants_in_wavenumbers
                        )
        return all_rotational_constants_in_wavenumbers[-1]

    @property
    def rotational_constants_in_MHz(self):
        all_rotational_constants_in_MHz = []
        for i, line_i in enumerate(self.contents):
            rotational_constants_in_MHz = []
            if line_i == "Rotational spectrum":
                for line_j in self.contents[i + 3 :]:
                    if len(line_j) == 0:
                        break
                    if "Rotational constants in MHz :" in line_j:
                        line_j_elements = line_j.split()
                        for line_j_element in line_j_elements:
                            if is_float(line_j_element):
                                rotational_constants_in_MHz.append(
                                    float(line_j_element)
                                )  # noqa: PERF401
                        all_rotational_constants_in_MHz.append(
                            rotational_constants_in_MHz
                        )
        return all_rotational_constants_in_MHz[-1]

    @property
    def vibrational_frequencies(self):
        vibrational_frequencies = []
        for i, line_i in enumerate(self.optimized_output_lines):
            if line_i == "VIBRATIONAL FREQUENCIES":
                for line_j in self.optimized_output_lines[i + 5 :]:
                    if "----------" in line_j:
                        break
                    # if 'Rotational constants in MHz :' in line_j:
                    line_j_elements = line_j.split()
                    if len(line_j_elements) != 0:
                        vibrational_frequencies.append(
                            float(line_j_elements[1])
                        )

        if self.molecule.is_monoatomic:
            # remove the first three frequencies (translations) for  monoatomic molecules
            vibrational_frequencies = vibrational_frequencies[3:]
        elif self.molecule.is_linear:
            # remove the first five frequencies (3 trans + 2 rot) for linear molecules
            vibrational_frequencies = vibrational_frequencies[5:]
        else:
            # remove the first six frequencies (3 trans + 3 rot) for non-linear molecules
            vibrational_frequencies = vibrational_frequencies[6:]

        return vibrational_frequencies

    @property
    def vib_freq_scale_factor(self):
        for i, line_i in enumerate(self.optimized_output_lines):
            if line_i == "VIBRATIONAL FREQUENCIES":
                for line_j in self.optimized_output_lines[i:]:
                    if "Scaling factor for frequencies =" in line_j:
                        line_j_elements = line_j.split()
                        return float(line_j_elements[-3])
        return None

    @property
    def molar_absorption_coefficients(self):
        """Eps  in units L/(mol*cm).

        The value under “eps” is the molar absorption coefficient, usually represented as ε.
        This number is directly proportional to the intensity of a given fundamental in an IR spectrum
        and is what is plotted by orca mapspc.
        """
        molar_absorption_coefficients = [
            0.0 for freq in self.vibrational_frequencies if freq == 0.0
        ]

        for i, line_i in enumerate(self.optimized_output_lines):
            if line_i == "IR SPECTRUM":
                for line_j in self.optimized_output_lines[i + 6 :]:
                    if (
                        "* The epsilon (eps) is given for a Dirac delta lineshape."
                        in line_j
                    ):
                        break
                    line_j_elements = line_j.split()
                    if len(line_j_elements) != 0:
                        molar_absorption_coefficients.append(
                            float(line_j_elements[2])
                        )
                return molar_absorption_coefficients
        return None

    @property
    def integrated_absorption_coefficients(self):
        """Units of km/mol.

        The values under “Int” are the integrated absorption coefficient [J. Comput. Chem., 2002, 23, 895.].
        """
        integrated_absorption_coefficients = [
            0.0 for freq in self.vibrational_frequencies if freq == 0.0
        ]

        for i, line_i in enumerate(self.optimized_output_lines):
            if line_i == "IR SPECTRUM":
                for line_j in self.optimized_output_lines[i + 6 :]:
                    if (
                        "* The epsilon (eps) is given for a Dirac delta lineshape."
                        in line_j
                    ):
                        break
                    line_j_elements = line_j.split()
                    if len(line_j_elements) != 0:
                        integrated_absorption_coefficients.append(
                            float(line_j_elements[3])
                        )
                return integrated_absorption_coefficients
        return None

    @property
    def transition_dipole_deriv_norm(self):
        """Units of a.u.

        “T**2” are the norm of the transition dipole derivatives,.
        """
        transition_dipole_deriv_norm = [
            0.0 for freq in self.vibrational_frequencies if freq == 0.0
        ]

        for i, line_i in enumerate(self.optimized_output_lines):
            if line_i == "IR SPECTRUM":
                for line_j in self.optimized_output_lines[i + 6 :]:
                    if (
                        "* The epsilon (eps) is given for a Dirac delta lineshape."
                        in line_j
                    ):
                        break
                    line_j_elements = line_j.split()
                    if len(line_j_elements) != 0:
                        transition_dipole_deriv_norm.append(
                            float(line_j_elements[4])
                        )
                return transition_dipole_deriv_norm
        return None

    @property
    def num_translation_and_rotation_modes(self):
        for i, line_i in enumerate(self.optimized_output_lines):
            if line_i == "IR SPECTRUM":
                for line_j in self.optimized_output_lines[i + 6 :]:
                    if (
                        "The first frequency considered to be a vibration is"
                        in line_j
                    ):
                        line_j_elements = line_j.split()
                        return int(line_j_elements[-1])
        return None

    @property
    def num_vibration_modes(self):
        for i, line_i in enumerate(self.contents):
            if line_i == "IR SPECTRUM":
                for line_j in self.contents[i + 6 :]:
                    if (
                        "The total number of vibrations considered is"
                        in line_j
                    ):
                        line_j_elements = line_j.split()
                        return int(line_j_elements[-1])
        return None

    # ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
    # *     THERMOCHEMISTRY      *
    # ** ** ** ** ** ** ** ** ** ** ** ** ** ** *
    @property
    def temperature_in_K(self):
        for i, line_i in enumerate(self.optimized_output_lines):
            if "THERMOCHEMISTRY" in line_i:
                for line_j in self.optimized_output_lines[i + 3 :]:
                    if "Temperature" in line_j:
                        line_j_elements = line_j.split()
                        return float(line_j_elements[-2])
        return None

    @property
    def pressure_in_atm(self):
        for i, line_i in enumerate(self.optimized_output_lines):
            if "THERMOCHEMISTRY" in line_i:
                for line_j in self.optimized_output_lines[i + 3 :]:
                    if "Pressure" in line_j:
                        line_j_elements = line_j.split()
                        return float(line_j_elements[-2])
        return None

    @property
    def total_mass_in_amu(self):
        """Total mass in amu."""
        for i, line_i in enumerate(self.contents):
            if "THERMOCHEMISTRY" in line_i:
                for line_j in self.contents[i + 3 :]:
                    if "Total Mass" in line_j:
                        line_j_elements = line_j.split()
                        return float(line_j_elements[-2])
        return None

    @property
    def mass(self):
        return self.total_mass_in_amu

    @property
    def moments_of_inertia(self):
        """Obtain moments of inertia from the output file of rotational
        constants in wavenumbers and convert to SI units (kg * m^2)."""
        all_moments_of_inertia = [
            units._hplanck
            / (8 * np.pi**2 * units._c * 1e2 * B)
            / (units._amu * (units.Ang / units.m) ** 2)
            for B in self.rotational_constants_in_wavenumbers
        ]
        return all_moments_of_inertia

    @property
    def internal_energy(self):
        """The inner energy is: U= E(el) + E(ZPE) + E(vib) + E(rot) + E(trans).

        E(el) = E(kin-el) + E(nuc-el) + E(el-el) + E(nuc-nuc)  is the total energy from the electronic structure
            calculation
        E(ZPE)  - the the zero temperature vibrational energy from the frequency calculation
        E(vib)  - the the finite temperature correction to E(ZPE) due to population of excited vibrational states
        E(rot)  - is the rotational thermal energy
        E(trans)- is the translational thermal energy.
        Default units are Hartree.
        """
        for i, line_i in enumerate(self.optimized_output_lines):
            if "INNER ENERGY" in line_i:
                for line_j in self.optimized_output_lines[i:]:
                    if "Total thermal energy" in line_j:
                        line_j_elements = line_j.split()
                        internal_energy_in_Hartree = float(line_j_elements[-2])
                        return internal_energy_in_Hartree
        return None

    @property
    def internal_energy_in_eV(self):
        return self.internal_energy * units.Hartree

    @property
    def electronic_energy(self):
        """E(el) = E(kin-el) + E(nuc-el) + E(el-el) + E(nuc-nuc).

        Total energy from the electronic structure calculation.
        Defaults to Hartree.
        """
        for i, line_i in enumerate(self.optimized_output_lines):
            if "INNER ENERGY" in line_i:
                for line_j in self.optimized_output_lines[i:]:
                    if "Electronic energy" in line_j:
                        line_j_elements = line_j.split()
                        electronic_energy_in_Hartree = float(
                            line_j_elements[-2]
                        )
                        return electronic_energy_in_Hartree
        return None

    @property
    def electronic_energy_in_eV(self):
        return self.electronic_energy * units.Hartree

    @property
    def zero_point_energy(self):
        """E(ZPE)  - the the zero temperature vibrational energy from the frequency calculation.
        Default units are Hartree."""
        for i, line_i in enumerate(self.optimized_output_lines):
            if "INNER ENERGY" in line_i:
                for line_j in self.optimized_output_lines[i:]:
                    if "Zero point energy" in line_j:
                        line_j_elements = line_j.split()
                        zpe_in_Hartree = float(line_j_elements[-4])
                        return zpe_in_Hartree
        return None

    @property
    def zero_point_energy_in_eV(self):
        return self.zero_point_energy * units.Hartree

    @property
    def thermal_vibration_correction(self):
        """E(vib)  - the the finite temperature correction to E(ZPE) due to population of excited vibrational states."""
        for i, line_i in enumerate(self.optimized_output_lines):
            if "INNER ENERGY" in line_i:
                for line_j in self.optimized_output_lines[i:]:
                    if "Thermal vibrational correction" in line_j:
                        line_j_elements = line_j.split()
                        thermal_vibration_correction_in_Hartree = float(
                            line_j_elements[-4]
                        )
                        return thermal_vibration_correction_in_Hartree
        return None

    @property
    def thermal_vibration_correction_in_eV(self):
        return self.thermal_vibration_correction * units.Hartree

    @property
    def thermal_rotation_correction(self):
        """E(rot)  - is the rotational thermal energy.
        Default units are Hartree."""
        for i, line_i in enumerate(self.optimized_output_lines):
            if "INNER ENERGY" in line_i:
                for line_j in self.optimized_output_lines[i:]:
                    if "Thermal rotational correction" in line_j:
                        line_j_elements = line_j.split()
                        thermal_rotation_correction_energy_in_Hartree = float(
                            line_j_elements[-4]
                        )
                        return thermal_rotation_correction_energy_in_Hartree
        return None

    @property
    def thermal_rotation_correction_in_eV(self):
        return self.thermal_rotation_correction * units.Hartree

    @property
    def thermal_translation_correction(self):
        """E(trans)- is the translational thermal energy."""
        for i, line_i in enumerate(self.optimized_output_lines):
            if "INNER ENERGY" in line_i:
                for line_j in self.optimized_output_lines[i:]:
                    if "Thermal translational correction" in line_j:
                        line_j_elements = line_j.split()
                        thermal_translation_correction_in_Hartree = float(
                            line_j_elements[-4]
                        )
                        return thermal_translation_correction_in_Hartree
        return None

    @property
    def thermal_translation_correction_in_eV(self):
        return self.thermal_translation_correction * units.Hartree

    @property
    def total_thermal_correction_due_to_trans_rot_vib(self):
        return (
            self.thermal_translation_correction
            + self.thermal_rotation_correction
            + self.thermal_vibration_correction
        )

    @property
    def total_correction(self):
        """Total correction due to Thermal (trans, rot, vib) + ZPE."""
        return (
            self.thermal_translation_correction
            + self.thermal_rotation_correction
            + self.thermal_vibration_correction
            + self.zero_point_energy
        )

    @property
    def enthalpy(self):
        """The enthalpy is H = U + kB*T.

        kB is Boltzmann's constant.
        Default units are Hartree.
        """
        for i, line_i in enumerate(self.optimized_output_lines):
            if line_i == "ENTHALPY":
                for line_j in self.optimized_output_lines[i:]:
                    if "Total Enthalpy" in line_j:
                        line_j_elements = line_j.split()
                        enthalpy_in_Hartree = float(line_j_elements[-2])
                        return enthalpy_in_Hartree
        return None

    @property
    def enthalpy_in_eV(self):
        return self.enthalpy * units.Hartree

    @property
    def thermal_enthalpy_correction(self):
        """kB*T term in the enthalpy is H = U + kB*T.

        kB is Boltzmann's constant.
        Default units are Hartree.
        """
        for i, line_i in enumerate(self.optimized_output_lines):
            if line_i == "ENTHALPY":
                for line_j in self.optimized_output_lines[i:]:
                    if "Thermal Enthalpy correction" in line_j:
                        line_j_elements = line_j.split()
                        thermal_enthalpy_correction_in_Hartree = float(
                            line_j_elements[-4]
                        )
                        return thermal_enthalpy_correction_in_Hartree
        return None

    @property
    def thermal_enthalpy_correction_in_eV(self):
        return self.thermal_enthalpy_correction * units.Hartree

    @property
    def electronic_entropy_no_temperature_in_SI(self):
        """Return electronic entropy in J/mol/K."""
        for i, line_i in enumerate(self.optimized_output_lines):
            if line_i == "ENTROPY":
                for line_j in self.optimized_output_lines[i + 10 :]:
                    if "Electronic entropy" in line_j:
                        line_j_elements = line_j.split()
                        electronic_entropy_hartree = float(line_j_elements[-4])
                        electronic_entropy_J_per_mol = (
                            electronic_entropy_hartree
                            * units.Hartree
                            * units.mol
                            / units.J
                        )
                        return (
                            electronic_entropy_J_per_mol
                            / self.temperature_in_K
                        )
        return None

    @property
    def vibrational_entropy_no_temperature_in_SI(self):
        """Return vibrational entropy in J/mol/K."""
        for i, line_i in enumerate(self.optimized_output_lines):
            if line_i == "ENTROPY":
                for line_j in self.optimized_output_lines[i + 10 :]:
                    if "Vibrational entropy" in line_j:
                        line_j_elements = line_j.split()
                        vibrational_entropy_hartree = float(
                            line_j_elements[-4]
                        )
                        vibrational_entropy_J_per_mol = (
                            vibrational_entropy_hartree
                            * units.Hartree
                            * units.mol
                            / units.J
                        )
                        return (
                            vibrational_entropy_J_per_mol
                            / self.temperature_in_K
                        )
        return None

    @property
    def rotational_entropy_no_temperature_in_SI(self):
        """Return rotational entropy in J/mol/K."""
        for i, line_i in enumerate(self.optimized_output_lines):
            if line_i == "ENTROPY":
                for line_j in self.optimized_output_lines[i + 10 :]:
                    if "Rotational entropy" in line_j:
                        line_j_elements = line_j.split()
                        rotational_entropy_hartree = float(line_j_elements[-4])
                        rotational_entropy_J_per_mol = (
                            rotational_entropy_hartree
                            * units.Hartree
                            * units.mol
                            / units.J
                        )
                        return (
                            rotational_entropy_J_per_mol
                            / self.temperature_in_K
                        )
        return None

    @property
    def translational_entropy_no_temperature_in_SI(self):
        """Return translational entropy in J/mol/K."""
        for i, line_i in enumerate(self.optimized_output_lines):
            if line_i == "ENTROPY":
                for line_j in self.optimized_output_lines[i + 10 :]:
                    if "Translational entropy" in line_j:
                        line_j_elements = line_j.split()
                        translational_entropy_hartree = float(
                            line_j_elements[-4]
                        )
                        translational_entropy_J_per_mol = (
                            translational_entropy_hartree
                            * units.Hartree
                            * units.mol
                            / units.J
                        )
                        return (
                            translational_entropy_J_per_mol
                            / self.temperature_in_K
                        )
        return None

    @property
    def entropy_in_J_per_mol_per_K(self):
        return (
            self.electronic_entropy_no_temperature_in_SI
            + self.translational_entropy_no_temperature_in_SI
            + self.rotational_entropy_no_temperature_in_SI
            + self.vibrational_entropy_no_temperature_in_SI
        )

    @property
    def entropy_TS(self):
        """The entropy contributions are T*S = T*(S(el)+S(vib)+S(rot)+S(trans)).

        ALREADY MULTIPLIED BY TEMPERATURE.
        The entropies will be listed as multiplied by the temperature to get units of energy.
        """
        for i, line_i in enumerate(self.optimized_output_lines):
            if line_i == "ENTROPY":
                for line_j in self.optimized_output_lines[i + 10 :]:
                    if "Final entropy term" in line_j:
                        line_j_elements = line_j.split()
                        entropy_hartree = float(line_j_elements[-4])
                        return entropy_hartree * units.Hartree
        return None

    @property
    def rotational_entropy_symmetry_correction_J_per_mol_per_K(self):
        """Return rotational entropy in J/mol/K for different symmetry numbers sn=1-12."""
        rotational_entropy_symmetry_correction_J_per_mol_per_K = {}
        for i, line_i in enumerate(self.optimized_output_lines):
            if "rotational entropy values for sn=1,12 :" in line_i:
                for line_j in self.optimized_output_lines[
                    i + 2 :
                ]:  # i+2 onwards
                    if "-------------------" in line_j:
                        break
                    new_line = line_j.replace("=", " ")
                    line_j_elements = new_line.split()
                    rotational_entropy_hartree = float(line_j_elements[-4])
                    rotational_entropy_J_per_mol = (
                        rotational_entropy_hartree
                        * units.Hartree
                        * units.mol
                        / units.J
                    )
                    rotational_entropy_J_per_mol_per_K = (
                        rotational_entropy_J_per_mol / self.temperature_in_K
                    )
                    rotational_entropy_J_per_mol_per_K = round(
                        rotational_entropy_J_per_mol_per_K, 6
                    )
                    rotational_entropy_symmetry_correction_J_per_mol_per_K[
                        int(line_j_elements[2])
                    ] = rotational_entropy_J_per_mol_per_K
                return rotational_entropy_symmetry_correction_J_per_mol_per_K
        return None

    @property
    def gibbs_free_energy(self):
        """The Gibbs free energy is G = H - T*S."""
        for i, line_i in enumerate(self.optimized_output_lines):
            if line_i == "GIBBS FREE ENERGY":
                for line_j in self.optimized_output_lines[i:]:
                    if "Final Gibbs free energy" in line_j:
                        line_j_elements = line_j.split()
                        entropy_hartree = float(line_j_elements[-2])
                        return entropy_hartree
        return None

    @property
    def gibbs_free_energy_in_eV(self):
        return self.gibbs_free_energy * units.Hartree

    # Below gives computing time/resources used by ORCA
    @cached_property
    def elapsed_walltime_by_jobs(self):
        elapsed_walltimes = []
        for line in self.contents:
            if line.startswith("TOTAL RUN TIME:"):
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

    @cached_property
    def cpu_runtime_by_jobs_core_hours(self):
        """CPU run time in core hours is total run time in hours times num
        of parallel processors."""
        # find the number of processors used
        for line in self.contents:
            if (
                "Program running with" in line
                and "parallel MPI-processes" in line
            ):
                match = re.search(orca_nproc_used_line_pattern, line)
                if match:
                    n_processors = int(match.group(1))
                    cpu_times = n_processors * self.total_elapsed_walltime
                    return cpu_times
            else:
                return self.total_elapsed_walltime

    @cached_property
    def service_units_by_jobs(self):
        """SUs defined as the JOB CPU time in hours."""
        return round(self.cpu_runtime_by_jobs_core_hours, 2)

    @cached_property
    def total_core_hours(self):
        return round(self.cpu_runtime_by_jobs_core_hours, 2)

    @cached_property
    def total_service_unit(self):
        return self.total_core_hours


class ORCAEngradFile(ORCAFileMixin):
    def __init__(self, filename):
        self.filename = filename
        """ Obtain energy and gradient of ORCA calculation"""
        self.natoms = self._get_natoms()
        self.energy = self._get_energy()
        self.gradient = self._get_gradient()
        self.molecule = self._get_molecule()

    def _get_natoms(self):
        for i, line in enumerate(self.contents):
            if "Number of atoms" in line:
                # check 3 lines following the match
                for content in self.contents[i + 1 : i + 4]:
                    try:
                        return int(content.split()[0])
                    except ValueError:
                        pass
        return None

    def _get_energy(self):
        """Get the total energy from the ORCA output file, in Hartree."""
        for i, line in enumerate(self.contents):
            if (
                "current total energy" in line
            ):  # in engrad file (all lower case)
                # check 3 lines following the match
                for content in self.contents[i + 1 : i + 4]:
                    try:
                        energy_in_hartree = float(content.split()[0])
                        return energy_in_hartree
                    except ValueError:
                        pass
        return None

    def _get_gradient(self):
        """Get the gradient from the ORCA output file, in Hartree/Bohr."""
        for i, line in enumerate(self.contents):
            if "current gradient" in line:
                # check 3N + 3 lines following the match, where N is number of atoms
                grad_data = []
                for content in self.contents[i + 1 : i + 3 * self.natoms + 4]:
                    try:
                        grad_value = float(
                            content.split()[0]
                        )  # in Hartree/Bohr
                        grad_data.append(grad_value)
                        # convert to eV/Angstrom
                        # grad_value * units.Hartree / units.Bohr

                    except ValueError:
                        pass
                return np.array(grad_data).reshape(self.natoms, 3)
        return None

    def _get_molecule(self):
        for i, line in enumerate(self.contents):
            if (
                "atomic numbers and current coordinates" in line
            ):  # in engrad file (all lower case)
                # read all the lines till the end
                symbols = []
                coords_tuple = []
                for content in self.contents[i + 1 :]:
                    if content.startswith("#") or len(content) == 0:
                        pass
                    try:
                        symbol = p.to_symbol(int(content.split()[0]))
                        symbols.append(symbol)
                        x_coord = (
                            float(content.split()[1]) * units.Bohr
                        )  # convert to Angstrom
                        y_coord = (
                            float(content.split()[2]) * units.Bohr
                        )  # convert to Angstrom
                        z_coord = (
                            float(content.split()[3]) * units.Bohr
                        )  # convert to Angstrom
                        coords_tuple.append((x_coord, y_coord, z_coord))
                    except ValueError:
                        pass
                return Molecule.from_symbols_and_positions_and_pbc_conditions(
                    list_of_symbols=symbols, positions=coords_tuple
                )
        return None


class ORCAQMMMOutput(ORCAOutput):
    """Class to parse the ORCA QMMM output files.

    This class is used to parse the QMMM output files generated by ORCA.
    It inherits from the ORCAOutput class and adds additional functionality
    specific to QMMM calculations.
    """

    def __init__(self, filename):
        super().__init__(filename)
        self.filename = filename

    @property
    def multiscale_model(self):
        pattern = re.compile(r"Multiscale model\s+\.{3}\s+(?P<model>\S+)")
        for line in self.contents:
            match = pattern.search(line)
            if match:
                return match.group("model")
        return None

    @property
    def qm2_method(self):
        pattern = re.compile(r"QM2 method\s+\.{3}\s+(?P<method>\S+)")
        for line in self.contents:
            match = pattern.search(line)
            if match:
                return match.group("method")
        return None

    @property
    def total_charge(self):
        pattern = re.compile(
            r"Charge of total system\s+\.{3}\s+(?P<charge>-?\d+)"
        )
        for line in self.contents:
            match = pattern.search(line)
            if match:
                return int(match.group("charge"))
        return None

    @property
    def scaling_factor_qm2(self):
        pattern = re.compile(
            r"Scaling factor for QM2 charges\s+\.{3}\s+(?P<scaling>\d+\.\d+)"
        )
        for line in self.contents:
            match = pattern.search(line)
            if match:
                return float(match.group("scaling"))
        return None

    @property
    def point_charges_in_qm_from_mm(self):
        point_charges_in_qm_from_mm, _ = self._get_point_charge_treatment()
        return point_charges_in_qm_from_mm

    @property
    def point_charges_in_qm_from_charge_shift(self):
        _, point_charges_in_qm_from_charge_shift = (
            self._get_point_charge_treatment()
        )
        return point_charges_in_qm_from_charge_shift

    @property
    def total_system_size(self):
        return self._get_partition_system_sizes()["QMMM"]

    @property
    def qm_system_size(self):
        return self.natoms

    @property
    def qm2_system_size(self):
        return self._get_partition_system_sizes()["QM2"]

    @property
    def number_of_link_atoms(self):
        pattern = re.compile(
            r"Number of link atoms\s+\.{3}\s+(?P<link_atoms>\d+)"
        )
        for line in self.contents:
            match = pattern.search(line)
            if match:
                return int(match.group("link_atoms"))
        return None

    @property
    def qm_plus_link_atoms_size(self):
        pattern = re.compile(
            r"Size of QM1 Subsystem plus link atoms\s+\.{3}\s+(?P<qm_plus_link_atoms>\d+)"
        )
        for line in self.contents:
            match = pattern.search(line)
            if match:
                return int(match.group("qm_plus_link_atoms"))
        return None

    @property
    def qm_region(self):
        pattern = re.compile(r"^QM1 Subsystem\s+\.{3}\s+(?P<qm_region>\d+)")
        qm_region = []
        reading_qm_lines = False
        for line in self.contents:
            match = pattern.search(line)
            if match:
                reading_qm_lines = True
                qm_atom_list = (
                    re.sub(r"^QM1 Subsystem\s+\.\.\.", "", line)
                    .strip()
                    .split()
                )
                for atom in qm_atom_list:
                    qm_region.append(int(atom))
            elif reading_qm_lines:
                if "..." in line.split():
                    break
                qm_atom_list = line.split()
                for atom in qm_atom_list:
                    qm_region.append(int(atom))
        # convert to 1-indexed
        qm_region = list(map(lambda x: x + 1, qm_region))
        qm_region = get_range_from_list(qm_region)
        return qm_region

    @property
    def qm2_energy_of_large_region(self):
        qm2_energy_of_large_region, _, _, _ = self._get_qmmm_energies()
        return qm2_energy_of_large_region

    @property
    def qm2_energy_of_small_region(self):
        _, qm2_energy_of_small_region, _, _ = self._get_qmmm_energies()
        return qm2_energy_of_small_region

    @property
    def qm_energy(self):
        _, _, _, qm_energy = self._get_qmmm_energies()
        return qm_energy

    @property
    def qm_qm2_energy(self):
        _, _, qm_qm2_energy, _ = self._get_qmmm_energies()
        return qm_qm2_energy

    def _get_qmmm_energies(self):
        """Get the QM/MM energies from the ORCA output file."""
        # todo: may need to refactor for 3-layer model
        qm2_energy_of_large_system = qm2_energy_of_small_system = (
            qm_qm2_energy
        ) = qm_energy = None
        pattern = re.compile(
            r"FINAL SINGLE POINT ENERGY(?:\s+\([^)]+\))?\s+(?P<energy>-?\d+\.\d+)"
        )
        # pattern = re.compile(r"FINAL SINGLE POINT ENERGY\s+\([^)]+\)\s+(?P<energy>-?\d+\.\d+)")
        for line in self.contents:
            match = pattern.search(line)
            if match:
                if "(L-QM2)" in line:
                    qm2_energy_of_large_system = (
                        float(match.group("energy")) * units.Hartree
                    )
                elif "(S-QM2)" in line:
                    qm2_energy_of_small_system = (
                        float(match.group("energy")) * units.Hartree
                    )
                elif "(QM/QM2)" in line:
                    qm_qm2_energy = (
                        float(match.group("energy")) * units.Hartree
                    )
                else:
                    qm_energy = float(match.group("energy")) * units.Hartree
            #
            # if "FINAL SINGLE POINT ENERGY" in line:
            #     match = re.search(
            #         r"FINAL SINGLE POINT ENERGY\s+\.{3}\s+(?P<energy>-?\d+\.\d+)",
            #         line,
            #     )
            #     if "(L-QM2)" in line:
            #         qm2_energy_of_large_system = (
            #             float(match.group("energy")) * units.Hartree
            #         )
            #     elif "(S-QM2)" in line:
            #         qm2_energy_of_small_system = (
            #             float(match.group("energy")) * units.Hartree
            #         )
            #     elif "(QM/QM2)" in line:
            #         qm_qm2_energy = (
            #             float(match.group("energy")) * units.Hartree
            #         )
            #     else:
            #         qm_energy = float(match.group("energy")) * units.Hartree
        return (
            qm2_energy_of_large_system,
            qm2_energy_of_small_system,
            qm_qm2_energy,
            qm_energy,
        )

    def _get_point_charge_treatment(self):
        pattern1 = re.compile(
            r"Point charges in QM calc\. from MM atoms\s*\.{3}\s*(\d+)"
        )
        pattern2 = re.compile(r"from charge shift scheme\s*\.{3}\s*(\d+)")
        point_charges_treatment = []
        for line in self.contents:
            match1 = pattern1.search(line)
            match2 = pattern2.search(line)
            if match1 is not None:
                point_charges_in_qm_from_mm = int(
                    re.sub(
                        r"Point charges in QM calc\. from MM atoms\s*\.{3}\s",
                        "",
                        line,
                    ).strip()
                )
                point_charges_treatment.append(point_charges_in_qm_from_mm)
            if match2 is not None:
                point_charges_in_qm_from_charge_shift = int(
                    re.sub(
                        r"from charge shift scheme\s*\.{3}\s", "", line
                    ).strip()
                )
                point_charges_treatment.append(
                    point_charges_in_qm_from_charge_shift
                )
        return point_charges_treatment[0], point_charges_treatment[1]

    def _get_partition_system_sizes(self):
        pattern1 = re.compile(
            r"Size of (?P<partition>\w+) System\s+\.{3}\s+(?P<atoms>\d+)"
        )
        pattern2 = re.compile(
            r"Size of (?P<partition>\w+) Subsystem\s+\.{3}\s+(?P<atoms>\d+)"
        )
        system_sizes = {}
        match = None
        for line in self.contents:
            if pattern1.search(line):
                match = pattern1.search(line)
            elif pattern2.search(line):
                match = pattern2.search(line)
            if match:
                partition = match.group("partition")
                atoms = int(match.group("atoms"))
                system_sizes[partition] = atoms
        return system_sizes
