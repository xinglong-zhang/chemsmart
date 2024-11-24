import logging
import math
import os
import re
from functools import cached_property
import numpy as np
from ase import units
from chemsmart.utils.mixins import FileMixin, ORCAFileMixin
from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.molecules.structure import CoordinateBlock
from chemsmart.utils.utils import is_float
from chemsmart.utils.repattern import (
    standard_coord_pattern,
    orca_input_coordinate_in_output,
)
from chemsmart.utils.periodictable import PeriodicTable


p = PeriodicTable()

logger = logging.getLogger(__name__)


class ORCAOutput(FileMixin, ORCAFileMixin):
    """ORCA output file with .out extension."""

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
            _line_contains_success_indicators(line) for line in self.contents[::-1]
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
        """
        optimized_geometry = {}
        for i, line_i in enumerate(self.optimized_output_lines):
            if "--- Optimized Parameters ---" in line_i:
                for line_j in self.optimized_output_lines[i + 5 :]:
                    if "---------------" in line_j:
                        break
                    line_elements = line_j.split()
                    if line_elements[1].lower().startswith("b"):  # bond
                        parameter = (
                            f"{line_elements[1]}{line_elements[2]}{line_elements[3]}"
                        )
                    if line_elements[1].lower().startswith("a"):  # angle
                        parameter = f"{line_elements[1]}{line_elements[2]}{line_elements[3]}{line_elements[4]}"
                    if line_elements[1].lower().startswith("d"):  # dihedral
                        parameter = (
                            f"{line_elements[1]}{line_elements[2]}{line_elements[3]}{line_elements[4]}"
                            f"{line_elements[5]}"
                        )
                    optimized_geometry[parameter] = float(line_elements[-1])
        ## the above result will return a dictionary storing the optimized parameters:
        ## optimized_geometry = { b(h1,o0) : 0.9627,  b(h2,o0) : 0.9627,  a(h1,o0, h2) : 103.35 }
        return optimized_geometry

    @property
    def final_structure(self):
        if self.optimized_output_lines is not None:
            return self._get_optimized_final_structure()
        return self._get_sp_structure()

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

    def _get_sp_structure(self):
        coordinate_lines = []
        for i, line in enumerate(self.contents):
            if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
                for line_j in self.contents[i:]:
                    pattern = re.compile(standard_coord_pattern)
                    if len(line_j) == 0:
                        break
                    if pattern.match(line_j):
                        coordinate_lines.append(line_j)
        cb = CoordinateBlock(coordinate_block=coordinate_lines)
        return cb.molecule

    @property
    def molecule(self):
        if self.final_structure is not None:
            return self.final_structure
        # Final structure is none as no "FINAL ENERGY EVALUATION AT THE STATIONARY POINT" line appear if the job
        # is not an optimisation job. In this case, get the atoms object from the SP output file.
        return self._get_atoms_from_sp_output_file()

    def _get_atoms_from_sp_output_file(self):
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
        # If atoms are not found, get them from the output file
        if molecule is None:
            molecule = self._get_input_structure_in_output()
        return molecule

    def _get_input_structure_in_output(self):
        """In ORCA output file, the input structure is rewritten and for single points, is same as the output structure.

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
                        return energy_in_hartree * units.Hartree

    @property
    def _get_sp_scf_energy(self):
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
                    return energy_in_hartree * units.Hartree

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
                return sp_energy_in_hartree * units.Hartree

    @property
    def final_nuclear_repulsion(self):
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
        return final_nuclear_repulsion_eV[-1]

    @property
    def final_electronic_energy(self):
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
        return final_electronic_energy_eV[-1]

    @property
    def one_electron_energy(self):
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
        return one_electron_energy_eV[-1]

    @property
    def two_electron_energy(self):
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
        return two_electron_energy_eV[-1]

    @property
    def max_cosx_asymmetry_energy(self):
        max_cosx_asymmetry_energy_hartree = []
        for line in self.contents:
            if "Max COSX asymmetry :" in line:
                energy_in_hartree = float(line.split()[-4])
                max_cosx_asymmetry_energy_hartree.append(energy_in_hartree)
        if len(max_cosx_asymmetry_energy_hartree) != 0:
            max_cosx_asymmetry_energy_eV = [
                value * units.Hartree for value in max_cosx_asymmetry_energy_hartree
            ]
            return max_cosx_asymmetry_energy_eV[-1]

    @property
    def potential_energy(self):
        potential_energy_hartree = []
        for line in self.contents:
            if "Potential Energy   :" in line:
                energy_in_hartree = float(line.split()[-4])
                potential_energy_hartree.append(energy_in_hartree)
        if len(potential_energy_hartree) != 0:
            potential_energy_eV = [
                value * units.Hartree for value in potential_energy_hartree
            ]
            return potential_energy_eV[-1]

    @property
    def kinetic_energy(self):
        kinetic_energy_hartree = []
        for line in self.contents:
            if "Kinetic Energy     :" in line:
                energy_in_hartree = float(line.split()[-4])
                kinetic_energy_hartree.append(energy_in_hartree)
        kinetic_energy_eV = [value * units.Hartree for value in kinetic_energy_hartree]
        if len(kinetic_energy_eV) != 0:
            return kinetic_energy_eV[-1]

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
        xc_energy_hartree = []
        for line in self.contents:
            if "E(XC)              :" in line:
                xc_energy_hartree.append(float(line.split()[-2]))
        if len(xc_energy_hartree) != 0:
            xc_energy_eV = [value * units.Hartree for value in xc_energy_hartree]
            return xc_energy_eV[-1]

    @property
    def dfet_embed_energy(self):
        dfet_embed_energy_hartree = []
        for line in self.contents:
            if "DFET-embed. en.    :" in line:
                dfet_embed_energy_hartree.append(float(line.split()[-2]))
        if len(dfet_embed_energy_hartree) != 0:
            dfet_embed_energy_eV = [
                value * units.Hartree for value in dfet_embed_energy_hartree
            ]
            return dfet_embed_energy_eV[-1]

    @property
    def orbital_occupancy(self):
        _, orbital_occupancy = self._get_orbital_energies_and_occupancy()
        return orbital_occupancy

    @property
    def orbital_energies(self):
        orbital_energies, _ = self._get_orbital_energies_and_occupancy()
        return orbital_energies

    def _get_orbital_energies_and_occupancy(self):
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
        occupied_energies = [energy for energy, occupancy in orbitals if occupancy == 2]
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
                    mulliken_atomic_charges[element_num] = float(line_j_elements[-1])
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
                        loewdin_atomic_charges[element_num] = float(line_j_elements[-1])
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
                        mayer_mulliken_gross_atomic_population[element_num] = float(
                            line_j_elements[2]
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
                all_mayer_total_nuclear_charge.append(mayer_total_nuclear_charge)
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
                        mayer_mulliken_gross_atomic_charge[element_num] = float(
                            line_j_elements[4]
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
                        mayer_total_valence[element_num] = float(line_j_elements[5])
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
                        mayer_bonded_valence[element_num] = float(line_j_elements[6])
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
                        mayer_free_valence[element_num] = float(line_j_elements[-1])
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
                        formatted_key = f"B({match[1]}{match[0]},{match[3]}{match[2]})"
                        formatted_value = float(match[4])
                        mayer_bond_orders_larger_than_zero_point_one[formatted_key] = (
                            formatted_value
                        )
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
                        all_hirshfeld_alpha_density.append(float(line_j.split()[-1]))
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
                        all_hirshfeld_beta_density.append(float(line_j.split()[-1]))
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
                        dict_key = f"{line_j_elements[1]}{int(line_j_elements[0])+1}"
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
    def hirshfeld_spins(self):
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
    def rotational_constants_in_wavenumbers(self):
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
                        vibrational_frequencies.append(float(line_j_elements[1]))
                return vibrational_frequencies
        return None

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
                        molar_absorption_coefficients.append(float(line_j_elements[2]))
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
                        transition_dipole_deriv_norm.append(float(line_j_elements[4]))
                return transition_dipole_deriv_norm
        return None

    @property
    def num_translation_and_rotation_modes(self):
        for i, line_i in enumerate(self.optimized_output_lines):
            if line_i == "IR SPECTRUM":
                for line_j in self.optimized_output_lines[i + 6 :]:
                    if "The first frequency considered to be a vibration is" in line_j:
                        line_j_elements = line_j.split()
                        return int(line_j_elements[-1])
        return None

    @property
    def num_vibration_modes(self):
        for i, line_i in enumerate(self.optimized_output_lines):
            if line_i == "IR SPECTRUM":
                for line_j in self.optimized_output_lines[i + 6 :]:
                    if "The total number of vibrations considered is" in line_j:
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
        for i, line_i in enumerate(self.optimized_output_lines):
            if "THERMOCHEMISTRY" in line_i:
                for line_j in self.optimized_output_lines[i + 3 :]:
                    if "Total Mass" in line_j:
                        line_j_elements = line_j.split()
                        return float(line_j_elements[-2])
        return None

    @property
    def internal_energy(self):
        """The inner energy is: U= E(el) + E(ZPE) + E(vib) + E(rot) + E(trans).

        E(el) = E(kin-el) + E(nuc-el) + E(el-el) + E(nuc-nuc)  is the total energy from the electronic structure
            calculation
        E(ZPE)  - the the zero temperature vibrational energy from the frequency calculation
        E(vib)  - the the finite temperature correction to E(ZPE) due to population of excited vibrational states
        E(rot)  - is the rotational thermal energy
        E(trans)- is the translational thermal energy.
        """
        for i, line_i in enumerate(self.optimized_output_lines):
            if "INNER ENERGY" in line_i:
                for line_j in self.optimized_output_lines[i:]:
                    if "Total thermal energy" in line_j:
                        line_j_elements = line_j.split()
                        internal_energy_in_Hartree = float(line_j_elements[-2])
                        return internal_energy_in_Hartree * units.Hartree
        return None

    @property
    def electronic_energy(self):
        """E(el) = E(kin-el) + E(nuc-el) + E(el-el) + E(nuc-nuc).

        Total energy from the electronic structure calculation.
        """
        for i, line_i in enumerate(self.optimized_output_lines):
            if "INNER ENERGY" in line_i:
                for line_j in self.optimized_output_lines[i:]:
                    if "Electronic energy" in line_j:
                        line_j_elements = line_j.split()
                        electronic_energy_in_Hartree = float(line_j_elements[-2])
                        return electronic_energy_in_Hartree * units.Hartree
        return None

    @property
    def zero_point_energy(self):
        """E(ZPE)  - the the zero temperature vibrational energy from the frequency calculation."""
        for i, line_i in enumerate(self.optimized_output_lines):
            if "INNER ENERGY" in line_i:
                for line_j in self.optimized_output_lines[i:]:
                    if "Zero point energy" in line_j:
                        line_j_elements = line_j.split()
                        zpe_in_Hartree = float(line_j_elements[-4])
                        return zpe_in_Hartree * units.Hartree
        return None

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
                        return thermal_vibration_correction_in_Hartree * units.Hartree
        return None

    @property
    def thermal_rotation_correction(self):
        """E(rot)  - is the rotational thermal energy."""
        for i, line_i in enumerate(self.optimized_output_lines):
            if "INNER ENERGY" in line_i:
                for line_j in self.optimized_output_lines[i:]:
                    if "Thermal rotational correction" in line_j:
                        line_j_elements = line_j.split()
                        thermal_rotation_correction_energy_in_Hartree = float(
                            line_j_elements[-4]
                        )
                        return (
                            thermal_rotation_correction_energy_in_Hartree
                            * units.Hartree
                        )
        return None

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
                        return thermal_translation_correction_in_Hartree * units.Hartree
        return None

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
        """
        for i, line_i in enumerate(self.optimized_output_lines):
            if line_i == "ENTHALPY":
                for line_j in self.optimized_output_lines[i:]:
                    if "Total Enthalpy" in line_j:
                        line_j_elements = line_j.split()
                        enthalpy_in_Hartree = float(line_j_elements[-2])
                        return enthalpy_in_Hartree * units.Hartree
        return None

    @property
    def thermal_enthalpy_correction(self):
        """kB*T term in the enthalpy is H = U + kB*T.

        kB is Boltzmann's constant.
        """
        for i, line_i in enumerate(self.optimized_output_lines):
            if line_i == "ENTHALPY":
                for line_j in self.optimized_output_lines[i:]:
                    if "Thermal Enthalpy correction" in line_j:
                        line_j_elements = line_j.split()
                        thermal_enthalpy_correction_in_Hartree = float(
                            line_j_elements[-4]
                        )
                        return thermal_enthalpy_correction_in_Hartree * units.Hartree
        return None

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
                        return electronic_entropy_J_per_mol / self.temperature_in_K
        return None

    @property
    def vibrational_entropy_no_temperature_in_SI(self):
        """Return vibrational entropy in J/mol/K."""
        for i, line_i in enumerate(self.optimized_output_lines):
            if line_i == "ENTROPY":
                for line_j in self.optimized_output_lines[i + 10 :]:
                    if "Vibrational entropy" in line_j:
                        line_j_elements = line_j.split()
                        vibrational_entropy_hartree = float(line_j_elements[-4])
                        vibrational_entropy_J_per_mol = (
                            vibrational_entropy_hartree
                            * units.Hartree
                            * units.mol
                            / units.J
                        )
                        return vibrational_entropy_J_per_mol / self.temperature_in_K
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
                        return rotational_entropy_J_per_mol / self.temperature_in_K
        return None

    @property
    def translational_entropy_no_temperature_in_SI(self):
        """Return translational entropy in J/mol/K."""
        for i, line_i in enumerate(self.optimized_output_lines):
            if line_i == "ENTROPY":
                for line_j in self.optimized_output_lines[i + 10 :]:
                    if "Translational entropy" in line_j:
                        line_j_elements = line_j.split()
                        translational_entropy_hartree = float(line_j_elements[-4])
                        translational_entropy_J_per_mol = (
                            translational_entropy_hartree
                            * units.Hartree
                            * units.mol
                            / units.J
                        )
                        return translational_entropy_J_per_mol / self.temperature_in_K
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
                for line_j in self.optimized_output_lines[i + 2 :]:  # i+2 onwards
                    if "-------------------" in line_j:
                        break
                    new_line = line_j.replace("=", " ")
                    line_j_elements = new_line.split()
                    rotational_entropy_hartree = float(line_j_elements[-4])
                    rotational_entropy_J_per_mol = (
                        rotational_entropy_hartree * units.Hartree * units.mol / units.J
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
                        return entropy_hartree * units.Hartree
        return None

    @property
    def total_run_time_hours(self):
        for line in self.contents:
            elapsed_walltime = []
            if line.startswith("TOTAL RUN TIME:"):
                n_days = float(line.lower().split("days")[0].strip().split()[-1])
                n_hours = float(line.lower().split("hours")[0].strip().split()[-1])
                n_minutes = float(line.lower().split("minutes")[0].strip().split()[-1])
                n_seconds = float(line.lower().split("seconds")[0].strip().split()[-1])
                total_seconds = (
                    n_days * 24 * 60 * 60
                    + n_hours * 60 * 60
                    + n_minutes * 60
                    + n_seconds
                )
                total_hours = round(total_seconds / 3600, 4)
                elapsed_walltime.append(total_hours)
        return sum(elapsed_walltime)

    # def read_settings(self):
    #     from chemsmart.jobs.orca.settings import ORCAJobSettings
    #
    #     dv = ORCAJobSettings.default()
    #     return ORCAJobSettings(
    #         ab_initio=self.ab_initio,
    #         functional=self.functional,
    #         dispersion=self.dispersion,
    #         basis=self.basis,
    #         aux_basis=self.aux_basis,
    #         extrapolation_basis=self.extrapolation_basis,
    #         defgrid=self.defgrid,
    #         scf_tol=self.scf_tol,
    #         scf_algorithm=self.scf_algorithm,
    #         scf_maxiter=self.scf_maxiter,
    #         scf_convergence=self.scf_convergence,
    #         charge=self.charge,
    #         multiplicity=self.multiplicity,
    #         gbw=dv.gbw,
    #         freq=self.freq,
    #         numfreq=self.numfreq,
    #         dipole=self.dipole,
    #         quadrupole=self.quadrupole,
    #         mdci_cutoff=self.mdci_cutoff,
    #         mdci_density=self.mdci_density,
    #         job_type=self.job_type,
    #         solvent_model=self.solvent_model,
    #         solvent_id=self.solvent_id,
    #         additional_route_parameters=dv.additional_route_parameters,
    #         route_to_be_written=dv.route_to_be_written,
    #         modred=dv.modred,
    #         gen_genecp=dv.gen_genecp,
    #         heavy_elements=dv.heavy_elements,
    #         heavy_elements_basis=dv.heavy_elements_basis,
    #         light_elements_basis=dv.light_elements_basis,
    #         custom_solvent=dv.custom_solvent,
    #         forces=dv.forces,
    #     )


class ORCAEngradFile(FileMixin):
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
        for i, line in enumerate(self.contents):
            if "current total energy" in line:  # in engrad file (all lower case)
                # check 3 lines following the match
                for content in self.contents[i + 1 : i + 4]:
                    try:
                        energy_in_hartree = float(content.split()[0])
                        # convert to eV
                        return energy_in_hartree * units.Hartree
                    except ValueError:
                        pass
        return None

    def _get_gradient(self):
        for i, line in enumerate(self.contents):
            if "current gradient" in line:
                # check 3N + 3 lines following the match, where N is number of atoms
                grad_data = []
                for content in self.contents[i + 1 : i + 3 * self.natoms + 4]:
                    try:
                        grad_value = float(content.split()[0])  # in Hartree/Bohr
                        # convert to eV/Angstrom
                        grad_data.append(grad_value * units.Hartree / units.Bohr)
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
