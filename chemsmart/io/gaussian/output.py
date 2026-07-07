import logging
import re
from functools import cached_property

import numpy as np
from ase import units

from chemsmart.io.molecules.structure import CoordinateBlock
from chemsmart.utils.constants import (
    cal_to_joules,
    energy_conversion,
    joule_per_mol_to_hartree,
    kcal_per_mol_to_hartree,
)
from chemsmart.utils.io import clean_duplicate_structure, create_molecule_list
from chemsmart.utils.mixins import GaussianFileMixin
from chemsmart.utils.periodictable import PeriodicTable
from chemsmart.utils.repattern import (
    basis_primitive_line_pattern,
    basis_shell_header_pattern,
    ecp_center_header_pattern,
    ecp_term_pattern,
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

        for line in reversed(contents):
            if not line:
                continue
            if "Normal termination of Gaussian" in line:
                logger.debug(f"File {self.filename} terminated normally.")
                return True
            break

        logger.debug(f"File {self.filename} has error termination.")
        return False

    @property
    def heavy_elements(self):
        """List of element symbols that use an explicitly defined basis set
        in the gen/genecp section.
        """
        return self._genecp_info["heavy_elements"]

    @property
    def heavy_elements_basis(self):
        """
        Structured basis specification for heavy elements.
        When the route string starts with #p, Gaussian prints the expanded
        orbital exponents and coefficients for elements that use an explicitly
        defined gen/genecp basis.  This property parses that information into
        a dict keyed by element symbol, where each value is an ordered list of
        shell dicts:
            {
                'Br': [
                    {'shell': 'S', 'primitives': [(exp, coef), ...]},
                    ...
                    {'shell': 'F', 'primitives': [(exp, coef)]},
                ],
            }
        Returns None when the information is unavailable.
        """
        return self._genecp_info["heavy_elements_basis"]

    @property
    def heavy_elements_ecp(self):
        """
        ECP (effective core potential) specification for heavy elements.
        Only available for genecp calculations with #p verbose output.
        Returns a dict keyed by element symbol:
            {
                'Ag': {
                    'n_valence_electrons': 19,
                    'channels': [
                        {
                            'name': 'F and up',
                            'terms': [(r_power, exponent, coefficient, so_coefficient), ...],
                        },
                        ...
                    ],
                }
            }
        Returns None when not a genecp calculation or information is unavailable.
        """
        return self._genecp_info["heavy_elements_ecp"]

    @property
    def light_elements(self):
        """List of element symbols that use a standard basis set
        in the gen/genecp section.
        """
        return self._genecp_info["light_elements"]

    @property
    def light_elements_basis(self):
        """Basis set used for light elements."""
        return self._genecp_info["light_elements_basis"]

    @cached_property
    def _genecp_info(self):
        """
        Parse gen/genecp basis set information from the output file.
        Only available when the route string starts with #p.
        """
        result = {
            "light_elements": None,
            "light_elements_basis": None,
            "heavy_elements": None,
            "heavy_elements_basis": None,
            "heavy_elements_ecp": None,
        }
        if self.gen_genecp is None:
            return result
        try:
            atom_symbols = self.symbols
        except Exception:
            return result
        if not atom_symbols:
            return result
        shell_header = re.compile(basis_shell_header_pattern)
        light_elements_set = set()
        light_elements_basis = None
        heavy_basis_by_element: dict[str, list[dict]] = {}

        for i, line in enumerate(self.contents):
            if not line.startswith("General basis read from cards:"):
                continue
            j = i + 1
            while j < len(self.contents):
                line_j = self.contents[j].strip()
                if not line_j:
                    j += 1
                    continue
                if not line_j.startswith("Centers:"):
                    break
                # Parse "Centers: N1 N2 ..." -> list of 1-based center indices
                center_nums = []
                for token in line_j.split()[1:]:
                    try:
                        center_nums.append(int(token))
                    except ValueError:
                        pass
                j += 1
                # Skip blank lines
                while j < len(self.contents) and not self.contents[j].strip():
                    j += 1
                if j >= len(self.contents):
                    break
                content_line = self.contents[j].strip()
                if shell_header.match(content_line):
                    # Explicit orbital specification -> heavy element
                    # Collect every line of this block up to "****"
                    block_lines = []
                    while (
                        j < len(self.contents)
                        and self.contents[j].strip() != "****"
                    ):
                        block_lines.append(self.contents[j].strip())
                        j += 1
                    shells = self._parse_explicit_basis_block(block_lines)
                    for cn in center_nums:
                        if 1 <= cn <= len(atom_symbols):
                            sym = atom_symbols[cn - 1]
                            heavy_basis_by_element[sym] = shells
                else:
                    # Named basis set -> light element
                    basis_name = content_line
                    if light_elements_basis is None:
                        light_elements_basis = basis_name
                    for cn in center_nums:
                        if 1 <= cn <= len(atom_symbols):
                            light_elements_set.add(atom_symbols[cn - 1])
                    # Advance to the next "****" separator
                    while (
                        j < len(self.contents)
                        and self.contents[j].strip() != "****"
                    ):
                        j += 1
                j += 1  # move past "****"
            break

        if light_elements_set:
            result["light_elements"] = p.sorted_periodic_table_list(
                list(light_elements_set)
            )
            result["light_elements_basis"] = light_elements_basis
        if heavy_basis_by_element:
            result["heavy_elements"] = p.sorted_periodic_table_list(
                list(heavy_basis_by_element.keys())
            )
            result["heavy_elements_basis"] = heavy_basis_by_element
        ecp = self._parse_pseudopotential_section()
        if ecp:
            result["heavy_elements_ecp"] = ecp

        return result

    @staticmethod
    def _parse_explicit_basis_block(block_lines):
        """
        Parse lines of an explicitly printed Gaussian basis block into a
        list of shell dictionaries.
        Each dict has the keys:
            shell – angular-momentum letter ('S', 'P', 'D', 'F', 'H')
            primitives – list of (exponent, coefficient) float tuples
        """
        shell_header = re.compile(basis_shell_header_pattern)
        primitive_line = re.compile(basis_primitive_line_pattern)
        shells = []
        current_shell = None
        for raw in block_lines:
            line = raw.strip()
            m_shell = shell_header.match(line)
            m_prim = primitive_line.search(line)
            if m_shell:
                if current_shell is not None:
                    shells.append(current_shell)
                current_shell = {"shell": m_shell.group(1), "primitives": []}
            elif m_prim and current_shell is not None:
                exp = float(m_prim.group(1).upper().replace("D", "E"))
                coef = float(m_prim.group(2).upper().replace("D", "E"))
                current_shell["primitives"].append((exp, coef))
        if current_shell is not None:
            shells.append(current_shell)
        return shells

    def _parse_pseudopotential_section(self):
        """
        Parse the Pseudopotential Parameters section into an element ->
        ECP information dictionary.
        Each ECP dict has the keys:
            n_valence_electrons – number of explicitly treated valence electrons
            channels – list of ECP channel dictionaries
        Each channel dict has the keys:
            name – channel label (e.g. 'F and up', 'S - F')
            terms – list of
                (r_power, exponent, coefficient, spin_orbit_coefficient)
                float tuples
        Elements with "No pseudopotential on this center." are omitted.
        Returns an empty dict if the section is not found.
        """
        try:
            atom_symbols = self.symbols
        except Exception:
            return {}
        center_re = re.compile(ecp_center_header_pattern)
        term_re = re.compile(ecp_term_pattern)
        result = {}
        eq_count = 0
        current = None

        def _flush_channel():
            if (
                current
                and current["_channel_name"]
                and current["_channel_terms"]
            ):
                current["channels"].append(
                    {
                        "name": current["_channel_name"],
                        "terms": current["_channel_terms"],
                    }
                )
            if current:
                current["_channel_name"] = None
                current["_channel_terms"] = []

        def _flush_center():
            nonlocal current
            if current and current["channels"]:
                sym = atom_symbols[current["center_num"] - 1]
                result[sym] = {
                    "n_valence_electrons": current["n_valence_electrons"],
                    "channels": current["channels"],
                }
            current = None

        for line in self.contents:
            if not line:
                continue
            if "Pseudopotential Parameters" in line:
                eq_count = 0
                continue
            if line.startswith("="):
                eq_count += 1
                if eq_count == 3:
                    _flush_channel()
                    _flush_center()
                    break
                continue
            if eq_count < 2:
                continue

            tokens = line.split()
            is_center = len(tokens) in (2, 3) and all(
                t.isdigit() for t in tokens
            )
            is_term = (
                len(tokens) == 4 and tokens[0].isdigit() and "." in tokens[1]
            )
            if is_center:
                _flush_channel()
                _flush_center()
                current = {
                    "center_num": int(tokens[0]),
                    "n_valence_electrons": (
                        int(tokens[2]) if len(tokens) == 3 else None
                    ),
                    "channels": [],
                    "_channel_name": None,
                    "_channel_terms": [],
                }
            elif "No pseudopotential on this center" in line:
                pass  # light element, skip
            elif is_term and current is not None:
                r_power = int(tokens[0])
                exp = float(tokens[1].upper().replace("D", "E"))
                coef = float(tokens[2].upper().replace("D", "E"))
                so_coef = float(tokens[3].upper().replace("D", "E"))
                current["_channel_terms"].append((r_power, exp, coef, so_coef))
            elif (
                current is not None
                and not center_re.match(line)
                and not term_re.match(line)
            ):
                _flush_channel()
                current["_channel_name"] = line
                current["_channel_terms"] = []
        return result

    @property
    def custom_solvent(self):
        """
        Parse custom/generic solvent parameters from the verbose Gaussian
        PCM/SMD section into a structured dictionary with standardized keys.
        Returns a stuctured representation of the custom/generic solvent
        parameters:
            {
                "SolventName": "1,1,1,3,3,3-hexafluoropropan-2-ol",
                "Eps": 16.7,
                "EpsInf": 1.625625,
                "HbondAcidity": 0.77,
                "HbondBasicity": 0.1,
                "SurfaceTensionAtInterface": 23.23,
                "CarbonAromaticity": 0.0,
                "ElectronegativeHalogenicity": 0.6,
            }
            or None if no custom solvent section is found.
        """
        non_standard_marker = "Using the following non-standard input for PCM:"
        if not any(
            line.startswith(non_standard_marker) for line in self.contents
        ):
            return None
        key_map = [
            ("Eps(infinity)", "EpsInf"),
            ("Eps ", "Eps"),
            ("Hydrogen bond acidity", "HbondAcidity"),
            ("Hydrogen bond basicity", "HbondBasicity"),
            ("Surface tension at interface", "SurfaceTensionAtInterface"),
            ("Carbon aromaticity", "CarbonAromaticity"),
            ("Electronegative halogenicity", "ElectronegativeHalogenicity"),
        ]
        params = {}
        inside = False
        for line in self.contents:
            if line.startswith("Solvent") and ":" in line:
                name = line.split(":", 1)[1].strip().rstrip(",").strip()
                params["SolventName"] = name
                inside = True
                continue
            if not inside:
                continue
            if line.startswith("---"):
                break
            for prefix, key in key_map:
                if line.startswith(prefix):
                    value_str = line.split("=", 1)[1].strip().split()[0]
                    params[key] = float(value_str)
                    break
        return params if params else None

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

        def _is_oldform_coordinate_line(line):
            elements = [element.strip() for element in line.split(",")]
            if len(elements) < 5:
                return False

            atom = elements[0]
            if atom in p.PERIODIC_TABLE:
                pass
            else:
                try:
                    atom_float = float(atom)
                    if not atom_float.is_integer():
                        return False
                except ValueError:
                    return False

            try:
                float(elements[2])
                float(elements[3])
                float(elements[4])
            except ValueError:
                return False

            return True

        def _normalize_oldform_coordinate_line(line):
            elements = [element.strip() for element in line.split(",")]
            return f"{elements[0]} {elements[2]} {elements[3]} {elements[4]}"

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
                if len(coordinates_block_lines_list) > 0:
                    break
            elif "Redundant internal coordinates found in file." in line:
                for j_line in self.contents[i + 1 :]:
                    if len(j_line) == 0:
                        if len(coordinates_block_lines_list) > 0:
                            break
                        continue
                    if j_line.startswith("Charge ="):
                        logger.debug(f"Skipping line: {j_line}")
                        continue
                    if _is_oldform_coordinate_line(j_line):
                        coordinates_block_lines_list.append(
                            _normalize_oldform_coordinate_line(j_line)
                        )
                    elif len(coordinates_block_lines_list) > 0:
                        break
                if len(coordinates_block_lines_list) > 0:
                    break
        cb = CoordinateBlock(coordinate_block=coordinates_block_lines_list)
        return cb

    @cached_property
    def all_structures(self):
        """
        Obtain all the structures from the output file.
        Use Standard orientations to get the
        structures; if not, use input orientations.
        Include their corresponding energy and forces if present.

        Each structure additionally carries these attributes when available:
          - point_group and rotational_constants — attached per step,
            aligned 1:1 with structures.
          - dipole_moment and dipole_moment_magnitude — attached to the
            final structure only.
          - mulliken_atomic_charges and mulliken_spin_densities — attached
            to the final structure only.
          - The following are attached to the final structure only, and only
            when a frequency section is present:
              vibrational_frequencies, vibrational_reduced_masses,
              vibrational_force_constants, vibrational_ir_intensities,
              vibrational_mode_symmetries, vibrational_modes,
              rotational_symmetry_number.
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
        rot_consts = list(self.all_rotational_constants) or None
        point_groups = list(self.all_point_groups) or None

        # Helper to drop the first item across all arrays (when present)
        def drop_first():
            nonlocal orientations, orientations_pbc, energies, forces, rot_consts, point_groups
            if orientations:
                orientations = orientations[1:]
            if orientations_pbc:
                orientations_pbc = orientations_pbc[1:]
            if energies:
                energies = energies[1:]
            if rot_consts:
                rot_consts = rot_consts[1:]
            if point_groups:
                point_groups = point_groups[1:]
            # Forces do not need to be dropped here because no
            # force computation occurs at the first link job.
            # This is intentional: the forces array is
            # already aligned with the relevant orientations.

        # Helper to keep only the last frame across all arrays
        def keep_last_only():
            nonlocal orientations, orientations_pbc, energies, forces, rot_consts, point_groups
            orientations = orientations[-1:] if orientations else []
            orientations_pbc = (
                orientations_pbc[-1:] if orientations_pbc else []
            )
            if energies:
                energies = energies[-1:]
            if forces:
                forces = forces[-1:]
            if rot_consts:
                rot_consts = rot_consts[-1:]
            if point_groups:
                point_groups = point_groups[-1:]

        # Right-trim auxiliaries to the number of
        # orientations (no data loss in orientations)
        def align_lengths_to_orientations():
            nonlocal orientations, orientations_pbc, energies, forces, rot_consts, point_groups
            n = len(orientations)
            if orientations_pbc and len(orientations_pbc) > n:
                orientations_pbc = orientations_pbc[:n]
            if energies and len(energies) > n:
                energies = energies[:n]
            if forces and len(forces) > n:
                forces = forces[:n]
            if rot_consts and len(rot_consts) > n:
                rot_consts = rot_consts[:n]
            if point_groups and len(point_groups) > n:
                point_groups = point_groups[:n]

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

        # Calculate is_optimized_structure_list
        is_optimized = [False] * num_structures_to_use
        if self.optimized_steps_indices and self.include_intermediate:
            for idx in self.optimized_steps_indices:
                if 0 <= idx < len(is_optimized):
                    is_optimized[idx] = True
        elif self.normal_termination:
            is_optimized[-1] = True

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
            is_optimized_structure_list=is_optimized,
            rotational_constants_list=rot_consts,
            point_groups_list=point_groups,
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
        if self.optimized_steps_indices and not self.include_intermediate:
            logger.debug(
                "Ignoring intermediate optimization steps (constrained opt)."
            )
            all_structures = [
                all_structures[i] for i in self.optimized_steps_indices
            ]
            # Since we filtered to only optimized steps, mark all as optimized
            for mol in all_structures:
                mol.is_optimized_structure = True

        logger.debug(
            "Attaching vibrational data to the final structure if available..."
        )

        last_mol = all_structures[-1]
        # Attach vibrational data to the final structure if available
        if self.num_vib_frequencies:
            all_structures[-1] = self._attach_vib_metadata(last_mol)
        # Attach Mulliken charges/spin densities and dipole moments to the final structure
        if self.mulliken_atomic_charges is not None:
            last_mol.mulliken_atomic_charges = self.mulliken_atomic_charges
        if self.mulliken_spin_densities is not None:
            last_mol.mulliken_spin_densities = self.mulliken_spin_densities
        if self.has_dipole_moment:
            last_mol.dipole_moment = self.all_dipole_moments[-1]
            last_mol.dipole_moment_magnitude = (
                self.all_dipole_moment_magnitudes[-1]
            )

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
        if self.all_structures:
            return self.all_structures[-1]
        return self.input_coordinates_block.molecule

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

        # Attach rotational symmetry number
        mol.rotational_symmetry_number = self.rotational_symmetry_number

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

    @cached_property
    def thermal_vibration_correction(self):
        """
        Thermal vibration correction in Hartree.
        """
        if self.zero_point_energy is not None:
            for i, line_i in enumerate(self.contents):
                if "E (Thermal)" in line_i and "CV" in line_i:
                    for line_j in self.contents[i + 1 :]:
                        if "Vibrational" in line_j:
                            return (
                                float(line_j.split()[1])
                                * kcal_per_mol_to_hartree
                                - self.zero_point_energy
                            )
        return None

    @cached_property
    def thermal_rotation_correction(self):
        """
        Thermal rotation correction in Hartree.
        """
        for i, line_i in enumerate(self.contents):
            if "E (Thermal)" in line_i and "CV" in line_i:
                for line_j in self.contents[i + 1 :]:
                    if "Rotational" in line_j:
                        return (
                            float(line_j.split()[1]) * kcal_per_mol_to_hartree
                        )
        return None

    @cached_property
    def thermal_translation_correction(self):
        """
        Thermal translation correction in Hartree.
        """
        for i, line_i in enumerate(self.contents):
            if "E (Thermal)" in line_i and "CV" in line_i:
                for line_j in self.contents[i + 1 :]:
                    if "Translational" in line_j:
                        return (
                            float(line_j.split()[1]) * kcal_per_mol_to_hartree
                        )
        return None

    @cached_property
    def thermal_energy_correction(self):
        """
        thermal correction to energy in Hartree.
        """
        for line in self.contents:
            if "Thermal correction to Energy=" in line:
                return float(line.split()[-1])
        return None

    @cached_property
    def thermal_enthalpy_correction(self):
        """
        thermal correction to enthalpy in Hartree.
        """
        for line in self.contents:
            if "Thermal correction to Enthalpy=" in line:
                return float(line.split()[-1])
        return None

    @cached_property
    def thermal_gibbs_free_energy_correction(self):
        """
        thermal correction to Gibbs free energy in Hartree.
        """
        for line in self.contents:
            if "Thermal correction to Gibbs Free Energy=" in line:
                return float(line.split()[-1])
        return None

    @cached_property
    def internal_energy(self):
        """
        Sum of electronic and thermal energies in Hartree.
        """
        for line in self.contents:
            if "Sum of electronic and thermal Energies=" in line:
                return float(line.split()[-1])
        return None

    @cached_property
    def enthalpy(self):
        """
        Sum of electronic and thermal enthalpies in Hartree.
        """
        for line in self.contents:
            if "Sum of electronic and thermal Enthalpies=" in line:
                return float(line.split()[-1])
        return None

    @cached_property
    def electronic_entropy_no_temperature_in_SI(self):
        """
        Electronic entropy in J/mol/K.
        """
        for i, line_i in enumerate(self.contents):
            if "E (Thermal)" in line_i and "CV" in line_i:
                for line_j in self.contents[i + 1 :]:
                    if "Electronic" in line_j:
                        return float(line_j.split()[-1]) * cal_to_joules
        return None

    @cached_property
    def electronic_entropy(self):
        """
        Electronic entropy in Hartree.
        """
        if self.electronic_entropy_no_temperature_in_SI is not None:
            electronic_entropy_hartree = (
                self.electronic_entropy_no_temperature_in_SI
                * joule_per_mol_to_hartree
            )
            return electronic_entropy_hartree
        return None

    @cached_property
    def vibrational_entropy_no_temperature_in_SI(self):
        """
        Vibrational entropy in J/mol/K.
        """
        for i, line_i in enumerate(self.contents):
            if "E (Thermal)" in line_i and "CV" in line_i:
                for line_j in self.contents[i + 1 :]:
                    if "Vibrational" in line_j:
                        return float(line_j.split()[-1]) * cal_to_joules
        return None

    @cached_property
    def vibrational_entropy(self):
        """
        Vibrational entropy in Hartree.
        """
        if self.vibrational_entropy_no_temperature_in_SI is not None:
            vibrational_entropy_hartree = (
                self.vibrational_entropy_no_temperature_in_SI
                * joule_per_mol_to_hartree
            )
            return vibrational_entropy_hartree
        return None

    @cached_property
    def rotational_entropy_no_temperature_in_SI(self):
        """
        Rotational entropy in J/mol/K.
        """
        for i, line_i in enumerate(self.contents):
            if "E (Thermal)" in line_i and "CV" in line_i:
                for line_j in self.contents[i + 1 :]:
                    if "Rotational" in line_j:
                        return float(line_j.split()[-1]) * cal_to_joules
        return None

    @cached_property
    def rotational_entropy(self):
        """
        Rotational entropy in Hartree.
        """
        if self.rotational_entropy_no_temperature_in_SI is not None:
            rotational_entropy_hartree = (
                self.rotational_entropy_no_temperature_in_SI
                * joule_per_mol_to_hartree
            )
            return rotational_entropy_hartree
        return None

    @cached_property
    def translational_entropy_no_temperature_in_SI(self):
        """
        Translational entropy in J/mol/K.
        """
        for i, line_i in enumerate(self.contents):
            if "E (Thermal)" in line_i and "CV" in line_i:
                for line_j in self.contents[i + 1 :]:
                    if "Translational" in line_j:
                        return float(line_j.split()[-1]) * cal_to_joules
        return None

    @cached_property
    def translational_entropy(self):
        """
        Translational entropy in Hartree.
        """
        if self.translational_entropy_no_temperature_in_SI is not None:
            translational_entropy_hartree = (
                self.translational_entropy_no_temperature_in_SI
                * joule_per_mol_to_hartree
            )
            return translational_entropy_hartree
        return None

    @cached_property
    def entropy_in_J_per_mol_per_K(self):
        """
        Total entropy in J/mol/K.
        """
        for i, line_i in enumerate(self.contents):
            if "E (Thermal)" in line_i and "CV" in line_i:
                for line_j in self.contents[i + 1 :]:
                    if "Total" in line_j:
                        return float(line_j.split()[-1]) * cal_to_joules
        return None

    @cached_property
    def entropy(self):
        """
        Total entropy in Hartree.
        """
        if self.entropy_in_J_per_mol_per_K is not None:
            total_entropy_hartree = (
                self.entropy_in_J_per_mol_per_K * joule_per_mol_to_hartree
            )
            return total_entropy_hartree
        return None

    @cached_property
    def entropy_times_temperature(self):
        """
        The entropy contributions are T*S = T*(S(el)+S(vib)+S(rot)+S(trans)).
        Return value in Hartree.
        """
        if (
            self.temperature_in_K
            and self.entropy_in_J_per_mol_per_K is not None
        ):
            entropy_ts_hartree = (
                self.entropy_in_J_per_mol_per_K
                * joule_per_mol_to_hartree
                * self.temperature_in_K
            )
            return entropy_ts_hartree
        return None

    @cached_property
    def gibbs_free_energy(self):
        """
        Sum of electronic and thermal free energies in Hartree.
        """
        for line in self.contents:
            if "Sum of electronic and thermal Free Energies=" in line:
                return float(line.split()[-1])
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
    def has_dipole_moment(self):
        """Check if the output file contains dipole moment calculations."""
        for line in self.contents:
            if "Dipole moment (field-independent basis, Debye):" in line:
                return True
        return False

    @cached_property
    def all_dipole_moments(self):
        """Obtain all dipole moments from the output file as [X, Y, Z] arrays in Debye."""
        list_of_all_dipole_moments, _ = (
            self._get_dipole_moments_for_molecules()
        )
        return list_of_all_dipole_moments

    @cached_property
    def all_dipole_moment_magnitudes(self):
        """Obtain all dipole moment magnitudes (total) from the output file in Debye."""
        _, list_of_all_dipole_moment_magnitudes = (
            self._get_dipole_moments_for_molecules()
        )
        return list_of_all_dipole_moment_magnitudes

    def _get_dipole_moments_for_molecules(self):
        """
        Obtain a list of dipole moments for all molecular geometries.
        Each dipole moment is stored as a np array of shape (3,) containing
        [X, Y, Z] components in Debye units.
        """
        all_dipole_moments = []
        all_dipole_moment_magnitudes = []
        for i, line in enumerate(self.contents):
            if "Dipole moment (field-independent basis, Debye):" in line:
                next_line = self.contents[i + 1]
                parts = next_line.split()
                x_idx = parts.index("X=") + 1
                y_idx = parts.index("Y=") + 1
                z_idx = parts.index("Z=") + 1
                tot_idx = parts.index("Tot=") + 1
                x_val = float(parts[x_idx])
                y_val = float(parts[y_idx])
                z_val = float(parts[z_idx])
                tot_val = float(parts[tot_idx])
                dipole_moment = np.array([x_val, y_val, z_val])
                all_dipole_moment_magnitudes.append(tot_val)
                all_dipole_moments.append(dipole_moment)
        if len(all_dipole_moments) == 0:
            return [], []
        return all_dipole_moments, all_dipole_moment_magnitudes

    @cached_property
    def num_dipole_moments(self):
        return len(self.all_dipole_moments)

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

    @cached_property
    def temperature_in_K(self):
        for line in self.contents:
            if "Temperature" in line and "Kelvin." in line:
                return float(line.split()[1])

    @cached_property
    def pressure_in_atm(self):
        for line in self.contents:
            if "Pressure" in line and "Atm." in line:
                return float(line.split()[-2])

    @property
    def mass(self):
        for line in self.contents:
            if "Molecular mass:" in line and "amu." in line:
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

    @cached_property
    def all_rotational_constants(self):
        """
        List of rotational constants (np.array in Hz) for each geometry step,
        in the order they appear in the file.
        """
        result = []
        for line in self.contents:
            if "Rotational constants (GHZ):" in line:
                vals = line.split("(GHZ):")[-1].split()
                result.append(np.array([float(v) * 1e9 for v in vals]))
        return result

    @cached_property
    def all_point_groups(self):
        """
        List of point group strings for each geometry step,
        in the order they appear in the file.
        """
        result = []
        for line in self.contents:
            if line.strip().startswith("Full point group"):
                parts = line.split()
                # format: Full point group  <PG>  NOp  <n>
                pg_idx = parts.index("group") + 1
                result.append(parts[pg_idx].upper())
        return result

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
        entropy_method="grimme",
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
            entropy_method (str): Entropy quasi-RRHO method ('grimme' or
                'truhlar'). Default 'grimme'.
            energy_units (str): Energy units for output values.
                Options: 'hartree', 'eV', 'kcal/mol', 'kJ/mol'. Default 'hartree'.
        """
        super().__init__(filename=filename)
        self.temperature = temperature
        self.concentration = concentration
        self.pressure = pressure
        self.cutoff_entropy_grimme = cutoff_entropy_grimme
        self.cutoff_enthalpy = cutoff_enthalpy
        self.entropy_method = entropy_method
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
                entropy_method=self.entropy_method,
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

    # =========================================================================
    # Multi-file pKa thermochemistry support
    # =========================================================================

    @staticmethod
    def compute_pka_thermochemistry(
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
        """Compute thermochemistry for pKa species (HA, A-, HRef, Ref-)."""
        from chemsmart.cli.pka import compute_pka_thermochemistry

        return compute_pka_thermochemistry(
            ha_file=ha_file,
            a_file=a_file,
            href_file=href_file,
            ref_file=ref_file,
            temperature=temperature,
            concentration=concentration,
            pressure=pressure,
            cutoff_entropy_grimme=cutoff_entropy_grimme,
            cutoff_enthalpy=cutoff_enthalpy,
            energy_units=energy_units,
        )

    @staticmethod
    def compute_pka(
        ha_gas_file,
        a_gas_file,
        href_gas_file=None,
        ref_gas_file=None,
        ha_solv_file=None,
        a_solv_file=None,
        href_solv_file=None,
        ref_solv_file=None,
        pka_reference=None,
        temperature=298.15,
        concentration=1.0,
        pressure=1.0,
        cutoff_entropy_grimme=100.0,
        cutoff_enthalpy=100.0,
        entropy_method="grimme",
        scheme="proton exchange",
        delta_G_proton=None,
    ):
        """
        Compute pKa using a dual-level thermodynamic cycle.

        **Proton exchange** (default): HA + Ref⁻ → A⁻ + HRef
            pKa = pKa_ref + ΔG_soln / (RT × ln10)

        **Direct dissociation** (``scheme='direct'``): HA → A⁻ + H⁺
            ΔG_diss = G_soln(A⁻) + G_soln(H⁺) - G_soln(HA)
            pKa = ΔG_diss / (2.303 × R × T)
        """
        from chemsmart.cli.pka import compute_pka

        return compute_pka(
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
            entropy_method=entropy_method,
            scheme=scheme,
            delta_G_proton=delta_G_proton,
        )

    @staticmethod
    def print_pka_summary(
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
        entropy_method: str = "grimme",
        scheme="proton exchange",
        delta_G_proton=None,
    ):
        """Print a formatted summary of a dual-level pKa calculation."""
        from chemsmart.cli.pka import print_pka_summary as _print_pka_summary

        return _print_pka_summary(
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
            entropy_method=entropy_method,
            scheme=scheme,
            delta_G_proton=delta_G_proton,
        )
