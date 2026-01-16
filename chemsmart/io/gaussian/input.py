from functools import cached_property

from chemsmart.io.gaussian.gengenecp import GenGenECPSection
from chemsmart.utils.mixins import GaussianFileMixin
from chemsmart.utils.utils import content_blocks_by_paragraph


class Gaussian16Input(GaussianFileMixin):
    """
    Parser for Gaussian 16 input files.
    """

    def __init__(self, filename):
        """
        Initialize Gaussian input file parser.
        """
        self.filename = filename

    @property
    def num_content_blocks(self):
        """
        Get the number of content blocks in the input file.
        """
        return len(self.content_groups)

    @cached_property
    def content_groups(self):
        """
        Parse Gaussian input file into content groups.

        The input file is organized into blocks separated by blank lines:
        - content_groups[0]: Header information and route string
        - content_groups[1]: Job title/description
        - content_groups[2]: Charge/multiplicity and molecular coordinates
        - content_groups[3:]: Additional sections (modred, gen/genecp, solvents)

        Note: Structure may vary depending on calculation type and may need
        updates for job-specific sections.
        """
        return content_blocks_by_paragraph(string_list=self.contents)

    @property
    def coordinate_block(self):
        """
        Get the molecular coordinate block.
        """
        from chemsmart.io.molecules.structure import CoordinateBlock

        cb = CoordinateBlock(coordinate_block=self.content_groups[2])
        return cb

    @property
    def modredundant_group(self):
        """
        Get the modredundant coordinates group if present.
        """
        if (
            "modred" in self.route_string or "modred" in self.route_string
        ) and self.num_content_blocks > 3:
            # Handle case where input file has opt=modred in route
            # but no modred section at the end
            return self.content_groups[3]
        return None

    @property
    def is_pbc(self):
        """
        Check if the calculation uses periodic boundary conditions.
        """
        for line in self.contents:
            line_elements = line.split()
            if len(line_elements) == 4 and line_elements[0].upper() == "TV":
                return True
        return False

    @property
    def translation_vectors(self):
        """
        Get translation vectors for periodic calculations.
        """
        return self.coordinate_block.translation_vectors

    @property
    def mem(self):
        """
        Get memory allocation specification.
        """
        return self._get_mem()

    @property
    def nproc(self):
        """
        Get number of processors specification.
        """
        return self._get_nproc()

    @property
    def charge(self):
        """
        Get molecular charge.
        """
        charge, _ = self._get_charge_and_multiplicity()
        return charge

    @property
    def multiplicity(self):
        """
        Get spin multiplicity.
        """
        _, multiplicity = self._get_charge_and_multiplicity()
        return multiplicity

    @property
    def route_string(self):
        """
        Get the complete route string.
        """
        return self._get_route()

    @property
    def gen_genecp_group(self):
        """
        Block of strings in the input file specifying the gen/genecp group.
        """
        return self._get_gen_genecp_group()

    @property
    def genecp_section(self):
        """
        Get parsed Gen/GenECP section object.
        """
        if self.gen_genecp_group:
            return GenGenECPSection.from_genecp_group(
                genecp_group=self.gen_genecp_group
            )
        return None

    @property
    def light_elements(self):
        """
        Get light elements from Gen/GenECP section.
        """
        if self.genecp_section:
            return self.genecp_section.light_elements
        return None

    @property
    def light_elements_basis(self):
        """
        Get basis set for light elements.
        """
        if self.genecp_section:
            return self.genecp_section.light_elements_basis
        return None

    @property
    def heavy_elements(self):
        """
        Get heavy elements from Gen/GenECP section.
        """
        if self.genecp_section:
            return self.genecp_section.heavy_elements
        return None

    @property
    def heavy_elements_basis(self):
        """
        Get basis set for heavy elements.
        """
        if self.genecp_section:
            return self.genecp_section.heavy_elements_basis
        return None

    @property
    def custom_solvent(self):
        """
        Get the custom solvent specification string.
        """
        return self._get_custom_solvent_string()

    @property
    def custom_solvent_group(self):
        """
        Get the custom solvent group.
        """
        return self._get_custom_solvent_group()

    @property
    def molecule(self):
        """
        Get the molecular structure object with charge and multiplicity.
        """
        molecule = self.coordinate_block.molecule
        # Update the charge and multiplicity of the molecule
        molecule.charge = self.charge
        molecule.multiplicity = self.multiplicity
        return molecule

    @property
    def num_atoms(self):
        """
        Get the number of atoms in the molecule.
        """
        return self.molecule.num_atoms

    @property
    def constrained_atoms(self):
        """
        Get atoms with coordinate constraints.
        """
        return self.coordinate_block.constrained_atoms

    @constrained_atoms.setter
    def constrained_atoms(self, value):
        """
        Set atoms with coordinate constraints.
        """
        self.constrained_atoms = value

    def _get_route(self):
        """
        Extract and format route strings from the input file.

        Route strings may span multiple lines starting with '#' and are
        terminated by blank lines or non-# lines. All route content is
        converted to lowercase for consistency.
        """
        route_strings = []
        current_block = []

        for line in self.contents:
            # Blank line ends current block if any
            if not line:
                if current_block:
                    route_strings.append(" ".join(current_block).lower())
                    current_block = []
                continue

            # Line starts or continues a '#' block
            if line.startswith("#"):
                # remove the leading '#' then trim again
                current_block.append(line)
            else:
                # Non-blank, non-# line also ends a block
                if current_block:
                    route_strings.append(" ".join(current_block).lower())
                    current_block = []
                # Otherwise ignore non-# lines

        # Flush last block if file ends with route content
        if current_block:
            route_strings.append(" ".join(current_block).lower())

        if self.is_link and route_strings:
            return route_strings[-1]

        return " ".join(route_strings)

    def _get_charge_and_multiplicity(self):
        """
        Extract molecular charge and spin multiplicity.
        """
        for line in self.contents:
            line_elements = line.split()
            if (
                len(line_elements) == 2
                and line_elements[0].replace("-", "").isdigit()
                and line_elements[1].isdigit()
            ):
                charge = int(line_elements[0])
                multiplicity = int(line_elements[1])
                return charge, multiplicity

    def _get_gen_genecp_group(self):
        """
        Extract Gen/GenECP basis set specification groups.
        """
        # Semi-empirical calculations don't use explicit basis sets
        if self.basis is None:
            return None

        # Only applicable for Gen/GenECP calculations
        if "gen" not in self.basis:
            return None

        # Determine group position based on additional input features
        if (
            "modred" in self.route_string
            and "solvent=generic" in self.route_string
        ):
            # Both modred and custom solvent present
            return self.content_groups[4:-1]
        if (
            "modred" in self.route_string
            and "solvent=generic" not in self.route_string
        ):
            # Only modred present
            return self.content_groups[4:]
        if (
            "modred" not in self.route_string
            and "solvent=generic" in self.route_string
        ):
            # Only custom solvent present
            return self.content_groups[3:-1]
        if (
            "modred" not in self.route_string
            and "solvent=generic" not in self.route_string
        ):
            # Neither modred nor custom solvent
            return self.content_groups[3:]

        return None

    def _get_custom_solvent_group(self):
        """
        Extract the custom solvent specification group.

        Custom solvent specifications appear as the last content group
        when 'solvent=generic' is specified in the route section.
        """
        if "solvent=generic" in self.route_string:
            return self.content_groups[-1]
        return None

    def _get_custom_solvent_string(self):
        """
        Convert custom solvent group to formatted string.
        """
        if self.custom_solvent_group is not None:
            return "\n".join(self.custom_solvent_group)
        return None
