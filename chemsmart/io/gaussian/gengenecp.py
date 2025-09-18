import logging
import os

from chemsmart.utils.periodictable import PeriodicTable
from chemsmart.utils.utils import (
    content_blocks_by_paragraph,
    write_list_of_lists_as_a_string_with_empty_line_between_lists,
)

pt = PeriodicTable()

logger = logging.getLogger(__name__)


class GenGenECPSection:
    """
    Parser and generator for Gaussian Gen/GenECP basis set sections.
    """

    def __init__(self, string):
        """
        Initialize Gen/GenECP section parser.
        """
        self.string = string

    @property
    def string_list(self):
        """
        Get section content as list of lines.
        """
        return self.string.split("\n")

    @property
    def light_elements_basis(self):
        """
        Get basis set name for light elements.

        The second line typically contains the basis set name for light
        elements (those not using ECPs).
        """
        try:
            return self.string_list[1]
        except IndexError:
            return None

    @property
    def heavy_elements(self):
        """
        Extract heavy elements that use ECPs.
        """
        elements = []
        for line in self.string_list[2:]:
            # Line should be in the format: '{element} 0'
            split_line = line.split()
            if len(split_line) != 2:
                continue

            element = split_line[0]
            if element in pt.PERIODIC_TABLE and split_line[1] == "0":
                elements.append(element)
        return pt.sorted_periodic_table_list(elements)

    @property
    def light_elements(self):
        """
        Extract light elements that use standard basis sets.

        Light elements are listed in the first line, space-separated,
        ending with '0'.
        """
        # First line contains light atoms
        line = self.string_list[0]
        if not line:
            return []

        if not line.endswith("0"):
            logger.warning(
                f"Line for light atoms should end with 0, but is "
                f"{line} instead."
            )

        split_line = line.split()
        elements = []
        for raw_element in split_line:
            element = raw_element.strip()
            if (
                element != "0"
                and element in pt.PERIODIC_TABLE
                and element not in elements
            ):
                elements.append(element)
        return pt.sorted_periodic_table_list(elements)

    @property
    def heavy_elements_basis(self):
        """
        Extract basis set name for heavy elements from comments.
        """
        for line in self.string.split("\n"):
            if "!   Basis set" in line:
                # From API, basis set info is in header line like:
                # '!   Basis set: def2-TZVPPD'
                return line.split(":")[-1].strip().lower()
            if "Basis Set Exchange Library" in line:
                # Directly obtained from Basis Set Exchange website
                return line.split()[1].lower()
        return None

    @property
    def _string_blocks(self):
        """
        Parse content into paragraph blocks.
        """
        return content_blocks_by_paragraph(string_list=self.string_list)

    @property
    def genecp_type(self):
        """
        Determine the type of Gen/GenECP section.
        """
        genecp_string_blocks = self._string_blocks
        genecp_type = ""
        if len(genecp_string_blocks) == 1:
            genecp_type = "gen"
        elif len(genecp_string_blocks) == 2:
            genecp_type = "genecp"
        return genecp_type

    @classmethod
    def from_genecp_path(cls, genecp_path):
        """
        Create GenGenECPSection from a file path.
        """
        if not os.path.exists(genecp_path):
            raise FileNotFoundError(
                f'Given gen/genecp path at "{genecp_path}" is not found!'
            )

        genecp_path = os.path.abspath(genecp_path)
        with open(genecp_path) as f:
            string = ""
            lines = f.readlines()
            if lines[-1] == "\n":
                lines = lines[:-1]
            for line in lines:
                # Concatenate lines preserving newlines
                string += line
        return cls(string)

    @classmethod
    def from_comfile(cls, comfile):
        """
        Create GenGenECPSection from a Gaussian input file.
        """
        from chemsmart.io.gaussian.input import Gaussian16Input

        comfile = Gaussian16Input(filename=comfile)
        return cls.from_genecp_group(comfile.gen_genecp_group)

    @classmethod
    def from_genecp_group(cls, genecp_group):
        """
        Create GenGenECPSection from the genecp_group string.
        """
        genecp_string = ""
        num_groups = len(genecp_group)
        for i in range(num_groups):
            # Add each string in the group with newlines
            for string in genecp_group[i]:
                genecp_string += string + "\n"
            genecp_string += "\n"
        return cls(genecp_string)

    @classmethod
    def from_bse_api(
        cls,
        light_elements,
        light_elements_basis,
        heavy_elements,
        heavy_elements_basis,
    ):
        """
        Create Gen/GenECP section using Basis Set Exchange API.

        This method generates a custom basis set section by combining
        light elements (standard basis sets) with heavy elements
        (basis sets with ECPs) using the Basis Set Exchange library.

        Args:
            light_elements (list): Light element symbols (no ECPs)
            light_elements_basis (str): Basis set name for light elements
            heavy_elements (list): Heavy element symbols (with ECPs)
            heavy_elements_basis (str): Basis set name for heavy elements
        """
        try:
            import basis_set_exchange as bse

            from chemsmart.io.gaussian import BSEMetadata

            bse_all_bases = BSEMetadata().all_bases_names()
        except ImportError as e:
            raise ImportError(
                "basis_set_exchange module needed.\n"
                "See https://github.com/MolSSI-BSE/basis_set_exchange "
                "for installation."
            ) from e

        # Sort element lists according to periodic table order
        heavy_elements = pt.sorted_periodic_table_list(
            list_of_elements=heavy_elements
        )
        heavy_elements_basis = heavy_elements_basis.lower()

        genecp_string = ""

        # Write light atoms and their basis set specification
        if len(light_elements) > 0:
            light_elements = pt.sorted_periodic_table_list(light_elements)
            light_elements_basis = light_elements_basis.lower()
            light_atoms_string = " ".join(light_elements) + " 0\n"
            genecp_string += light_atoms_string

            # Handle def2- basis set naming convention
            if "def2-" in light_elements_basis:
                light_elements_basis = light_elements_basis.replace(
                    "def2-", "def2"
                )

            genecp_string += light_elements_basis + "\n"
            # separate light atoms basis from beginning of heavy atoms gen/genecp basis
            genecp_string += "****\n"
        # Generate heavy atom basis set content from BSE API
        # Handle def2 basis set naming convention
        if "def2" in heavy_elements_basis and "-" not in heavy_elements_basis:
            heavy_elements_basis = heavy_elements_basis.replace(
                "def2", "def2-"
            )

        assert heavy_elements_basis in bse_all_bases, (
            f"BSE basis for {heavy_elements} given is "
            f"{heavy_elements_basis}.\n"
            f"This is not in BSE available bases: {bse_all_bases}."
        )

        # Retrieve basis set data from BSE in Gaussian format
        heavy_atoms_gengenecp_basis = bse.get_basis(
            name=heavy_elements_basis,
            elements=heavy_elements,
            fmt="gaussian94",
            header=True,
        )

        # Parse the basis set data into blocks
        heavy_atoms_gengenecp_basis_list = heavy_atoms_gengenecp_basis.split(
            "\n"
        )
        heavy_atoms_gengenecp_basis_blocks = content_blocks_by_paragraph(
            string_list=heavy_atoms_gengenecp_basis_list
        )

        # Write header block containing basis set information
        header_block = heavy_atoms_gengenecp_basis_blocks[0]
        for line in header_block:
            if line:
                genecp_string += line + "\n"

        # Write the actual basis set data blocks
        heavy_atoms_gengenecp_basis_string = (
            write_list_of_lists_as_a_string_with_empty_line_between_lists(
                heavy_atoms_gengenecp_basis_blocks[1:]
            )
        )
        genecp_string += heavy_atoms_gengenecp_basis_string
        return cls(string=genecp_string)
