import logging
import os

from chemsmart.io.utils import (
    content_blocks_by_paragraph,
    write_list_of_lists_as_a_string_with_empty_line_between_lists,
)
from chemsmart.utils.periodictable import PeriodicTable as pt

logger = logging.getLogger(__name__)


class GenGenECPSection:
    def __init__(self, string):
        self.string = string

    @property
    def string_list(self):
        return self.string.split("\n")

    @property
    def light_elements_basis(self):
        """Second line is the light atoms basis."""
        try:
            return self.string_list[1]
        except IndexError:
            return None

    @property
    def heavy_elements(self):
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
        # First line for light atoms
        line = self.string_list[0]
        if not line:
            return []

        if not line.endswith("0"):
            logger.warning(
                f"Line for light atoms should end with 0, but is {line} instead."
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
        for line in self.string.split("\n"):
            if "!   Basis set" in line:
                # from api, the info of basis set is in header with line
                # '!   Basis set: def2-TZVPPD'
                return line.split(":")[-1].strip().lower()
            if "Basis Set Exchange Library" in line:
                # directly obtained from basissetexchange website
                return line.split()[1].lower()
        return None

    @property
    def _string_blocks(self):
        return content_blocks_by_paragraph(string_list=self.string_list)

    @property
    def genecp_type(self):
        genecp_string_blocks = self._string_blocks
        genecp_type = ""
        if len(genecp_string_blocks) == 1:
            genecp_type = "gen"
        elif len(genecp_string_blocks) == 2:
            genecp_type = "genecp"
        return genecp_type

    @classmethod
    def from_genecp_path(cls, genecp_path):
        if not os.path.exists(genecp_path):
            raise FileNotFoundError(
                f'Given gen/genecp path at "{genecp_path}" is not found!'
            )

        string = ""
        genecp_path = os.path.abspath(genecp_path)
        with open(genecp_path) as f:
            for line in f.readlines():
                string += line
        return cls(string)

    @classmethod
    def from_comfile(cls, comfile):
        from pyatoms.io.gaussian.inputs import Gaussian16Input

        comfile = Gaussian16Input(comfile=comfile)
        return cls.from_genecp_group(comfile.gen_genecp)

    @classmethod
    def from_genecp_group(cls, genecp_group):
        # genecp group from gen_genecp attribute of inputs.py
        genecp_string = ""
        num_groups = len(genecp_group)
        for i in range(num_groups):
            for string in genecp_group[
                i
            ]:  # for each string in the list; so need to add end of line '\n' below
                genecp_string += string + "\n"
            genecp_string += "\n"
        return cls(genecp_string)

    @classmethod
    def from_bse_api(
        cls, light_elements, light_elements_basis, heavy_elements, heavy_elements_basis
    ):
        """Create ECP from basis set exchange api.

        :param light_elements: list of light atoms as elements in string
        :param light_elements_basis: basis set for light atoms
        :param heavy_elements: list of heavy atoms as elements in string
        :param heavy_elements_basis: basis set for heavy atoms; to be obtained from api
        """
        try:
            import basis_set_exchange as bse

            from pyatoms.io.gaussian import BSEMetadata

            bse_all_bases = BSEMetadata().all_bases_names()
        except ImportError as e:
            raise ImportError(
                "basis_set_exchange module needed.\n"
                "see https://github.com/MolSSI-BSE/basis_set_exchange for installation."
            ) from e

        heavy_elements = pt.sorted_periodic_table_list(heavy_elements)
        heavy_elements_basis = heavy_elements_basis.lower()

        genecp_string = ""

        # write light atoms and light atoms basis
        # check if light_elements is not empty

        light_atoms_string = ""
        if len(light_elements) == 0:
            pass
        else:
            light_elements = pt.sorted_periodic_table_list(light_elements)
            light_elements_basis = light_elements_basis.lower()
            light_atoms_string = " ".join(light_elements) + " 0\n"
            genecp_string += light_atoms_string
            if "def2-" in light_elements_basis:
                light_elements_basis = light_elements_basis.replace("def2-", "def2")

            genecp_string += light_elements_basis + "\n"
            genecp_string += "****\n"  # separate light atoms basis from beginning of heavy atoms gen/genecp basis

        # write heavy atom basis (from api)
        if "def2" in heavy_elements_basis and "-" not in heavy_elements_basis:
            heavy_elements_basis = heavy_elements_basis.replace("def2", "def2-")

        assert heavy_elements_basis in bse_all_bases, (
            f"BSE basis for {heavy_elements} given is {heavy_elements_basis}.\n"
            f"This is not in BSE available bases: {bse_all_bases}. "
        )

        heavy_atoms_gengenecp_basis = bse.get_basis(
            name=heavy_elements_basis,
            elements=heavy_elements,
            fmt="gaussian94",
            header=True,
        )

        heavy_atoms_gengenecp_basis_list = heavy_atoms_gengenecp_basis.split("\n")
        heavy_atoms_gengenecp_basis_blocks = content_blocks_by_paragraph(
            string_list=heavy_atoms_gengenecp_basis_list
        )

        # first block is header; write header info, which contains basis set name
        header_block = heavy_atoms_gengenecp_basis_blocks[0]
        for line in header_block:
            if line:
                genecp_string += line + "\n"

        heavy_atoms_gengenecp_basis_string = (
            write_list_of_lists_as_a_string_with_empty_line_between_lists(
                heavy_atoms_gengenecp_basis_blocks[1:]
            )
        )
        genecp_string += heavy_atoms_gengenecp_basis_string
        return cls(string=genecp_string)
