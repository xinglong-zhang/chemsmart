from functools import cached_property
from chemsmart.utils.mixins import GaussianFileMixin
from chemsmart.utils.utils import content_blocks_by_paragraph
from chemsmart.io.gaussian.gengenecp import GenGenECPSection


class XTBInput(XTBFileMixin):
    def __init__(self, filename):
        self.filename = filename

    @property
    def num_content_blocks(self):
        return len(self.content_groups)

    @cached_property
    def content_groups(self):
        """XTB input file with content groups.
        # content_groups[0] gives the number of atoms
        # content_groups[1] gives the xyz coordinates
        """
        return content_blocks_by_paragraph(string_list=self.contents)

    @property
    def natoms(self):
        return self.molecule.natoms

    @property
    def coordinate_block(self):
        from chemsmart.io.molecules.structure import CoordinateBlock

        cb = CoordinateBlock(coordinate_block=self.content_groups[1])
        return cb

    @property
    def molecule(self):
        return self.coordinate_block.molecule