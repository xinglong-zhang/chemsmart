import numpy as np
from chemsmart.utils.utils import is_float
from chemsmart.io.gaussian.input import Gaussian16Input
from chemsmart.io.molecules.structure import CoordinateBlock
from chemsmart.utils.utils import content_blocks_by_paragraph


class TestUtils:
    def test_is_float(self):
        assert is_float("-1.0")
        assert is_float("1.9")
        assert is_float("-0.1")
        assert not is_float("-1")
        assert not is_float("1")

    def test_content_blocking(self, gaussian_opt_inputfile):
        g16_input = Gaussian16Input(filename=gaussian_opt_inputfile)
        content_blocks = content_blocks_by_paragraph(g16_input.contents)
        assert len(content_blocks) == 3
        cb_string = "\n".join(content_blocks[2])
        cb = CoordinateBlock(coordinate_block=cb_string)
        assert cb.molecule.empirical_formula == "C7H5ClO"
        assert cb.molecule.translation_vectors is None
        assert all(
            np.isclose(
                cb.molecule.positions[0],
                [-0.544821, -1.169457, 0.000127],
                atol=1e-4,
            )
        )
