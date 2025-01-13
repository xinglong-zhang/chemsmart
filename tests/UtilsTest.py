import numpy as np
from chemsmart.utils.utils import is_float
from chemsmart.io.gaussian.input import Gaussian16Input
from chemsmart.io.molecules.structure import CoordinateBlock
from chemsmart.utils.utils import content_blocks_by_paragraph
from chemsmart.utils.utils import get_list_from_string_range


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

    def test_get_list_from_string_range(self):
        s1 = "1-3"
        s2 = "1,3"
        s3 = "1,2,3"
        s4 = "1-3,5"
        s5 = "1-3,5-7"
        s6 = "1-3,5-7,10"
        s7 = "1,2,3,5-7,10"
        s8 = "[1-3,28-31,34-41]"
        s9 = "1-3,28-31,34-41"
        assert get_list_from_string_range(s1) == [1, 2, 3]
        assert get_list_from_string_range(s2) == [1, 3]
        assert get_list_from_string_range(s3) == [1, 2, 3]
        assert get_list_from_string_range(s4) == [1, 2, 3, 5]
        assert get_list_from_string_range(s5) == [1, 2, 3, 5, 6, 7]
        assert get_list_from_string_range(s6) == [1, 2, 3, 5, 6, 7, 10]
        assert get_list_from_string_range(s7) == [1, 2, 3, 5, 6, 7, 10]
        assert get_list_from_string_range(s8) == [
            1,
            2,
            3,
            28,
            29,
            30,
            31,
            34,
            35,
            36,
            37,
            38,
            39,
            40,
            41,
        ]
        assert get_list_from_string_range(s9) == [
            1,
            2,
            3,
            28,
            29,
            30,
            31,
            34,
            35,
            36,
            37,
            38,
            39,
            40,
            41,
        ]
