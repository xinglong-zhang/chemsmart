import numpy as np
import pytest

from chemsmart.io.gaussian.input import Gaussian16Input
from chemsmart.io.molecules.structure import CoordinateBlock
from chemsmart.utils.utils import (
    cmp_with_ignore,
    content_blocks_by_paragraph,
    get_list_from_string_range,
    get_range_from_list,
    is_float,
    str_indices_range_to_list,
    string2index_1based,
)


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

    def test_cmp_with_ignore_string(
        self,
        gaussian_written_sp_from_nhc_singlet_log_with_custom_basis_from_api_file,
        gaussian_written_sp_from_nhc_singlet_log_with_custom_basis_from_api_file_v2,
    ):
        assert cmp_with_ignore(
            gaussian_written_sp_from_nhc_singlet_log_with_custom_basis_from_api_file,
            gaussian_written_sp_from_nhc_singlet_log_with_custom_basis_from_api_file_v2,
            ignore_string="Version",
        )

    def test_cmp_with_ignore_list(
        self,
        gaussian_written_opt_file,
        gaussian_written_opt_file_with_route,
    ):
        assert cmp_with_ignore(
            gaussian_written_opt_file,
            gaussian_written_opt_file_with_route,
            ignore_string=["#", "job"],
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

    def test_get_range_from_list(self):
        s1 = [1, 2, 3, 5, 6, 7]
        range = get_range_from_list(s1)
        assert range == ["1-3", "5-7"]


class TestGetListFromStringRange:
    def test_get_list_from_string_range(self):
        s1 = "[1-3,28-31,34-41]"
        s1_list = get_list_from_string_range(string_of_range=s1)
        assert s1_list == [
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

        s2 = "1-3,28-31,34-41"
        s2_list = get_list_from_string_range(string_of_range=s2)
        assert s2_list == [
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

        s3 = "1,3,33,37,42,43,44,45"
        s3_list = get_list_from_string_range(string_of_range=s3)
        assert s3_list == [1, 3, 33, 37, 42, 43, 44, 45]

    def test_get_list_from_string(self):
        s1 = "1:9"
        s1_list = str_indices_range_to_list(str_indices=s1)
        assert s1_list == [1, 2, 3, 4, 5, 6, 7, 8]

        s2 = "1,2,4"
        s2_list = str_indices_range_to_list(str_indices=s2)
        assert s2_list == [1, 2, 4]

        s3 = "1-9"
        s3_list = str_indices_range_to_list(str_indices=s3)
        assert s3_list == [1, 2, 3, 4, 5, 6, 7, 8]

        s4 = "[1-9]"
        s4_list = str_indices_range_to_list(str_indices=s4)
        assert s4_list == [1, 2, 3, 4, 5, 6, 7, 8]

        s6 = "2:3"
        s6_list = str_indices_range_to_list(str_indices=s6)
        assert s6_list == [2]


class TestString2Index1Based:
    def test_single_integer(self):
        assert string2index_1based("1") == 0  # 1-based -> 0-based
        assert string2index_1based("5") == 4  # 1-based -> 0-based
        assert string2index_1based("10") == 9  # 1-based -> 0-based

    def test_slice(self):
        result = string2index_1based("1:5")
        assert isinstance(result, slice)
        list1 = list(range(5))
        assert list1[result] == [0, 1, 2, 3]
        assert result.start == 0  # 1-based start -> 0-based
        assert result.stop == 4  # 1-based stop remains same
        assert result.step is None

        result = string2index_1based("3:10")
        assert isinstance(result, slice)
        assert result.start == 2  # 1-based start -> 0-based
        assert result.stop == 9
        assert result.step is None

    def test_slice_with_step(self):
        result = string2index_1based("1:10:2")
        assert isinstance(result, slice)
        assert result.start == 0  # 1-based -> 0-based
        assert result.stop == 9
        assert result.step == 2

        result = string2index_1based("2:8:3")
        assert isinstance(result, slice)
        assert result.start == 1  # 1-based -> 0-based
        assert result.stop == 7
        assert result.step == 3

    def test_open_ended_slice(self):
        result = string2index_1based("5:")
        assert isinstance(result, slice)
        assert result.start == 4  # 1-based -> 0-based
        assert result.stop is None
        assert result.step is None

        result = string2index_1based(":5")
        assert isinstance(result, slice)
        assert result.start is None
        assert result.stop == 4  # 1-based stop remains same
        assert result.step is None

        result = string2index_1based(":")
        assert isinstance(result, slice)
        assert result.start is None
        assert result.stop is None
        assert result.step is None

    def test_invalid_inputs(self):
        # Invalid integer
        with pytest.raises(ValueError):
            string2index_1based("invalid")

        # Slice with non-integer values
        with pytest.raises(ValueError):
            string2index_1based("a:b")

        # Mixed invalid formats
        with pytest.raises(ValueError):
            string2index_1based("1:x:2")
