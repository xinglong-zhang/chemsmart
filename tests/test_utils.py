import subprocess

import numpy as np
import pytest

from chemsmart.io.gaussian.input import Gaussian16Input
from chemsmart.io.molecules.structure import CoordinateBlock, Molecule
from chemsmart.utils.io import (
    clean_duplicate_structure,
    clean_label,
    convert_string_indices_to_pymol_id_indices,
    create_molecule_list,
    line_of_all_integers,
    line_of_integer_followed_by_floats,
)
from chemsmart.utils.utils import (
    cmp_with_ignore,
    content_blocks_by_paragraph,
    convert_string_index_from_1_based_to_0_based,
    get_list_from_string_range,
    get_range_from_list,
    is_float,
    iterative_compare,
    naturally_sorted,
    return_objects_and_indices_from_string_index,
    run_command,
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
        assert not is_float("abc")

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

    def test_get_indices_from_string(self):
        """Test the conversion of string indices
        to a list of integers; 1-based indices."""
        objects = ["a", "b", "c", "d", "e", "f", "g", "h"]
        s1 = "1:3"  # standard python slicing
        s2 = "1,2,4"
        s3 = "1-3"  # user-defined slicing, 1-3 inclusive
        s4 = "[1-3]"  # user-defined slicing, 1-3 inclusive
        s5 = "2:3"  # standard python slicing
        s6 = "1"  # single string index
        s7 = "-1"  # single python last index
        s8 = "0"  # this will raise an error, as 1-based indices are expected
        assert objects[convert_string_index_from_1_based_to_0_based(s1)] == [
            "a",
            "b",
        ]
        assert [
            objects[i]
            for i in convert_string_index_from_1_based_to_0_based(s2)
        ] == ["a", "b", "d"]
        assert [
            objects[i]
            for i in convert_string_index_from_1_based_to_0_based(s3)
        ] == ["a", "b", "c"]
        assert [
            objects[i]
            for i in convert_string_index_from_1_based_to_0_based(s4)
        ] == ["a", "b", "c"]
        assert objects[convert_string_index_from_1_based_to_0_based(s5)] == [
            "b"
        ]
        assert [objects[convert_string_index_from_1_based_to_0_based(s6)]] == [
            "a"
        ]
        assert [objects[convert_string_index_from_1_based_to_0_based(s7)]] == [
            "h"
        ]
        with pytest.raises(ValueError):
            convert_string_index_from_1_based_to_0_based(s8)

    def test_iterative_compare_list_of_elements(self):
        list1 = [1, 2, 3, 4, 5]
        unique_list1 = iterative_compare(list1)
        assert unique_list1 == list1

        list2 = [1, 2, 3, 4, 5, 1, 2, 3, 4, 5]
        unique_list2 = iterative_compare(list2)
        assert unique_list2 == list1

        list3 = [1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5]
        unique_list3 = iterative_compare(list3)
        assert unique_list3 == list1

    def test_iterative_compare_list_of_lists(self):
        list1 = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        unique_list1 = iterative_compare(list1)
        assert unique_list1 == list1

        list2 = [
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9],
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9],
        ]
        unique_list2 = iterative_compare(list2)
        assert unique_list2 == list1

        list3 = [
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9],
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9],
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9],
        ]
        unique_list3 = iterative_compare(list3)
        assert unique_list3 == list1

    def test_iterative_compare_list_of_tuples(self):
        list1 = [(1, 2, 3), (4, 5, 6), (7, 8, 9)]
        unique_list1 = iterative_compare(list1)
        assert unique_list1 == list1

        list2 = [
            (1, 2, 3),
            (4, 5, 6),
            (7, 8, 9),
            (1, 2, 3),
            (4, 5, 6),
            (7, 8, 9),
        ]
        unique_list2 = iterative_compare(list2)
        assert unique_list2 == list1

        list3 = [
            (1, 2, 3),
            (4, 5, 6),
            (7, 8, 9),
            (1, 2, 3),
            (4, 5, 6),
            (7, 8, 9),
            (1, 2, 3),
            (4, 5, 6),
            (7, 8, 9),
        ]
        unique_list3 = iterative_compare(list3)
        assert unique_list3 == list1

    def test_iterative_compare_list_of_string(self):
        list1 = ["a", "b", "c", "d", "e"]
        unique_list1 = iterative_compare(list1)
        assert unique_list1 == list1

        list2 = ["a", "b", "c", "d", "e", "a", "b", "c", "d", "e"]
        unique_list2 = iterative_compare(list2)
        assert unique_list2 == list1

        list3 = [
            "a",
            "b",
            "c",
            "d",
            "e",
            "a",
            "b",
            "c",
            "d",
            "e",
            "a",
            "b",
            "c",
            "d",
            "e",
        ]
        unique_list3 = iterative_compare(list3)
        assert unique_list3 == list1

    def test_iterative_compare_list_of_dicts(self):
        dict1 = {"a": 1, "b": 2, "c": 3}
        dict2 = {"d": 4, "e": 5, "f": 6}
        dict3 = {"g": 7, "h": 8, "i": 9}
        list1 = [dict1, dict2, dict3]
        unique_list1 = iterative_compare(list1)
        assert unique_list1 == list1

        dict4 = {"a": 1, "b": 2, "c": 3}
        dict5 = {"d": 4, "e": 5, "f": 6}
        dict6 = {"g": 7, "h": 8, "i": 9}
        list2 = [dict4, dict5, dict6, dict1, dict2, dict3]
        unique_list2 = iterative_compare(list2)
        assert unique_list2 == list1

        dict7 = {"a": 1, "b": 2, "c": 3}
        dict8 = {"d": 4, "e": 5, "f": 6}
        dict9 = {"g": 7, "h": 8, "i": 9}
        list3 = [dict7, dict8, dict9, dict1, dict2, dict3, dict4, dict5, dict6]
        unique_list3 = iterative_compare(list3)
        assert unique_list3 == list1

        dict11 = {"a": 11, "b": 12, "c": 13}
        dict12 = {"d": 14, "e": 15, "f": 16}
        dict13 = {"g": 17, "h": 18, "i": 19}
        list4 = [dict11, dict12, dict13, dict1, dict2, dict3]
        unique_list4 = iterative_compare(list4)
        assert len(unique_list4) == 6

        list5 = [dict1, dict11]
        unique_list5 = iterative_compare(list5)
        assert len(unique_list5) == 2

        list6 = [dict1, dict11, dict1]
        unique_list6 = iterative_compare(list6)
        assert len(unique_list6) == 2

    def test_get_range_from_list(self):
        s1 = [1, 2, 3, 5, 6, 7]
        range = get_range_from_list(s1)
        assert range == ["1-3", "5-7"]

        s2 = [1, 34, 45, 46, 48, 50]
        range = get_range_from_list(s2)
        assert range == ["1", "34", "45-46", "48", "50"]

        s3 = [28, 45, 60, 89]
        range = get_range_from_list(s3)
        assert range == ["28", "45", "60", "89"]

        s4 = [
            18,
            19,
            20,
            21,
            23,
            25,
            27,
            29,
            30,
            31,
            33,
            35,
            37,
            39,
            41,
            43,
            46,
            47,
            48,
            49,
            51,
            53,
            55,
            57,
            61,
            62,
            63,
            64,
            66,
            68,
            70,
            72,
        ]
        range = get_range_from_list(s4)
        assert range == [
            "18-21",
            "23",
            "25",
            "27",
            "29-31",
            "33",
            "35",
            "37",
            "39",
            "41",
            "43",
            "46-49",
            "51",
            "53",
            "55",
            "57",
            "61-64",
            "66",
            "68",
            "70",
            "72",
        ]


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
        assert s3_list == [1, 2, 3, 4, 5, 6, 7, 8, 9]

        s4 = "[1-9]"
        s4_list = str_indices_range_to_list(str_indices=s4)
        assert s4_list == [1, 2, 3, 4, 5, 6, 7, 8, 9]

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


class TestParseIndexSpecification:
    """Tests for the new unified parse_index_specification function."""

    def test_ase_style_single_indices(self):
        """Test ASE-style single index specifications."""
        from chemsmart.utils.utils import parse_index_specification

        assert parse_index_specification("1") == 0
        assert parse_index_specification("5") == 4
        assert parse_index_specification("-1") == -1
        assert parse_index_specification("-2") == -2

    def test_ase_style_slices(self):
        """Test ASE-style slice specifications."""
        from chemsmart.utils.utils import parse_index_specification

        # Basic slices
        result = parse_index_specification("1:5")
        assert isinstance(result, slice)
        assert result == slice(0, 4)

        # All items
        result = parse_index_specification(":")
        assert isinstance(result, slice)
        assert result == slice(None, None)

        # Open-ended slices
        result = parse_index_specification("5:")
        assert result == slice(4, None)

        result = parse_index_specification(":5")
        assert result == slice(None, 4)

        # With step
        result = parse_index_specification("::2")
        assert result == slice(None, None, 2)

        result = parse_index_specification("1:10:2")
        assert result == slice(0, 9, 2)

    def test_free_format_comma_separated(self):
        """Test free-format comma-separated specifications."""
        from chemsmart.utils.utils import parse_index_specification

        assert parse_index_specification("1,3,5") == [0, 2, 4]
        assert parse_index_specification("1,2,4") == [0, 1, 3]

    def test_free_format_with_negative_indices(self):
        """Test free-format with negative indices."""
        from chemsmart.utils.utils import parse_index_specification

        assert parse_index_specification("1,-1") == [0, -1]
        assert parse_index_specification("1,3,-1") == [0, 2, -1]
        assert parse_index_specification("2,-2") == [1, -2]
        assert parse_index_specification("-1,-2") == [-1, -2]

    def test_free_format_hyphen_ranges(self):
        """Test free-format hyphen-based range specifications."""
        from chemsmart.utils.utils import parse_index_specification

        # Simple range (inclusive)
        assert parse_index_specification("1-5") == [0, 1, 2, 3, 4]
        assert parse_index_specification("2-4") == [1, 2, 3]

        # Range with brackets
        assert parse_index_specification("[1-5]") == [0, 1, 2, 3, 4]

    def test_free_format_mixed(self):
        """Test free-format mixed specifications."""
        from chemsmart.utils.utils import parse_index_specification

        # Mix of ranges and individual indices
        assert parse_index_specification("1-3,5") == [0, 1, 2, 4]
        assert parse_index_specification("1-3,5,7-9") == [0, 1, 2, 4, 6, 7, 8]

        # Mix with negative indices
        assert parse_index_specification("1-2,-1") == [0, 1, -1]
        assert parse_index_specification("1,3-5,-1") == [0, 2, 3, 4, -1]

    def test_invalid_inputs(self):
        """Test that invalid inputs raise ValueError."""
        from chemsmart.utils.utils import parse_index_specification

        # Index 0 is not allowed (1-based indexing)
        with pytest.raises(ValueError):
            parse_index_specification("0")

        with pytest.raises(ValueError):
            parse_index_specification("1,0,3")

    def test_with_actual_lists(self):
        """Test parse_index_specification with actual list indexing."""
        from chemsmart.utils.utils import parse_index_specification

        objects = ["a", "b", "c", "d", "e", "f", "g", "h"]

        # Single index
        idx = parse_index_specification("1")
        assert objects[idx] == "a"

        # Negative index
        idx = parse_index_specification("-1")
        assert objects[idx] == "h"

        # Slice
        idx = parse_index_specification("1:4")
        assert objects[idx] == ["a", "b", "c"]

        idx = parse_index_specification("1:7:2")
        assert objects[idx] == ["a", "c", "e"]

        idx = parse_index_specification("1:8:2")
        assert objects[idx] == ["a", "c", "e", "g"]

        idx = parse_index_specification("::2")
        assert objects[idx] == ["a", "c", "e", "g"]

        # All
        idx = parse_index_specification(":")
        assert objects[idx] == objects

        # Comma-separated
        idx = parse_index_specification("1,3,5")
        assert [objects[i] for i in idx] == ["a", "c", "e"]

        # Range
        idx = parse_index_specification("1-3")
        assert [objects[i] for i in idx] == ["a", "b", "c"]

        # Mixed with negative
        idx = parse_index_specification("1,-1")
        assert [objects[i] for i in idx] == ["a", "h"]

        idx = parse_index_specification("1,3,-1")
        assert [objects[i] for i in idx] == ["a", "c", "h"]

    def test_duplicate_detection_enabled(self):
        """Test duplicate detection when allow_duplicates=False."""
        from chemsmart.utils.utils import parse_index_specification

        # Test explicit duplicates should fail
        with pytest.raises(ValueError, match="Index overlap detected"):
            parse_index_specification(
                "5,-1", total_count=5, allow_duplicates=False
            )

    def test_boundary_checking_enabled(self):
        """Test boundary checking when allow_out_of_range=False."""
        from chemsmart.utils.utils import parse_index_specification

        # Test out of range should fail (10
        # structures, index 11 is out of range)
        with pytest.raises(
            ValueError, match="Index 11 is out of range.*10 structures"
        ):
            parse_index_specification(
                "11", total_count=10, allow_out_of_range=False
            )

        # Test negative out of range should fail (-11 with 10 structures)
        with pytest.raises(
            ValueError,
            match="Negative index -11 is out of range.*10 structures",
        ):
            parse_index_specification(
                "-11", total_count=10, allow_out_of_range=False
            )

        # Test range extending beyond bounds (8-11 with 10 structures)
        with pytest.raises(
            ValueError, match="Index 11 is out of range.*10 structures"
        ):
            parse_index_specification(
                "8-11", total_count=10, allow_out_of_range=False
            )

    def test_parse_index_duplicate_detection_disabled(self):
        """Test duplicate detection when allow_duplicates=False."""
        from chemsmart.utils.utils import parse_index_specification

        # Test explicit duplicates should fail
        with pytest.raises(ValueError, match="Index overlap detected"):
            parse_index_specification(
                "1,4,-2", total_count=5, allow_duplicates=False
            )

        # Test negative and positive indices pointing to same structure
        with pytest.raises(ValueError, match="Index overlap detected"):
            parse_index_specification(
                "1,-5", total_count=5, allow_duplicates=False
            )

    def test_parse_index_duplicate_detection_enabled(self):
        """Test duplicate detection when allow_duplicates=True."""
        from chemsmart.utils.utils import parse_index_specification

        # Test duplicates are allowed - should return all indices normalized
        result = parse_index_specification(
            "1,4,-2", total_count=5, allow_duplicates=True
        )
        # After normalization: [1-1=0, 4-1=3, 5+(-2)=3]
        assert result == [0, 3, 3]

        # Test negative and positive indices
        # pointing to same structure are allowed
        result = parse_index_specification(
            "1,-5", total_count=5, allow_duplicates=True
        )
        # After normalization: [1-1=0, 5+(-5)=0]
        assert result == [0, 0]

    def test_parse_index_boundary_detection_disabled(self):
        """Test boundary detection when allow_out_of_range=False."""
        from chemsmart.utils.utils import parse_index_specification

        # Test out-of-range positive index should fail
        with pytest.raises(ValueError, match="out of range"):
            parse_index_specification(
                "8", total_count=5, allow_out_of_range=False
            )

        # Test out-of-range negative index should fail
        with pytest.raises(ValueError, match="out of range"):
            parse_index_specification(
                "-6", total_count=5, allow_out_of_range=False
            )

        # Test range with out-of-bounds indices should fail
        with pytest.raises(ValueError, match="out of range"):
            parse_index_specification(
                "3-8", total_count=5, allow_out_of_range=False
            )

    def test_parse_index_boundary_detection_enabled(self):
        """Test boundary detection when allow_out_of_range=True."""
        from chemsmart.utils.utils import parse_index_specification

        # Test out-of-range indices are filtered out, valid ones remain
        result = parse_index_specification(
            "3,8,2", total_count=5, allow_out_of_range=True
        )
        assert result == [
            2,
            1,
        ]  # Only indices 3 and 2 (0-based: 2, 1) are valid, 8 is filtered

        # Test all out-of-range should raise error
        with pytest.raises(
            ValueError, match="All specified indices are out of range"
        ):
            parse_index_specification(
                "8,9,10", total_count=5, allow_out_of_range=True
            )


class TestIOUtilities:
    def test_clean_duplicate_structure(self):
        orientations = [
            np.array([1, 2, 3]),
            np.array([4, 5, 6]),
            np.array([4, 5, 6]),
        ]
        clean_duplicate_structure(orientations)
        assert len(orientations) == 2  # Should remove the duplicate

    def test_create_molecule_list(self):
        orientations = [np.array([[0, 0, 0]]), np.array([[1, 1, 1]])]
        orientations_pbc = [None, None]
        energies = [1.0, 2.0]
        forces = [[np.array([0, 0, 0])], [np.array([0, 0, 0])]]
        symbols = ["H"]
        charge = 0
        multiplicity = 1
        frozen_atoms = None
        pbc_conditions = [False]
        molecules = create_molecule_list(
            orientations,
            orientations_pbc,
            energies,
            forces,
            symbols,
            charge,
            multiplicity,
            frozen_atoms,
            pbc_conditions,
        )

        assert len(molecules) == 2
        assert isinstance(molecules[0], Molecule)
        assert isinstance(molecules[1], Molecule)
        assert molecules[0].energy == 1.0
        assert molecules[1].energy == 2.0

    @pytest.mark.parametrize(
        "line,allow_sign,expected",
        [
            ("0", True, True),
            ("1 2 3", True, True),
            ("+1 -2 +003 0", True, True),
            ("   10    20   30   ", True, True),
            ("+0 -0 0", True, True),
            ("+1 -2", False, False),  # signs not allowed
            ("1 2 3", False, True),
            ("001 0002 3", False, True),
            ("", True, False),  # empty
            ("   ", True, False),  # whitespace only
            ("1.0 2 3", True, False),  # float present
            ("1e3 2 3", True, False),  # scientific notation is not int()
            ("1 two 3", True, False),  # non-numeric
        ],
    )
    def test_line_of_all_integers(self, line, allow_sign, expected):
        assert line_of_all_integers(line, allow_sign=allow_sign) is expected

    @pytest.mark.parametrize(
        "line,expected",
        [
            # Valid: first token int; remaining are
            # proper floats (decimal or exponent)
            ("3 1.0 -2.3 4e-2", True),
            ("0 .5 5. 5.0 -0.3E+2", True),
            ("-1 +.3 -0.5e2", True),
            ("+4  .7", True),
            # Invalid: remaining tokens are plain integers (assuming
            # your float pattern requires decimal/exponent)
            ("3 1 2 3", False),
            # Invalid: not enough floats (only an
            # integer). Recommended behavior = False.
            ("+4", False),
            ("7   ", False),
            # Invalid: bad first token or malformed floats
            ("3.0 1.0 2.0", False),  # first token is not an integer
            ("x 1.0 2.0", False),  # first token non-numeric
            ("2 1.0 nope", False),  # invalid float token
            ("2 1.0 2.0e", False),  # malformed exponent
            # Whitespace / empty
            ("   ", False),
            ("", False),
        ],
    )
    def test_line_of_integer_followed_by_floats(self, line, expected):
        assert line_of_integer_followed_by_floats(line) is expected

    def test_header_like_then_data_like(self):
        # Typical ORCA header/data pattern
        header = "0 1 2 3 4 5"
        data = "0   0.123   -0.456   7.89   1e-2   .3"
        assert line_of_all_integers(header) is True
        assert line_of_integer_followed_by_floats(data) is True

    def test_trailing_and_leading_spaces(self):
        assert line_of_all_integers("   1   2   3   ") is True
        assert (
            line_of_integer_followed_by_floats("   5    1.0   2e0   .3   ")
            is True
        )

    def test_reject_plain_ints_as_floats(self):
        # Ensures your float regex isn't too permissive
        assert line_of_integer_followed_by_floats("2 3") is False
        assert (
            line_of_integer_followed_by_floats("2 3.") is True
        )  # decimal present
        assert (
            line_of_integer_followed_by_floats("2 3e0") is True
        )  # exponent present

    def test_clean_label(self):
        # spaces -> "_"
        assert clean_label("label with space") == "label_with_space"

        # commas -> "_"
        assert clean_label("label,with,comma") == "label_with_comma"

        # periods and parentheses -> "_"
        assert clean_label("Fig. 1(a)") == "Fig_1_a"

        # apostrophe -> "_prime_"
        assert clean_label("O'Hara") == "O_prime_Hara"

        # asterisk -> "_star_"
        assert clean_label("label*") == "label_star"

        # combination of several special characters
        assert (
            clean_label("O'Hara* test, v1.0") == "O_prime_Hara_star_test_v1_0"
        )

        # --- edge cases around underscore collapsing/stripping ---
        # 1) Empty string input
        assert clean_label("") == ""

        # 2) String with only special characters
        # "***" -> "_star__star__star_" -> collapse + strip -> "star_star_star"
        assert clean_label("***") == "star_star_star"

        # 3) Leading/trailing underscores after conversion
        # "*label*" -> "_star_label_star_" ->
        # collapse + strip -> "star_label_star"
        assert clean_label("*label*") == "star_label_star"

        # 4) Multiple consecutive special characters
        # "label...test" -> "label___test" -> collapse -> "label_test"
        assert clean_label("label...test") == "label_test"

    @pytest.mark.parametrize(
        "input_str, expected",
        [
            ("1-10", "id 1-10"),
            ("11", "id 11"),
            ("1-10,11", "id 1-10 or id 11"),
            ("1-10,11,14,19-30", "id 1-10 or id 11 or id 14 or id 19-30"),
        ],
    )
    def test_basic_conversion(self, input_str, expected):
        assert (
            convert_string_indices_to_pymol_id_indices(input_str) == expected
        )

    def test_conversion_strips_whitespace(self):
        input_str = " 1-10,  11 ,14 ,  19-30  "
        expected = "id 1-10 or id 11 or id 14 or id 19-30"
        assert (
            convert_string_indices_to_pymol_id_indices(input_str) == expected
        )

    def test_trailing_comma_is_ignored(self):
        input_str = "1-10,"
        expected = "id 1-10"
        assert (
            convert_string_indices_to_pymol_id_indices(input_str) == expected
        )

    @pytest.mark.parametrize("bad_input", ["", "   ", ",,,", ",  ,"])
    def test_raises_value_error_on_empty_or_invalid_input(self, bad_input):
        with pytest.raises(ValueError):
            convert_string_indices_to_pymol_id_indices(bad_input)


class TestNaturallySorted:
    def test_empty_list(self):
        """Test sorting an empty list."""
        assert naturally_sorted([]) == []

    def test_single_item(self):
        """Test sorting a list with one item."""
        assert naturally_sorted(["item1"]) == ["item1"]

    def test_numeric_order(self):
        """Test sorting strings with numbers in natural order."""
        input_list = ["z10", "z2", "z1"]
        expected = ["z1", "z2", "z10"]
        assert naturally_sorted(input_list) == expected

    def test_mixed_case(self):
        """Test sorting with mixed case letters."""
        input_list = ["Z1", "z2", "Z10", "z1"]
        expected = ["Z1", "z1", "z2", "Z10"]
        assert naturally_sorted(input_list) == expected

    def test_alphanumeric(self):
        """Test sorting alphanumeric strings."""
        input_list = ["a11", "a1", "b2", "b10"]
        expected = ["a1", "a11", "b2", "b10"]
        assert naturally_sorted(input_list) == expected

    def test_file_names(self):
        """Test sorting typical file names."""
        input_list = ["file10.txt", "file2.txt", "file1.txt"]
        expected = ["file1.txt", "file2.txt", "file10.txt"]
        assert naturally_sorted(input_list) == expected

    def test_special_characters(self):
        """Test sorting with special characters."""
        input_list = ["item-2", "item_10", "item_1"]
        expected = ["item-2", "item_1", "item_10"]
        assert naturally_sorted(input_list) == expected

    def test_mixed_types(self):
        """Test sorting with mixed formats
        (numbers, letters, and empty strings)."""
        input_list = ["100", "2", "abc", "", "Z", "z1"]
        expected = ["", "2", "100", "abc", "Z", "z1"]
        assert naturally_sorted(input_list) == expected

    def test_large_numbers(self):
        """Test sorting with large numbers."""
        input_list = ["item1000", "item999", "item10000"]
        expected = ["item999", "item1000", "item10000"]
        assert naturally_sorted(input_list) == expected

    def test_no_numbers(self):
        """Test sorting strings without numbers."""
        input_list = ["zebra", "Apple", "banana"]
        expected = ["Apple", "banana", "zebra"]
        assert naturally_sorted(input_list) == expected


class TestRunCommand:
    """Tests for the run_command utility function."""

    def test_list_command_success(self, mock_popen):
        """Test running a command provided as
        a list with successful execution."""
        mock_process = mock_popen.return_value
        mock_process.communicate.return_value = ("dir contents\n", "")
        mock_process.returncode = 0

        result = run_command(["ls", "-l"])
        assert result == "dir contents"
        mock_popen.assert_called_once_with(
            ["ls", "-l"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

    def test_string_command_success(self, mock_popen):
        """Test running a command provided as
        a string with successful execution."""
        mock_process = mock_popen.return_value
        mock_process.communicate.return_value = ("dir contents\n", "")
        mock_process.returncode = 0

        result = run_command("ls -l")
        assert result == "dir contents"
        mock_popen.assert_called_once_with(
            ["ls", "-l"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

    def test_string_command_with_quotes(self, mock_popen):
        """Test running a string command with quoted arguments."""
        mock_process = mock_popen.return_value
        mock_process.communicate.return_value = ("committed\n", "")
        mock_process.returncode = 0

        result = run_command("git commit -m 'initial commit'")
        assert result == "committed"
        mock_popen.assert_called_once_with(
            ["git", "commit", "-m", "initial commit"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

    def test_command_failure(self, mock_popen, capture_log):
        """Test running a command that fails with non-zero return code."""
        mock_process = mock_popen.return_value
        mock_process.communicate.return_value = ("", "command not found")
        mock_process.returncode = 1

        result = run_command(["invalid_cmd"])
        assert result is None
        assert (
            "Error running ['invalid_cmd']: command not found"
            in capture_log.text
        )

    def test_command_exception(self, mock_popen, capture_log):
        """Test handling an exception during command execution."""
        mock_popen.side_effect = OSError("Permission denied")

        result = run_command(["ls", "-l"])
        assert result is None
        assert (
            "Exception while running ['ls', '-l']: Permission denied"
            in capture_log.text
        )

    def test_invalid_input_type(self, capture_log):
        """Test handling invalid input type (neither string nor list)."""
        result = run_command(123)
        assert result is None
        assert (
            "Invalid command type: <class 'int'>. Expected str or list."
            in capture_log.text
        )


class TestReturnObjectsAndIndicesFromStringIndex:
    """Tests for the return_objects_and_indices_from_string_index
    utility function."""

    def test_single_index_string(self):
        """Test single index as a string (1-based)."""
        objects = ["a", "b", "c", "d", "e"]
        result_objects, result_indices = (
            return_objects_and_indices_from_string_index(objects, "1")
        )
        assert result_objects == "a"
        assert result_indices == 1

    def test_single_index_middle(self):
        """Test single index in middle of list."""
        objects = ["a", "b", "c", "d", "e"]
        result_objects, result_indices = (
            return_objects_and_indices_from_string_index(objects, "3")
        )
        assert result_objects == "c"
        assert result_indices == 3

    def test_single_index_last(self):
        """Test single index at end of list."""
        objects = ["a", "b", "c", "d", "e"]
        result_objects, result_indices = (
            return_objects_and_indices_from_string_index(objects, "5")
        )
        assert result_objects == "e"
        assert result_indices == 5

    def test_single_negative_index(self):
        """Test negative index (last item)."""
        objects = ["a", "b", "c", "d", "e"]
        result_objects, result_indices = (
            return_objects_and_indices_from_string_index(objects, "-1")
        )
        assert result_objects == "e"
        assert result_indices == 5

    def test_slice_range(self):
        """Test slice with start and stop (1-based, exclusive stop)."""
        objects = ["a", "b", "c", "d", "e"]
        result_objects, result_indices = (
            return_objects_and_indices_from_string_index(objects, "2:4")
        )
        assert result_objects == ["b", "c"]
        assert result_indices == [2, 3]

    def test_slice_from_start(self):
        """Test slice from beginning."""
        objects = ["a", "b", "c", "d", "e"]
        result_objects, result_indices = (
            return_objects_and_indices_from_string_index(objects, "1:3")
        )
        assert result_objects == ["a", "b"]
        assert result_indices == [1, 2]

    def test_slice_to_end(self):
        """Test slice to end using open-ended slice."""
        objects = ["a", "b", "c", "d", "e"]
        result_objects, result_indices = (
            return_objects_and_indices_from_string_index(objects, "3:")
        )
        assert result_objects == ["c", "d", "e"]
        assert result_indices == [3, 4, 5]

    def test_slice_from_beginning(self):
        """Test slice from beginning to index."""
        objects = ["a", "b", "c", "d", "e"]
        result_objects, result_indices = (
            return_objects_and_indices_from_string_index(objects, ":3")
        )
        assert result_objects == ["a", "b"]
        assert result_indices == [1, 2]

    def test_slice_all(self):
        """Test slice selecting all elements using ':'."""
        objects = ["a", "b", "c", "d", "e"]
        result_objects, result_indices = (
            return_objects_and_indices_from_string_index(objects, ":")
        )
        assert result_objects == ["a", "b", "c", "d", "e"]
        assert result_indices == [1, 2, 3, 4, 5]

    def test_slice_with_step(self):
        """Test slice with step parameter."""
        objects = ["a", "b", "c", "d", "e", "f", "g", "h"]
        result_objects, result_indices = (
            return_objects_and_indices_from_string_index(objects, "1:8:2")
        )
        assert result_objects == ["a", "c", "e", "g"]
        assert result_indices == [1, 3, 5, 7]

    def test_user_defined_range(self):
        """Test user-defined range format (comma-separated)."""
        objects = ["a", "b", "c", "d", "e"]
        result_objects, result_indices = (
            return_objects_and_indices_from_string_index(objects, "1,3,5")
        )
        assert result_objects == ["a", "c", "e"]
        assert result_indices == [1, 3, 5]

    def test_user_defined_range_with_hyphen(self):
        """Test user-defined range with hyphen notation."""
        objects = ["a", "b", "c", "d", "e", "f", "g", "h"]
        result_objects, result_indices = (
            return_objects_and_indices_from_string_index(objects, "1-3")
        )
        assert result_objects == ["a", "b", "c"]
        assert result_indices == [1, 2, 3]

    def test_user_defined_range_complex(self):
        """Test complex user-defined range."""
        objects = ["a", "b", "c", "d", "e", "f", "g", "h"]
        result_objects, result_indices = (
            return_objects_and_indices_from_string_index(objects, "1-3,5,7-8")
        )
        assert result_objects == ["a", "b", "c", "e", "g", "h"]
        assert result_indices == [1, 2, 3, 5, 7, 8]

    def test_range_with_brackets(self):
        """Test user-defined range with brackets."""
        objects = ["a", "b", "c", "d", "e"]
        result_objects, result_indices = (
            return_objects_and_indices_from_string_index(objects, "[1-3]")
        )
        assert result_objects == ["a", "b", "c"]
        assert result_indices == [1, 2, 3]

    def test_with_integer_objects(self):
        """Test with list of integers as objects."""
        objects = [10, 20, 30, 40, 50]
        result_objects, result_indices = (
            return_objects_and_indices_from_string_index(objects, "2:4")
        )
        assert result_objects == [20, 30]
        assert result_indices == [2, 3]

    def test_with_mixed_objects(self):
        """Test with list of mixed types as objects."""
        objects = [1, "two", 3.0, [4], {"five": 5}]
        result_objects, result_indices = (
            return_objects_and_indices_from_string_index(objects, "1,3,5")
        )
        assert result_objects == [1, 3.0, {"five": 5}]
        assert result_indices == [1, 3, 5]

    def test_specified_indices_5_to_8(self):
        """Test that specified indices are preserved
        (e.g., 5:8 gives indices 5,6,7)."""
        objects = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"]
        result_objects, result_indices = (
            return_objects_and_indices_from_string_index(objects, "5:8")
        )
        assert result_objects == ["e", "f", "g"]
        assert result_indices == [5, 6, 7]

    def test_empty_list_raises_index_error(self):
        """Test that accessing empty list raises IndexError."""
        objects = []
        with pytest.raises(IndexError):
            return_objects_and_indices_from_string_index(objects, "1")

    def test_index_zero_raises_value_error(self):
        """Test that index 0 raises ValueError (1-based indexing required)."""
        objects = ["a", "b", "c"]
        with pytest.raises(ValueError):
            return_objects_and_indices_from_string_index(objects, "0")

    def test_out_of_range_raises_index_error(self):
        """Test that out-of-range index raises IndexError."""
        objects = ["a", "b", "c"]
        with pytest.raises(IndexError):
            return_objects_and_indices_from_string_index(objects, "10")


class TestPKaTableParsing:
    """Tests for the pKa table parsing utility functions."""

    def test_parse_pka_table_txt(self, tmp_path):
        """Test parsing a whitespace-delimited .txt table file."""
        from chemsmart.utils.utils import parse_pka_table

        table_file = tmp_path / "molecules.txt"
        table_file.write_text(
            "# filepath proton_index charge multiplicity\n"
            "filepath proton_index charge multiplicity\n"  # header row
            "mol1.xyz 10 0 1\n"
            "mol2.xyz 15 1 2\n"
            "mol3.xyz 8 -1 1\n"
        )

        entries = parse_pka_table(str(table_file))
        assert len(entries) == 3
        assert entries[0].filepath == "mol1.xyz"
        assert entries[0].proton_index == 10
        assert entries[0].charge == 0
        assert entries[0].multiplicity == 1
        assert entries[1].filepath == "mol2.xyz"
        assert entries[1].proton_index == 15
        assert entries[1].charge == 1
        assert entries[1].multiplicity == 2
        assert entries[2].filepath == "mol3.xyz"
        assert entries[2].proton_index == 8
        assert entries[2].charge == -1
        assert entries[2].multiplicity == 1

    def test_parse_pka_table_csv(self, tmp_path):
        """Test parsing a comma-delimited .csv table file."""
        from chemsmart.utils.utils import parse_pka_table

        table_file = tmp_path / "molecules.csv"
        table_file.write_text(
            "filepath,proton_index,charge,multiplicity\n"
            "path/to/mol1.xyz,5,0,1\n"
            "path/to/mol2.xyz,12,-2,3\n"
        )

        entries = parse_pka_table(str(table_file))
        assert len(entries) == 2
        assert entries[0].filepath == "path/to/mol1.xyz"
        assert entries[0].proton_index == 5
        assert entries[0].charge == 0
        assert entries[0].multiplicity == 1
        assert entries[1].filepath == "path/to/mol2.xyz"
        assert entries[1].proton_index == 12
        assert entries[1].charge == -2
        assert entries[1].multiplicity == 3

    def test_parse_pka_table_skip_comments_and_empty_lines(self, tmp_path):
        """Test that comments and empty lines are skipped."""
        from chemsmart.utils.utils import parse_pka_table

        table_file = tmp_path / "molecules.txt"
        table_file.write_text(
            "# This is a comment\n"
            "filepath proton_index charge multiplicity\n"
            "\n"
            "mol1.xyz 1 0 1\n"
            "# Another comment\n"
            "\n"
            "mol2.xyz 2 0 1\n"
        )

        entries = parse_pka_table(str(table_file))
        assert len(entries) == 2

    def test_parse_pka_table_invalid_column_count(self, tmp_path):
        """Test that invalid column count raises ValueError."""
        from chemsmart.utils.utils import parse_pka_table

        table_file = tmp_path / "bad.txt"
        table_file.write_text(
            "filepath proton_index charge\n"  # header
            "mol1.xyz 10 0\n"  # missing multiplicity
        )

        with pytest.raises(ValueError, match="expected 4 columns"):
            parse_pka_table(str(table_file))

    def test_parse_pka_table_invalid_integer(self, tmp_path):
        """Test that invalid integer values raise ValueError."""
        from chemsmart.utils.utils import parse_pka_table

        table_file = tmp_path / "bad.txt"
        table_file.write_text(
            "filepath proton_index charge multiplicity\n"
            "mol1.xyz abc 0 1\n"  # invalid proton_index
        )

        with pytest.raises(ValueError, match="not an integer"):
            parse_pka_table(str(table_file))

    def test_parse_pka_table_empty_raises(self, tmp_path):
        """Test that empty table raises ValueError."""
        from chemsmart.utils.utils import parse_pka_table

        table_file = tmp_path / "empty.txt"
        table_file.write_text("# Only comments\n")

        with pytest.raises(ValueError, match="No valid entries"):
            parse_pka_table(str(table_file))

    def test_parse_pka_table_file_not_found(self):
        """Test that missing file raises FileNotFoundError."""
        from chemsmart.utils.utils import parse_pka_table

        with pytest.raises(FileNotFoundError):
            parse_pka_table("/nonexistent/path/molecules.txt")

    def test_pka_table_entry_validate_missing_file(self, tmp_path):
        """Test PKaTableEntry validation catches missing files."""
        from chemsmart.utils.utils import PKaTableEntry

        entry = PKaTableEntry(
            filepath="/nonexistent/file.xyz",
            proton_index=10,
            charge=0,
            multiplicity=1,
            row_number=1,
        )

        with pytest.raises(ValueError, match="File not found"):
            entry.validate()

    def test_pka_table_entry_validate_invalid_proton_index(self, tmp_path):
        """Test PKaTableEntry validation catches invalid proton_index."""
        from chemsmart.utils.utils import PKaTableEntry

        # Create a real file for this test
        test_file = tmp_path / "test.xyz"
        test_file.write_text("1\n\nH 0 0 0\n")

        entry = PKaTableEntry(
            filepath=str(test_file),
            proton_index=0,  # Invalid: must be >= 1
            charge=0,
            multiplicity=1,
        )

        with pytest.raises(ValueError, match="proton_index must be >= 1"):
            entry.validate()

    def test_pka_table_entry_validate_invalid_multiplicity(self, tmp_path):
        """Test PKaTableEntry validation catches invalid multiplicity."""
        from chemsmart.utils.utils import PKaTableEntry

        test_file = tmp_path / "test.xyz"
        test_file.write_text("1\n\nH 0 0 0\n")

        entry = PKaTableEntry(
            filepath=str(test_file),
            proton_index=1,
            charge=0,
            multiplicity=0,  # Invalid: must be >= 1
        )

        with pytest.raises(ValueError, match="multiplicity must be >= 1"):
            entry.validate()

    def test_validate_pka_table_entries(self, tmp_path):
        """Test batch validation of PKaTableEntry list."""
        from chemsmart.utils.utils import (
            PKaTableEntry,
            validate_pka_table_entries,
        )

        # Create test files
        file1 = tmp_path / "mol1.xyz"
        file1.write_text("1\n\nH 0 0 0\n")
        file2 = tmp_path / "mol2.xyz"
        file2.write_text("1\n\nH 0 0 0\n")

        entries = [
            PKaTableEntry(str(file1), 1, 0, 1),
            PKaTableEntry(str(file2), 2, -1, 2),
        ]

        # Should not raise
        result = validate_pka_table_entries(entries, check_file_exists=True)
        assert result == entries

    def test_pka_table_entry_repr(self):
        """Test PKaTableEntry string representation."""
        from chemsmart.utils.utils import PKaTableEntry

        entry = PKaTableEntry(
            filepath="mol.xyz",
            proton_index=10,
            charge=0,
            multiplicity=1,
        )

        repr_str = repr(entry)
        assert "PKaTableEntry" in repr_str
        assert "mol.xyz" in repr_str
        assert "'proton_index': 10" in repr_str
        assert "'charge': 0" in repr_str
        assert "'multiplicity': 1" in repr_str

    def test_pka_table_entry_from_headers_and_row_dynamic(self):
        from chemsmart.utils.utils import PKaTableEntry

        headers = ["ha_file", "proton_index", "charge", "mult", "route"]
        row = ["mol.xyz", "10", "0", "1", "opt freq"]

        entry = PKaTableEntry.from_headers_and_row(
            headers, row, row_number=2
        )._data

        assert entry["filepath"] == "mol.xyz"
        assert entry["proton_index"] == "10"  # raw row value preserved
        assert entry["charge"] == "0"
        assert entry["multiplicity"] == "1"
        assert entry["route"] == "opt freq"

    def test_pka_table_entry_dict_and_kwargs_helpers(self):
        """to_dict/to_kwargs should support forwarding to downstream settings."""
        from chemsmart.utils.utils import PKaTableEntry

        entry = PKaTableEntry(
            {
                "filepath": "mol.xyz",
                "charge": 0,
                "multiplicity": 1,
                "extra": None,
            }
        )

        data = entry.to_dict()
        assert data["filepath"] == "mol.xyz"
        assert data["charge"] == 0

        forwarded = entry.to_kwargs(rename_map={"filepath": "filename"})
        assert "filename" in forwarded
        assert forwarded["filename"] == "mol.xyz"

        forwarded_no_none = entry.to_kwargs(drop_none=True)
        assert "extra" not in forwarded_no_none

    def test_pka_table_entry_alias_resolution(self, tmp_path):
        """Alias resolution should keep backward-compatible attribute access."""
        from chemsmart.utils.utils import PKaTableEntry

        test_file = tmp_path / "test.xyz"
        test_file.write_text("1\n\nH 0 0 0\n")

        entry = PKaTableEntry(
            {
                "ha_file": str(test_file),
                "pi": 1,
                "q": 0,
                "mult": 1,
            }
        )

        # canonical access via aliases
        assert entry.filepath == str(test_file)
        assert entry.proton_index == 1
        assert entry.charge == 0
        assert entry.multiplicity == 1

        # and validation still passes
        entry.validate()

    # Tests for the pKa output-table parsing utilities (--output-table).

    def test_parse_pka_output_table_csv(self, tmp_path):
        """Test parsing a CSV output table with canonical column names."""
        from chemsmart.utils.utils import parse_pka_output_table

        # Create dummy output files
        for name in [
            "acid_opt.log",
            "base_opt.log",
            "ref_opt.log",
            "refbase_opt.log",
            "acid_sp.log",
            "base_sp.log",
            "ref_sp.log",
            "refbase_sp.log",
        ]:
            (tmp_path / name).write_text("dummy")

        table_file = tmp_path / "outputs.csv"
        table_file.write_text(
            "basename,ha_gas,a_gas,hb_gas,b_gas,ha_sp,a_sp,hb_sp,b_sp,pka_ref\n"
            f"system1,{tmp_path}/acid_opt.log,{tmp_path}/base_opt.log,"
            f"{tmp_path}/ref_opt.log,{tmp_path}/refbase_opt.log,"
            f"{tmp_path}/acid_sp.log,{tmp_path}/base_sp.log,"
            f"{tmp_path}/ref_sp.log,{tmp_path}/refbase_sp.log,6.75\n"
        )

        entries = parse_pka_output_table(str(table_file))
        assert len(entries) == 1
        assert entries[0].basename == "system1"
        assert entries[0]["ha_gas"] == f"{tmp_path}/acid_opt.log"
        assert entries[0]["pka_ref"] == 6.75

    def test_parse_pka_output_table_alias_columns(self, tmp_path):
        """Test that aliased column names (e.g., HA_optimization_output) work."""
        from chemsmart.utils.utils import parse_pka_output_table

        for name in [
            "a.log",
            "b.log",
            "c.log",
            "d.log",
            "e.log",
            "f.log",
            "g.log",
            "h.log",
        ]:
            (tmp_path / name).write_text("dummy")

        table_file = tmp_path / "outputs.csv"
        table_file.write_text(
            "basename,HA_optimization_output,A_optimization_output,"
            "HB_optimization_output,B_optimization_output,"
            "HA_single_point_output,A_single_point_output,"
            "HB_single_point_output,B_single_point_output,reference_pka\n"
            f"sys1,{tmp_path}/a.log,{tmp_path}/b.log,"
            f"{tmp_path}/c.log,{tmp_path}/d.log,"
            f"{tmp_path}/e.log,{tmp_path}/f.log,"
            f"{tmp_path}/g.log,{tmp_path}/h.log,14.0\n"
        )

        entries = parse_pka_output_table(str(table_file))
        assert len(entries) == 1
        assert entries[0].basename == "sys1"
        assert entries[0]["ha_gas"] == f"{tmp_path}/a.log"
        assert entries[0]["pka_ref"] == 14.0

    def test_parse_pka_output_table_missing_required_columns(self, tmp_path):
        """Test that missing required columns raise ValueError."""
        from chemsmart.utils.utils import parse_pka_output_table

        table_file = tmp_path / "bad.csv"
        table_file.write_text("basename,ha_gas\n" "sys1,file.log\n")

        with pytest.raises(ValueError, match="missing required columns"):
            parse_pka_output_table(str(table_file))

    def test_parse_pka_output_table_empty(self, tmp_path):
        """Test that an empty output table raises ValueError."""
        from chemsmart.utils.utils import parse_pka_output_table

        table_file = tmp_path / "empty.csv"
        table_file.write_text("# comment only\n")

        with pytest.raises(ValueError, match="No valid entries"):
            parse_pka_output_table(str(table_file))

    def test_parse_pka_output_table_file_not_found(self):
        """Test that a missing table file raises FileNotFoundError."""
        from chemsmart.utils.utils import parse_pka_output_table

        with pytest.raises(FileNotFoundError):
            parse_pka_output_table("/nonexistent/path.csv")

    def test_resolve_pka_output_references_carry_forward(self, tmp_path):
        """Test that blank reference cells are filled from the previous row."""
        from chemsmart.utils.utils import (
            PKaOutputTableEntry,
            resolve_pka_output_references,
        )

        entries = [
            PKaOutputTableEntry(
                {
                    "basename": "sys1",
                    "ha_gas": "a1.log",
                    "a_gas": "b1.log",
                    "href_gas": "ref.log",
                    "ref_gas": "refbase.log",
                    "ha_sp": "a1_sp.log",
                    "a_sp": "b1_sp.log",
                    "href_sp": "ref_sp.log",
                    "ref_sp": "refbase_sp.log",
                    "pka_ref": 6.75,
                },
                row_number=2,
            ),
            PKaOutputTableEntry(
                {
                    "basename": "sys2",
                    "ha_gas": "a2.log",
                    "a_gas": "b2.log",
                    "href_gas": None,
                    "ref_gas": None,
                    "ha_sp": "a2_sp.log",
                    "a_sp": "b2_sp.log",
                    "href_sp": None,
                    "ref_sp": None,
                    "pka_ref": None,
                },
                row_number=3,
            ),
        ]

        resolve_pka_output_references(entries)

        # Second row should have inherited reference from first row
        assert entries[1]["href_gas"] == "ref.log"
        assert entries[1]["ref_gas"] == "refbase.log"
        assert entries[1]["href_sp"] == "ref_sp.log"
        assert entries[1]["ref_sp"] == "refbase_sp.log"
        assert entries[1]["pka_ref"] == 6.75

    def test_resolve_pka_output_references_first_row_blank_raises(self):
        """Test that blank reference in first row raises ValueError."""
        from chemsmart.utils.utils import (
            PKaOutputTableEntry,
            resolve_pka_output_references,
        )

        entries = [
            PKaOutputTableEntry(
                {
                    "basename": "sys1",
                    "ha_gas": "a.log",
                    "a_gas": "b.log",
                    "href_gas": None,
                    "ref_gas": None,
                    "ha_sp": "a_sp.log",
                    "a_sp": "b_sp.log",
                    "href_sp": None,
                    "ref_sp": None,
                    "pka_ref": None,
                },
                row_number=2,
            ),
        ]

        with pytest.raises(ValueError, match="blank in row 2"):
            resolve_pka_output_references(entries)

    def test_resolve_pka_output_references_partial_carry_forward(self):
        """Test carry-forward when only some ref columns are blank."""
        from chemsmart.utils.utils import (
            PKaOutputTableEntry,
            resolve_pka_output_references,
        )

        entries = [
            PKaOutputTableEntry(
                {
                    "basename": "sys1",
                    "ha_gas": "a.log",
                    "a_gas": "b.log",
                    "href_gas": "ref1.log",
                    "ref_gas": "refbase1.log",
                    "ha_sp": "a_sp.log",
                    "a_sp": "b_sp.log",
                    "href_sp": "ref1_sp.log",
                    "ref_sp": "refbase1_sp.log",
                    "pka_ref": 6.75,
                },
                row_number=2,
            ),
            PKaOutputTableEntry(
                {
                    "basename": "sys2",
                    "ha_gas": "a2.log",
                    "a_gas": "b2.log",
                    "href_gas": "ref2.log",
                    "ref_gas": None,  # partially blank
                    "ha_sp": "a2_sp.log",
                    "a_sp": "b2_sp.log",
                    "href_sp": None,
                    "ref_sp": None,
                    "pka_ref": 4.5,  # new pka_ref
                },
                row_number=3,
            ),
        ]

        resolve_pka_output_references(entries)

        assert entries[1]["href_gas"] == "ref2.log"  # kept own value
        assert entries[1]["ref_gas"] == "refbase1.log"  # carried forward
        assert entries[1]["href_sp"] == "ref1_sp.log"  # carried forward
        assert entries[1]["ref_sp"] == "refbase1_sp.log"  # carried forward
        assert entries[1]["pka_ref"] == 4.5  # kept own value

    def test_pka_output_table_entry_validate_valid(self, tmp_path):
        """Test validation passes for a complete, valid entry."""
        from chemsmart.utils.utils import PKaOutputTableEntry

        for name in [
            "a.log",
            "b.log",
            "c.log",
            "d.log",
            "e.log",
            "f.log",
            "g.log",
            "h.log",
        ]:
            (tmp_path / name).write_text("dummy")

        entry = PKaOutputTableEntry(
            {
                "basename": "test",
                "ha_gas": str(tmp_path / "a.log"),
                "a_gas": str(tmp_path / "b.log"),
                "href_gas": str(tmp_path / "c.log"),
                "ref_gas": str(tmp_path / "d.log"),
                "ha_sp": str(tmp_path / "e.log"),
                "a_sp": str(tmp_path / "f.log"),
                "href_sp": str(tmp_path / "g.log"),
                "ref_sp": str(tmp_path / "h.log"),
                "pka_ref": 6.75,
            },
            row_number=2,
        )

        # Should not raise
        entry.validate(check_file_exists=True)

    def test_pka_output_table_entry_validate_missing_file(self, tmp_path):
        """Test validation catches missing files."""
        from chemsmart.utils.utils import PKaOutputTableEntry

        entry = PKaOutputTableEntry(
            {
                "basename": "test",
                "ha_gas": "/nonexistent/file.log",
                "a_gas": "/nonexistent/file2.log",
                "href_gas": "/nonexistent/file3.log",
                "ref_gas": "/nonexistent/file4.log",
                "ha_sp": "/nonexistent/file5.log",
                "a_sp": "/nonexistent/file6.log",
                "href_sp": "/nonexistent/file7.log",
                "ref_sp": "/nonexistent/file8.log",
                "pka_ref": 6.75,
            },
            row_number=2,
        )

        with pytest.raises(ValueError, match="File not found"):
            entry.validate(check_file_exists=True)

    def test_pka_output_table_entry_validate_missing_basename(self):
        """Test validation catches missing basename."""
        from chemsmart.utils.utils import PKaOutputTableEntry

        entry = PKaOutputTableEntry(
            {
                "basename": "",
                "ha_gas": "a.log",
                "a_gas": "b.log",
                "href_gas": "c.log",
                "ref_gas": "d.log",
                "ha_sp": "e.log",
                "a_sp": "f.log",
                "href_sp": "g.log",
                "ref_sp": "h.log",
                "pka_ref": 6.75,
            },
            row_number=1,
        )

        with pytest.raises(ValueError, match="Missing basename"):
            entry.validate(check_file_exists=False)

    def test_pka_output_table_entry_repr(self):
        """Test PKaOutputTableEntry string representation."""
        from chemsmart.utils.utils import PKaOutputTableEntry

        entry = PKaOutputTableEntry(
            {"basename": "my_system", "ha_gas": "a.log"},
            row_number=5,
        )

        repr_str = repr(entry)
        assert "PKaOutputTableEntry" in repr_str
        assert "my_system" in repr_str
        assert "5" in repr_str

    def test_pka_output_table_entry_to_dict(self):
        """Test to_dict returns all stored data."""
        from chemsmart.utils.utils import PKaOutputTableEntry

        data = {"basename": "sys", "ha_gas": "a.log", "pka_ref": 6.75}
        entry = PKaOutputTableEntry(data)
        d = entry.to_dict()
        assert d["basename"] == "sys"
        assert d["ha_gas"] == "a.log"
        assert d["pka_ref"] == 6.75

    def test_export_pka_results_table_csv(self, tmp_path):
        """Test that export writes CSV with pka column appended."""
        import pandas as pd

        from chemsmart.utils.utils import (
            PKaOutputTableEntry,
            export_pka_results_table,
        )

        entries = [
            PKaOutputTableEntry(
                {
                    "basename": "sys1",
                    "ha_gas": "a.log",
                    "a_gas": "b.log",
                    "hb_gas": "c.log",
                    "b_gas": "d.log",
                    "ha_sp": "e.log",
                    "a_sp": "f.log",
                    "hb_sp": "g.log",
                    "b_sp": "h.log",
                    "pka_ref": 6.75,
                }
            ),
        ]
        results = [
            {"pKa": 12.34, "delta_G_soln_kcal_mol": 7.89, "basename": "sys1"},
        ]

        out_path = tmp_path / "results.csv"
        export_pka_results_table(entries, results, str(out_path))

        df = pd.read_csv(str(out_path))
        assert "pka" in df.columns
        assert "delta_G_soln_kcal_mol" in df.columns
        assert len(df) == 1
        assert abs(df["pka"].iloc[0] - 12.34) < 0.01

    def test_export_pka_results_table_txt(self, tmp_path):
        """Test that export writes tab-delimited TXT for non-.csv extension."""
        import pandas as pd

        from chemsmart.utils.utils import (
            PKaOutputTableEntry,
            export_pka_results_table,
        )

        entries = [
            PKaOutputTableEntry(
                {
                    "basename": "sys1",
                    "ha_gas": "a.log",
                    "a_gas": "b.log",
                    "hb_gas": "c.log",
                    "b_gas": "d.log",
                    "ha_sp": "e.log",
                    "a_sp": "f.log",
                    "hb_sp": "g.log",
                    "b_sp": "h.log",
                    "pka_ref": 6.75,
                }
            ),
        ]
        results = [
            {"pKa": 5.0, "delta_G_soln_kcal_mol": -2.1, "basename": "sys1"},
        ]

        out_path = tmp_path / "results.txt"
        export_pka_results_table(entries, results, str(out_path))

        df = pd.read_csv(str(out_path), sep="\t")
        assert "pka" in df.columns
        assert len(df) == 1

    def test_parse_and_resolve_multi_row_table(self, tmp_path):
        """End-to-end test: parse → resolve → validate on a multi-row table."""
        from chemsmart.utils.utils import (
            parse_pka_output_table,
            resolve_pka_output_references,
        )

        # Create dummy files
        for name in [
            "a1.log",
            "b1.log",
            "a2.log",
            "b2.log",
            "a1_sp.log",
            "b1_sp.log",
            "a2_sp.log",
            "b2_sp.log",
            "ref.log",
            "refbase.log",
            "ref_sp.log",
            "refbase_sp.log",
        ]:
            (tmp_path / name).write_text("dummy")

        table_file = tmp_path / "outputs.csv"
        table_file.write_text(
            "basename,ha_gas,a_gas,href_gas,ref_gas,ha_sp,a_sp,href_sp,ref_sp,pka_ref\n"
            f"sys1,{tmp_path}/a1.log,{tmp_path}/b1.log,"
            f"{tmp_path}/ref.log,{tmp_path}/refbase.log,"
            f"{tmp_path}/a1_sp.log,{tmp_path}/b1_sp.log,"
            f"{tmp_path}/ref_sp.log,{tmp_path}/refbase_sp.log,6.75\n"
            f"sys2,{tmp_path}/a2.log,{tmp_path}/b2.log,"
            f",,"
            f"{tmp_path}/a2_sp.log,{tmp_path}/b2_sp.log,"
            f",,\n"
        )

        entries = parse_pka_output_table(str(table_file))
        assert len(entries) == 2

        resolve_pka_output_references(entries)

        # Second row should have inherited reference from first
        assert entries[1]["href_gas"] == f"{tmp_path}/ref.log"
        assert entries[1]["ref_gas"] == f"{tmp_path}/refbase.log"
        assert entries[1]["href_sp"] == f"{tmp_path}/ref_sp.log"
        assert entries[1]["ref_sp"] == f"{tmp_path}/refbase_sp.log"
        assert entries[1]["pka_ref"] == 6.75

        # Validate all entries
        for entry in entries:
            entry.validate(check_file_exists=True)
