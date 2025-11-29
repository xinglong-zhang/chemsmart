import subprocess

import numpy as np
import pytest

from chemsmart.io.gaussian.input import Gaussian16Input
from chemsmart.io.molecules.structure import CoordinateBlock, Molecule
from chemsmart.utils.io import (
    clean_duplicate_structure,
    clean_label,
    create_molecule_list,
    line_of_all_integers,
    line_of_integer_followed_by_floats,
)
from chemsmart.utils.utils import (
    cmp_with_ignore,
    content_blocks_by_paragraph,
    convert_string_index_from_1_based_to_0_based,
    get_list_from_string_range,
    is_float,
    iterative_compare,
    naturally_sorted,
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
        """Test the conversion of string indices to a list of integers; 1-based indices."""
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
            # Valid: first token int; remaining are proper floats (decimal or exponent)
            ("3 1.0 -2.3 4e-2", True),
            ("0 .5 5. 5.0 -0.3E+2", True),
            ("-1 +.3 -0.5e2", True),
            ("+4  .7", True),
            # Invalid: remaining tokens are plain integers (assuming your float pattern requires decimal/exponent)
            ("3 1 2 3", False),
            # Invalid: not enough floats (only an integer). Recommended behavior = False.
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
        # "*label*" -> "_star_label_star_" -> collapse + strip -> "star_label_star"
        assert clean_label("*label*") == "star_label_star"

        # 4) Multiple consecutive special characters
        # "label...test" -> "label___test" -> collapse -> "label_test"
        assert clean_label("label...test") == "label_test"


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
        """Test sorting with mixed formats (numbers, letters, and empty strings)."""
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
        """Test running a command provided as a list with successful execution."""
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
        """Test running a command provided as a string with successful execution."""
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
