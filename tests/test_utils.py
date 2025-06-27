import subprocess

import numpy as np
import pytest

from chemsmart.io.gaussian.input import Gaussian16Input
from chemsmart.io.molecules.structure import CoordinateBlock, Molecule
from chemsmart.utils.io import clean_duplicate_structure, create_molecule_list
from chemsmart.utils.utils import (
    cmp_with_ignore,
    content_blocks_by_paragraph,
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
