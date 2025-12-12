"""Tests for general utility functions and classes."""

import os
import tempfile

import numpy as np
import pytest

from chemsmart.utils.utils import (
    OrderedSet,
    content_blocks_by_paragraph,
    convert_modred_list_to_string,
    convert_string_index_from_1_based_to_0_based,
    extract_number,
    get_list_from_string_range,
    get_prepend_string_for_modred,
    get_prepend_string_list_from_modred_free_format,
    is_float,
    iterative_compare,
    kabsch_align,
    kabsch_align2,
    return_objects_from_string_index,
    string2index_1based,
    strip_out_comments,
    two_files_have_similar_contents,
    two_lists_have_similar_contents,
    update_dict_with_existing_keys,
    write_list_of_lists_as_a_string_with_empty_line_between_lists,
)


class TestOrderedSet:
    """Tests for the OrderedSet class."""

    def test_empty_initialization(self):
        """Test that empty OrderedSet initializes correctly."""
        os_empty = OrderedSet()
        assert len(os_empty) == 0
        assert list(os_empty) == []

    def test_initialization_with_iterable(self):
        """Test initialization with an iterable."""
        os_init = OrderedSet([1, 2, 3])
        assert len(os_init) == 3
        assert list(os_init) == [1, 2, 3]

    def test_initialization_removes_duplicates(self):
        """Test that initialization removes duplicates."""
        os_dup = OrderedSet([1, 2, 2, 3, 1])
        assert len(os_dup) == 3
        assert list(os_dup) == [1, 2, 3]

    def test_add_new_item(self):
        """Test adding a new item."""
        os_add = OrderedSet()
        os_add.add(1)
        os_add.add(2)
        assert list(os_add) == [1, 2]

    def test_add_duplicate_item(self):
        """Test adding a duplicate item."""
        os_add = OrderedSet([1, 2])
        os_add.add(1)
        assert list(os_add) == [1, 2]

    def test_remove_existing_item(self):
        """Test removing an existing item."""
        os_rem = OrderedSet([1, 2, 3])
        os_rem.remove(2)
        assert list(os_rem) == [1, 3]

    def test_remove_non_existing_item(self):
        """Test removing a non-existing item does nothing."""
        os_rem = OrderedSet([1, 2, 3])
        os_rem.remove(5)  # Should not raise error
        assert list(os_rem) == [1, 2, 3]

    def test_contains(self):
        """Test __contains__ method."""
        os_con = OrderedSet([1, 2, 3])
        assert 1 in os_con
        assert 5 not in os_con

    def test_len(self):
        """Test __len__ method."""
        os_len = OrderedSet([1, 2, 3, 4, 5])
        assert len(os_len) == 5

    def test_iteration(self):
        """Test iteration preserves order."""
        items = [3, 1, 4, 1, 5, 9, 2, 6]
        os_iter = OrderedSet(items)
        expected = [3, 1, 4, 5, 9, 2, 6]  # Duplicates removed
        assert list(os_iter) == expected


class TestIsFloat:
    """Tests for the is_float function."""

    def test_float_string(self):
        """Test that float strings are recognized."""
        assert is_float("3.14") is True
        assert is_float("-2.5") is True
        assert is_float("1e-10") is True
        assert is_float("1.5E+3") is True

    def test_integer_string_not_float(self):
        """Test that integer strings are not recognized as floats."""
        assert is_float("123") is False
        assert is_float("-456") is False
        assert is_float("0") is False

    def test_non_numeric_string(self):
        """Test that non-numeric strings return False."""
        assert is_float("abc") is False
        assert is_float("1.2.3") is False
        assert is_float("") is False


class TestStripOutComments:
    """Tests for the strip_out_comments function."""

    def test_single_line_with_comment(self):
        """Test stripping comment from single line."""
        result = strip_out_comments("code # comment")
        assert result == "code"

    def test_multiple_lines_with_comments(self):
        """Test stripping comments from multiple lines."""
        input_str = "line1 # comment1\nline2 # comment2\nline3"
        result = strip_out_comments(input_str)
        assert result == "line1\nline2\nline3"

    def test_no_comments(self):
        """Test string without comments."""
        input_str = "no comments here"
        result = strip_out_comments(input_str)
        assert result == "no comments here"

    def test_empty_string(self):
        """Test empty string."""
        assert strip_out_comments("") == ""

    def test_only_comment(self):
        """Test line with only a comment."""
        result = strip_out_comments("# only comment")
        assert result == ""


class TestContentBlocksByParagraph:
    """Tests for the content_blocks_by_paragraph function."""

    def test_single_block(self):
        """Test single block without empty lines."""
        string_list = ["line1", "line2", "line3"]
        result = content_blocks_by_paragraph(string_list)
        assert result == [["line1", "line2", "line3"]]

    def test_two_blocks(self):
        """Test two blocks separated by empty line."""
        string_list = ["line1", "line2", "", "line3", "line4"]
        result = content_blocks_by_paragraph(string_list)
        assert result == [["line1", "line2"], ["line3", "line4"]]

    def test_multiple_empty_lines(self):
        """Test multiple consecutive empty lines."""
        string_list = ["line1", "", "", "line2"]
        result = content_blocks_by_paragraph(string_list)
        assert result == [["line1"], ["line2"]]

    def test_empty_list(self):
        """Test empty input list."""
        result = content_blocks_by_paragraph([])
        assert result == []


class TestWriteListOfListsAsString:
    """Tests for write_list_of_lists_as_a_string_with_empty_line_between_lists."""

    def test_single_list(self):
        """Test single list."""
        list_of_lists = [["a", "b"]]
        result = write_list_of_lists_as_a_string_with_empty_line_between_lists(
            list_of_lists
        )
        assert result == "a\nb\n"

    def test_two_lists(self):
        """Test two lists with empty line between."""
        list_of_lists = [["a", "b"], ["c", "d"]]
        result = write_list_of_lists_as_a_string_with_empty_line_between_lists(
            list_of_lists
        )
        assert result == "a\nb\n\nc\nd\n"

    def test_empty_list_of_lists(self):
        """Test empty list of lists."""
        result = write_list_of_lists_as_a_string_with_empty_line_between_lists(
            []
        )
        assert result == ""


class TestGetListFromStringRange:
    """Tests for the get_list_from_string_range function."""

    def test_simple_range(self):
        """Test simple range string."""
        result = get_list_from_string_range("1-3")
        assert result == [1, 2, 3]

    def test_comma_separated(self):
        """Test comma-separated values."""
        result = get_list_from_string_range("1,3,5")
        assert result == [1, 3, 5]

    def test_mixed_range_and_values(self):
        """Test mixed ranges and values."""
        result = get_list_from_string_range("1-3,5,7-9")
        assert result == [1, 2, 3, 5, 7, 8, 9]

    def test_with_brackets(self):
        """Test range with brackets."""
        result = get_list_from_string_range("[1-3,5]")
        assert result == [1, 2, 3, 5]


class TestString2Index1Based:
    """Tests for the string2index_1based function."""

    def test_single_positive_index(self):
        """Test single positive index."""
        result = string2index_1based("5")
        assert result == 4  # Converted to 0-based

    def test_negative_index(self):
        """Test negative index."""
        result = string2index_1based("-1")
        assert result == -1  # Negative stays as-is

    def test_slice_notation(self):
        """Test slice notation."""
        result = string2index_1based("1:5")
        assert result == slice(0, 4, None)

    def test_zero_index_converts_to_negative_one(self):
        """Test that zero index converts to -1 (0-1=-1)."""
        # The function converts 0 to -1 (0-based conversion)
        result = string2index_1based("0")
        assert result == -1

    def test_slice_with_step(self):
        """Test slice with step."""
        result = string2index_1based("1:10:2")
        assert result == slice(0, 9, 2)


class TestConvertStringIndex:
    """Tests for convert_string_index_from_1_based_to_0_based."""

    def test_positive_integer(self):
        """Test positive integer conversion."""
        result = convert_string_index_from_1_based_to_0_based("5")
        assert result == 4

    def test_negative_integer(self):
        """Test negative integer stays the same."""
        result = convert_string_index_from_1_based_to_0_based("-1")
        assert result == -1

    def test_zero_raises_error(self):
        """Test that zero raises error."""
        with pytest.raises(ValueError):
            convert_string_index_from_1_based_to_0_based("0")


class TestReturnObjectsFromStringIndex:
    """Tests for the return_objects_from_string_index function."""

    def test_single_index(self):
        """Test single index."""
        objects = ["a", "b", "c", "d"]
        result = return_objects_from_string_index(objects, "2")
        assert result == "b"  # 2 -> index 1 (0-based)

    def test_negative_index(self):
        """Test negative index."""
        objects = ["a", "b", "c", "d"]
        result = return_objects_from_string_index(objects, "-1")
        assert result == "d"

    def test_slice_index(self):
        """Test slice index."""
        objects = ["a", "b", "c", "d"]
        result = return_objects_from_string_index(objects, "1:3")
        assert result == ["a", "b"]


class TestGetPrependStringForModred:
    """Tests for the get_prepend_string_for_modred function."""

    def test_bond_two_atoms(self):
        """Test bond (2 atoms) returns 'B'."""
        result = get_prepend_string_for_modred([1, 2])
        assert result == "B"

    def test_angle_three_atoms(self):
        """Test angle (3 atoms) returns 'A'."""
        result = get_prepend_string_for_modred([1, 2, 3])
        assert result == "A"

    def test_dihedral_four_atoms(self):
        """Test dihedral (4 atoms) returns 'D'."""
        result = get_prepend_string_for_modred([1, 2, 3, 4])
        assert result == "D"

    def test_invalid_length_raises_error(self):
        """Test invalid length raises error."""
        with pytest.raises(ValueError):
            get_prepend_string_for_modred([1])
        with pytest.raises(ValueError):
            get_prepend_string_for_modred([1, 2, 3, 4, 5])


class TestConvertModredListToString:
    """Tests for the convert_modred_list_to_string function."""

    def test_simple_list(self):
        """Test simple list conversion."""
        result = convert_modred_list_to_string([1, 2, 3])
        assert result == "1 2 3"

    def test_single_item(self):
        """Test single item list."""
        result = convert_modred_list_to_string([5])
        assert result == "5"


class TestGetPrependStringListFromModredFreeFormat:
    """Tests for get_prepend_string_list_from_modred_free_format."""

    def test_list_of_lists(self):
        """Test list of lists input."""
        result = get_prepend_string_list_from_modred_free_format(
            [[1, 2], [3, 4, 5]], program="gaussian"
        )
        assert result == ["B 1 2", "A 3 4 5"]

    def test_single_list(self):
        """Test single list input."""
        result = get_prepend_string_list_from_modred_free_format(
            [1, 2, 3, 4], program="gaussian"
        )
        assert result == ["D 1 2 3 4"]

    def test_orca_program_0_indexed(self):
        """Test ORCA uses 0-indexed."""
        result = get_prepend_string_list_from_modred_free_format(
            [1, 2], program="orca"
        )
        assert result == ["B 0 1"]

    def test_invalid_input_raises_error(self):
        """Test invalid input raises error."""
        with pytest.raises(ValueError):
            get_prepend_string_list_from_modred_free_format("not a list")


class TestTwoFilesHaveSimilarContents:
    """Tests for the two_files_have_similar_contents function."""

    def test_identical_files(self):
        """Test identical files return True."""
        with (
            tempfile.NamedTemporaryFile(mode="w", delete=False) as f1,
            tempfile.NamedTemporaryFile(mode="w", delete=False) as f2,
        ):
            f1.write("line1\nline2\n")
            f2.write("line1\nline2\n")
            f1.flush()
            f2.flush()
            temp_name1 = f1.name
            temp_name2 = f2.name

        result = two_files_have_similar_contents(temp_name1, temp_name2)
        assert result is True

        os.unlink(temp_name1)
        os.unlink(temp_name2)

    def test_different_files(self):
        """Test different files return False."""
        with (
            tempfile.NamedTemporaryFile(mode="w", delete=False) as f1,
            tempfile.NamedTemporaryFile(mode="w", delete=False) as f2,
        ):
            f1.write("line1\nline2\n")
            f2.write("line1\nline3\n")
            f1.flush()
            f2.flush()
            temp_name1 = f1.name
            temp_name2 = f2.name

        result = two_files_have_similar_contents(temp_name1, temp_name2)
        assert result is False

        os.unlink(temp_name1)
        os.unlink(temp_name2)

    def test_with_ignored_string(self):
        """Test with ignored string."""
        with (
            tempfile.NamedTemporaryFile(mode="w", delete=False) as f1,
            tempfile.NamedTemporaryFile(mode="w", delete=False) as f2,
        ):
            f1.write("line1\ntimestamp: 12345\n")
            f2.write("line1\ntimestamp: 67890\n")
            f1.flush()
            f2.flush()
            temp_name1 = f1.name
            temp_name2 = f2.name

        result = two_files_have_similar_contents(
            temp_name1, temp_name2, ignored_string="timestamp"
        )
        assert result is True

        os.unlink(temp_name1)
        os.unlink(temp_name2)


class TestTwoListsHaveSimilarContents:
    """Tests for the two_lists_have_similar_contents function."""

    def test_identical_lists(self):
        """Test identical lists return True."""
        list1 = ["a", "b", "c"]
        list2 = ["a", "b", "c"]
        assert two_lists_have_similar_contents(list1, list2) is True

    def test_different_lists(self):
        """Test different lists return False."""
        list1 = ["a", "b", "c"]
        list2 = ["a", "b", "d"]
        assert two_lists_have_similar_contents(list1, list2) is False

    def test_different_lengths(self):
        """Test lists with different lengths return False."""
        list1 = ["a", "b"]
        list2 = ["a", "b", "c"]
        assert two_lists_have_similar_contents(list1, list2) is False

    def test_with_ignore_string(self):
        """Test with ignored string."""
        list1 = ["a", "timestamp: 123", "c"]
        list2 = ["a", "timestamp: 456", "c"]
        assert (
            two_lists_have_similar_contents(
                list1, list2, ignore_string="timestamp"
            )
            is True
        )


class TestUpdateDictWithExistingKeys:
    """Tests for the update_dict_with_existing_keys function."""

    def test_update_existing_keys(self):
        """Test updating existing keys."""
        dict1 = {"a": 1, "b": 2}
        dict2 = {"a": 10}
        result = update_dict_with_existing_keys(dict1, dict2)
        assert result == {"a": 10, "b": 2}

    def test_non_existing_key_raises_error(self):
        """Test non-existing key raises error."""
        dict1 = {"a": 1}
        dict2 = {"b": 2}
        with pytest.raises(ValueError):
            update_dict_with_existing_keys(dict1, dict2)


class TestIterativeCompare:
    """Tests for the iterative_compare function."""

    def test_removes_duplicates(self):
        """Test that duplicates are removed."""
        result = iterative_compare([1, 2, 2, 3, 1])
        assert result == [1, 2, 3]

    def test_preserves_order(self):
        """Test that order is preserved."""
        result = iterative_compare([3, 1, 4, 1, 5])
        assert result == [3, 1, 4, 5]

    def test_empty_list(self):
        """Test empty list returns empty."""
        result = iterative_compare([])
        assert result == []


class TestKabschAlign:
    """Tests for the kabsch_align function."""

    def test_identical_structures(self):
        """Test alignment of identical structures."""
        P = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        Q = P.copy()
        p_aligned, q_aligned, R, t, rmsd = kabsch_align(P, Q)
        assert np.isclose(rmsd, 0, atol=1e-10)

    def test_translated_structure(self):
        """Test alignment of translated structure."""
        P = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        Q = P + np.array([5, 5, 5])
        p_aligned, q_aligned, R, t, rmsd = kabsch_align(P, Q)
        assert np.isclose(rmsd, 0, atol=1e-10)

    def test_rotated_structure(self):
        """Test alignment of rotated structure."""
        P = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        # 90 degree rotation around z-axis
        R_rot = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        Q = np.dot(P, R_rot.T)
        p_aligned, q_aligned, R, t, rmsd = kabsch_align(P, Q)
        assert np.isclose(rmsd, 0, atol=1e-5)

    def test_dimension_mismatch_raises_error(self):
        """Test that dimension mismatch raises error."""
        P = np.array([[0, 0, 0], [1, 0, 0]])
        Q = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        with pytest.raises(AssertionError):
            kabsch_align(P, Q)


class TestKabschAlign2:
    """Tests for the kabsch_align2 function."""

    def test_identical_structures(self):
        """Test alignment of identical structures."""
        P = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
        Q = P.copy()
        aligned_P, Q_out, R, t, rmsd = kabsch_align2(P, Q)
        assert np.isclose(rmsd, 0, atol=1e-10)


class TestExtractNumber:
    """Tests for the extract_number function."""

    def test_extract_simple_number(self):
        """Test extracting simple number."""
        assert extract_number("c123") == 123
        assert extract_number("c1") == 1

    def test_extract_number_with_asterisk(self):
        """Test extracting number with asterisk."""
        assert extract_number("c456*") == 456

    def test_no_number_returns_inf(self):
        """Test no number returns infinity."""
        assert extract_number("abc") == float("inf")
