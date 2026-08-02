"""Tests for general utility functions and classes."""

import os
import tempfile
import time
from types import SimpleNamespace

import numpy as np
import pytest

from chemsmart.utils.utils import (
    OrderedSet,
    content_blocks_by_paragraph,
    convert_list_to_gaussian_frozen_list,
    convert_list_to_orca_frozen_list,
    convert_modred_list_to_string,
    convert_string_index_from_1_based_to_0_based,
    extract_number,
    file_cache,
    get_key_by_value_and_number,
    get_list_from_string_range,
    get_prepend_string_for_modred,
    get_prepend_string_list_from_modred_free_format,
    get_value_by_number,
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


class TestFileCache:
    """file_cache had no existing coverage anywhere."""

    def test_no_file_arguments_bypasses_caching(self):
        calls = []

        @file_cache()
        def add(a, b):
            calls.append((a, b))
            return a + b

        assert add(1, 2) == 3
        assert add(1, 2) == 3
        # every call actually invokes the function since there are no
        # file-path arguments to key the cache on
        assert len(calls) == 2

    def test_caches_result_for_unchanged_file(self, tmp_path):
        calls = []
        path = tmp_path / "f.txt"
        path.write_text("hello")

        @file_cache()
        def read_it(filepath):
            calls.append(filepath)
            with open(filepath) as f:
                return f.read()

        assert read_it(str(path)) == "hello"
        assert read_it(str(path)) == "hello"
        assert len(calls) == 1

    def test_copy_result_true_prevents_cache_pollution(self, tmp_path):
        path = tmp_path / "f.txt"
        path.write_text("x")

        @file_cache(copy_result=True)
        def get_list(filepath):
            return [1, 2, 3]

        first = get_list(str(path))
        first.append(999)
        second = get_list(str(path))
        assert second == [1, 2, 3]

    def test_copy_result_false_leaks_mutations_into_cache(self, tmp_path):
        path = tmp_path / "f.txt"
        path.write_text("x")

        @file_cache(copy_result=False)
        def get_list(filepath):
            return [1, 2, 3]

        first = get_list(str(path))
        first.append(999)
        second = get_list(str(path))
        assert second == [1, 2, 3, 999]

    def test_unwraps_staticmethod(self, tmp_path):
        path = tmp_path / "f.txt"
        path.write_text("x")
        calls = []

        class Foo:
            @file_cache()
            @staticmethod
            def bar(filepath):
                calls.append(filepath)
                return "result"

        assert Foo.bar(str(path)) == "result"
        assert Foo.bar(str(path)) == "result"
        assert len(calls) == 1

    def test_rejects_classmethod(self):
        with pytest.raises(
            ValueError, match="Unable to use this with classmethod"
        ):

            class Baz:
                @file_cache()
                @classmethod
                def qux(cls, x):
                    return x

    def test_recently_modified_integer_mtime_uses_content_hash(self, tmp_path):
        """When a file's mtime is a whole number of seconds and was
        modified recently, the cache key falls back to a content hash
        instead of the raw mtime (guards against low-resolution
        filesystem timestamps within the same second)."""
        path = tmp_path / "f.txt"
        path.write_text("hello")
        # force an integer mtime so the hash-based branch is taken
        now_int = int(time.time())
        os.utime(str(path), (now_int, now_int))
        calls = []

        @file_cache()
        def read_it(filepath):
            calls.append(filepath)
            with open(filepath) as f:
                return f.read()

        assert read_it(str(path)) == "hello"
        assert read_it(str(path)) == "hello"
        assert len(calls) == 1

        # changing the content (same integer mtime) invalidates the cache
        os.utime(str(path), (now_int, now_int))
        path.write_text("changed")
        os.utime(str(path), (now_int, now_int))
        assert read_it(str(path)) == "changed"
        assert len(calls) == 2


class TestConvertListToFrozenList:
    def test_gaussian_frozen_list_is_1_indexed(self):
        molecule = SimpleNamespace(chemical_symbols=["C", "H", "H", "H"])
        assert convert_list_to_gaussian_frozen_list([1, 3], molecule) == [
            -1,
            0,
            -1,
            0,
        ]

    def test_orca_frozen_list_is_0_indexed(self):
        molecule = SimpleNamespace(chemical_symbols=["C", "H", "H", "H"])
        assert convert_list_to_orca_frozen_list([0, 2], molecule) == [
            -1,
            0,
            -1,
            0,
        ]


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
    """Tests for
    write_list_of_lists_as_a_string_with_empty_line_between_lists."""

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

    def test_actual_int_positive_converted(self):
        """Passing a real int (not a numeric string) takes the direct
        int branch rather than the string parser."""
        assert convert_string_index_from_1_based_to_0_based(5) == 4

    def test_actual_int_negative_returned_as_is(self):
        assert convert_string_index_from_1_based_to_0_based(-2) == -2

    def test_actual_int_zero_raises_error(self):
        with pytest.raises(ValueError, match="out of range"):
            convert_string_index_from_1_based_to_0_based(0)

    def test_slice_passed_through_unchanged(self):
        s = slice(1, 5, 2)
        assert convert_string_index_from_1_based_to_0_based(s) is s

    def test_invalid_type_raises_value_error(self):
        with pytest.raises(ValueError, match="Invalid index type"):
            convert_string_index_from_1_based_to_0_based(1.5)


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

    def test_list_index(self):
        """A comma-separated index string resolves to a list of
        indices, selecting multiple objects at once."""
        objects = ["a", "b", "c", "d"]
        result = return_objects_from_string_index(objects, "1,3")
        assert result == ["a", "c"]


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


class TestGetValueByNumber:
    def test_returns_value_for_matching_key_number(self):
        data = {"atom1": "C", "atom2": "H", "atom10": "O"}
        assert get_value_by_number(2, data) == "H"
        assert get_value_by_number(10, data) == "O"

    def test_returns_none_when_no_key_matches(self):
        data = {"atom1": "C", "atom2": "H"}
        assert get_value_by_number(99, data) is None


class TestGetKeyByValueAndNumber:
    def test_returns_key_matching_both_value_and_number(self):
        data = {"charge1": 0, "charge2": 1, "charge3": 0}
        assert get_key_by_value_and_number(1, 2, data) == "charge2"

    def test_returns_none_when_value_does_not_match(self):
        data = {"charge1": 0, "charge2": 1}
        assert get_key_by_value_and_number(99, 2, data) is None

    def test_returns_none_when_number_does_not_match(self):
        data = {"charge1": 0, "charge2": 1}
        assert get_key_by_value_and_number(1, 99, data) is None

    def test_skips_keys_with_no_trailing_number(self):
        data = {"nonumberkey": 0, "charge2": 1}
        assert get_key_by_value_and_number(1, 2, data) == "charge2"


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
