"""Tests for IO utility functions."""

import os
import tempfile

import numpy as np
import pytest

from chemsmart.utils.io import (
    clean_duplicate_structure,
    clean_label,
    convert_string_indices_to_pymol_id_indices,
    get_program_type_from_file,
    increment_numbers,
    load_molecules_from_paths,
    match_outfile_pattern,
    remove_keyword,
)


class TestIncrementNumbers:
    """Tests for the increment_numbers function."""

    def test_single_number(self):
        """Test incrementing a single number."""
        result = increment_numbers("atom 5", increment=1)
        assert result == "atom 6"

    def test_multiple_numbers(self):
        """Test incrementing multiple numbers."""
        result = increment_numbers("bond 1 2", increment=1)
        assert result == "bond 2 3"

    def test_custom_increment(self):
        """Test custom increment value."""
        result = increment_numbers("index 10", increment=5)
        assert result == "index 15"

    def test_no_numbers(self):
        """Test string with no numbers."""
        result = increment_numbers("no numbers here")
        assert result == "no numbers here"

    def test_zero_padding_removed(self):
        """Test that zero padding is handled."""
        result = increment_numbers("item 001")
        assert result == "item 2"


class TestRemoveKeyword:
    """Tests for the remove_keyword function."""

    def test_remove_keyword(self):
        """Test removing a keyword."""
        result = remove_keyword("the quick brown fox", "quick")
        assert result == "the  brown fox"

    def test_case_insensitive(self):
        """Test case-insensitive removal."""
        result = remove_keyword("The QUICK brown fox", "quick")
        assert result == "The  brown fox"

    def test_no_match(self):
        """Test when keyword is not present."""
        result = remove_keyword("hello world", "missing")
        assert result == "hello world"

    def test_word_boundary(self):
        """Test that word boundaries are respected."""
        # Should not match "quickly" when searching for "quick"
        result = remove_keyword("he quickly ran", "quick")
        assert result == "he quickly ran"

    def test_multiple_occurrences(self):
        """Test removing multiple occurrences."""
        result = remove_keyword("opt opt freq", "opt")
        assert result == "  freq"


class TestMatchOutfilePattern:
    """Tests for the match_outfile_pattern function."""

    def test_gaussian_pattern(self):
        """Test Gaussian output pattern detection."""
        assert match_outfile_pattern("Entering Gaussian System") == "gaussian"
        assert match_outfile_pattern("Gaussian, Inc.") == "gaussian"
        assert match_outfile_pattern("Gaussian(R)") == "gaussian"

    def test_orca_pattern(self):
        """Test ORCA output pattern detection."""
        assert match_outfile_pattern("* O   R   C   A *") == "orca"
        assert match_outfile_pattern("Your ORCA version") == "orca"

    def test_xtb_pattern(self):
        """Test xTB output pattern detection."""
        assert match_outfile_pattern("x T B") == "xtb"
        assert match_outfile_pattern("xtb version 6.5.0") == "xtb"

    def test_crest_pattern(self):
        """Test CREST output pattern detection."""
        assert match_outfile_pattern("C R E S T") == "crest"
        assert match_outfile_pattern("$ crest") == "crest"

    def test_no_match(self):
        """Test when no pattern matches."""
        assert match_outfile_pattern("random text") is None
        assert match_outfile_pattern("") is None


class TestGetOutfileFormat:
    """Tests for the get_outfile_format function."""

    def test_gaussian_output_detection(self, gaussian_singlet_opt_outfile):
        """Test detection of Gaussian output file."""
        result = get_program_type_from_file(gaussian_singlet_opt_outfile)
        assert result == "gaussian"

    def test_orca_output_detection(self, water_output_gas_path):
        """Test detection of ORCA output file."""
        result = get_program_type_from_file(water_output_gas_path)
        assert result == "orca"

    def test_unknown_format(self):
        """Test unknown file format."""
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as f:
            f.write("This is just random text\n" * 10)
            f.flush()
            temp_name = f.name
        result = get_program_type_from_file(temp_name)
        assert result == "unknown"
        os.unlink(f.name)


class TestCleanDuplicateStructure:
    """Tests for the clean_duplicate_structure function."""

    def test_removes_duplicate_last_structure(self):
        """Test that duplicate last structure is removed."""
        orientations = [
            np.array([[0, 0, 0], [1, 0, 0]]),
            np.array([[0, 0, 0], [1.5, 0, 0]]),
            np.array([[0, 0, 0], [1.5, 0, 0]]),  # Duplicate
        ]
        clean_duplicate_structure(orientations)
        assert len(orientations) == 2

    def test_no_duplicate(self):
        """Test that non-duplicate structures are preserved."""
        orientations = [
            np.array([[0, 0, 0], [1, 0, 0]]),
            np.array([[0, 0, 0], [1.5, 0, 0]]),
            np.array([[0, 0, 0], [2, 0, 0]]),
        ]
        clean_duplicate_structure(orientations)
        assert len(orientations) == 3

    def test_single_structure(self):
        """Test single structure is preserved."""
        orientations = [np.array([[0, 0, 0], [1, 0, 0]])]
        clean_duplicate_structure(orientations)
        assert len(orientations) == 1

    def test_empty_list(self):
        """Test empty list handling."""
        orientations = []
        clean_duplicate_structure(orientations)
        assert len(orientations) == 0


class TestCleanLabel:
    """Tests for the clean_label function."""

    def test_simple_label(self):
        """Test simple label cleaning."""
        assert clean_label("methane") == "methane"

    def test_label_with_spaces(self):
        """Test label with spaces."""
        result = clean_label("methyl chloride")
        assert " " not in result
        assert result == "methyl_chloride"

    def test_label_with_special_chars(self):
        """Test label with special characters."""
        result = clean_label("compound (A)")
        assert "(" not in result
        assert ")" not in result

    def test_label_with_prime(self):
        """Test label with prime character."""
        result = clean_label("C2'")
        assert result == "C2_prime"

    def test_label_with_star(self):
        """Test label with star character."""
        result = clean_label("pi*")
        assert result == "pi_star"

    def test_no_multiple_underscores(self):
        """Test that multiple underscores are collapsed."""
        result = clean_label("a   b")
        assert "__" not in result

    def test_no_leading_trailing_underscores(self):
        """Test no leading/trailing underscores."""
        result = clean_label(" test ")
        assert not result.startswith("_")
        assert not result.endswith("_")


class TestConvertStringIndicesToPymolIdIndices:
    """Tests for the convert_string_indices_to_pymol_id_indices function."""

    def test_single_index(self):
        """Test single index conversion."""
        result = convert_string_indices_to_pymol_id_indices("11")
        assert result == "id 11"

    def test_range_index(self):
        """Test range index conversion."""
        result = convert_string_indices_to_pymol_id_indices("1-10")
        assert result == "id 1-10"

    def test_multiple_indices(self):
        """Test multiple indices conversion."""
        result = convert_string_indices_to_pymol_id_indices("1-10,11")
        assert result == "id 1-10 or id 11"

    def test_complex_indices(self):
        """Test complex index string."""
        result = convert_string_indices_to_pymol_id_indices("1-10,11,14,19-30")
        assert result == "id 1-10 or id 11 or id 14 or id 19-30"

    def test_empty_raises_error(self):
        """Test empty string raises error."""
        with pytest.raises(ValueError):
            convert_string_indices_to_pymol_id_indices("")

    def test_whitespace_only_raises_error(self):
        """Test whitespace-only string raises error."""
        with pytest.raises(ValueError):
            convert_string_indices_to_pymol_id_indices("   ")


class TestLoadMoleculesFromPaths:
    """Tests for the load_molecules_from_paths function."""

    def test_load_with_none_index(self, gaussian_singlet_opt_outfile):
        """Test loading molecules with None index defaults to '-1'."""
        # This test ensures that index=None doesn't cause a TypeError
        molecules = load_molecules_from_paths(
            [gaussian_singlet_opt_outfile],
            index=None,
            add_index_suffix_for_single=False,
            check_exists=True,
        )
        assert len(molecules) > 0
        assert all(mol is not None for mol in molecules)

    def test_load_with_explicit_index(self, gaussian_singlet_opt_outfile):
        """Test loading molecules with explicit index."""
        molecules = load_molecules_from_paths(
            [gaussian_singlet_opt_outfile],
            index="-1",
            add_index_suffix_for_single=False,
            check_exists=True,
        )
        assert len(molecules) > 0
        assert all(mol is not None for mol in molecules)


class TestSelectItemsByIndex:
    """Tests for the select_items_by_index function."""

    def setup_method(self):
        """Set up test fixtures."""
        from chemsmart.utils.io import select_items_by_index

        self.select_fn = select_items_by_index
        self.test_list = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"]

    # ---- Tests for None or ":" (return all items) ----
    def test_none_returns_all(self):
        """Test that None index_spec returns all items."""
        result = self.select_fn(self.test_list, None)
        assert result == self.test_list

    def test_colon_returns_all(self):
        """Test that ':' index_spec returns all items."""
        result = self.select_fn(self.test_list, ":")
        assert result == self.test_list

    # ---- Tests for int type (0-based Python indexing) ----
    def test_int_positive_index(self):
        """Test direct int index (0-based)."""
        result = self.select_fn(self.test_list, 0)
        assert result == ["a"]

    def test_int_middle_index(self):
        """Test direct int index in middle (0-based)."""
        result = self.select_fn(self.test_list, 4)
        assert result == ["e"]

    def test_int_negative_index(self):
        """Test direct int negative index (0-based)."""
        result = self.select_fn(self.test_list, -1)
        assert result == ["j"]

    def test_int_negative_index_second_last(self):
        """Test direct int negative index for second last (0-based)."""
        result = self.select_fn(self.test_list, -2)
        assert result == ["i"]

    # ---- Tests for slice type (0-based Python indexing) ----
    def test_slice_from_start(self):
        """Test direct slice from start."""
        result = self.select_fn(self.test_list, slice(None, 5))
        assert result == ["a", "b", "c", "d", "e"]

    def test_slice_from_index_to_end(self):
        """Test direct slice from index to end."""
        result = self.select_fn(self.test_list, slice(4, None))
        assert result == ["e", "f", "g", "h", "i", "j"]

    def test_slice_with_step(self):
        """Test direct slice with step."""
        result = self.select_fn(self.test_list, slice(None, None, 2))
        assert result == ["a", "c", "e", "g", "i"]

    def test_slice_range(self):
        """Test direct slice range."""
        result = self.select_fn(self.test_list, slice(2, 6))
        assert result == ["c", "d", "e", "f"]

    def test_slice_with_negative_indices(self):
        """Test direct slice with negative indices."""
        result = self.select_fn(self.test_list, slice(-3, None))
        assert result == ["h", "i", "j"]

    def test_int_out_of_range_raises(self):
        """Test that out-of-range int index raises IndexError."""
        with pytest.raises(IndexError):
            self.select_fn(self.test_list, 100)

    def test_int_negative_out_of_range_raises(self):
        """Test that out-of-range negative int index raises IndexError."""
        with pytest.raises(IndexError):
            self.select_fn(self.test_list, -100)

    def test_slice_empty_result(self):
        """Test that slice returning empty list works correctly."""
        result = self.select_fn(self.test_list, slice(5, 5))
        assert result == []

    def test_slice_beyond_range_returns_partial(self):
        """Test that slice beyond range returns partial results."""
        result = self.select_fn(self.test_list, slice(8, 20))
        assert result == ["i", "j"]

    def test_slice_completely_out_of_range_returns_empty(self):
        """Test that slice completely out of range returns empty list."""
        result = self.select_fn(self.test_list, slice(20, 30))
        assert result == []

    # ---- Tests for single index ----
    def test_single_positive_index_middle(self):
        """Test single positive index in the middle."""
        result = self.select_fn(self.test_list, "5")
        assert result == ["e"]

    def test_single_negative_index(self):
        """Test single negative index (-1 = last)."""
        result = self.select_fn(self.test_list, "-1")
        assert result == ["j"]

    def test_single_negative_index_second_last(self):
        """Test single negative index (-2 = second last)."""
        result = self.select_fn(self.test_list, "-2")
        assert result == ["i"]

    # ---- Tests for comma-separated indices ----
    def test_comma_separated_indices(self):
        """Test comma-separated specific indices (1-based)."""
        result = self.select_fn(self.test_list, "1,3,5")
        assert result == ["a", "c", "e"]

    def test_comma_separated_with_negative(self):
        """Test comma-separated indices including negative."""
        result = self.select_fn(self.test_list, "1,-1")
        assert result == ["a", "j"]

    def test_comma_separated_multiple_negatives(self):
        """Test comma-separated multiple negative indices."""
        result = self.select_fn(self.test_list, "-3,-2,-1")
        assert result == ["h", "i", "j"]

    # ---- Tests for hyphen-based ranges ----
    def test_hyphen_range_middle(self):
        """Test hyphen range in the middle of list."""
        result = self.select_fn(self.test_list, "3-6")
        assert result == ["c", "d", "e", "f"]

    # ---- Tests for mixed comma and hyphen ----
    def test_mixed_ranges_and_singles(self):
        """Test mixed comma-separated singles and ranges."""
        result = self.select_fn(self.test_list, "1-3,5,7-9")
        assert result == ["a", "b", "c", "e", "g", "h", "i"]

    def test_mixed_with_negative(self):
        """Test mixed ranges with negative indices."""
        result = self.select_fn(self.test_list, "1,3,-1")
        assert result == ["a", "c", "j"]

    # ---- Tests for colon-based slicing (ASE-style) ----
    def test_slice_from_start(self):
        """Test slice from start to specific index."""
        result = self.select_fn(self.test_list, ":5")
        assert result == ["a", "b", "c", "d"]

    def test_slice_from_index_to_end(self):
        """Test slice from specific index to end."""
        result = self.select_fn(self.test_list, "5:")
        assert result == ["e", "f", "g", "h", "i", "j"]

    def test_slice_with_step(self):
        """Test slice with step."""
        result = self.select_fn(self.test_list, "::2")
        assert result == ["a", "c", "e", "g", "i"]

    def test_slice_from3to9_with_step(self):
        """Test slice with step."""
        result = self.select_fn(self.test_list, "3:9:2")
        assert result == ["c", "e", "g"]

    def test_slice_range_exclusive_end(self):
        """Test slice range has exclusive end (1:5 -> items 1-4)."""
        result = self.select_fn(self.test_list, "1:5")
        assert result == ["a", "b", "c", "d"]

    # ---- Tests for allow_duplicates ----
    def test_allow_duplicates_true(self):
        """Test that duplicates are allowed when allow_duplicates=True."""
        result = self.select_fn(self.test_list, "1,1,2", allow_duplicates=True)
        assert result == ["a", "a", "b"]

    def test_allow_duplicates_false_raises(self):
        """Test that allow_duplicates=False raises on duplicates."""
        with pytest.raises(ValueError):
            self.select_fn(self.test_list, "1,1,2", allow_duplicates=False)

    def test_allow_duplicates_false_negative_positive_same(self):
        """Test that -10 and 1 (same item) raises when duplicates not allowed."""
        with pytest.raises(ValueError):
            self.select_fn(self.test_list, "1,-10", allow_duplicates=False)

    # ---- Tests for allow_out_of_range ----
    def test_allow_out_of_range_true_filters(self):
        """Test that out of range indices are filtered when allowed."""
        result = self.select_fn(
            self.test_list, "1,15", allow_out_of_range=True
        )
        # 15 is out of range for 10 items, should be filtered
        assert result == ["a"]

    def test_allow_out_of_range_false_raises(self):
        """Test that allow_out_of_range=False raises on out of bounds."""
        with pytest.raises(ValueError):
            self.select_fn(self.test_list, "1,9-11", allow_out_of_range=False)

    # ---- Tests for invalid inputs ----
    def test_zero_index_raises(self):
        """Test that index 0 raises ValueError (1-based indexing)."""
        with pytest.raises(ValueError):
            self.select_fn(self.test_list, "0")

    def test_zero_in_range_raises(self):
        """Test that 0 in range raises ValueError."""
        with pytest.raises(ValueError):
            self.select_fn(self.test_list, "0-5")

    # ---- Tests for empty list ----
    def test_empty_list_with_none(self):
        """Test empty list with None returns empty."""
        result = self.select_fn([], None)
        assert result == []

    def test_empty_list_with_colon(self):
        """Test empty list with ':' returns empty."""
        result = self.select_fn([], ":")
        assert result == []

    # ---- Tests for single item list ----
    def test_single_item_list_select_first(self):
        """Test selecting from single item list."""
        result = self.select_fn(["only"], "1")
        assert result == ["only"]

    def test_single_item_list_select_last(self):
        """Test selecting last from single item list."""
        result = self.select_fn(["only"], "-1")
        assert result == ["only"]

    # ---- Test for bracket format ----
    def test_bracket_format(self):
        """Test that bracket format [1-3] works."""
        result = self.select_fn(self.test_list, "[1-3]")
        assert result == ["a", "b", "c"]

    # ---- Test return type is always list ----
    def test_return_type_is_list(self):
        """Test that return type is always a list."""
        assert isinstance(self.select_fn(self.test_list, None), list)
        assert isinstance(self.select_fn(self.test_list, "1"), list)
        assert isinstance(self.select_fn(self.test_list, "1-3"), list)
        assert isinstance(self.select_fn(self.test_list, ":"), list)
        assert isinstance(self.select_fn(self.test_list, "::2"), list)
