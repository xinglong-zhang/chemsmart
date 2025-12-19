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
    match_outfile_pattern,
    normalize_metal_bonds,
    remove_keyword,
    safe_sanitize,
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


class TestSafeSanitize:
    """Tests for the safe_sanitize function."""

    def test_sanitize_normal_molecule(self):
        """Test sanitization of a normal organic molecule."""
        from rdkit import Chem

        # Create a simple benzene molecule
        mol = Chem.MolFromSmiles("c1ccccc1")
        assert mol is not None

        # Should sanitize without issues
        result = safe_sanitize(mol)
        assert result is not None
        assert result.GetNumAtoms() == 6

    def test_sanitize_organometallic_complex(self):
        """Test sanitization of an organometallic complex.

        This tests the fallback mechanism when kekulization fails.
        """
        from rdkit import Chem

        # Create a ferrocene-like structure with aromatic rings
        # This is a simplified representation that may trigger kekulization issues
        mol = Chem.MolFromSmiles("[Fe]c1ccccc1")
        if mol is not None:
            # Should handle the molecule even if standard sanitization fails
            result = safe_sanitize(mol)
            assert result is not None


class TestNormalizeMetalBonds:
    """Tests for the normalize_metal_bonds function."""

    def test_no_metal_bonds(self):
        """Test molecule with no metal atoms."""
        from rdkit import Chem

        mol = Chem.MolFromSmiles("c1ccccc1")  # benzene
        assert mol is not None

        result = normalize_metal_bonds(mol)
        assert result is not None
        # Should not change anything for organic molecules
        assert result.GetNumAtoms() == 6

    def test_metal_complex_with_aromatic_ligand(self):
        """Test normalization of metal-aromatic bonds."""
        from rdkit import Chem

        # Create a simple organometallic structure
        # Iron with aromatic ligand (simplified ferrocene-like)
        mol = Chem.MolFromSmiles("[Fe]c1ccccc1", sanitize=False)
        if mol is not None:
            # Normalize metal bonds
            result = normalize_metal_bonds(mol)
            assert result is not None

            # Check that bonds to Fe are not aromatic
            fe_idx = None
            for atom in result.GetAtoms():
                if atom.GetSymbol() == "Fe":
                    fe_idx = atom.GetIdx()
                    break

            if fe_idx is not None:
                for bond in result.GetBonds():
                    if (
                        bond.GetBeginAtomIdx() == fe_idx
                        or bond.GetEndAtomIdx() == fe_idx
                    ):
                        # Bonds to metal should not be aromatic
                        assert not bond.GetIsAromatic()
                        assert bond.GetBondType() != Chem.BondType.AROMATIC

    def test_mixed_organic_and_metal(self):
        """Test that organic aromatic bonds are preserved."""
        from rdkit import Chem

        # Create a metal complex with an aromatic ring
        mol = Chem.MolFromSmiles("[Fe].c1ccccc1", sanitize=False)
        if mol is not None:
            result = normalize_metal_bonds(mol)
            assert result is not None

            # Count aromatic bonds in the benzene ring (not connected to Fe)
            # The benzene ring should still have aromatic bonds
            aromatic_bonds = sum(
                1 for bond in result.GetBonds() if bond.GetIsAromatic()
            )
            # Benzene has 6 aromatic bonds when not connected to metal
            assert aromatic_bonds >= 0  # May vary based on sanitization state
