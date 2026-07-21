"""
Direct unit tests for :class:`BasePreprocessor` in
``chemsmart.jobs.iterate.iterate`` — the shared substituent-detection
and atom-extraction logic used by ``SkeletonPreprocessor`` and
``SubstituentPreprocessor``.

Uses the ``methanol_molecule`` fixture (C-O-H with 3 more H on carbon)
already defined in conftest.py, whose real 3D geometry produces a
correctly bonded connectivity graph via ``Molecule.to_graph()``.
"""

import pytest

from chemsmart.jobs.iterate.iterate import BasePreprocessor


class TestMaxBondingCapacity:
    def test_carbon_max_valence(self):
        assert BasePreprocessor._get_max_bonding_capacity("C") == 4

    def test_oxygen_max_valence(self):
        assert BasePreprocessor._get_max_bonding_capacity("O") == 2

    def test_hydrogen_max_valence(self):
        assert BasePreprocessor._get_max_bonding_capacity("H") == 1

    def test_nitrogen_max_valence(self):
        assert BasePreprocessor._get_max_bonding_capacity("N") == 3


class TestHasAvailableBondingPosition:
    def test_saturated_carbon_has_no_capacity(self, methanol_molecule):
        # Carbon (index 0, 1-based link_index=1) already has 4 bonds
        # (O + 3H) in methanol.
        pre = BasePreprocessor(molecule=methanol_molecule, link_index=1)
        assert not pre._has_available_bonding_position()

    def test_oxygen_is_saturated_too(self, methanol_molecule):
        # Oxygen (index 1, 1-based link_index=2) has 2 bonds: C and
        # hydroxyl H, matching its max valence of 2.
        pre = BasePreprocessor(molecule=methanol_molecule, link_index=2)
        assert not pre._has_available_bonding_position()


class TestDetectSubstituent:
    def test_smallest_component_is_a_lone_hydrogen(self, methanol_molecule):
        # Breaking any single C-H bond isolates one H atom (the
        # smallest possible component), smaller than breaking the C-O
        # bond (which isolates the 2-atom O-H hydroxyl group).
        pre = BasePreprocessor(molecule=methanol_molecule, link_index=1)
        substituent = pre.detect_substituent()
        assert len(substituent) == 1
        removed_index = substituent[0]
        assert methanol_molecule.chemical_symbols[removed_index] == "H"

    def test_no_neighbors_returns_empty_list(self):
        from chemsmart.io.molecules.structure import Molecule

        # A single isolated atom has no bonds at all.
        lone_atom = Molecule(symbols=["Ar"], positions=[[0.0, 0.0, 0.0]])
        pre = BasePreprocessor(molecule=lone_atom, link_index=1)
        assert pre.detect_substituent() == []


class TestComplementAndExtraction:
    def test_get_complement_indices(self, methanol_molecule):
        pre = BasePreprocessor(molecule=methanol_molecule, link_index=1)
        complement = pre._get_complement_indices([0, 2])
        assert complement == [1, 3, 4, 5]

    def test_extract_by_indices_raises_for_empty(self, methanol_molecule):
        pre = BasePreprocessor(molecule=methanol_molecule, link_index=1)
        with pytest.raises(ValueError, match="no atoms"):
            pre._extract_by_indices([])

    def test_extract_by_indices_builds_subset_molecule(
        self, methanol_molecule
    ):
        pre = BasePreprocessor(molecule=methanol_molecule, link_index=1)
        # keep only C (0) and O (1)
        extracted = pre._extract_by_indices([1, 0])
        assert extracted.chemical_symbols == ["C", "O"]
        assert extracted.charge == methanol_molecule.charge
        assert extracted.multiplicity == methanol_molecule.multiplicity


class TestRunAutoDetectAndNewLinkIndex:
    def test_run_auto_detect_removes_smallest_component(
        self, methanol_molecule
    ):
        pre = BasePreprocessor(molecule=methanol_molecule, link_index=1)
        result = pre._run_auto_detect()
        # One H atom removed from the original 6-atom methanol.
        assert len(result) == 5
        assert pre._final_indices is not None

    def test_get_new_link_index_after_run(self, methanol_molecule):
        pre = BasePreprocessor(molecule=methanol_molecule, link_index=1)
        pre._run_auto_detect()
        # Carbon (original index 0) should still be present and be the
        # first atom (index 0) in the reduced/sorted molecule -> 1-based 1.
        assert pre.get_new_link_index() == 1

    def test_get_new_link_index_uses_fallback_when_run_not_called(
        self, methanol_molecule
    ):
        pre = BasePreprocessor(molecule=methanol_molecule, link_index=1)
        assert pre._final_indices is None
        # Should not raise, using _get_fallback_indices() internally.
        assert pre.get_new_link_index() == 1

    def test_get_new_link_index_raises_if_link_atom_removed(
        self, methanol_molecule
    ):
        # Use the hydroxyl H (index 2, 1-based link_index=3) as the link
        # atom, then manually simulate it being removed from the final
        # extracted indices.
        pre = BasePreprocessor(molecule=methanol_molecule, link_index=3)
        pre._final_indices = [0, 1, 3, 4, 5]  # excludes original index 2
        with pytest.raises(ValueError, match="was removed"):
            pre.get_new_link_index()
