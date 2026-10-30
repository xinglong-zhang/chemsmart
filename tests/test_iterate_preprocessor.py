"""Tests for the iterate preprocessors in ``chemsmart.jobs.iterate.iterate``:
:class:`BasePreprocessor` and its :class:`SkeletonPreprocessor` and
:class:`SubstituentPreprocessor` subclasses.
"""

import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.iterate.iterate import (
    BasePreprocessor,
    SkeletonPreprocessor,
    SubstituentPreprocessor,
)


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

    def test_ring_link_atom_has_no_severable_component(self):
        import numpy as np

        from chemsmart.io.molecules.structure import Molecule

        # A 3-membered carbon ring: removing any single bond to the
        # link atom never disconnects its neighbor (the other path
        # around the ring keeps everything in one component), so the
        # inner "which component contains the neighbor" search never
        # breaks for any neighbor.
        ring = Molecule(
            symbols=["C", "C", "C"],
            positions=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.51, 0.0, 0.0],
                    [0.755, 1.308, 0.0],
                ]
            ),
        )
        pre = BasePreprocessor(molecule=ring, link_index=3)
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

    def test_extract_by_indices_subsets_frozen_atoms(self):
        import numpy as np

        from chemsmart.io.molecules.structure import Molecule

        mol = Molecule(
            symbols=["C", "O", "H"],
            positions=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.4, 0.0, 0.0],
                    [-1.0, 0.0, 0.0],
                ]
            ),
            frozen_atoms=[1, 0, 1],
        )
        pre = BasePreprocessor(molecule=mol, link_index=1)
        extracted = pre._extract_by_indices([2, 0])
        assert extracted.chemical_symbols == ["C", "H"]
        assert extracted.frozen_atoms == [1, 1]


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


class TestSkeletonPreprocessorEarlyReturn:
    def test_unsaturated_link_atom_returns_molecule_unchanged(self):
        # A CH2 fragment: carbon has only 2 bonds (valence capacity 4),
        # so it has an available bonding position and no removal
        # should happen at all.
        ch2 = Molecule(
            symbols=["C", "H", "H"],
            positions=[
                [0.0, 0.0, 0.0],
                [1.09, 0.0, 0.0],
                [-0.36, 1.03, 0.0],
            ],
        )
        pre = SkeletonPreprocessor(molecule=ch2, link_index=1)
        result = pre.run()
        assert result is ch2
        assert pre._final_indices == [0, 1, 2]


class TestSkeletonPreprocessorAutoDetect:
    def test_no_skeleton_indices_falls_back_to_auto_detect(
        self, methanol_molecule
    ):
        pre = SkeletonPreprocessor(molecule=methanol_molecule, link_index=1)
        result = pre.run()
        # Auto-detect removes the smallest component (a single H).
        assert len(result) == 5


class TestSkeletonPreprocessorWithSkeletonIndices:
    def test_link_atom_not_in_skeleton_indices_raises(self, ethanol_molecule):
        # link_index=2 (C1, 1-based) but skeleton_indices only mentions
        # atom 1 (O) and 3 (C2), excluding C1 itself.
        pre = SkeletonPreprocessor(
            molecule=ethanol_molecule,
            link_index=2,
            skeleton_indices=[1, 3],
        )
        with pytest.raises(ValueError, match="must be included"):
            pre.run()

    def test_all_branches_are_skeleton_returns_unchanged(
        self, methanol_molecule
    ):
        # Every atom is declared part of the skeleton, so there are no
        # non-skeleton branches to remove.
        pre = SkeletonPreprocessor(
            molecule=methanol_molecule,
            link_index=1,
            skeleton_indices=[1, 2, 3, 4, 5, 6],
        )
        result = pre.run()
        assert result is methanol_molecule
        assert pre._final_indices == list(range(6))

    def test_single_non_skeleton_branch_removed(self, ethanol_molecule):
        # link atom = O (1-based 1); skeleton = O-C1-C2 backbone
        # (1-based [1, 2, 3]). The only non-skeleton branch reachable
        # from O (excluding the skeleton-bound C1 branch) is the lone
        # hydroxyl hydrogen (atom 9, 1-based / index 8, 0-based).
        pre = SkeletonPreprocessor(
            molecule=ethanol_molecule,
            link_index=1,
            skeleton_indices=[1, 2, 3],
        )
        result = pre.run()
        assert len(result) == 8
        assert pre._final_indices == [0, 1, 2, 3, 4, 5, 6, 7]

    def test_multiple_non_skeleton_branches_falls_back_to_auto_detect(
        self, ethanol_molecule
    ):
        # link atom = C1 (1-based 2); skeleton = C1-C2 only (1-based
        # [2, 3]). From C1, both the O-branch and each individual
        # methylene hydrogen (H3, H4) are separate non-skeleton
        # branches, forcing the ambiguous multi-branch fallback.
        pre = SkeletonPreprocessor(
            molecule=ethanol_molecule,
            link_index=2,
            skeleton_indices=[2, 3],
        )
        result = pre.run()
        # Auto-detect fallback removes the smallest single component,
        # which is one of the lone methylene hydrogens (H3 or H4).
        assert len(result) == 8


class TestFindNonSkeletonBranchesEdgeCases:
    def test_link_atom_with_no_neighbors_returns_empty_list(self):
        # An isolated carbon (far from everything else, no bonds).
        mol = Molecule(
            symbols=["C", "H", "C"],
            positions=[
                [0.0, 0.0, 0.0],
                [1.09, 0.0, 0.0],
                [20.0, 20.0, 20.0],
            ],
        )
        pre = SkeletonPreprocessor(
            molecule=mol, link_index=3, skeleton_indices=[3]
        )
        assert pre._find_non_skeleton_branches() == []

    def test_ring_branch_revisits_link_atom_mid_traversal(self):
        # A 3-membered carbon ring: DFS starting from either neighbor
        # of the link atom (atom 2, the ring's second carbon) traverses
        # around the ring and encounters the link atom again from the
        # opposite side, exercising the `node == excluded` skip for a
        # non-start node in _dfs_collect_branch.
        mol = Molecule(
            symbols=["C", "C", "C"],
            positions=[
                [0.0, 0.0, 0.0],
                [1.51, 0.0, 0.0],
                [0.755, 1.308, 0.0],
            ],
        )
        pre = SkeletonPreprocessor(
            molecule=mol, link_index=2, skeleton_indices=[2]
        )
        branches = pre._find_non_skeleton_branches()
        # Both other ring atoms are reachable from each other's branch,
        # since it's a ring; the link atom itself never appears.
        for branch in branches:
            assert 1 not in branch  # 0-based index of the link atom


class TestGetFallbackIndices:
    def test_returns_sorted_skeleton_indices_when_set(self, ethanol_molecule):
        pre = SkeletonPreprocessor(
            molecule=ethanol_molecule,
            link_index=1,
            skeleton_indices=[3, 1, 2],
        )
        assert pre._get_fallback_indices() == [0, 1, 2]

    def test_falls_back_to_auto_detect_when_not_set(self, methanol_molecule):
        pre = SkeletonPreprocessor(molecule=methanol_molecule, link_index=1)
        # No skeleton_indices provided, so this delegates to the base
        # class's auto-detect implementation.
        fallback = pre._get_fallback_indices()
        assert isinstance(fallback, list)
        assert fallback == sorted(fallback)


class TestSubstituentPreprocessorEarlyReturn:
    def test_unsaturated_link_atom_returns_unchanged(self):
        ch2 = Molecule(
            symbols=["C", "H", "H"],
            positions=[
                [0.0, 0.0, 0.0],
                [1.09, 0.0, 0.0],
                [-0.36, 1.03, 0.0],
            ],
        )
        pre = SubstituentPreprocessor(molecule=ch2, link_index=1)
        result = pre.run()
        assert result is ch2
        assert pre._final_indices == [0, 1, 2]

    def test_fallback_indices_when_unsaturated(self):
        ch2 = Molecule(
            symbols=["C", "H", "H"],
            positions=[
                [0.0, 0.0, 0.0],
                [1.09, 0.0, 0.0],
                [-0.36, 1.03, 0.0],
            ],
        )
        pre = SubstituentPreprocessor(molecule=ch2, link_index=1)
        assert pre._final_indices is None
        assert pre.get_new_link_index() == 1


class TestSubstituentPreprocessorAutoDetect:
    def test_saturated_link_atom_removes_smallest_group(
        self, methanol_molecule
    ):
        pre = SubstituentPreprocessor(molecule=methanol_molecule, link_index=1)
        result = pre.run()
        assert len(result) == 5
        assert pre._final_indices is not None

    def test_fallback_indices_when_saturated_and_run_not_called(
        self, methanol_molecule
    ):
        # link atom is saturated (no available bonding position), so
        # _get_fallback_indices() must delegate to the base class's
        # auto-detect implementation instead of the "unchanged" shortcut.
        pre = SubstituentPreprocessor(molecule=methanol_molecule, link_index=1)
        assert pre._final_indices is None
        fallback = pre._get_fallback_indices()
        assert isinstance(fallback, list)
        assert fallback == sorted(fallback)
        assert pre.get_new_link_index() == 1
