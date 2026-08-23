"""
Direct unit tests for :class:`SkeletonPreprocessor` in
``chemsmart.jobs.iterate.iterate``.

Uses the ``ethanol_molecule`` fixture (O-C-C backbone, index 0=O,
1=C(bonded to O), 2=C(methyl), 3-4=H on C1, 5-7=H on C2, 8=hydroxyl H)
to exercise the "user-provided skeleton_indices" branches, and a small
synthetic unsaturated CH2 fragment for the early-return
(already-has-bonding-capacity) path.
"""

import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.iterate.iterate import SkeletonPreprocessor


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
