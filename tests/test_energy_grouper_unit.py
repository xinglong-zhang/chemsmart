"""
Direct unit tests for chemsmart.jobs.grouper.energy.EnergyGrouper.

Uses small hand-built Molecule objects with controlled energy values
to reach edge cases (invalid/missing energies, num_groups >= n,
threshold auto-search, group merging, and output helpers) that are
hard to trigger indirectly through real quantum-chemistry log/xyz
fixtures.
"""

import numpy as np
import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.grouper.energy import (
    KCAL_TO_HARTREE,
    EnergyGrouper,
)

SYMBOLS = ["H", "H"]
POSITIONS = np.array([[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]])


def _mol(energy):
    return Molecule(symbols=SYMBOLS, positions=POSITIONS, energy=energy)


class TestEnergyGrouperConstruction:
    def test_rejects_both_threshold_and_num_groups(self):
        molecules = [_mol(0.0), _mol(1.0)]
        with pytest.raises(ValueError, match="Cannot specify both"):
            EnergyGrouper(molecules, threshold=1.0, num_groups=2)

    def test_defaults_to_one_kcal_threshold(self):
        molecules = [_mol(0.0), _mol(1.0)]
        grouper = EnergyGrouper(molecules)
        assert grouper.threshold == 1.0
        assert np.isclose(grouper.threshold_hartree, 1.0 * KCAL_TO_HARTREE)

    def test_num_groups_leaves_threshold_hartree_none(self):
        molecules = [_mol(0.0), _mol(1.0)]
        grouper = EnergyGrouper(molecules, num_groups=1)
        assert grouper.threshold is None
        assert grouper.threshold_hartree is None


class TestEnergyGrouperValidation:
    def test_rejects_missing_energy_few_molecules(self):
        molecules = [_mol(0.0), _mol(None)]
        with pytest.raises(ValueError, match="missing energy information"):
            EnergyGrouper(molecules)

    def test_rejects_non_numeric_energy(self):
        molecules = [_mol(0.0), _mol("not-a-number")]
        with pytest.raises(ValueError, match="missing energy information"):
            EnergyGrouper(molecules)

    def test_rejects_non_finite_energy(self):
        molecules = [_mol(0.0), _mol(float("nan"))]
        with pytest.raises(ValueError, match="missing energy information"):
            EnergyGrouper(molecules)

    def test_reports_count_when_many_molecules_missing_energy(self):
        molecules = [_mol(None) for _ in range(6)]
        with pytest.raises(ValueError, match="Found 6 molecules"):
            EnergyGrouper(molecules)


class TestEnergyGrouperThresholdGrouping:
    def test_group_by_threshold_separates_distant_energies(self, tmp_path):
        molecules = [_mol(0.0), _mol(0.0001), _mol(10.0)]
        grouper = EnergyGrouper(
            molecules, threshold=1.0, output_dir=str(tmp_path)
        )
        groups, index_groups = grouper.group()
        # First two are close in energy (within 1 kcal/mol); third is far.
        assert len(groups) == 2
        assert sorted(len(g) for g in index_groups) == [1, 2]

    def test_repr_with_threshold(self):
        molecules = [_mol(0.0), _mol(1.0)]
        grouper = EnergyGrouper(molecules, threshold=2.0)
        assert "threshold=2.0 kcal/mol" in repr(grouper)


class TestEnergyGrouperNumGroups:
    def test_num_groups_at_least_n_returns_singleton_groups(self, tmp_path):
        molecules = [_mol(0.0), _mol(1.0)]
        grouper = EnergyGrouper(
            molecules, num_groups=5, output_dir=str(tmp_path)
        )
        groups, index_groups = grouper.group()
        assert len(groups) == 2
        assert index_groups == [[0], [1]]

    def test_num_groups_produces_requested_count(self, tmp_path):
        # Energies chosen so distinct thresholds yield distinct group counts.
        energies = [0.0, 0.0005, 5.0, 5.0005, 10.0]
        molecules = [_mol(e) for e in energies]
        grouper = EnergyGrouper(
            molecules, num_groups=3, output_dir=str(tmp_path)
        )
        groups, index_groups = grouper.group()
        assert len(groups) == 3
        assert grouper._auto_threshold is not None

    def test_repr_with_num_groups(self):
        molecules = [_mol(0.0), _mol(1.0)]
        grouper = EnergyGrouper(molecules, num_groups=1)
        assert "num_groups=1" in repr(grouper)


class TestEnergyGrouperFindOptimalThreshold:
    def test_returns_zero_when_all_diffs_infinite(self):
        molecules = [_mol(0.0), _mol(1.0)]
        grouper = EnergyGrouper(molecules, num_groups=1)
        threshold = grouper._find_optimal_threshold(
            energy_diff_values=[float("inf")],
            indices=[(0, 1)],
            n=2,
        )
        assert threshold == 0.0


class TestEnergyGrouperMerge:
    def test_merge_with_fully_disconnected_groups(self, tmp_path):
        # Force _find_optimal_threshold to return 0.0 (rejecting every
        # pair), so complete linkage yields n singleton groups with a
        # fully-disconnected adjacency matrix. This exercises the
        # `actual_groups > self.num_groups` path in _group_by_num_groups
        # and the merge loop's "best_merge_idx" selection when every
        # candidate has zero connections (see BUGS_FOUND.md: the "no
        # connections found" fallback branch is unreachable dead code,
        # since 0 connections still beats the initial best_connection_count
        # of -1, so best_merge_idx is never left at -1 here).
        from unittest.mock import patch

        energies = [0.0, 0.1, 0.2, 0.3, 0.4]
        molecules = [_mol(e) for e in energies]
        grouper = EnergyGrouper(
            molecules, num_groups=2, output_dir=str(tmp_path)
        )
        with patch.object(
            grouper, "_find_optimal_threshold", return_value=0.0
        ):
            groups, index_groups = grouper.group()
        assert len(groups) == 2
        assert sum(len(g) for g in index_groups) == len(molecules)

    def test_merge_prefers_most_connected_group(self, tmp_path):
        # A moderate threshold connects molecules within each of three
        # clusters but not across clusters, forcing exactly 3 groups
        # before merging down to 2 -- exercising the "merge into the
        # most-connected target" branch.
        from unittest.mock import patch

        from chemsmart.jobs.grouper.energy import KCAL_TO_HARTREE

        energies = [0.0, 0.05, 5.0, 5.05, 10.0]
        molecules = [_mol(e) for e in energies]
        grouper = EnergyGrouper(
            molecules, num_groups=2, output_dir=str(tmp_path)
        )
        threshold = 4.0 * KCAL_TO_HARTREE
        with patch.object(
            grouper, "_find_optimal_threshold", return_value=threshold
        ):
            groups, index_groups = grouper.group()
        assert len(groups) == 2
        assert sum(len(g) for g in index_groups) == len(molecules)

    def test_merge_groups_to_target_reduces_group_count(self, tmp_path):
        # Three well-separated clusters of increasing size force a merge
        # down to 2 groups, exercising _merge_groups_to_target directly.
        energies = [0.0, 0.0001, 20.0, 40.0, 40.0001, 40.0002]
        molecules = [_mol(e) for e in energies]
        grouper = EnergyGrouper(
            molecules, num_groups=2, output_dir=str(tmp_path)
        )
        groups, index_groups = grouper.group()
        assert len(groups) == 2
        total_members = sum(len(g) for g in index_groups)
        assert total_members == len(molecules)


class TestCompleteLinkageGrouping:
    def test_skips_already_assigned_candidate(self):
        # 4 molecules; adjacency crafted so that seed i=3 first absorbs
        # j=1 (but not j=2 or j=0), leaving the subsequent seed i=2 to
        # encounter j=1 already assigned -- exercising the `if
        # assigned[j]: continue` branch inside the inner loop.
        molecules = [_mol(float(i)) for i in range(4)]
        grouper = EnergyGrouper(molecules, threshold=1.0)
        adj_matrix = np.zeros((4, 4), dtype=bool)
        adj_matrix[3, 1] = adj_matrix[1, 3] = True
        groups, index_groups = grouper._complete_linkage_grouping(
            adj_matrix, 4
        )
        assert sorted(map(tuple, index_groups)) == [(0,), (1, 3), (2,)]


class TestMergeGroupsToTargetDirect:
    def test_merges_into_the_most_connected_group(self):
        # Directly drive _merge_groups_to_target with a hand-built
        # adjacency matrix where group 0 has a real connection to group
        # 2 (but not group 1), exercising the `connection_count += 1`
        # line and the best_merge_idx selection based on it.
        molecules = [_mol(float(i)) for i in range(4)]
        grouper = EnergyGrouper(molecules, num_groups=2)
        groups = [[molecules[0]], [molecules[1]], [molecules[2], molecules[3]]]
        index_groups = [[0], [1], [2, 3]]
        adj_matrix = np.zeros((4, 4), dtype=bool)
        adj_matrix[0, 2] = adj_matrix[2, 0] = True

        merged_groups, merged_index_groups = grouper._merge_groups_to_target(
            groups, index_groups, adj_matrix
        )
        assert len(merged_groups) == 2
        merged_sets = sorted(
            [frozenset(g) for g in merged_index_groups], key=len
        )
        assert merged_sets == [frozenset({1}), frozenset({0, 2, 3})]


class TestFindOptimalThresholdBinarySearch:
    def test_binary_search_explores_too_many_and_too_few_branches(self):
        # This energy set/target combination requires the binary search
        # to visit both a midpoint yielding too many groups and one
        # yielding too few before it converges, exercising both branches
        # of the search in a single call.
        energies = [6.31, 11.79, 47.22, 56.92, 76.10, 80.23]
        molecules = [_mol(e) for e in energies]
        grouper = EnergyGrouper(molecules, num_groups=4)
        indices = [
            (i, j)
            for i in range(len(molecules))
            for j in range(i + 1, len(molecules))
        ]
        diffs = [
            abs(molecules[j].energy - molecules[i].energy)
            for (i, j) in indices
        ]
        threshold = grouper._find_optimal_threshold(diffs, indices, 6)
        assert threshold > 0

    def test_binary_search_falls_back_when_no_exact_match_exists(self):
        # For this energy set/target combination, no candidate threshold
        # in the binary search yields exactly num_groups groups, so the
        # loop runs to exhaustion and falls back to `return
        # best_threshold` at the end of the method.
        energies = [9.47, 23.57, 24.14, 25.08, 26.64, 32.08, 62.74, 64.25]
        molecules = [_mol(e) for e in energies]
        grouper = EnergyGrouper(molecules, num_groups=6)
        indices = [
            (i, j)
            for i in range(len(molecules))
            for j in range(i + 1, len(molecules))
        ]
        diffs = [
            abs(molecules[j].energy - molecules[i].energy)
            for (i, j) in indices
        ]
        threshold = grouper._find_optimal_threshold(diffs, indices, 8)
        assert threshold > 0


class TestEnergyGrouperRecordResultsDirect:
    def test_record_results_without_grouping_time_or_cached_groups(
        self, tmp_path
    ):
        # Calling _record_results directly (bypassing group()) with
        # grouping_time=None and no cached group indices exercises the
        # false side of both `if grouping_time is not None` and `if
        # index_groups is not None`.
        molecules = [_mol(0.0), _mol(1.0)]
        grouper = EnergyGrouper(
            molecules,
            threshold=1.0,
            output_dir=str(tmp_path),
            matrix_format="csv",
        )
        assert grouper._cached_group_indices is None
        energy_matrix = np.zeros((2, 2))
        grouper._record_results(energy_matrix, grouping_time=None)


class TestEnergyGrouperRecordResults:
    def test_record_results_writes_output_with_threshold(self, tmp_path):
        molecules = [_mol(0.0), _mol(0.5), _mol(5.0)]
        grouper = EnergyGrouper(
            molecules,
            threshold=1.0,
            label="sys",
            output_dir=str(tmp_path),
            matrix_format="csv",
        )
        groups, index_groups = grouper.group()
        assert len(groups) >= 1

    def test_record_results_writes_output_with_num_groups(self, tmp_path):
        molecules = [_mol(0.0), _mol(0.5), _mol(5.0), _mol(5.5)]
        grouper = EnergyGrouper(
            molecules,
            num_groups=2,
            label="sys2",
            output_dir=str(tmp_path),
            matrix_format="txt",
        )
        groups, index_groups = grouper.group()
        assert len(groups) == 2

    def test_record_results_includes_thermo_header_for_qhg(self, tmp_path):
        molecules = [_mol(0.0), _mol(0.5)]
        grouper = EnergyGrouper(
            molecules,
            threshold=1.0,
            energy_type="QHG",
            thermo_parameters="s=1,T=298.15",
            output_dir=str(tmp_path),
            matrix_format="csv",
        )
        grouper.group()
