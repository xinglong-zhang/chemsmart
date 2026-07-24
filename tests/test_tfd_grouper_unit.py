"""
Direct unit tests for chemsmart.jobs.grouper.tfd.TorsionFingerprintGrouper.

Complements tests/test_groupers.py::Test_TorsionFingerprint_grouper (which
exercises realistic threshold-based grouping) by targeting construction
edge cases, conformer-preparation failure paths, num_groups-based
grouping/merging, and output-recording branches that are hard to reach
through typical threshold-based usage.
"""

from unittest.mock import patch

import numpy as np
import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.xyz.xyzfile import XYZFile
from chemsmart.jobs.grouper.tfd import TorsionFingerprintGrouper

SYMBOLS = ["H", "H"]
POSITIONS = np.array([[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]])


def _mol():
    return Molecule(symbols=SYMBOLS, positions=POSITIONS)


class TestTorsionFingerprintGrouperConstruction:
    def test_rejects_both_threshold_and_num_groups(self):
        with pytest.raises(ValueError, match="Cannot specify both"):
            TorsionFingerprintGrouper(
                [_mol(), _mol()], threshold=0.1, num_groups=2
            )

    def test_rejects_invalid_max_dev(self):
        with pytest.raises(ValueError, match="max_dev must be either"):
            TorsionFingerprintGrouper([_mol(), _mol()], max_dev="bogus")

    def test_defaults_threshold_when_neither_specified(self):
        grouper = TorsionFingerprintGrouper([_mol(), _mol()])
        assert grouper.threshold == 0.1
        assert grouper.num_groups is None


class TestPrepareConformerMoleculeEdgeCases:
    def test_empty_molecules_list(self):
        grouper = TorsionFingerprintGrouper.__new__(TorsionFingerprintGrouper)
        grouper.molecules = []
        grouper.ignore_hydrogens = False
        grouper._prepare_conformer_molecule()
        assert grouper.rdkit_mol is None
        assert grouper.valid_conformer_ids == []

    def test_rdkit_mol_none_short_circuits(self):
        grouper = TorsionFingerprintGrouper.__new__(TorsionFingerprintGrouper)
        grouper.molecules = [_mol()]
        grouper.ignore_hydrogens = False
        with patch.object(Molecule, "to_rdkit", return_value=None):
            grouper._prepare_conformer_molecule()
        assert grouper.rdkit_mol is None
        assert grouper.valid_conformer_ids == []

    def test_add_conformer_failure_is_logged_and_skipped(self):
        grouper = TorsionFingerprintGrouper([_mol(), _mol()])
        # rdkit_mol was already prepared successfully in __init__; now
        # force a failure while re-preparing to hit the except branch
        # inside the per-molecule conformer-adding loop.
        with patch(
            "chemsmart.jobs.grouper.tfd.Chem.Conformer",
            side_effect=RuntimeError("boom"),
        ):
            grouper._prepare_conformer_molecule()
        assert grouper.valid_conformer_ids == []

    def test_ignore_hydrogens_removes_hs_successfully(self):
        grouper = TorsionFingerprintGrouper(
            [_mol(), _mol()], ignore_hydrogens=True
        )
        assert grouper.rdkit_mol is not None

    def test_ignore_hydrogens_falls_back_on_failure(self):
        with patch(
            "chemsmart.jobs.grouper.tfd.Chem.RemoveHs",
            side_effect=RuntimeError("cannot remove"),
        ):
            grouper = TorsionFingerprintGrouper(
                [_mol(), _mol()], ignore_hydrogens=True
            )
        assert grouper.rdkit_mol is not None


class TestCalculateTFDEdgeCases:
    def test_returns_inf_when_no_rdkit_mol(self):
        grouper = TorsionFingerprintGrouper([_mol(), _mol()])
        grouper.rdkit_mol = None
        assert grouper._calculate_tfd((0, 1)) == float("inf")

    def test_returns_inf_on_calculation_exception(self):
        grouper = TorsionFingerprintGrouper([_mol(), _mol()])
        with patch(
            "chemsmart.jobs.grouper.tfd.TorsionFingerprints"
            ".GetTFDBetweenConformers",
            side_effect=RuntimeError("rdkit failure"),
        ):
            assert grouper._calculate_tfd((0, 1)) == float("inf")


class TestGroupEdgeCases:
    def test_empty_molecules_returns_empty_groups(self):
        grouper = TorsionFingerprintGrouper.__new__(TorsionFingerprintGrouper)
        grouper.molecules = []
        groups, index_groups = grouper.group()
        assert groups == []
        assert index_groups == []

    def test_single_molecule_returns_single_group(self):
        grouper = TorsionFingerprintGrouper.__new__(TorsionFingerprintGrouper)
        grouper.molecules = [_mol()]
        groups, index_groups = grouper.group()
        assert groups == [[grouper.molecules[0]]]
        assert index_groups == [[0]]

    def test_no_valid_conformers_falls_back_to_singleton_groups(self):
        molecules = [_mol(), _mol(), _mol()]
        grouper = TorsionFingerprintGrouper(molecules)
        grouper.valid_conformer_ids = []
        groups, index_groups = grouper.group()
        assert index_groups == [[0], [1], [2]]


class TestNumGroupsGrouping:
    NUM_PROCS = 1

    def test_num_groups_at_least_n_returns_singletons(
        self, multiple_molecules_xyz_file, tmp_path
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)[:3]
        grouper = TorsionFingerprintGrouper(
            molecules,
            num_groups=10,
            num_procs=self.NUM_PROCS,
            output_dir=str(tmp_path),
        )
        groups, index_groups = grouper.group()
        assert index_groups == [[0], [1], [2]]

    def test_num_groups_produces_requested_count(
        self, multiple_molecules_xyz_file, tmp_path
    ):
        xyz_file = XYZFile(filename=multiple_molecules_xyz_file)
        molecules = xyz_file.get_molecules(index=":", return_list=True)
        grouper = TorsionFingerprintGrouper(
            molecules,
            num_groups=5,
            num_procs=self.NUM_PROCS,
            output_dir=str(tmp_path),
        )
        groups, index_groups = grouper.group()
        assert len(groups) == 5
        assert grouper._auto_threshold is not None

    def test_merge_groups_to_target_reduces_group_count(self):
        grouper = TorsionFingerprintGrouper([_mol(), _mol()])
        grouper.num_groups = 2
        groups = [
            [_mol()],
            [_mol()],
            [_mol(), _mol()],
        ]
        index_groups = [[0], [1], [2, 3]]
        merged_groups, merged_index_groups = grouper._merge_groups_to_target(
            groups, index_groups
        )
        assert len(merged_groups) == 2
        sizes = sorted(len(g) for g in merged_index_groups)
        assert sizes == [1, 3]

    def test_find_optimal_threshold_returns_zero_for_all_infinite(self):
        grouper = TorsionFingerprintGrouper([_mol(), _mol()])
        grouper.num_groups = 1
        threshold = grouper._find_optimal_tfd_threshold(
            tfd_values=[float("inf")], indices=[(0, 1)], n=2
        )
        assert threshold == 0.0

    def test_group_by_num_groups_triggers_merge_when_search_undershoots(
        self,
    ):
        # Force _find_optimal_tfd_threshold to return 0.0 (rejecting
        # every pair), so complete linkage yields n singleton groups --
        # more than the requested num_groups -- forcing
        # _group_by_num_groups to call _merge_groups_to_target.
        molecules = [_mol() for _ in range(5)]
        grouper = TorsionFingerprintGrouper.__new__(TorsionFingerprintGrouper)
        grouper.molecules = molecules
        grouper.num_groups = 2
        indices = [(i, j) for i in range(5) for j in range(i + 1, 5)]
        tfd_values = [0.1 + 0.01 * k for k in range(len(indices))]
        with patch.object(
            grouper, "_find_optimal_tfd_threshold", return_value=0.0
        ):
            groups, index_groups = grouper._group_by_num_groups(
                tfd_values, indices, 5
            )
        assert len(groups) == 2
        assert sum(len(g) for g in index_groups) == 5

    def test_find_optimal_threshold_falls_back_when_no_exact_match(self):
        # This TFD set/target combination has no candidate threshold in
        # the binary search that yields exactly num_groups groups, so
        # the search exhausts and falls back to `return best_threshold`.
        molecules = [object() for _ in range(8)]
        grouper = TorsionFingerprintGrouper.__new__(TorsionFingerprintGrouper)
        grouper.molecules = molecules
        grouper.num_groups = 6
        values = [9.47, 23.57, 24.14, 25.08, 26.64, 32.08, 62.74, 64.25]
        indices = [
            (i, j)
            for i in range(len(values))
            for j in range(i + 1, len(values))
        ]
        tfd_values = [abs(values[j] - values[i]) for (i, j) in indices]
        threshold = grouper._find_optimal_tfd_threshold(tfd_values, indices, 8)
        assert threshold > 0


class TestRecordResultsBranches:
    def test_record_results_with_num_groups_and_auto_threshold(self, tmp_path):
        molecules = [_mol(), _mol(), _mol()]
        grouper = TorsionFingerprintGrouper(
            molecules,
            num_groups=2,
            output_dir=str(tmp_path),
            matrix_format="csv",
        )
        grouper._auto_threshold = 0.05
        tfd_matrix = np.zeros((3, 3))
        grouper._record_results(
            tfd_matrix=tfd_matrix,
            grouping_time=1.23,
            groups=[[molecules[0], molecules[1]], [molecules[2]]],
            index_groups=[[0, 1], [2]],
        )

    def test_record_results_without_groups_or_grouping_time(self, tmp_path):
        molecules = [_mol(), _mol()]
        grouper = TorsionFingerprintGrouper(
            molecules,
            threshold=0.1,
            output_dir=str(tmp_path),
            matrix_format="csv",
        )
        tfd_matrix = np.zeros((2, 2))
        grouper._record_results(
            tfd_matrix=tfd_matrix,
            grouping_time=None,
            groups=None,
            index_groups=None,
        )


class TestRepr:
    def test_repr_with_threshold(self):
        grouper = TorsionFingerprintGrouper([_mol(), _mol()], threshold=0.2)
        text = repr(grouper)
        assert "threshold=0.2" in text
        assert "TorsionFingerprintGrouper" in text

    def test_repr_with_num_groups(self):
        grouper = TorsionFingerprintGrouper([_mol(), _mol()], num_groups=1)
        text = repr(grouper)
        assert "num_groups=1" in text
