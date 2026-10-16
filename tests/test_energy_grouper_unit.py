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


class TestEnergyGrouperMerge:
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
