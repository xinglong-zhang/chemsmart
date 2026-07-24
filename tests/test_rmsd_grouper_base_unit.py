"""
Direct unit tests for the shared RMSDGrouper base logic in
chemsmart.jobs.grouper.rmsd, exercised through BasicRMSDGrouper (the
simplest concrete subclass).

Complements tests/test_groupers.py (which exercises realistic
threshold-based grouping across all RMSD grouper variants) by targeting
num_groups-based grouping/merging, calculate_full_rmsd_matrix,
_record_results branches, and __repr__ -- edge cases that are hard to
reach through typical threshold-based CLI usage.
"""

import os
from unittest.mock import patch

import numpy as np
import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.grouper.rmsd import (
    BasicRMSDGrouper,
    HungarianRMSDGrouper,
    IRMSDGrouper,
    PymolRMSDGrouper,
    RMSDGrouperSharedMemory,
    SpyRMSDGrouper,
)

SYMBOLS = ["H", "H"]


def _mol(x_offset=0.0):
    positions = np.array(
        [[0.0 + x_offset, 0.0, 0.0], [0.74 + x_offset, 0.0, 0.0]]
    )
    return Molecule(symbols=SYMBOLS, positions=positions)


def _formaldehyde(x_offset=0.0, y_rotate=False):
    # One C, one O, two H -- exercises the Hungarian assignment's
    # single-atom-per-element shortcut (C and O) alongside the real
    # linear_sum_assignment path (the two H atoms).
    positions = np.array(
        [
            [0.0 + x_offset, 0.0, 0.0],
            [1.2 + x_offset, 0.0, 0.0],
            [-0.6 + x_offset, 0.9, 0.0],
            [-0.6 + x_offset, -0.9, 0.0],
        ]
    )
    if y_rotate:
        # Swap the two symmetric H atoms' order to force the Hungarian
        # algorithm to actually search for the optimal assignment.
        positions = positions[[0, 1, 3, 2]]
    return Molecule(symbols=["C", "O", "H", "H"], positions=positions)


class TestRMSDGrouperConstruction:
    def test_rejects_both_threshold_and_num_groups(self):
        with pytest.raises(ValueError, match="Cannot specify both"):
            BasicRMSDGrouper([_mol(), _mol()], threshold=0.5, num_groups=2)


class TestGetHeavyAtoms:
    def test_returns_all_atoms_when_not_ignoring_hydrogens(self):
        molecules = [_mol(), _mol()]
        grouper = BasicRMSDGrouper(molecules, ignore_hydrogens=False)
        positions, symbols = grouper._get_heavy_atoms(molecules[0])
        assert symbols == ["H", "H"]
        assert positions.shape == (2, 3)

    def test_filters_hydrogens_when_ignoring(self):
        mol = Molecule(
            symbols=["C", "H", "H"],
            positions=np.array(
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]]
            ),
        )
        grouper = BasicRMSDGrouper([mol, mol], ignore_hydrogens=True)
        positions, symbols = grouper._get_heavy_atoms(mol)
        assert symbols == ["C"]
        assert positions.shape == (1, 3)


class TestRMSDGrouperNumGroups:
    def test_num_groups_at_least_n_returns_singletons(self, tmp_path):
        molecules = [_mol(), _mol(1.0)]
        grouper = BasicRMSDGrouper(
            molecules, num_groups=5, output_dir=str(tmp_path)
        )
        groups, index_groups = grouper.group()
        assert index_groups == [[0], [1]]

    def test_num_groups_produces_requested_count(self, tmp_path):
        molecules = [
            _mol(0.0),
            _mol(0.001),
            _mol(5.0),
            _mol(5.001),
            _mol(10.0),
        ]
        grouper = BasicRMSDGrouper(
            molecules,
            num_groups=3,
            align_molecules=False,
            output_dir=str(tmp_path),
        )
        groups, index_groups = grouper.group()
        assert len(groups) == 3
        assert grouper._auto_threshold is not None

    def test_repr_with_num_groups(self):
        grouper = BasicRMSDGrouper([_mol(), _mol()], num_groups=1)
        assert "num_groups=1" in repr(grouper)

    def test_repr_with_threshold(self):
        grouper = BasicRMSDGrouper([_mol(), _mol()], threshold=0.3)
        assert "threshold=0.3" in repr(grouper)


class TestFindOptimalThreshold:
    def test_returns_zero_when_all_values_infinite(self):
        grouper = BasicRMSDGrouper([_mol(), _mol()], num_groups=1)
        threshold = grouper._find_optimal_threshold(
            rmsd_values=[float("inf")], indices=[(0, 1)], n=2
        )
        assert threshold == 0.0

    def test_falls_back_when_no_exact_match_exists(self):
        # Same crafted energy-like separations used for the analogous
        # EnergyGrouper/TorsionFingerprintGrouper tests: no candidate
        # threshold yields exactly num_groups groups, forcing the binary
        # search to exhaust and fall back to `return best_threshold`.
        molecules = [object() for _ in range(8)]
        grouper = BasicRMSDGrouper.__new__(BasicRMSDGrouper)
        grouper.molecules = molecules
        grouper.num_groups = 6
        values = [9.47, 23.57, 24.14, 25.08, 26.64, 32.08, 62.74, 64.25]
        indices = [
            (i, j)
            for i in range(len(values))
            for j in range(i + 1, len(values))
        ]
        rmsd_values = [abs(values[j] - values[i]) for (i, j) in indices]
        threshold = grouper._find_optimal_threshold(rmsd_values, indices, 8)
        assert threshold > 0


class TestMergeGroupsToTarget:
    def test_merge_reduces_group_count_and_prefers_connected_group(self):
        grouper = BasicRMSDGrouper.__new__(BasicRMSDGrouper)
        grouper.num_groups = 2
        groups = [["a"], ["b"], ["c", "d"]]
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

    def test_group_by_num_groups_triggers_merge(self, tmp_path):
        molecules = [_mol(float(i) * 0.001) for i in range(5)]
        grouper = BasicRMSDGrouper(
            molecules, num_groups=2, output_dir=str(tmp_path)
        )
        with patch.object(
            grouper, "_find_optimal_threshold", return_value=0.0
        ):
            groups, index_groups = grouper.group()
        assert len(groups) == 2
        assert sum(len(g) for g in index_groups) == 5


class TestCalculateFullRmsdMatrix:
    def test_returns_symmetric_matrix_without_output_file(self):
        molecules = [_mol(0.0), _mol(0.001), _mol(5.0)]
        grouper = BasicRMSDGrouper(molecules, num_procs=1)
        matrix = grouper.calculate_full_rmsd_matrix()
        assert matrix.shape == (3, 3)
        assert np.allclose(matrix, matrix.T)
        assert np.allclose(np.diag(matrix), 0.0)

    def test_writes_output_file_when_requested(self, tmp_path):
        molecules = [_mol(0.0), _mol(0.001)]
        grouper = BasicRMSDGrouper(
            molecules, num_procs=1, output_dir=str(tmp_path)
        )
        output_file = tmp_path / "nested" / "matrix_out"
        grouper.calculate_full_rmsd_matrix(output_file=str(output_file))
        # record() writes to the grouper's own output_dir/group_result
        # folder, not directly to output_file -- just confirm the parent
        # directory of output_file was created as instructed.
        assert output_file.parent.is_dir()


class TestRecordResultsBranches:
    def test_record_results_with_num_groups_and_grouping_time(self, tmp_path):
        molecules = [_mol(0.0), _mol(0.001), _mol(5.0)]
        grouper = BasicRMSDGrouper(
            molecules,
            num_groups=2,
            output_dir=str(tmp_path),
            matrix_format="csv",
        )
        groups, index_groups = grouper.group()
        assert len(groups) == 2

    def test_record_results_with_threshold_and_no_grouping_time(
        self, tmp_path
    ):
        molecules = [_mol(0.0), _mol(0.001)]
        grouper = BasicRMSDGrouper(
            molecules,
            threshold=0.5,
            output_dir=str(tmp_path),
            matrix_format="txt",
        )
        rmsd_matrix = np.zeros((2, 2))
        grouper._record_results(rmsd_matrix, grouping_time=None)


class TestHungarianRMSDGrouperSpecifics:
    def test_incompatible_symbols_return_inf(self):
        mol1 = _mol()
        mol2 = Molecule(
            symbols=["C", "H"],
            positions=np.array([[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]),
        )
        grouper = HungarianRMSDGrouper([mol1, mol2])
        assert grouper._calculate_rmsd((0, 1)) == np.inf

    def test_single_atom_per_element_shortcut_and_alignment(self):
        mol1 = _formaldehyde()
        mol2 = _formaldehyde(y_rotate=True)
        grouper = HungarianRMSDGrouper([mol1, mol2], align_molecules=True)
        rmsd = grouper._calculate_rmsd((0, 1))
        # C and O each appear once (shortcut path); the two equivalent
        # H atoms are optimally matched regardless of input order, so
        # RMSD after Kabsch alignment should be ~0.
        assert rmsd < 1e-6

    def test_without_alignment_skips_kabsch(self):
        mol1 = _formaldehyde()
        mol2 = _formaldehyde(x_offset=5.0)
        grouper = HungarianRMSDGrouper([mol1, mol2], align_molecules=False)
        rmsd = grouper._calculate_rmsd((0, 1))
        # Without alignment, a pure translation still contributes to RMSD.
        assert rmsd > 1.0


class TestRMSDGrouperSharedMemory:
    def test_group_with_alignment(self, tmp_path):
        molecules = [_mol(0.0), _mol(0.001), _mol(5.0)]
        grouper = RMSDGrouperSharedMemory(
            molecules,
            threshold=0.5,
            num_procs=1,
            align_molecules=True,
            output_dir=str(tmp_path),
        )
        groups, index_groups = grouper.group()
        assert sum(len(g) for g in index_groups) == 3
        assert grouper._cached_groups is not None

    def test_group_without_alignment(self, tmp_path):
        molecules = [_mol(0.0), _mol(0.001), _mol(5.0)]
        grouper = RMSDGrouperSharedMemory(
            molecules,
            threshold=0.5,
            num_procs=1,
            align_molecules=False,
            output_dir=str(tmp_path),
        )
        groups, index_groups = grouper.group()
        assert sum(len(g) for g in index_groups) == 3

    def test_calculate_rmsd_returns_inf_for_mismatched_shapes(self):
        # Call _init_worker/_calculate_rmsd directly in-process (rather
        # than through multiprocessing.Pool) so coverage.py can actually
        # trace them -- code that only runs inside pool worker
        # subprocesses isn't recorded without special coverage
        # configuration for subprocess tracing.
        from multiprocessing import RawArray

        import numpy as _np

        molecules = [_mol(), _mol()]
        grouper = RMSDGrouperSharedMemory(molecules, threshold=0.5)
        shared_pos = RawArray("d", 2 * 2 * 3)
        RMSDGrouperSharedMemory._init_worker(shared_pos, (2, 2, 3))

        import chemsmart.jobs.grouper.rmsd as rmsd_module

        rmsd_module.shared_positions = [
            _np.zeros((2, 3)),
            _np.zeros((3, 3)),
        ]
        assert grouper._calculate_rmsd((0, 1)) == np.inf

    def test_calculate_rmsd_with_alignment_in_process(self):
        from multiprocessing import RawArray

        molecules = [_mol(0.0), _mol(0.001)]
        grouper = RMSDGrouperSharedMemory(molecules, align_molecules=True)
        n, num_atoms = 2, 2
        shared_pos = RawArray("d", n * num_atoms * 3)
        pos_np = np.frombuffer(shared_pos, dtype=np.float64).reshape(
            n, num_atoms, 3
        )
        for i, mol in enumerate(molecules):
            pos_np[i] = mol.positions
        RMSDGrouperSharedMemory._init_worker(shared_pos, (n, num_atoms, 3))
        rmsd = grouper._calculate_rmsd((0, 1))
        assert rmsd < 0.01

    def test_calculate_rmsd_without_alignment_in_process(self):
        from multiprocessing import RawArray

        molecules = [_mol(0.0), _mol(5.0)]
        grouper = RMSDGrouperSharedMemory(molecules, align_molecules=False)
        n, num_atoms = 2, 2
        shared_pos = RawArray("d", n * num_atoms * 3)
        pos_np = np.frombuffer(shared_pos, dtype=np.float64).reshape(
            n, num_atoms, 3
        )
        for i, mol in enumerate(molecules):
            pos_np[i] = mol.positions
        RMSDGrouperSharedMemory._init_worker(shared_pos, (n, num_atoms, 3))
        rmsd = grouper._calculate_rmsd((0, 1))
        # Without alignment, the translation directly contributes to RMSD.
        assert rmsd > 1.0

    def test_complete_linkage_grouping_skips_already_assigned(self):
        molecules = [_mol() for _ in range(4)]
        grouper = RMSDGrouperSharedMemory(molecules, threshold=0.5)
        adj_matrix = np.zeros((4, 4), dtype=bool)
        adj_matrix[3, 1] = adj_matrix[1, 3] = True
        groups, index_groups = grouper._complete_linkage_grouping(
            adj_matrix, 4
        )
        assert sorted(map(tuple, index_groups)) == [(0,), (1, 3), (2,)]


class TestSpyRMSDGrouperSpecifics:
    def test_mismatched_coord_length_returns_inf(self):
        mol1 = _mol()
        mol2 = Molecule(
            symbols=["H", "H", "H"],
            positions=np.array(
                [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0], [1.5, 0.0, 0.0]]
            ),
        )
        grouper = SpyRMSDGrouper([mol1, mol2])
        assert grouper._calculate_rmsd((0, 1)) == np.inf
        assert grouper.best_isomorphisms[(0, 1)] is None

    def test_mismatched_atomic_numbers_returns_inf(self):
        mol1 = _mol()
        mol2 = Molecule(
            symbols=["C", "H"],
            positions=np.array([[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]),
        )
        grouper = SpyRMSDGrouper([mol1, mol2])
        assert grouper._calculate_rmsd((0, 1)) == np.inf
        assert grouper.best_isomorphisms[(0, 1)] is None

    def test_non_tuple_result_from_symmrmsd(self):
        grouper = SpyRMSDGrouper([_mol(), _mol(0.001)])
        with patch("spyrmsd.rmsd.symmrmsd", return_value=0.123):
            rmsd = grouper._calculate_rmsd((0, 1))
        assert rmsd == pytest.approx(0.123)
        assert grouper.best_isomorphisms[(0, 1)] is None

    def test_exception_during_symmrmsd_returns_inf(self):
        grouper = SpyRMSDGrouper([_mol(), _mol(0.001)])
        with patch(
            "spyrmsd.rmsd.symmrmsd",
            side_effect=RuntimeError("boom"),
        ):
            rmsd = grouper._calculate_rmsd((0, 1))
        assert rmsd == np.inf
        assert grouper.best_isomorphisms[(0, 1)] is None

    def test_get_best_isomorphism_direct_hit(self):
        grouper = SpyRMSDGrouper([_mol(), _mol(0.001)])
        grouper.best_isomorphisms[(0, 1)] = ([0, 1], [1, 0])
        assert grouper.get_best_isomorphism(0, 1) == ([0, 1], [1, 0])

    def test_get_best_isomorphism_reversed_hit_with_mapping(self):
        grouper = SpyRMSDGrouper([_mol(), _mol(0.001)])
        grouper.best_isomorphisms[(1, 0)] = ([0, 1], [1, 0])
        assert grouper.get_best_isomorphism(0, 1) == ([1, 0], [0, 1])

    def test_get_best_isomorphism_reversed_hit_without_mapping(self):
        grouper = SpyRMSDGrouper([_mol(), _mol(0.001)])
        grouper.best_isomorphisms[(1, 0)] = None
        assert grouper.get_best_isomorphism(0, 1) is None

    def test_get_best_isomorphism_not_found(self):
        grouper = SpyRMSDGrouper([_mol(), _mol(0.001)])
        assert grouper.get_best_isomorphism(5, 6) is None


def _cleanup_pymol_grouper(grouper):
    """Avoid __del__ triggering a real cmd.quit() during teardown,
    mirroring the pattern used in tests/test_groupers.py."""
    import shutil

    if getattr(grouper, "_temp_dir", None):
        shutil.rmtree(grouper._temp_dir, ignore_errors=True)
        grouper._temp_dir = None
    grouper.cmd = None


@pytest.fixture()
def shared_pymol_grouper():
    # A single real PymolRMSDGrouper reused across several tests below.
    # Constructing PyMOL's embedded, process-wide session is expensive
    # and (empirically) prone to hangs if done many times back-to-back
    # in the same process alongside tests/test_groupers.py's own PyMOL
    # tests -- so we build it once per test instead of once per test
    # method.
    grouper = PymolRMSDGrouper([_mol(), _mol(0.001)])
    yield grouper
    _cleanup_pymol_grouper(grouper)


class TestPymolRMSDGrouperSpecifics:
    def test_calculate_rmsd_negative_result_becomes_inf(
        self, shared_pymol_grouper
    ):
        with patch.object(
            shared_pymol_grouper.cmd, "align", return_value=[-1.0]
        ):
            assert shared_pymol_grouper._calculate_rmsd((0, 1)) == float("inf")

    def test_calculate_rmsd_non_finite_result_becomes_inf(
        self, shared_pymol_grouper
    ):
        with patch.object(
            shared_pymol_grouper.cmd, "align", return_value=[float("nan")]
        ):
            assert shared_pymol_grouper._calculate_rmsd((0, 1)) == float("inf")

    def test_calculate_rmsd_none_result_becomes_inf(
        self, shared_pymol_grouper
    ):
        with patch.object(
            shared_pymol_grouper.cmd, "align", return_value=None
        ):
            assert shared_pymol_grouper._calculate_rmsd((0, 1)) == float("inf")

    def test_calculate_rmsd_exception_becomes_inf(self, shared_pymol_grouper):
        with patch.object(
            shared_pymol_grouper.cmd,
            "align",
            side_effect=RuntimeError("boom"),
        ):
            assert shared_pymol_grouper._calculate_rmsd((0, 1)) == float("inf")

    def test_calculate_rmsd_uses_cache_on_second_call(
        self, shared_pymol_grouper
    ):
        with patch.object(
            shared_pymol_grouper.cmd, "align", return_value=[0.5]
        ) as mock_align:
            first = shared_pymol_grouper._calculate_rmsd((0, 1))
            second = shared_pymol_grouper._calculate_rmsd((1, 0))
        assert first == second == 0.5
        mock_align.assert_called_once()

    def test_repr_with_threshold(self, shared_pymol_grouper):
        assert "threshold=0.5" in repr(shared_pymol_grouper)

    def test_repr_with_num_groups(self):
        # This is the one repr variant that needs a distinct instance
        # (num_groups vs. threshold is fixed at construction time).
        grouper = PymolRMSDGrouper(
            [_mol(), _mol(0.001)], threshold=None, num_groups=1
        )
        try:
            assert "num_groups=1" in repr(grouper)
        finally:
            _cleanup_pymol_grouper(grouper)

    def test_del_quits_pymol_and_removes_temp_dir(self):
        # NOTE: __del__ is called explicitly here for coverage, but the
        # real `cmd` object is a process-wide singleton PyMOL session
        # shared by every PymolRMSDGrouper in this process. If we left
        # `grouper.cmd` pointing at the real session, Python's normal
        # garbage collection would call __del__ again later (unmocked),
        # issuing a *real* cmd.quit() that kills the shared session for
        # every subsequent test. So we always null out cmd/_temp_dir
        # afterward, same as _cleanup_pymol_grouper does elsewhere.
        grouper = PymolRMSDGrouper([_mol(), _mol(0.001)])
        temp_dir = grouper._temp_dir
        assert os.path.isdir(temp_dir)
        try:
            with patch.object(grouper.cmd, "quit") as mock_quit:
                grouper.__del__()
                mock_quit.assert_called_once()
            assert not os.path.isdir(temp_dir)
        finally:
            _cleanup_pymol_grouper(grouper)

    def test_del_swallows_quit_exception(self):
        grouper = PymolRMSDGrouper([_mol(), _mol(0.001)])
        try:
            with patch.object(
                grouper.cmd, "quit", side_effect=RuntimeError("boom")
            ):
                grouper.__del__()  # Should not raise
        finally:
            _cleanup_pymol_grouper(grouper)

    def test_del_swallows_rmtree_exception(self):
        grouper = PymolRMSDGrouper([_mol(), _mol(0.001)])
        try:
            # Also mock quit() here: __del__ calls cmd.quit() before the
            # rmtree cleanup, and a real quit() tears down PyMOL's
            # process-wide embedded session for good, hanging every
            # later test (including in test_groupers.py) that tries to
            # launch a new PymolRMSDGrouper.
            with (
                patch.object(grouper.cmd, "quit"),
                patch("shutil.rmtree", side_effect=RuntimeError("boom")),
            ):
                grouper.__del__()  # Should not raise
        finally:
            _cleanup_pymol_grouper(grouper)


def _make_irmsd_grouper(**kwargs):
    # find_irmsd_command() is called at __init__ time and raises
    # RuntimeError if no real irmsd binary is found. Mocking it lets us
    # exercise the rest of IRMSDGrouper's logic (XYZ writing, output
    # parsing, subprocess error handling, and _record_results) without
    # needing the actual external tool installed.
    with patch(
        "chemsmart.jobs.grouper.rmsd.find_irmsd_command",
        return_value="/fake/bin/irmsd",
    ):
        return IRMSDGrouper([_mol(), _mol(0.001)], **kwargs)


class TestIRMSDGrouperConstruction:
    def test_raises_when_irmsd_command_not_found(self):
        with patch(
            "chemsmart.jobs.grouper.rmsd.find_irmsd_command",
            return_value=None,
        ):
            with pytest.raises(RuntimeError, match="irmsd command not found"):
                IRMSDGrouper([_mol(), _mol(0.001)])

    def test_explicit_threshold_is_not_overridden(self):
        grouper = _make_irmsd_grouper(threshold=0.3)
        assert grouper.threshold == 0.3


class TestIRMSDGrouperSpecifics:
    def test_incompatible_symbols_return_inf(self):
        grouper = _make_irmsd_grouper()
        mol2 = Molecule(
            symbols=["C", "H"],
            positions=np.array([[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]]),
        )
        grouper.molecules.append(mol2)
        assert grouper._calculate_rmsd((0, 2)) == np.inf

    def test_successful_run_parses_rmsd_and_inversion(self):
        grouper = _make_irmsd_grouper(inversion="auto")
        fake_stdout = "some header\niRMSD: 0.4567\nInversion check: off\n"
        with patch(
            "subprocess.run",
            return_value=type(
                "Result",
                (),
                {"returncode": 0, "stdout": fake_stdout, "stderr": ""},
            )(),
        ):
            rmsd = grouper._calculate_rmsd((0, 1))
        assert rmsd == pytest.approx(0.4567)
        assert grouper._actual_inversion == "off"

    def test_nonzero_returncode_returns_inf(self):
        grouper = _make_irmsd_grouper()
        with patch(
            "subprocess.run",
            return_value=type(
                "Result",
                (),
                {"returncode": 1, "stdout": "", "stderr": "boom"},
            )(),
        ):
            assert grouper._calculate_rmsd((0, 1)) == np.inf

    def test_timeout_returns_inf(self):
        import subprocess as subprocess_module

        grouper = _make_irmsd_grouper()
        with patch(
            "subprocess.run",
            side_effect=subprocess_module.TimeoutExpired(
                cmd="irmsd", timeout=60
            ),
        ):
            assert grouper._calculate_rmsd((0, 1)) == np.inf

    def test_unexpected_exception_returns_inf(self):
        grouper = _make_irmsd_grouper()
        with patch("subprocess.run", side_effect=RuntimeError("boom")):
            assert grouper._calculate_rmsd((0, 1)) == np.inf

    def test_temp_file_cleanup_swallows_unlink_error(self):
        grouper = _make_irmsd_grouper()
        with (
            patch("subprocess.run", side_effect=RuntimeError("boom")),
            patch("os.unlink", side_effect=OSError("already gone")),
        ):
            # Should not raise despite os.unlink failing during cleanup.
            assert grouper._calculate_rmsd((0, 1)) == np.inf

    def test_ignore_hydrogens_adds_heavy_flag(self):
        grouper = _make_irmsd_grouper(ignore_hydrogens=True)
        captured_cmd = {}

        def fake_run(cmd, **kwargs):
            captured_cmd["cmd"] = cmd
            return type(
                "Result",
                (),
                {"returncode": 0, "stdout": "iRMSD: 0.1\n", "stderr": ""},
            )()

        with patch("subprocess.run", side_effect=fake_run):
            grouper._calculate_rmsd((0, 1))
        assert "--heavy" in captured_cmd["cmd"]

    def test_non_standard_inversion_omits_inversion_flag(self):
        # inversion is only normalized with .lower() in __init__, not
        # validated against {"on", "off", "auto"} -- so a nonstandard
        # value is kept as-is and skips the --inversion flag entirely.
        grouper = _make_irmsd_grouper(inversion="weird")
        captured_cmd = {}

        def fake_run(cmd, **kwargs):
            captured_cmd["cmd"] = cmd
            return type(
                "Result",
                (),
                {"returncode": 0, "stdout": "iRMSD: 0.1\n", "stderr": ""},
            )()

        with patch("subprocess.run", side_effect=fake_run):
            grouper._calculate_rmsd((0, 1))
        assert "--inversion" not in captured_cmd["cmd"]

    def test_parse_irmsd_output_handles_unparsable_rmsd_line(self):
        grouper = _make_irmsd_grouper()
        rmsd, inversion = grouper._parse_irmsd_output(
            "iRMSD: not-a-number\n", parse_inversion=True
        )
        assert rmsd == np.inf
        assert inversion is None

    def test_parse_irmsd_output_skips_inversion_when_not_requested(self):
        grouper = _make_irmsd_grouper()
        rmsd, inversion = grouper._parse_irmsd_output(
            "iRMSD: 0.2\nInversion check: on\n", parse_inversion=False
        )
        assert rmsd == pytest.approx(0.2)
        assert inversion is None

    def test_record_results_includes_inversion_header(self, tmp_path):
        grouper = _make_irmsd_grouper(
            output_dir=str(tmp_path), matrix_format="csv"
        )
        grouper._actual_inversion = "auto"
        rmsd_matrix = np.zeros((2, 2))
        grouper._record_results(rmsd_matrix, grouping_time=1.0)
