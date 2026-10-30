"""Tests for :class:`IterateJobRunner` in ``chemsmart.jobs.iterate.runner``,
including combination generation and the multiprocessing ``run`` path.
"""

import os
import queue
from unittest.mock import MagicMock, patch

from chemsmart.jobs.iterate.job import IterateJob
from chemsmart.jobs.iterate.runner import (
    IterateCombination,
    IterateJobRunner,
    IterateMoleculePool,
    _run_combination_task,
    _run_combination_worker,
)
from chemsmart.jobs.iterate.settings import IterateJobSettings


class FakeQueue:
    """Stand-in for a ``multiprocessing.Manager().Queue()``.

    ``get_nowait_plan`` is a list consumed one entry per call to
    ``get_nowait``: ``None`` means "return the next queued item (or
    raise ``queue.Empty`` if none)", and any exception instance/class
    means "raise this instead". Once the plan is exhausted, falls back
    to the "return next item or raise Empty" behavior.
    """

    def __init__(self, get_nowait_plan=None):
        self._items = []
        self._plan = list(get_nowait_plan or [])

    def put(self, item):
        self._items.append(item)

    def get_nowait(self):
        if self._plan:
            effect = self._plan.pop(0)
            if effect is not None:
                raise effect
        if self._items:
            return self._items.pop(0)
        raise queue.Empty()


class FakeProcess:
    """Stand-in for a ``multiprocessing.Process``.

    ``alive`` controls every ``is_alive()`` call for this process
    (kept simple: a process is either always "still running" or
    already "finished", which is all the watchdog logic needs).
    """

    def __init__(self, target=None, args=(), daemon=False, alive=False):
        self.target = target
        self.args = args
        self.daemon = daemon
        self.pid = 4242
        self._alive = alive
        self.terminated = False
        self.killed = False
        self.joined_timeouts = []

    def start(self):
        pass

    def is_alive(self):
        return self._alive

    def join(self, timeout=None):
        self.joined_timeouts.append(timeout)

    def terminate(self):
        self.terminated = True

    def kill(self):
        self.killed = True


def _make_pool_and_combo(label_suffix="a"):
    skeleton = MagicMock()
    substituent = MagicMock()
    pool = IterateMoleculePool(
        skeletons=[skeleton], substituents=[substituent]
    )
    combo = IterateCombination(
        skeleton_idx=0,
        skeleton_label=f"skel{label_suffix}",
        skeleton_link_index=1,
        skeleton_indices=None,
        substituent_idx=0,
        substituent_label=f"sub{label_suffix}",
        substituent_link_index=1,
    )
    return pool, combo


def _patched_manager(fake_queue):
    fake_manager = MagicMock()
    fake_manager.Queue.return_value = fake_queue
    fake_manager.shutdown.return_value = None
    return fake_manager


class TestRunCombinationTask:
    """_run_combination_task is a plain module-level function (no
    multiprocessing needed to call it directly); its three
    collaborators are mocked so this exercises only the task's own
    success/failure/exception branches."""

    def _make_pool_and_combo(self):
        skeleton = MagicMock()
        skeleton.copy.return_value = skeleton
        substituent = MagicMock()
        substituent.copy.return_value = substituent
        pool = IterateMoleculePool(
            skeletons=[skeleton], substituents=[substituent]
        )
        combo = IterateCombination(
            skeleton_idx=0,
            skeleton_label="s",
            skeleton_link_index=1,
            skeleton_indices=None,
            substituent_idx=0,
            substituent_label="sub",
            substituent_link_index=1,
        )
        return pool, combo

    def test_success_returns_generated_molecule(self):
        pool, combo = self._make_pool_and_combo()
        with (
            patch(
                "chemsmart.jobs.iterate.runner.SkeletonPreprocessor"
            ) as mock_skel_pre,
            patch(
                "chemsmart.jobs.iterate.runner.SubstituentPreprocessor"
            ) as mock_sub_pre,
            patch(
                "chemsmart.jobs.iterate.runner.IterateAnalyzer"
            ) as mock_analyzer,
        ):
            mock_skel_pre.return_value.run.return_value = "processed_skel"
            mock_skel_pre.return_value.get_new_link_index.return_value = 1
            mock_sub_pre.return_value.run.return_value = "processed_sub"
            mock_sub_pre.return_value.get_new_link_index.return_value = 1
            mock_analyzer.return_value.run.return_value = "result-molecule"

            label, result = _run_combination_task(combo, pool)

        assert label == combo.label
        assert result == "result-molecule"

    def test_analyzer_returning_none_is_reported_as_failure(self):
        pool, combo = self._make_pool_and_combo()
        with (
            patch("chemsmart.jobs.iterate.runner.SkeletonPreprocessor"),
            patch("chemsmart.jobs.iterate.runner.SubstituentPreprocessor"),
            patch(
                "chemsmart.jobs.iterate.runner.IterateAnalyzer"
            ) as mock_analyzer,
        ):
            mock_analyzer.return_value.run.return_value = None
            label, result = _run_combination_task(combo, pool)

        assert label == combo.label
        assert result is None

    def test_exception_is_caught_and_returns_none(self):
        pool, combo = self._make_pool_and_combo()
        with patch(
            "chemsmart.jobs.iterate.runner.SkeletonPreprocessor",
            side_effect=RuntimeError("boom"),
        ):
            label, result = _run_combination_task(combo, pool)

        assert label == combo.label
        assert result is None


class TestRunCombinationWorker:
    def test_puts_task_result_on_queue(self):
        combo = IterateCombination(
            skeleton_idx=0,
            skeleton_label="s",
            skeleton_link_index=1,
            skeleton_indices=None,
            substituent_idx=0,
            substituent_label="sub",
            substituent_link_index=1,
        )
        pool = IterateMoleculePool()
        result_queue = MagicMock()

        with patch(
            "chemsmart.jobs.iterate.runner._run_combination_task",
            return_value=(combo.label, "mol"),
        ):
            _run_combination_worker(combo, pool, result_queue)

        result_queue.put.assert_called_once_with((combo.label, "mol"))

    def test_worker_panic_puts_none_result_on_queue(self):
        combo = IterateCombination(
            skeleton_idx=0,
            skeleton_label="s",
            skeleton_link_index=1,
            skeleton_indices=None,
            substituent_idx=0,
            substituent_label="sub",
            substituent_link_index=1,
        )
        pool = IterateMoleculePool()
        result_queue = MagicMock()

        with patch(
            "chemsmart.jobs.iterate.runner._run_combination_task",
            side_effect=RuntimeError("panic"),
        ):
            _run_combination_worker(combo, pool, result_queue)

        result_queue.put.assert_called_once_with((combo.label, None))


class TestIterateCombinationLabel:
    def test_label_format(self):
        combo = IterateCombination(
            skeleton_idx=0,
            skeleton_label="benzene",
            skeleton_link_index=5,
            skeleton_indices=None,
            substituent_idx=0,
            substituent_label="methyl",
            substituent_link_index=1,
        )
        assert combo.label == "benzene_5_methyl_1"


class TestIterateJobRunnerBasics:
    def test_executable_is_none(self):
        runner = IterateJobRunner()
        assert runner.executable is None

    def test_get_command_is_none(self):
        runner = IterateJobRunner()
        assert runner._get_command(job=MagicMock()) is None

    def test_defaults_scratch_to_false(self):
        assert IterateJobRunner.SCRATCH is False
        runner = IterateJobRunner()
        assert runner.scratch is False


class TestRunSingleFakeMode:
    def test_fake_mode_returns_none_without_running(self):
        runner = IterateJobRunner(fake=True)
        combo = IterateCombination(
            skeleton_idx=0,
            skeleton_label="s",
            skeleton_link_index=1,
            skeleton_indices=None,
            substituent_idx=0,
            substituent_label="sub",
            substituent_link_index=1,
        )
        pool = IterateMoleculePool()
        label, result = runner.run_single(combo, pool)
        assert label == combo.label
        assert result is None

    def test_non_fake_mode_delegates_to_task_runner(self):
        runner = IterateJobRunner(fake=False)
        combo = IterateCombination(
            skeleton_idx=0,
            skeleton_label="s",
            skeleton_link_index=1,
            skeleton_indices=None,
            substituent_idx=0,
            substituent_label="sub",
            substituent_link_index=1,
        )
        pool = IterateMoleculePool()
        with patch(
            "chemsmart.jobs.iterate.runner._run_combination_task",
            return_value=("s_1_sub_1", "mol-sentinel"),
        ) as mock_task:
            result = runner.run_single(combo, pool)

        mock_task.assert_called_once_with(combo, pool)
        assert result == ("s_1_sub_1", "mol-sentinel")


class TestRunCombinationsEarlyReturns:
    """run_combinations' fake-mode and empty-combinations early
    returns don't touch multiprocessing at all, unlike the rest of the
    method (intentionally not exercised here, see module docstring)."""

    def test_fake_mode_returns_none_results_without_running(self):
        runner = IterateJobRunner(fake=True)
        combo = IterateCombination(
            skeleton_idx=0,
            skeleton_label="s",
            skeleton_link_index=1,
            skeleton_indices=None,
            substituent_idx=0,
            substituent_label="sub",
            substituent_link_index=1,
        )
        result = runner.run_combinations(IterateMoleculePool(), [combo])
        assert result == [(combo.label, None)]

    def test_empty_combinations_returns_empty_list(self):
        runner = IterateJobRunner(fake=False)
        assert runner.run_combinations(IterateMoleculePool(), []) == []


class TestLoadMolecule:
    def test_returns_none_without_file_path(self):
        runner = IterateJobRunner()
        molecule, label = runner._load_molecule(
            {"label": "mymol"}, "skeleton", 0
        )
        assert molecule is None
        assert label == "mymol"

    def test_default_label_when_missing(self):
        runner = IterateJobRunner()
        molecule, label = runner._load_molecule({}, "skeleton", 2)
        assert molecule is None
        assert label == "skeleton3"

    def test_loads_molecule_from_file(self, single_molecule_xyz_file):
        runner = IterateJobRunner()
        molecule, label = runner._load_molecule(
            {"file_path": single_molecule_xyz_file, "label": "myskel"},
            "skeleton",
            0,
        )
        assert molecule is not None
        assert label == "myskel"

    def test_link_index_out_of_bounds_returns_none(
        self, single_molecule_xyz_file
    ):
        runner = IterateJobRunner()
        molecule, label = runner._load_molecule(
            {
                "file_path": single_molecule_xyz_file,
                "label": "myskel",
                "link_index": [99999],
            },
            "skeleton",
            0,
        )
        assert molecule is None
        assert label == "myskel"

    def test_skeleton_indices_out_of_bounds_returns_none(
        self, single_molecule_xyz_file
    ):
        runner = IterateJobRunner()
        molecule, label = runner._load_molecule(
            {
                "file_path": single_molecule_xyz_file,
                "label": "myskel",
                "skeleton_indices": [99999],
            },
            "skeleton",
            0,
        )
        assert molecule is None
        assert label == "myskel"

    def test_load_failure_returns_none(self):
        runner = IterateJobRunner()
        molecule, label = runner._load_molecule(
            {"file_path": "/no/such/file.xyz", "label": "bad"},
            "skeleton",
            0,
        )
        assert molecule is None
        assert label == "bad"

    def test_link_index_as_plain_int_is_normalized_to_list(
        self, single_molecule_xyz_file
    ):
        """link_index is normally a list (as set by the CLI), but the
        runner defensively wraps a bare int in a list too."""
        runner = IterateJobRunner()
        molecule, label = runner._load_molecule(
            {
                "file_path": single_molecule_xyz_file,
                "label": "myskel",
                "link_index": 1,
            },
            "skeleton",
            0,
        )
        assert molecule is not None
        assert label == "myskel"

    def test_skeleton_indices_as_plain_int_is_normalized_to_list(
        self, single_molecule_xyz_file
    ):
        runner = IterateJobRunner()
        molecule, label = runner._load_molecule(
            {
                "file_path": single_molecule_xyz_file,
                "label": "myskel",
                "skeleton_indices": 1,
            },
            "skeleton",
            0,
        )
        assert molecule is not None
        assert label == "myskel"

    def test_link_index_as_non_int_non_list_falsy_value_is_left_alone(
        self, single_molecule_xyz_file
    ):
        """A link_index that's neither a list nor an int (e.g. an empty
        string, if the CLI's own normalization were somehow skipped)
        isn't wrapped in a list; falsy values just skip the
        out-of-bounds check below."""
        runner = IterateJobRunner()
        molecule, label = runner._load_molecule(
            {
                "file_path": single_molecule_xyz_file,
                "label": "myskel",
                "link_index": "",
            },
            "skeleton",
            0,
        )
        assert molecule is not None
        assert label == "myskel"

    def test_skeleton_indices_as_non_int_non_list_falsy_value_is_left_alone(
        self, single_molecule_xyz_file
    ):
        runner = IterateJobRunner()
        molecule, label = runner._load_molecule(
            {
                "file_path": single_molecule_xyz_file,
                "label": "myskel",
                "skeleton_indices": "",
            },
            "skeleton",
            0,
        )
        assert molecule is not None
        assert label == "myskel"


class TestGenerateCombinations:
    def test_generates_cross_product_of_valid_configs(
        self, single_molecule_xyz_file
    ):
        settings = IterateJobSettings()
        settings.skeleton_list = [
            {
                "file_path": single_molecule_xyz_file,
                "label": "skel1",
                "link_index": [1, 2],
            }
        ]
        settings.substituent_list = [
            {
                "file_path": single_molecule_xyz_file,
                "label": "sub1",
                "link_index": [1],
            }
        ]
        job = IterateJob(settings=settings)
        runner = IterateJobRunner()

        pool, combinations = runner._generate_combinations(job)

        assert len(pool.skeletons) == 1
        assert len(pool.substituents) == 1
        # 2 skeleton link indices x 1 substituent -> 2 combinations
        assert len(combinations) == 2
        labels = {c.label for c in combinations}
        assert labels == {"skel1_1_sub1_1", "skel1_2_sub1_1"}

    def test_skips_skeleton_without_link_index(self, single_molecule_xyz_file):
        settings = IterateJobSettings()
        settings.skeleton_list = [
            {"file_path": single_molecule_xyz_file, "label": "skel1"}
        ]
        settings.substituent_list = [
            {
                "file_path": single_molecule_xyz_file,
                "label": "sub1",
                "link_index": [1],
            }
        ]
        job = IterateJob(settings=settings)
        runner = IterateJobRunner()

        pool, combinations = runner._generate_combinations(job)
        assert combinations == []

    def test_skips_substituent_without_link_index(
        self, single_molecule_xyz_file
    ):
        settings = IterateJobSettings()
        settings.skeleton_list = [
            {
                "file_path": single_molecule_xyz_file,
                "label": "skel1",
                "link_index": [1],
            }
        ]
        settings.substituent_list = [
            {"file_path": single_molecule_xyz_file, "label": "sub1"}
        ]
        job = IterateJob(settings=settings)
        runner = IterateJobRunner()

        pool, combinations = runner._generate_combinations(job)
        assert combinations == []

    def test_empty_lists_produce_no_combinations(self):
        settings = IterateJobSettings()
        job = IterateJob(settings=settings)
        runner = IterateJobRunner()

        pool, combinations = runner._generate_combinations(job)
        assert pool.skeletons == []
        assert pool.substituents == []
        assert combinations == []

    def test_skips_skeleton_and_substituent_that_fail_to_load(self):
        """A missing file makes _load_molecule return None, which
        _generate_combinations must skip rather than pool."""
        settings = IterateJobSettings()
        settings.skeleton_list = [
            {
                "file_path": "/no/such/file.xyz",
                "label": "skel1",
                "link_index": [1],
            }
        ]
        settings.substituent_list = [
            {
                "file_path": "/no/such/file2.xyz",
                "label": "sub1",
                "link_index": [1],
            }
        ]
        job = IterateJob(settings=settings)
        runner = IterateJobRunner()

        pool, combinations = runner._generate_combinations(job)
        assert pool.skeletons == []
        assert pool.substituents == []
        assert combinations == []


class TestWriteOutputs:
    def _make_molecule(self):
        from chemsmart.io.molecules.structure import Molecule

        return Molecule(
            symbols=["H", "H"],
            positions=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        )

    def test_merged_output_writes_all_successful_molecules(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = IterateJob(outputfile="results")
        runner = IterateJobRunner()
        mol = self._make_molecule()
        results = [("combo1", mol), ("combo2", None)]

        runner._write_outputs(results, job)

        assert os.path.exists("results.xyz")
        with open("results.xyz") as f:
            content = f.read()
        assert "combo1" in content
        assert content.count("2\n") == 1  # only one successful molecule

    def test_separate_outputs_writes_one_file_per_result(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = IterateJob(
            separate_outputs=True, output_directory=str(tmp_path / "out")
        )
        runner = IterateJobRunner()
        mol = self._make_molecule()
        results = [("combo1", mol), ("combo2", None)]

        runner._write_outputs(results, job)

        assert os.path.exists(str(tmp_path / "out" / "combo1.xyz"))
        assert not os.path.exists(str(tmp_path / "out" / "combo2.xyz"))

    def test_separate_outputs_falls_back_to_cwd_without_directory(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = IterateJob(separate_outputs=True, output_directory=None)
        runner = IterateJobRunner()
        mol = self._make_molecule()
        results = [("combo1", mol)]

        runner._write_outputs(results, job)

        assert os.path.exists("combo1.xyz")

    def test_separate_outputs_logs_and_continues_on_write_failure(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = IterateJob(separate_outputs=True, output_directory=str(tmp_path))
        runner = IterateJobRunner()
        mol = self._make_molecule()

        with patch("builtins.open", side_effect=OSError("disk full")):
            runner._write_outputs([("combo1", mol)], job)  # should not raise

        assert not os.path.exists("combo1.xyz")


class TestIterateJobRunnerRun:
    def test_fake_mode_no_op(self):
        runner = IterateJobRunner(fake=True)
        job = MagicMock()
        runner.run(job)
        # No assertions needed beyond "does not raise"; fake mode
        # short-circuits before touching job.settings at all.

    def test_no_valid_combinations_warns_and_returns(self, caplog):
        runner = IterateJobRunner(fake=False)
        job = IterateJob(settings=IterateJobSettings())
        with caplog.at_level("WARNING"):
            runner.run(job)
        assert "No valid combinations" in caplog.text

    def test_full_run_generates_writes_and_uses_run_combinations(
        self, single_molecule_xyz_file, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        settings = IterateJobSettings()
        settings.skeleton_list = [
            {
                "file_path": single_molecule_xyz_file,
                "label": "skel1",
                "link_index": [1],
            }
        ]
        settings.substituent_list = [
            {
                "file_path": single_molecule_xyz_file,
                "label": "sub1",
                "link_index": [1],
            }
        ]
        job = IterateJob(settings=settings, outputfile="results")
        runner = IterateJobRunner(fake=False)

        with (
            patch.object(
                runner,
                "run_combinations",
                return_value=[("skel1_1_sub1_1", None)],
            ) as mock_run_combinations,
            patch.object(runner, "_write_outputs") as mock_write,
        ):
            runner.run(job)

        mock_run_combinations.assert_called_once()
        mock_write.assert_called_once()


class TestRunCombinationsQueueDrainExceptions:
    def test_generic_exception_from_get_nowait_breaks_drain_loop(self):
        """A non-Empty exception raised by the queue must be swallowed
        the same way queue.Empty is (line 341-342), not propagated."""
        runner = IterateJobRunner()
        pool, combo = _make_pool_and_combo()
        result_mol = MagicMock()

        fake_queue = FakeQueue(get_nowait_plan=[RuntimeError("transient")])
        fake_queue.put((combo.label, result_mol))

        fake_process = FakeProcess(alive=False)

        with (
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Manager",
                return_value=_patched_manager(fake_queue),
            ),
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Process",
                return_value=fake_process,
            ),
        ):
            results = runner.run_combinations(
                pool, [combo], nprocs=1, timeout=120
            )

        # The RuntimeError on the first get_nowait() call just breaks the
        # drain loop for that pass; the result is still picked up on a
        # later drain (main loop or final double-check).
        assert dict(results)[combo.label] is result_mol

    def test_generic_exception_in_final_double_check_is_swallowed(self):
        """Same as above, but for the final "double check queue one
        last time" drain after the main loop exits (line 393-394)."""
        runner = IterateJobRunner()
        pool, combo = _make_pool_and_combo()

        # First call (inside the main loop's single pass): queue.Empty
        # (nothing arrived yet). Process is already "finished", so the
        # main loop exits after one pass. Second call (the final
        # double-check): a generic exception, which must be swallowed
        # rather than propagated.
        fake_queue = FakeQueue(
            get_nowait_plan=[queue.Empty(), RuntimeError("boom")]
        )
        fake_process = FakeProcess(alive=False)

        with (
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Manager",
                return_value=_patched_manager(fake_queue),
            ),
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Process",
                return_value=fake_process,
            ),
        ):
            results = runner.run_combinations(
                pool, [combo], nprocs=1, timeout=120
            )

        # Never arrived on the queue at all -> counted as a crashed/
        # missing worker.
        assert dict(results)[combo.label] is None

    def test_final_double_check_picks_up_late_arriving_result(self):
        """A result that wasn't yet on the queue during the main loop's
        last drain, but is present after the loop exits, is still
        picked up by the final double-check (line 390)."""
        runner = IterateJobRunner()
        pool, combo = _make_pool_and_combo()
        result_mol = MagicMock()

        # First get_nowait() call (main loop's only pass): Empty.
        # Second call (final double-check): the real result.
        fake_queue = FakeQueue(get_nowait_plan=[queue.Empty()])
        fake_process = FakeProcess(alive=False)

        def _late_put():
            fake_queue.put((combo.label, result_mol))

        # Simulate the result "arriving" between the main loop's drain
        # and the final double-check by seeding it right after the
        # first (Empty) call has been consumed from the plan.
        original_get_nowait = fake_queue.get_nowait
        call_count = {"n": 0}

        def _get_nowait_with_late_arrival():
            call_count["n"] += 1
            if call_count["n"] == 2:
                _late_put()
            return original_get_nowait()

        fake_queue.get_nowait = _get_nowait_with_late_arrival

        with (
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Manager",
                return_value=_patched_manager(fake_queue),
            ),
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Process",
                return_value=fake_process,
            ),
        ):
            results = runner.run_combinations(
                pool, [combo], nprocs=1, timeout=120
            )

        assert dict(results)[combo.label] is result_mol


class TestRunCombinationsTimeoutWatchdog:
    def test_stuck_process_is_terminated_then_killed(self):
        """A process still alive after both the timeout window and a
        terminate() attempt is force-killed (line 365-370)."""
        runner = IterateJobRunner()
        pool, combo = _make_pool_and_combo()

        fake_queue = FakeQueue()
        # Always "alive": guarantees both is_alive() checks in the
        # watchdog branch (before terminate, and after terminate+join)
        # see a still-running process.
        fake_process = FakeProcess(alive=True)

        with (
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Manager",
                return_value=_patched_manager(fake_queue),
            ),
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Process",
                return_value=fake_process,
            ),
        ):
            # timeout=0 -> immediately considered overdue on the very
            # first watchdog pass.
            results = runner.run_combinations(
                pool, [combo], nprocs=1, timeout=0
            )

        assert fake_process.terminated is True
        assert fake_process.killed is True
        assert dict(results)[combo.label] is None

    def test_duplicate_label_explicit_none_not_appended_twice(self):
        """Two combinations that happen to share the same label, both
        resolving to an explicit (non-timed-out) None result, must
        only add that label to failed_labels once (line 413->401)."""
        runner = IterateJobRunner()
        pool, combo_a = _make_pool_and_combo()
        combo_b = IterateCombination(
            skeleton_idx=0,
            skeleton_label=combo_a.skeleton_label,
            skeleton_link_index=combo_a.skeleton_link_index,
            skeleton_indices=None,
            substituent_idx=0,
            substituent_label=combo_a.substituent_label,
            substituent_link_index=combo_a.substituent_link_index,
        )
        assert combo_a.label == combo_b.label

        fake_queue = FakeQueue()
        fake_queue.put((combo_a.label, None))
        fake_process = FakeProcess(alive=False)

        with (
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Manager",
                return_value=_patched_manager(fake_queue),
            ),
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Process",
                return_value=fake_process,
            ),
        ):
            results = runner.run_combinations(
                pool, [combo_a, combo_b], nprocs=1, timeout=120
            )

        assert [mol for _, mol in results] == [None, None]


class TestRunCombinationsResultBookkeeping:
    def test_missing_result_counts_as_crashed_worker(self):
        """A process that exits without ever writing to the queue is
        reported as a crashed worker (line 403-407)."""
        runner = IterateJobRunner()
        pool, combo = _make_pool_and_combo()

        fake_queue = FakeQueue()  # never gets anything put on it
        fake_process = FakeProcess(alive=False)

        with (
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Manager",
                return_value=_patched_manager(fake_queue),
            ),
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Process",
                return_value=fake_process,
            ),
        ):
            results = runner.run_combinations(
                pool, [combo], nprocs=1, timeout=120
            )

        assert dict(results)[combo.label] is None

    def test_explicit_none_result_counts_as_failed_not_timed_out(self):
        """A worker that legitimately finishes and reports None (e.g.
        analyzer found no valid geometry) is counted as a failure, not
        a timeout (line 413-414)."""
        runner = IterateJobRunner()
        pool, combo = _make_pool_and_combo()

        fake_queue = FakeQueue()
        fake_queue.put((combo.label, None))
        fake_process = FakeProcess(alive=False)

        with (
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Manager",
                return_value=_patched_manager(fake_queue),
            ),
            patch(
                "chemsmart.jobs.iterate.runner.multiprocessing.Process",
                return_value=fake_process,
            ),
        ):
            results = runner.run_combinations(
                pool, [combo], nprocs=1, timeout=120
            )

        assert dict(results)[combo.label] is None
