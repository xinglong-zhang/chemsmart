"""
Direct unit tests for the abstract :class:`chemsmart.jobs.runner.JobRunner`
base class and its module-level helper functions, using a minimal
concrete subclass to exercise shared behavior (scratch resolution,
``run()`` orchestration, ``from_job`` dispatch, cleanup, and the
scratch-directory-deletion safety checks) that every real runner
(Gaussian, ORCA, thermochemistry, grouper, NCIPLOT, mol) inherits.
"""

import os
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.jobs.runner import (
    JobRunner,
    decide_phase_transition,
    get_configured_max_submitters,
    get_submitter_worker_count,
)


class ConcreteJobRunner(JobRunner):
    JOBTYPES = ["concrete_runner_type"]
    PROGRAM = "concrete"
    FAKE = False
    SCRATCH = False

    @property
    def executable(self):
        return None

    def _get_command(self, job):
        return "echo hi"

    def _create_process(self, job, command, env):
        return MagicMock()


@pytest.fixture()
def runner(pbs_server):
    return ConcreteJobRunner(server=pbs_server, scratch=False)


class TestJobRunnerConstruction:
    def test_rejects_invalid_server_type(self):
        with pytest.raises(ValueError, match="server must be instance"):
            ConcreteJobRunner(server=12345, scratch=False)

    def test_accepts_server_instance(self, pbs_server):
        runner = ConcreteJobRunner(server=pbs_server, scratch=False)
        assert runner.server is pbs_server

    def test_num_cores_defaults_from_server(self, pbs_server):
        runner = ConcreteJobRunner(server=pbs_server, scratch=False)
        assert runner.num_cores == pbs_server.num_cores

    def test_num_cores_explicit_override(self, pbs_server):
        runner = ConcreteJobRunner(
            server=pbs_server, scratch=False, num_cores=99
        )
        assert runner.num_cores == 99

    def test_repr_contains_class_and_server(self, runner):
        text = repr(runner)
        assert "ConcreteJobRunner" in text


class TestJobRunnerScratchDirSetter:
    def test_setter_rejects_nonexistent_path(self, runner):
        with pytest.raises(FileNotFoundError):
            runner.scratch_dir = "/path/does/not/exist/anywhere"

    def test_setter_accepts_existing_path(self, runner, tmp_path):
        runner.scratch_dir = str(tmp_path)
        assert runner._scratch_dir == str(tmp_path)

    def test_setter_expands_user_home(self, runner):
        # Should not raise for a bare "~" (always exists) and should be
        # expanded to an absolute path.
        runner.scratch_dir = "~"
        assert runner._scratch_dir == os.path.expanduser("~")


class TestJobRunnerSetScratch:
    def test_uses_explicit_scratch_dir_when_set(self, pbs_server, tmp_path):
        runner = ConcreteJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(tmp_path)
        )
        assert runner._scratch_dir == str(tmp_path)

    def test_disables_scratch_when_no_source_available(self, pbs_server):
        with (
            patch.object(type(pbs_server), "scratch_dir", new=None),
            patch("chemsmart.jobs.runner.user_settings") as mock_user_settings,
        ):
            mock_user_settings.scratch = None
            runner = ConcreteJobRunner(server=pbs_server, scratch=True)
        assert runner.scratch is False

    def test_raises_when_resolved_scratch_dir_missing(self, pbs_server):
        with patch.object(
            type(pbs_server), "scratch_dir", new="/does/not/exist/at/all"
        ):
            with pytest.raises(FileNotFoundError):
                ConcreteJobRunner(server=pbs_server, scratch=True)


class TestJobRunnerRunOrchestration:
    def test_run_calls_hooks_in_order(self, runner):
        job = MagicMock()
        calls = []
        runner._prerun = MagicMock(
            side_effect=lambda j: calls.append("prerun")
        )
        runner._write_input = MagicMock(
            side_effect=lambda j: calls.append("write_input")
        )
        runner._get_command = MagicMock(
            side_effect=lambda j: calls.append("get_command") or "cmd"
        )
        runner._update_os_environ = MagicMock(
            side_effect=lambda j: calls.append("update_env") or {}
        )
        runner._create_process = MagicMock(
            side_effect=lambda j, command, env: calls.append("create_process")
            or MagicMock()
        )
        runner._run = MagicMock(
            side_effect=lambda p, **kw: calls.append("run")
        )
        runner._postrun = MagicMock(
            side_effect=lambda j: calls.append("postrun")
        )
        runner._postrun_cleanup = MagicMock(
            side_effect=lambda j: calls.append("postrun_cleanup")
        )

        runner.run(job)

        assert calls == [
            "prerun",
            "write_input",
            "get_command",
            "update_env",
            "create_process",
            "run",
            "postrun",
            "postrun_cleanup",
        ]

    def test_default_run_communicates_and_polls_process(self, runner):
        process = MagicMock()
        process.poll.return_value = 0
        result = runner._run(process)
        process.communicate.assert_called_once()
        process.poll.assert_called_once()
        assert result == 0


class TestJobRunnerCopy:
    def test_copy_returns_shallow_copy(self, runner):
        copied = runner.copy()
        assert copied is not runner
        assert copied.server is runner.server


class TestJobRunnerFromJob:
    def test_from_job_dispatches_to_registered_runner(self, pbs_server):
        from chemsmart.jobs.thermochemistry.job import ThermochemistryJob
        from chemsmart.jobs.thermochemistry.runner import (
            ThermochemistryJobRunner,
        )

        job = MagicMock(spec=ThermochemistryJob)
        job.TYPE = "thermochemistry"

        result = JobRunner.from_job(job, server=pbs_server, fake=True)
        assert isinstance(result, ThermochemistryJobRunner)

    def test_from_job_raises_for_unknown_jobtype(self, pbs_server):
        job = MagicMock()
        job.TYPE = "totally_unregistered_jobtype_xyz"
        with pytest.raises(ValueError, match="Could not find any runners"):
            JobRunner.from_job(job, server=pbs_server)


class TestJobRunnerErrFileCleanup:
    def test_remove_err_files_removes_existing_suffixes(
        self, runner, tmp_path
    ):
        job = MagicMock()
        job.folder = str(tmp_path)
        job.label = "myjob"
        for suffix in [".err", ".pbserr", ".slurmerr"]:
            (tmp_path / f"myjob{suffix}").write_text("x")

        runner._remove_err_files(job)

        for suffix in [".err", ".pbserr", ".slurmerr"]:
            assert not (tmp_path / f"myjob{suffix}").exists()

    def test_remove_err_files_no_op_when_missing(self, runner, tmp_path):
        job = MagicMock()
        job.folder = str(tmp_path)
        job.label = "myjob"
        # Should not raise even though no err files exist.
        runner._remove_err_files(job)

    def test_append_suffix_to_job_label_appends_once(self, runner):
        job = MagicMock()
        job.label = "myjob"
        runner._append_suffix_to_job_label(job, "_fake")
        assert job.label == "myjob_fake"
        runner._append_suffix_to_job_label(job, "_fake")
        assert job.label == "myjob_fake"

    def test_append_suffix_no_op_for_falsy_suffix(self, runner):
        job = MagicMock()
        job.label = "myjob"
        runner._append_suffix_to_job_label(job, "")
        assert job.label == "myjob"


class TestJobRunnerPostrunCleanup:
    def test_cleanup_removes_err_files_when_complete(self, runner):
        job = MagicMock()
        job.is_complete.return_value = True
        runner._remove_err_files = MagicMock()
        runner._delete_scratch_directory = MagicMock()
        runner.scratch = False

        runner._postrun_cleanup(job)

        runner._remove_err_files.assert_called_once_with(job)
        runner._delete_scratch_directory.assert_not_called()

    def test_cleanup_no_op_when_incomplete(self, runner):
        job = MagicMock()
        job.is_complete.return_value = False
        runner._remove_err_files = MagicMock()

        runner._postrun_cleanup(job)

        runner._remove_err_files.assert_not_called()

    def test_cleanup_deletes_scratch_when_enabled(self, runner):
        job = MagicMock()
        job.is_complete.return_value = True
        runner._remove_err_files = MagicMock()
        runner._delete_scratch_directory = MagicMock()
        runner.scratch = True
        runner.delete_scratch = True

        runner._postrun_cleanup(job)

        runner._delete_scratch_directory.assert_called_once()


class TestJobRunnerDeleteScratchDirectory:
    def test_no_op_without_running_directory_attr(self, runner):
        # runner has no `running_directory` set at all.
        runner._delete_scratch_directory()  # should not raise

    def test_refuses_when_running_directory_equals_scratch_root(
        self, runner, tmp_path
    ):
        runner.running_directory = str(tmp_path)
        runner._scratch_dir = str(tmp_path)
        with patch("chemsmart.jobs.runner.rmtree") as mock_rmtree:
            runner._delete_scratch_directory()
        mock_rmtree.assert_not_called()

    def test_refuses_when_running_directory_outside_scratch(
        self, runner, tmp_path
    ):
        scratch_root = tmp_path / "scratch"
        scratch_root.mkdir()
        outside_dir = tmp_path / "outside"
        outside_dir.mkdir()

        runner.running_directory = str(outside_dir)
        runner._scratch_dir = str(scratch_root)
        with patch("chemsmart.jobs.runner.rmtree") as mock_rmtree:
            runner._delete_scratch_directory()
        mock_rmtree.assert_not_called()
        assert outside_dir.exists()

    def test_deletes_when_running_directory_inside_scratch(
        self, runner, tmp_path
    ):
        scratch_root = tmp_path / "scratch"
        job_dir = scratch_root / "myjob"
        job_dir.mkdir(parents=True)

        runner.running_directory = str(job_dir)
        runner._scratch_dir = str(scratch_root)

        runner._delete_scratch_directory()

        assert not job_dir.exists()
        assert scratch_root.exists()


class TestModuleLevelHelpers:
    def test_get_configured_max_submitters_from_env(self, monkeypatch):
        monkeypatch.setenv("CHEMSMART_MAX_SUBMITTERS", "7")
        assert get_configured_max_submitters(jobrunner=None) == 7

    def test_get_configured_max_submitters_from_runner_num_cores(
        self, monkeypatch
    ):
        monkeypatch.delenv("CHEMSMART_MAX_SUBMITTERS", raising=False)
        runner = MagicMock(spec=[])
        runner.num_cores = 4
        assert get_configured_max_submitters(jobrunner=runner) == 4

    def test_get_configured_max_submitters_falls_back_to_cpu_count(
        self, monkeypatch
    ):
        monkeypatch.delenv("CHEMSMART_MAX_SUBMITTERS", raising=False)
        with patch("os.cpu_count", return_value=16):
            assert get_configured_max_submitters(jobrunner=None) == 16

    def test_get_submitter_worker_count_bounds_by_num_jobs(self, monkeypatch):
        monkeypatch.setenv("CHEMSMART_MAX_SUBMITTERS", "10")
        assert get_submitter_worker_count(None, num_jobs=3) == 3

    def test_get_submitter_worker_count_bounds_by_max_submitters(
        self, monkeypatch
    ):
        monkeypatch.setenv("CHEMSMART_MAX_SUBMITTERS", "2")
        assert get_submitter_worker_count(None, num_jobs=10) == 2

    def test_get_submitter_worker_count_at_least_one_for_zero_jobs(self):
        assert get_submitter_worker_count(None, num_jobs=0) == 1


class TestDecidePhaseTransition:
    def test_proceeds_with_no_failures(self):
        decision = decide_phase_transition(phase_name="opt")
        assert decision.proceed is True
        assert decision.should_raise is False

    def test_raises_on_failures(self):
        decision = decide_phase_transition(
            phase_name="opt", failures=["worker1 failed"]
        )
        assert decision.proceed is False
        assert decision.should_raise is True
        assert "worker1 failed" in decision.message

    def test_halts_without_raising_when_incomplete_required(self):
        decision = decide_phase_transition(
            phase_name="opt", require_complete=True, is_complete=False
        )
        assert decision.proceed is False
        assert decision.should_raise is False
        assert "incomplete" in decision.message
