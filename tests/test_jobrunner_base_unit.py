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
    _executable_class_for_program,
    decide_phase_transition,
    get_configured_max_submitters,
    get_submitter_worker_count,
    run_phase_jobs,
)
from chemsmart.settings.server import Server


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

    def test_server_as_string_resolves_via_from_servername(self):
        fake_server = MagicMock(spec=Server)
        fake_server.num_cores = 8
        fake_server.num_gpus = 0
        fake_server.mem_gb = 32
        with patch.object(
            Server, "from_servername", return_value=fake_server
        ) as mock_from_servername:
            runner = ConcreteJobRunner(server="myserver", scratch=False)
        mock_from_servername.assert_called_once_with("myserver")
        assert runner.server is fake_server

    def test_num_gpus_explicit_override(self, pbs_server):
        runner = ConcreteJobRunner(
            server=pbs_server, scratch=False, num_gpus=3
        )
        assert runner.num_gpus == 3

    def test_mem_gb_explicit_override(self, pbs_server):
        runner = ConcreteJobRunner(
            server=pbs_server, scratch=False, mem_gb=256
        )
        assert runner.mem_gb == 256

    def test_servername_property_delegates_to_server(self, runner):
        assert runner.servername == runner.server.name


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

    def test_setter_accepts_none_without_validation(self, runner):
        """Covers the setter's `if value is not None:` False arm: a
        None value skips expanduser/exists-check entirely."""
        runner._scratch_dir = "/some/previous/value"
        runner.scratch_dir = None
        assert runner._scratch_dir is None

    def test_getter_computes_and_caches_when_unset(self, pbs_server, tmp_path):
        """Covers the getter's `if self._scratch_dir is None:` branch:
        the first access to .scratch_dir must resolve and cache it via
        _set_scratch(), rather than staying None forever."""
        with patch.object(type(pbs_server), "scratch_dir", new=str(tmp_path)):
            runner = ConcreteJobRunner(server=pbs_server, scratch=True)
            assert runner._scratch_dir is None
            resolved = runner.scratch_dir
        assert resolved == str(tmp_path)
        assert runner._scratch_dir == str(tmp_path)


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

    def test_prerun_write_input_postrun_are_no_ops_by_default(self, runner):
        """ConcreteJobRunner doesn't override these optional hooks, so
        the base class's own pass-through bodies run for real."""
        job = MagicMock()
        assert runner._prerun(job) is None
        assert runner._write_input(job) is None
        assert runner._postrun(job) is None


class TestUpdateOsEnviron:
    def test_no_executable_returns_plain_environ_copy(self, runner):
        job = MagicMock()
        env = runner._update_os_environ(job)
        assert env == dict(os.environ)

    def test_empty_executable_env_returns_plain_environ_copy(self, runner):
        job = MagicMock()
        with patch.object(
            type(runner),
            "executable",
            new=property(lambda self: MagicMock(env={})),
        ):
            env = runner._update_os_environ(job)
        assert env == dict(os.environ)

    def test_executable_env_vars_are_applied_with_user_expansion(self, runner):
        job = MagicMock()
        with patch.object(
            type(runner),
            "executable",
            new=property(
                lambda self: MagicMock(
                    env={"MY_STR_VAR": "~/mydir", "MY_INT_VAR": 5}
                )
            ),
        ):
            env = runner._update_os_environ(job)
        assert env["MY_STR_VAR"] == os.path.expanduser("~/mydir")
        assert env["MY_INT_VAR"] == "5"


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

    def test_from_job_treats_notimplemented_jobtypes_as_empty_list(
        self, pbs_server
    ):
        """A registered runner subclass that never overrides JOBTYPES
        (still NotImplemented) must be treated as supporting no
        jobtypes, not crash the `jobtype in runner_jobtypes` check."""

        class _RunnerWithoutJobtypes(JobRunner):
            PROGRAM = "concrete_no_jobtypes"
            FAKE = False
            SCRATCH = False

            @property
            def executable(self):
                return None

            def _get_command(self, job):
                return "echo hi"

            def _create_process(self, job, command, env):
                return MagicMock()

        job = MagicMock()
        job.TYPE = "concrete_runner_type"
        result = JobRunner.from_job(job, server=pbs_server, fake=True)
        assert isinstance(result, ConcreteJobRunner)


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

    def test_logs_error_when_scratch_dir_missing(
        self, runner, tmp_path, caplog
    ):
        """Covers the `if not sd.exists() or not sd.is_dir():` branch:
        running_directory exists, but the resolved scratch_dir itself
        doesn't -- refuse to proceed instead of raising."""
        running_dir = tmp_path / "running"
        running_dir.mkdir()

        runner.running_directory = str(running_dir)
        runner._scratch_dir = str(tmp_path / "does_not_exist_scratch")

        with (
            patch("chemsmart.jobs.runner.rmtree") as mock_rmtree,
            caplog.at_level("ERROR"),
        ):
            runner._delete_scratch_directory()

        mock_rmtree.assert_not_called()
        assert "doesn't exist or is not a directory" in caplog.text

    def test_rmtree_failure_is_logged_not_raised(self, runner, tmp_path):
        scratch_root = tmp_path / "scratch"
        job_dir = scratch_root / "myjob"
        job_dir.mkdir(parents=True)

        runner.running_directory = str(job_dir)
        runner._scratch_dir = str(scratch_root)

        with patch(
            "chemsmart.jobs.runner.rmtree",
            side_effect=OSError("simulated rmtree failure"),
        ):
            runner._delete_scratch_directory()  # should not raise


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

    def test_get_configured_max_submitters_zero_env_falls_through(
        self, monkeypatch
    ):
        """A non-positive env value is treated as unset (covers
        _positive_int_or_none's `parsed <= 0` branch), falling through
        to the next resolution step instead of returning 0."""
        monkeypatch.setenv("CHEMSMART_MAX_SUBMITTERS", "0")
        runner = MagicMock(spec=[])
        runner.num_cores = 5
        assert get_configured_max_submitters(jobrunner=runner) == 5

    def test_get_configured_max_submitters_from_runner_max_submitters(
        self, monkeypatch
    ):
        monkeypatch.delenv("CHEMSMART_MAX_SUBMITTERS", raising=False)
        runner = MagicMock(spec=["max_submitters"])
        runner.max_submitters = 9
        assert get_configured_max_submitters(jobrunner=runner) == 9

    def test_get_configured_max_submitters_from_server_max_submitters(
        self, monkeypatch
    ):
        monkeypatch.delenv("CHEMSMART_MAX_SUBMITTERS", raising=False)
        runner = MagicMock(spec=["server"])
        runner.server = MagicMock(spec=["max_submitters"])
        runner.server.max_submitters = 6
        assert get_configured_max_submitters(jobrunner=runner) == 6

    def test_get_configured_max_submitters_from_server_num_cores(
        self, monkeypatch
    ):
        monkeypatch.delenv("CHEMSMART_MAX_SUBMITTERS", raising=False)
        runner = MagicMock(spec=["server"])
        runner.server = MagicMock(spec=["num_cores"])
        runner.server.num_cores = 12
        assert get_configured_max_submitters(jobrunner=runner) == 12

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


class TestExecutableClassForProgram:
    def test_falsy_program_returns_none(self):
        assert _executable_class_for_program(None) is None
        assert _executable_class_for_program("") is None

    def test_notimplemented_program_returns_none(self):
        assert _executable_class_for_program(NotImplemented) is None


class TestRunPhaseJobs:
    def test_delegates_to_job_execute_phase_jobs(self):
        with patch(
            "chemsmart.jobs.runner.Job._execute_phase_jobs"
        ) as mock_execute:
            run_phase_jobs(
                parent_runner="runner-sentinel",
                jobs=["job1"],
                stop_on_incomplete=True,
                phase_label="opt",
            )
        mock_execute.assert_called_once_with(
            parent_runner="runner-sentinel",
            jobs=["job1"],
            jobs_factory=None,
            stop_on_incomplete=True,
            before_run=None,
            logger_obj=None,
            phase_label="opt",
        )
