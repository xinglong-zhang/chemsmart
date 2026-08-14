"""
Direct unit tests for the abstract :class:`chemsmart.jobs.job.Job` base
class, using a minimal concrete subclass to exercise shared behavior
(folder handling, run/skip-completed logic, runner propagation, phase
job orchestration, and backup helpers) that every job subclass inherits.
"""

import os
from unittest.mock import MagicMock

import pytest

from chemsmart.jobs.job import Job


class ConcreteJob(Job):
    """Minimal concrete subclass for exercising base ``Job`` behavior."""

    TYPE = "concrete"

    def __init__(self, complete=False, output=None, **kwargs):
        super().__init__(**kwargs)
        self._complete = complete
        self._output_obj = output
        self.ran = False

    def _run(self, **kwargs):
        self.ran = True

    def _job_is_complete(self):
        return self._complete

    def _output(self):
        return self._output_obj

    def _backup_files(self, **kwargs):
        pass


@pytest.fixture()
def concrete_job(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    job = ConcreteJob(
        molecule=None, label="myjob", jobrunner=None, skip_completed=True
    )
    return job


class TestJobFolderHandling:
    def test_folder_defaults_to_cwd(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        job = ConcreteJob(molecule=None, label="x", jobrunner=None)
        assert job.folder == os.path.abspath(str(tmp_path))

    def test_folder_setter(self, concrete_job, tmp_path):
        new_folder = str(tmp_path / "elsewhere")
        concrete_job.folder = new_folder
        assert concrete_job.folder == new_folder

    def test_set_folder_method(self, concrete_job, tmp_path):
        new_folder = str(tmp_path / "elsewhere2")
        concrete_job.set_folder(new_folder)
        assert concrete_job.folder == new_folder

    def test_base_folder_returns_basename(self, concrete_job):
        concrete_job.folder = "/a/b/myfolder"
        assert concrete_job.base_folder() == "myfolder"

    def test_base_folder_none_when_folder_none(self, concrete_job):
        concrete_job.folder = None
        assert concrete_job.base_folder() is None

    def test_joblog_path(self, concrete_job):
        expected = os.path.join(concrete_job.folder, "myjob.joblog")
        assert concrete_job.joblog == expected

    def test_repr_contains_folder_label_jobrunner(self, concrete_job):
        text = repr(concrete_job)
        assert "ConcreteJob" in text
        assert concrete_job.folder in text
        assert "myjob" in text


class TestJobRunAndComplete:
    def test_run_executes_when_not_complete(self, concrete_job):
        concrete_job._complete = False
        concrete_job.run()
        assert concrete_job.ran is True

    def test_run_skips_when_complete_and_skip_completed(self, concrete_job):
        concrete_job._complete = True
        concrete_job.skip_completed = True
        concrete_job.run()
        assert concrete_job.ran is False

    def test_run_executes_when_complete_but_skip_completed_false(
        self, concrete_job
    ):
        concrete_job._complete = True
        concrete_job.skip_completed = False
        concrete_job.run()
        assert concrete_job.ran is True

    def test_is_complete_delegates_to_job_is_complete(self, concrete_job):
        concrete_job._complete = True
        assert concrete_job.is_complete() is True
        concrete_job._complete = False
        assert concrete_job.is_complete() is False


class TestJobOptimizedStructure:
    def test_optimized_structure_none_when_no_output(self, concrete_job):
        concrete_job._output_obj = None
        assert concrete_job.optimized_structure() is None

    def test_optimized_structure_none_when_not_normal_termination(
        self, concrete_job
    ):
        output = MagicMock()
        output.normal_termination = False
        concrete_job._output_obj = output
        assert concrete_job.optimized_structure() is None

    def test_optimized_structure_returned_on_normal_termination(
        self, concrete_job
    ):
        output = MagicMock()
        output.normal_termination = True
        output.optimized_structure = "the_structure"
        concrete_job._output_obj = output
        assert concrete_job.optimized_structure() == "the_structure"

    def test_intermediate_optimization_points_empty_without_output(
        self, concrete_job
    ):
        concrete_job._output_obj = None
        assert concrete_job._intermediate_optimization_points() == []

    def test_intermediate_optimization_points_from_output(self, concrete_job):
        output = MagicMock()
        output.all_structures = ["s1", "s2"]
        concrete_job._output_obj = output
        assert concrete_job._intermediate_optimization_points() == [
            "s1",
            "s2",
        ]


class TestJobRunnerPropagation:
    def test_propagate_runner_returns_none_for_falsy_runner(self):
        job = ConcreteJob(molecule=None, label="x", jobrunner=None)
        result = Job._propagate_runner(None, job)
        assert result is None
        assert job.jobrunner is None

    def test_propagate_runner_copies_and_assigns(self):
        job = ConcreteJob(molecule=None, label="x", jobrunner=None)
        parent_runner = MagicMock()
        child_copy = MagicMock()
        parent_runner.copy.return_value = child_copy

        result = Job._propagate_runner(parent_runner, job)

        parent_runner.copy.assert_called_once()
        assert result is child_copy
        assert job.jobrunner is child_copy

    def test_propagate_runner_overrides_resources(self):
        job = ConcreteJob(molecule=None, label="x", jobrunner=None)
        parent_runner = MagicMock()
        child_copy = MagicMock()
        parent_runner.copy.return_value = child_copy

        Job._propagate_runner(
            parent_runner, job, num_cores=8, mem_gb=16, num_gpus=1
        )

        assert child_copy.num_cores == 8
        assert child_copy.mem_gb == 16
        assert child_copy.num_gpus == 1


class TestJobExecutePhaseJobs:
    def test_no_op_when_no_jobs_and_no_factory(self):
        # Should not raise.
        Job._execute_phase_jobs(parent_runner=None, jobs=None)

    def test_runs_each_job_and_propagates_runner(self):
        parent_runner = MagicMock()
        child_copy = MagicMock()
        parent_runner.copy.return_value = child_copy

        job1 = ConcreteJob(molecule=None, label="j1", jobrunner=None)
        job2 = ConcreteJob(molecule=None, label="j2", jobrunner=None)

        Job._execute_phase_jobs(parent_runner=parent_runner, jobs=[job1, job2])

        assert job1.ran is True
        assert job2.ran is True
        assert job1.jobrunner is child_copy
        assert job2.jobrunner is child_copy

    def test_jobs_factory_used_when_provided(self):
        job1 = ConcreteJob(molecule=None, label="j1", jobrunner=None)
        factory = MagicMock(return_value=[job1])

        Job._execute_phase_jobs(
            parent_runner=None, jobs=None, jobs_factory=factory
        )

        factory.assert_called_once()
        assert job1.ran is True

    def test_before_run_called_before_jobs(self):
        calls = []
        job1 = ConcreteJob(molecule=None, label="j1", jobrunner=None)

        def before_run():
            calls.append("before")

        Job._execute_phase_jobs(
            parent_runner=None,
            jobs=[job1],
            before_run=before_run,
        )
        assert calls == ["before"]
        assert job1.ran is True

    def test_stop_on_incomplete_breaks_loop(self):
        job1 = ConcreteJob(
            molecule=None, label="j1", jobrunner=None, complete=False
        )
        job2 = ConcreteJob(
            molecule=None, label="j2", jobrunner=None, complete=True
        )

        Job._execute_phase_jobs(
            parent_runner=None,
            jobs=[job1, job2],
            stop_on_incomplete=True,
        )

        assert job1.ran is True
        # job2 should never run since job1 stayed incomplete and broke the loop.
        assert job2.ran is False


class TestJobBackup:
    def test_previous_backup_folders_empty_initially(self, concrete_job):
        assert concrete_job._previous_backup_folders() == []

    def test_create_backup_folder_name_uses_prefix_and_folder(
        self, concrete_job
    ):
        name = concrete_job._create_backup_folder_name(prefix="bk")
        assert name.startswith(os.path.join(concrete_job.folder, "bk."))

    def test_backup_file_noop_when_file_missing(self, concrete_job, tmp_path):
        missing = str(tmp_path / "does_not_exist.txt")
        # Should not raise.
        concrete_job.backup_file(missing)

    def test_backup_file_copies_file(self, concrete_job, tmp_path):
        src = tmp_path / "data.txt"
        src.write_text("hello")
        backup_folder = tmp_path / "backup"

        concrete_job.backup_file(str(src), folder=str(backup_folder))

        copied = backup_folder / "data.txt"
        assert copied.exists()
        assert copied.read_text() == "hello"
        assert src.exists()  # original not removed

    def test_backup_file_moves_when_remove_true(self, concrete_job, tmp_path):
        src = tmp_path / "data.txt"
        src.write_text("hello")
        backup_folder = tmp_path / "backup"

        concrete_job.backup_file(
            str(src), folder=str(backup_folder), remove=True
        )

        assert not src.exists()
        assert (backup_folder / "data.txt").exists()

    def test_backup_calls_backup_files(self, concrete_job):
        concrete_job._backup_files = MagicMock()
        concrete_job.backup(extra="kwarg")
        concrete_job._backup_files.assert_called_once_with(extra="kwarg")


class TestJobFromMolecule:
    def test_from_molecule_constructs_job(self):
        job = ConcreteJob.from_molecule(
            molecule=None, label="from_mol", jobrunner=None
        )
        assert isinstance(job, ConcreteJob)
        assert job.label == "from_mol"
