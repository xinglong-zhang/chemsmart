"""Tests for JobRunner class and no_run_in_parallel flag."""

from unittest.mock import Mock

import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.runner import JobRunner, get_submitter_worker_count


class MockMolecule(Molecule):
    def __init__(self, symbols=None, **kwargs):
        self.symbols = symbols or ["C"]
        self.coordinates = [[0.0, 0.0, 0.0]]
        self.charge = 0
        self.multiplicity = 1
        self.num_atoms = len(self.symbols)

    @property
    def has_vibrations(self):
        return True

    def vibrationally_displaced(self, *args, **kwargs):
        return MockMolecule()


class TestJobRunner:
    """Test JobRunner functionality."""

    def test_jobrunner_default_no_run_in_parallel(self, pbs_server):
        """Test that no_run_in_parallel defaults to False."""
        runner = JobRunner(server=pbs_server, fake=True)
        assert hasattr(runner, "no_run_in_parallel")
        assert runner.no_run_in_parallel is False

    def test_jobrunner_no_run_in_parallel_true(self, pbs_server):
        """Test that no_run_in_parallel can be set to True."""
        runner = JobRunner(
            server=pbs_server, fake=True, no_run_in_parallel=True
        )
        assert runner.no_run_in_parallel is True

    def test_jobrunner_no_run_in_parallel_false(self, pbs_server):
        """Test that no_run_in_parallel can be set to False."""
        runner = JobRunner(
            server=pbs_server, fake=True, no_run_in_parallel=False
        )
        assert runner.no_run_in_parallel is False

    def test_submitter_worker_count_uses_policy_cap(self, pbs_server):
        runner = JobRunner(server=pbs_server, fake=True)
        assert get_submitter_worker_count(runner, 1000) == runner.num_cores

    def test_submitter_worker_count_honors_env_override(
        self, pbs_server, monkeypatch
    ):
        monkeypatch.setenv("CHEMSMART_MAX_SUBMITTERS", "3")
        runner = JobRunner(server=pbs_server, fake=True)
        assert get_submitter_worker_count(runner, 1000) == 3


class TestSerialExecution:
    """Test serial execution enforcement for jobs with subjobs."""

    def test_crest_job_serial_execution_with_completion(
        self, pbs_server, gaussian_jobrunner_no_scratch, mocker
    ):
        """Test that GaussianCrestJob runs all jobs when they complete successfully."""
        from chemsmart.jobs.gaussian.crest import GaussianCrestJob
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings

        # Create mock molecules
        mock_molecules = [MockMolecule(), MockMolecule(), MockMolecule()]

        # Create settings
        settings = GaussianJobSettings()

        # Create job with no_run_in_parallel=True
        gaussian_jobrunner_no_scratch.no_run_in_parallel = True
        job = GaussianCrestJob(
            molecules=mock_molecules,
            settings=settings,
            label="test_crest",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        # Mock conformer jobs that all complete successfully
        mock_jobs = []
        for i in range(3):
            mock_job = Mock()
            mock_job.label = f"conf_{i}"
            mock_job.is_complete.return_value = True
            mock_jobs.append(mock_job)

        # Patch _prepare_all_jobs instead of the property
        mocker.patch.object(job, "_prepare_all_jobs", return_value=mock_jobs)

        # Run the job
        job._run()

        # Verify all jobs were run
        for mock_job in mock_jobs:
            mock_job.run.assert_called_once()


class TestBatchJobRefactor:
    """Regression tests for the shared BatchJob orchestration layer."""

    @staticmethod
    def _dummy_batch_cls():
        from chemsmart.jobs.batch import BatchJob

        class DummyBatchJob(BatchJob):
            PROGRAM = "dummy"

            def _configure_runner_for_node(self, runner, node, job):
                return runner

        return DummyBatchJob

    def test_split_jobs_across_nodes_balances_chunks(self):
        from chemsmart.jobs.batch import BatchJob

        jobs = [1, 2, 3, 4, 5]
        chunks = BatchJob._split_jobs_across_nodes(jobs, 3)

        assert chunks == [[1, 2], [3, 4], [5]]

    def test_orca_and_gaussian_batch_pin_runner_to_slurm_node(
        self, pbs_server, monkeypatch
    ):
        """Engine wrappers should pin child commands via srun --nodelist."""
        from chemsmart.jobs.gaussian.batch import GaussianBatchJob
        from chemsmart.jobs.orca.batch import OrcaBatchJob

        monkeypatch.setenv("SLURM_JOB_NODELIST", "nodeA")

        for batch_cls in (GaussianBatchJob, OrcaBatchJob):
            runner = JobRunner(server=pbs_server, fake=True)
            runner._get_command = Mock(return_value="run.exe job.inp")
            batch = batch_cls(jobs=[], jobrunner=runner)
            pinned = batch._configure_runner_for_node(
                runner=runner, node="nodeA", job=Mock()
            )
            command = pinned._get_command(Mock())
            assert command.startswith(
                "srun --nodelist=nodeA --exclusive -N1 -n1 "
            )
            assert command.endswith("run.exe job.inp")

    def test_batch_serial_mode_when_unset(self, pbs_server):
        """Test that BatchJob runs all jobs when they are unset."""
        dummy_batch_cls = self._dummy_batch_cls()
        runner = JobRunner(
            server=pbs_server, fake=True, no_run_in_parallel=True
        )

        batch = dummy_batch_cls(
            jobs=[],
            jobrunner=runner,
        )

        assert batch.no_run_in_parallel is True

    def test_batch_serial_mode_keeps_explicit_false(self, pbs_server):
        """Test that BatchJob runs all jobs when they are unset."""
        dummy_batch_cls = self._dummy_batch_cls()
        runner = JobRunner(
            server=pbs_server, fake=True, no_run_in_parallel=True
        )

        batch = dummy_batch_cls(
            jobs=[],
            no_run_in_parallel=False,
            jobrunner=runner,
        )

        assert batch.no_run_in_parallel is False

    def test_batch_serial_mode_keeps_explicit_true(self, pbs_server):
        dummy_batch_cls = self._dummy_batch_cls()
        runner = JobRunner(
            server=pbs_server, fake=True, no_run_in_parallel=False
        )

        batch = dummy_batch_cls(
            jobs=[],
            no_run_in_parallel=True,
            jobrunner=runner,
        )

        assert batch.no_run_in_parallel is True

    def test_batch_writes_success_and_failed_logs(self, pbs_server, tmp_path):
        from chemsmart.jobs.batch import BatchExecutionError

        dummy_batch_cls = self._dummy_batch_cls()
        runner = JobRunner(
            server=pbs_server, fake=True, no_run_in_parallel=False
        )

        success_job = Mock()
        success_job.label = "success_job"
        success_job.run.return_value = None
        success_job.is_complete.return_value = True

        fail_job = Mock()
        fail_job.label = "fail_job"
        fail_job.run.side_effect = RuntimeError("fail")
        fail_job.is_complete.return_value = False

        batch = dummy_batch_cls(
            jobs=[success_job, fail_job],
            no_run_in_parallel=False,
            write_outcome_logs=True,
            jobrunner=runner,
            label="batch_logs_test",
        )
        batch.folder = str(tmp_path)

        with pytest.raises(BatchExecutionError):
            batch.run()

        success_lines = (tmp_path / "success.log").read_text().splitlines()
        failed_lines = (tmp_path / "failed.log").read_text().splitlines()

        assert "success_job" in success_lines
        assert any(line.startswith("fail_job\t") for line in failed_lines)

    def test_batch_run_raises_with_failed_job_summary(
        self, pbs_server, tmp_path
    ):
        from chemsmart.jobs.batch import BatchExecutionError

        dummy_batch_cls = self._dummy_batch_cls()
        runner = JobRunner(
            server=pbs_server, fake=True, no_run_in_parallel=True
        )

        fail_job = Mock()
        fail_job.label = "fail_job"
        fail_job.run.side_effect = RuntimeError("fail")
        fail_job.is_complete.return_value = False

        batch = dummy_batch_cls(
            jobs=[fail_job],
            no_run_in_parallel=True,
            jobrunner=runner,
            label="batch_raise_test",
        )
        batch.folder = str(tmp_path)

        with pytest.raises(BatchExecutionError, match="fail_job"):
            batch.run()

    def test_serial_fail_fast_false_attempts_all_then_raises(self, pbs_server):
        from chemsmart.jobs.batch import BatchExecutionError

        dummy_batch_cls = self._dummy_batch_cls()
        runner = JobRunner(
            server=pbs_server, fake=True, no_run_in_parallel=True
        )

        fail_job = Mock(label="fail_job")
        fail_job.run.side_effect = RuntimeError("fail")
        fail_job.is_complete.return_value = False

        ok_job = Mock(label="ok_job")
        ok_job.run.return_value = None
        ok_job.is_complete.return_value = True

        batch = dummy_batch_cls(
            jobs=[fail_job, ok_job],
            no_run_in_parallel=True,
            fail_fast=False,
            jobrunner=runner,
        )

        with pytest.raises(BatchExecutionError, match="1 of 2") as exc_info:
            batch.run()

        assert "not started" not in str(exc_info.value)
        fail_job.run.assert_called_once()
        ok_job.run.assert_called_once()

    def test_serial_fail_fast_true_omits_unstarted_from_report(
        self, pbs_server
    ):
        from chemsmart.jobs.batch import BatchExecutionError

        dummy_batch_cls = self._dummy_batch_cls()
        runner = JobRunner(
            server=pbs_server, fake=True, no_run_in_parallel=True
        )

        fail_job = Mock(label="fail_job")
        fail_job.run.side_effect = RuntimeError("fail")
        fail_job.is_complete.return_value = False

        later_job = Mock(label="later_job")
        later_job.run.return_value = None
        later_job.is_complete.return_value = True

        batch = dummy_batch_cls(
            jobs=[fail_job, later_job],
            no_run_in_parallel=True,
            fail_fast=True,
            jobrunner=runner,
        )

        with pytest.raises(
            BatchExecutionError, match="1 attempted, 1 failed, 1 not started"
        ) as exc_info:
            batch.run()

        assert "fail_job" in str(exc_info.value)
        assert "later_job" not in str(exc_info.value)
        fail_job.run.assert_called_once()
        later_job.run.assert_not_called()
        assert batch._jobs_not_started == 1

    def test_parallel_fail_fast_still_submits_all(self, pbs_server):
        from chemsmart.jobs.batch import BatchExecutionError

        dummy_batch_cls = self._dummy_batch_cls()
        runner = JobRunner(
            server=pbs_server, fake=True, no_run_in_parallel=False
        )

        fail_job = Mock(label="fail_job")
        fail_job.run.side_effect = RuntimeError("fail")
        fail_job.is_complete.return_value = False

        ok_job = Mock(label="ok_job")
        ok_job.run.return_value = None
        ok_job.is_complete.return_value = True

        batch = dummy_batch_cls(
            jobs=[fail_job, ok_job],
            no_run_in_parallel=False,
            fail_fast=True,
            jobrunner=runner,
        )

        with pytest.raises(BatchExecutionError, match="1 of 2"):
            batch.run()

        fail_job.run.assert_called_once()
        ok_job.run.assert_called_once()

    def test_run_child_jobs_as_batch_forwards_serial_mode(self, pbs_server):
        from chemsmart.jobs.batch import run_child_jobs_as_batch

        dummy_batch_cls = self._dummy_batch_cls()
        runner = JobRunner(
            server=pbs_server, fake=True, no_run_in_parallel=True
        )
        parent = Mock()
        parent.label = "parent"
        parent.jobrunner = runner

        ok_job = Mock(label="child")
        ok_job.run.return_value = None
        ok_job.is_complete.return_value = True

        batch = run_child_jobs_as_batch(
            batch_cls=dummy_batch_cls,
            jobs=[ok_job],
            parent=parent,
            label_suffix="_batch",
            fail_fast=False,
        )

        assert batch.no_run_in_parallel is True
        assert batch.fail_fast is False
        assert batch.label == "parent_batch"
        ok_job.run.assert_called_once()

    def test_batch_run_multi_node_records_node_future_exception(
        self, pbs_server, mocker
    ):
        dummy_batch_cls = self._dummy_batch_cls()
        runner = JobRunner(
            server=pbs_server, fake=True, no_run_in_parallel=False
        )

        job1 = Mock()
        job1.label = "job1"
        job2 = Mock()
        job2.label = "job2"

        batch = dummy_batch_cls(
            jobs=[job1, job2],
            no_run_in_parallel=False,
            jobrunner=runner,
            label="batch_node_future_failure",
        )

        def _chunk_side_effect(jobs_chunk, node, **kwargs):
            if node == "n2":
                raise RuntimeError("node future failure")
            return [
                {
                    "label": jobs_chunk[0].label,
                    "success": True,
                    "error": "",
                    "node": node,
                }
            ]

        mocker.patch.object(
            batch, "_run_chunk_on_node", side_effect=_chunk_side_effect
        )

        outcomes = batch._run_multi_node(nodes=["n1", "n2"])

        assert any(item["success"] for item in outcomes)
        node_failure = [
            item
            for item in outcomes
            if (not item["success"]) and item["node"] == "n2"
        ]
        assert len(node_failure) == 1
        assert node_failure[0]["label"].startswith("node:n2:")
        assert "node future failure" in node_failure[0]["error"]


class TestGaussianBatchDelegation:
    """Tests for Gaussian multi-subjob workflows using GaussianBatchJob."""

    def test_crest_job_uses_batch_parallel(
        self, pbs_server, gaussian_jobrunner_no_scratch, mocker
    ):
        from chemsmart.jobs.gaussian.crest import GaussianCrestJob
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings

        mock_molecules = [MockMolecule(), MockMolecule(), MockMolecule()]
        settings = GaussianJobSettings()

        gaussian_jobrunner_no_scratch.no_run_in_parallel = False
        job = GaussianCrestJob(
            molecules=mock_molecules,
            settings=settings,
            label="test_crest",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        mock_jobs = [Mock(label="conf_1"), Mock(label="conf_2")]
        mocker.patch.object(job, "_prepare_all_jobs", return_value=mock_jobs)

        mock_batch_cls = mocker.patch(
            "chemsmart.jobs.gaussian.crest.GaussianBatchJob"
        )
        mock_batch = mock_batch_cls.return_value

        job._run_all_jobs()

        mock_batch_cls.assert_called_once()
        call_kwargs = mock_batch_cls.call_args.kwargs
        assert call_kwargs["jobs"] == mock_jobs
        assert call_kwargs["no_run_in_parallel"] is False
        assert call_kwargs["fail_fast"] is False
        assert call_kwargs["label"] == "test_crest_batch"
        assert call_kwargs["jobrunner"] == gaussian_jobrunner_no_scratch
        mock_batch.run.assert_called_once()

    def test_traj_job_uses_batch_parallel(
        self, pbs_server, gaussian_jobrunner_no_scratch, mocker
    ):
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings
        from chemsmart.jobs.gaussian.traj import GaussianTrajJob

        mock_molecules = [MockMolecule(), MockMolecule(), MockMolecule()]
        settings = GaussianJobSettings()

        gaussian_jobrunner_no_scratch.no_run_in_parallel = False
        job = GaussianTrajJob(
            molecules=mock_molecules,
            settings=settings,
            label="test_traj",
            jobrunner=gaussian_jobrunner_no_scratch,
            proportion_structures_to_use=1.0,
        )

        mock_jobs = [Mock(label="traj_1"), Mock(label="traj_2")]
        mocker.patch.object(job, "_prepare_all_jobs", return_value=mock_jobs)

        mock_batch_cls = mocker.patch(
            "chemsmart.jobs.gaussian.traj.GaussianBatchJob"
        )
        mock_batch = mock_batch_cls.return_value

        job._run_all_jobs()

        mock_batch_cls.assert_called_once()
        call_kwargs = mock_batch_cls.call_args.kwargs
        assert call_kwargs["jobs"] == mock_jobs
        assert call_kwargs["no_run_in_parallel"] is False
        assert call_kwargs["fail_fast"] is False
        assert call_kwargs["label"] == "test_traj_batch"
        assert call_kwargs["jobrunner"] == gaussian_jobrunner_no_scratch
        mock_batch.run.assert_called_once()

    def test_traj_job_serial_mode_uses_batch_without_fail_fast(
        self, pbs_server, gaussian_jobrunner_no_scratch, mocker
    ):
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings
        from chemsmart.jobs.gaussian.traj import GaussianTrajJob

        mock_molecules = [MockMolecule(), MockMolecule(), MockMolecule()]
        settings = GaussianJobSettings()

        gaussian_jobrunner_no_scratch.no_run_in_parallel = True
        job = GaussianTrajJob(
            molecules=mock_molecules,
            settings=settings,
            label="test_traj_serial",
            jobrunner=gaussian_jobrunner_no_scratch,
            proportion_structures_to_use=1.0,
        )

        mock_jobs = [Mock(label="traj_1"), Mock(label="traj_2")]
        mocker.patch.object(job, "_prepare_all_jobs", return_value=mock_jobs)

        mock_batch_cls = mocker.patch(
            "chemsmart.jobs.gaussian.traj.GaussianBatchJob"
        )
        mock_batch = mock_batch_cls.return_value

        job._run_all_jobs()

        mock_batch_cls.assert_called_once()
        call_kwargs = mock_batch_cls.call_args.kwargs
        assert call_kwargs["jobs"] == mock_jobs
        assert call_kwargs["no_run_in_parallel"] is True
        assert call_kwargs["fail_fast"] is False
        assert call_kwargs["label"] == "test_traj_serial_batch"
        mock_batch.run.assert_called_once()

    def test_crest_job_serial_mode_uses_batch_without_fail_fast(
        self, pbs_server, gaussian_jobrunner_no_scratch, mocker
    ):
        from chemsmart.jobs.gaussian.crest import GaussianCrestJob
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings

        mock_molecules = [MockMolecule(), MockMolecule()]
        settings = GaussianJobSettings()

        gaussian_jobrunner_no_scratch.no_run_in_parallel = True
        job = GaussianCrestJob(
            molecules=mock_molecules,
            settings=settings,
            label="test_crest_serial",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        mock_jobs = [Mock(label="conf_1"), Mock(label="conf_2")]
        mocker.patch.object(job, "_prepare_all_jobs", return_value=mock_jobs)

        mock_batch_cls = mocker.patch(
            "chemsmart.jobs.gaussian.crest.GaussianBatchJob"
        )
        mock_batch = mock_batch_cls.return_value

        job._run_all_jobs()

        mock_batch_cls.assert_called_once()
        call_kwargs = mock_batch_cls.call_args.kwargs
        assert call_kwargs["no_run_in_parallel"] is True
        assert call_kwargs["fail_fast"] is False
        mock_batch.run.assert_called_once()


class TestRunListFailureAggregation:
    """Legacy list path used by pKa should fail like BatchJob."""

    def test_parallel_list_execution_raises_aggregated_failures(
        self, pbs_server, mocker
    ):
        import click

        from chemsmart.cli.run import process_pipeline, run
        from chemsmart.jobs.batch import BatchExecutionError

        jobrunner = JobRunner(
            server=pbs_server, fake=True, no_run_in_parallel=False
        )
        ctx = click.Context(run)
        ctx.ensure_object(dict)
        ctx.obj["jobrunner"] = jobrunner

        ok_job = Mock(label="ok_job", TYPE="g16opt")
        ok_job.run.return_value = None
        fail_job = Mock(label="fail_job", TYPE="g16opt")
        fail_job.run.side_effect = RuntimeError("boom")

        mocker.patch.object(JobRunner, "from_job", return_value=jobrunner)

        with pytest.raises(BatchExecutionError, match="fail_job") as exc_info:
            process_pipeline.__wrapped__(ctx, [ok_job, fail_job])

        assert "1 of 2 list job(s) failed" in str(exc_info.value)
        ok_job.run.assert_called_once()
        fail_job.run.assert_called_once()
