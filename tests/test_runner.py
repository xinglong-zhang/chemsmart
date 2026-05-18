"""Tests for JobRunner class and run_in_serial flag."""

from unittest.mock import Mock, PropertyMock

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

    def test_jobrunner_default_run_in_serial(self, pbs_server):
        """Test that run_in_serial defaults to False."""
        runner = JobRunner(server=pbs_server, fake=True)
        assert hasattr(runner, "run_in_serial")
        assert runner.run_in_serial is False

    def test_jobrunner_run_in_serial_true(self, pbs_server):
        """Test that run_in_serial can be set to True."""
        runner = JobRunner(server=pbs_server, fake=True, run_in_serial=True)
        assert runner.run_in_serial is True

    def test_jobrunner_run_in_serial_false(self, pbs_server):
        """Test that run_in_serial can be set to False."""
        runner = JobRunner(server=pbs_server, fake=True, run_in_serial=False)
        assert runner.run_in_serial is False

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

        # Create job with run_in_serial=True
        gaussian_jobrunner_no_scratch.run_in_serial = True
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

    def test_batch_serial_mode_when_unset(self, pbs_server):
        """Test that BatchJob runs all jobs when they are unset."""
        dummy_batch_cls = self._dummy_batch_cls()
        runner = JobRunner(server=pbs_server, fake=True, run_in_serial=True)

        batch = dummy_batch_cls(
            jobs=[],
            jobrunner=runner,
        )

        assert batch.run_in_serial is True

    def test_batch_serial_mode_keeps_explicit_false(self, pbs_server):
        """Test that BatchJob runs all jobs when they are unset."""
        dummy_batch_cls = self._dummy_batch_cls()
        runner = JobRunner(server=pbs_server, fake=True, run_in_serial=True)

        batch = dummy_batch_cls(
            jobs=[],
            run_in_serial=False,
            jobrunner=runner,
        )

        assert batch.run_in_serial is False

    def test_batch_serial_mode_keeps_explicit_true(self, pbs_server):
        dummy_batch_cls = self._dummy_batch_cls()
        runner = JobRunner(server=pbs_server, fake=True, run_in_serial=False)

        batch = dummy_batch_cls(
            jobs=[],
            run_in_serial=True,
            jobrunner=runner,
        )

        assert batch.run_in_serial is True

    def test_batch_writes_success_and_failed_logs(self, pbs_server, tmp_path):
        from chemsmart.jobs.batch import BatchExecutionError

        dummy_batch_cls = self._dummy_batch_cls()
        runner = JobRunner(server=pbs_server, fake=True, run_in_serial=False)

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
            run_in_serial=False,
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
        runner = JobRunner(server=pbs_server, fake=True, run_in_serial=True)

        fail_job = Mock()
        fail_job.label = "fail_job"
        fail_job.run.side_effect = RuntimeError("fail")
        fail_job.is_complete.return_value = False

        batch = dummy_batch_cls(
            jobs=[fail_job],
            run_in_serial=True,
            jobrunner=runner,
            label="batch_raise_test",
        )
        batch.folder = str(tmp_path)

        with pytest.raises(BatchExecutionError, match="fail_job"):
            batch.run()

    def test_batch_run_multi_node_records_node_future_exception(
        self, pbs_server, mocker
    ):
        dummy_batch_cls = self._dummy_batch_cls()
        runner = JobRunner(server=pbs_server, fake=True, run_in_serial=False)

        job1 = Mock()
        job1.label = "job1"
        job2 = Mock()
        job2.label = "job2"

        batch = dummy_batch_cls(
            jobs=[job1, job2],
            run_in_serial=False,
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

        gaussian_jobrunner_no_scratch.run_in_serial = False
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
        call_kwargs = mock_batch_cls.call_args[1]
        assert call_kwargs["jobs"] == mock_jobs
        assert call_kwargs["run_in_serial"] is False
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

        gaussian_jobrunner_no_scratch.run_in_serial = False
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
        call_kwargs = mock_batch_cls.call_args[1]
        assert call_kwargs["jobs"] == mock_jobs
        assert call_kwargs["run_in_serial"] is False
        assert call_kwargs["label"] == "test_traj_batch"
        assert call_kwargs["jobrunner"] == gaussian_jobrunner_no_scratch
        mock_batch.run.assert_called_once()

    def test_traj_job_serial_execution_stops_on_incomplete(
        self, pbs_server, gaussian_jobrunner_no_scratch, mocker
    ):
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings
        from chemsmart.jobs.gaussian.traj import GaussianTrajJob

        mock_molecules = [MockMolecule(), MockMolecule(), MockMolecule()]
        settings = GaussianJobSettings()

        gaussian_jobrunner_no_scratch.run_in_serial = True
        job = GaussianTrajJob(
            molecules=mock_molecules,
            settings=settings,
            label="test_traj_serial",
            jobrunner=gaussian_jobrunner_no_scratch,
            proportion_structures_to_use=1.0,
        )

        mock_job1 = Mock()
        mock_job1.label = "traj_1"
        mock_job1.is_complete.return_value = False

        mock_job2 = Mock()
        mock_job2.label = "traj_2"
        mock_job2.is_complete.return_value = True

        mocker.patch.object(
            job,
            "_prepare_all_jobs",
            return_value=[mock_job1, mock_job2],
        )

        job._run_all_jobs()

        mock_job1.run.assert_called_once()
        mock_job2.run.assert_not_called()
