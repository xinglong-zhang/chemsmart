"""Tests for JobRunner class and run_in_serial flag."""

from unittest.mock import Mock

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.runner import JobRunner


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


class TestSerialExecution:
    """Test serial execution enforcement for jobs with subjobs."""

    def test_link_job_serial_execution_stops_on_incomplete(
        self, pbs_server, gaussian_jobrunner_no_scratch, mocker
    ):
        """Test that GaussianLinkJob stops serial execution if a job doesn't complete."""
        from chemsmart.jobs.gaussian.link import GaussianLinkJob
        from chemsmart.jobs.gaussian.settings import GaussianLinkJobSettings

        # Create a mock molecule
        mock_molecule = MockMolecule()

        # Create settings for IRC job
        settings = GaussianLinkJobSettings()
        settings.jobtype = "irc"
        settings.direction = "both"

        # Create job with run_in_serial=True
        gaussian_jobrunner_no_scratch.run_in_serial = True
        job = GaussianLinkJob(
            molecule=mock_molecule,
            settings=settings,
            label="test_irc",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        # Mock the IRC jobs and their behavior
        mock_irc_job1 = Mock()
        mock_irc_job1.settings.jobtype = "ircf"
        mock_irc_job1.is_complete.return_value = False  # First job fails

        mock_irc_job2 = Mock()
        mock_irc_job2.settings.jobtype = "ircr"
        mock_irc_job2.is_complete.return_value = True

        # Patch _get_irc_jobs to return our mocked jobs
        mocker.patch.object(
            job, "_get_irc_jobs", return_value=[mock_irc_job1, mock_irc_job2]
        )

        # Run the job
        job._run()

        # Verify first job was run
        mock_irc_job1.run.assert_called_once()

        # Verify second job was NOT run (stopped after first failure)
        mock_irc_job2.run.assert_not_called()

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

    def test_qrc_job_default_behavior_runs_all(
        self, pbs_server, gaussian_jobrunner_no_scratch, mocker
    ):
        """Test that GaussianQRCJob runs all jobs when run_in_serial is False."""
        from chemsmart.jobs.gaussian.qrc import GaussianQRCJob
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings

        # Create mock molecule
        mock_molecule = MockMolecule()

        # Create settings
        settings = GaussianJobSettings()

        # Create job with run_in_serial=False (default)
        gaussian_jobrunner_no_scratch.run_in_serial = False
        job = GaussianQRCJob(
            molecule=mock_molecule,
            settings=settings,
            label="test_qrc",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        # Mock QRC jobs - even if one fails, all should run in default mode
        mock_job1 = Mock()
        mock_job1.label = "qrc_f"
        mock_job1.is_complete.return_value = False

        mock_job2 = Mock()
        mock_job2.label = "qrc_r"
        mock_job2.is_complete.return_value = True

        # Patch _prepare_both_qrc_jobs instead of the property
        mocker.patch.object(
            job, "_prepare_both_qrc_jobs", return_value=[mock_job1, mock_job2]
        )

        # Run the job
        job._run()

        # Verify both jobs were run despite first one not completing
        mock_job1.run.assert_called_once()
        mock_job2.run.assert_called_once()


class TestBatchJobRefactor:
    """Regression tests for the shared BatchJob orchestration layer."""

    def test_split_jobs_across_nodes_balances_chunks(self):
        from chemsmart.jobs.batch import BatchJob

        jobs = [1, 2, 3, 4, 5]
        chunks = BatchJob._split_jobs_across_nodes(jobs, 3)

        assert chunks == [[1, 2], [3, 4], [5]]


class TestOrcaQRCBatchExecution:
    """Tests for ORCA QRC execution via the shared batch layer."""

    def test_orca_qrc_job_uses_batch_parallel(
        self, pbs_server, orca_jobrunner_no_scratch, mocker
    ):
        from chemsmart.jobs.orca.qrc import ORCAQRCJob
        from chemsmart.jobs.orca.settings import ORCAJobSettings

        mock_molecule = MockMolecule()
        settings = ORCAJobSettings(functional="B3LYP", basis="def2-SVP")
        settings.jobtype = "qrc"

        orca_jobrunner_no_scratch.run_in_serial = False
        job = ORCAQRCJob(
            molecule=mock_molecule,
            settings=settings,
            label="test_orca_qrc",
            jobrunner=orca_jobrunner_no_scratch,
        )

        mock_job1 = Mock()
        mock_job1.label = "orca_qrc_f"
        mock_job2 = Mock()
        mock_job2.label = "orca_qrc_r"

        mocker.patch.object(
            job,
            "_prepare_both_qrc_jobs",
            return_value=[mock_job1, mock_job2],
        )

        mock_batch_cls = mocker.patch("chemsmart.jobs.orca.qrc.OrcaBatchJob")
        mock_batch = mock_batch_cls.return_value

        job._run_both_jobs()

        mock_batch_cls.assert_called_once()
        call_kwargs = mock_batch_cls.call_args[1]
        assert call_kwargs["run_in_serial"] is False
        assert len(call_kwargs["jobs"]) == 2
        assert call_kwargs["label"] == "test_orca_qrc_batch"
        assert call_kwargs["jobrunner"] == orca_jobrunner_no_scratch
        mock_batch.run.assert_called_once()

    def test_orca_qrc_job_serial_execution_stops_on_incomplete(
        self, pbs_server, orca_jobrunner_no_scratch, mocker
    ):
        from chemsmart.jobs.orca.qrc import ORCAQRCJob
        from chemsmart.jobs.orca.settings import ORCAJobSettings

        mock_molecule = MockMolecule()
        settings = ORCAJobSettings(functional="B3LYP", basis="def2-SVP")
        settings.jobtype = "qrc"

        orca_jobrunner_no_scratch.run_in_serial = True
        job = ORCAQRCJob(
            molecule=mock_molecule,
            settings=settings,
            label="test_orca_qrc_serial",
            jobrunner=orca_jobrunner_no_scratch,
        )

        mock_job1 = Mock()
        mock_job1.label = "orca_qrc_f"
        mock_job1.is_complete.return_value = False

        mock_job2 = Mock()
        mock_job2.label = "orca_qrc_r"
        mock_job2.is_complete.return_value = True

        mocker.patch.object(
            job,
            "_prepare_both_qrc_jobs",
            return_value=[mock_job1, mock_job2],
        )

        job._run_both_jobs()

        mock_job1.run.assert_called_once()
        mock_job2.run.assert_not_called()
