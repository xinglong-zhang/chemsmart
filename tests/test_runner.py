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


class TestGaussianPkaReferenceExecution:
    """Regression tests for Gaussian pKa reference workflow support."""

    @staticmethod
    def _write_reference_xyz(path):
        path.write_text("2\nreference\nC 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")

    @staticmethod
    def _build_gaussian_pka_settings(reference_file):
        from chemsmart.jobs.gaussian.settings import GaussianpKaJobSettings

        return GaussianpKaJobSettings(
            proton_index=2,
            scheme="proton exchange",
            reference_file=str(reference_file),
            reference_proton_index=2,
            reference_charge=0,
            reference_multiplicity=1,
            charge=0,
            multiplicity=1,
            functional="B3LYP",
            basis="6-31G(d)",
        )

    def test_gaussian_pka_job_prepares_reference_opt_jobs(
        self, tmp_path, pbs_server, gaussian_jobrunner_no_scratch
    ):
        from chemsmart.jobs.gaussian.pka import GaussianpKaJob

        reference_file = tmp_path / "reference.xyz"
        self._write_reference_xyz(reference_file)
        settings = self._build_gaussian_pka_settings(reference_file)
        target_molecule = Molecule(
            symbols=["C", "H"],
            positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.1]],
            charge=0,
            multiplicity=1,
        )

        job = GaussianpKaJob(
            molecule=target_molecule,
            settings=settings,
            label="test_gaussian_pka_ref",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        assert job.has_reference_jobs is True
        assert len(job.ref_opt_jobs) == 2
        assert job.ref_acid_job is job.ref_opt_jobs[0]
        assert job.ref_conjugate_base_job is job.ref_opt_jobs[1]
        assert job.ref_acid_job.label == "test_gaussian_pka_ref_HRef_opt"
        assert (
            job.ref_conjugate_base_job.label == "test_gaussian_pka_ref_Ref_opt"
        )
        assert list(job.ref_acid_job.molecule.symbols) == ["C", "H"]
        assert list(job.ref_conjugate_base_job.molecule.symbols) == ["C"]

    def test_gaussian_pka_run_ref_sp_jobs_creates_reference_sp_jobs(
        self, tmp_path, pbs_server, gaussian_jobrunner_no_scratch, mocker
    ):
        from chemsmart.jobs.gaussian.pka import GaussianpKaJob

        reference_file = tmp_path / "reference.xyz"
        self._write_reference_xyz(reference_file)
        settings = self._build_gaussian_pka_settings(reference_file)
        target_molecule = Molecule(
            symbols=["C", "H"],
            positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.1]],
            charge=0,
            multiplicity=1,
        )

        job = GaussianpKaJob(
            molecule=target_molecule,
            settings=settings,
            label="test_gaussian_pka_ref_sp",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        href_optimized = Molecule(
            symbols=["N", "H"],
            positions=[[1.0, 0.0, 0.0], [1.0, 0.0, 1.0]],
            charge=0,
            multiplicity=1,
        )
        ref_optimized = Molecule(
            symbols=["N"],
            positions=[[2.0, 0.0, 0.0]],
            charge=-1,
            multiplicity=1,
        )

        mocker.patch.object(
            job, "_ref_opt_jobs_are_complete", return_value=True
        )
        mocker.patch.object(
            job.ref_acid_job,
            "_output",
            return_value=Mock(
                normal_termination=True,
                molecule=href_optimized,
            ),
        )
        mocker.patch.object(
            job.ref_conjugate_base_job,
            "_output",
            return_value=Mock(
                normal_termination=True,
                molecule=ref_optimized,
            ),
        )
        run_spy = mocker.patch(
            "chemsmart.jobs.gaussian.pka.GaussianSinglePointJob.run"
        )

        job._run_ref_sp_jobs()

        assert len(job.ref_sp_jobs) == 2
        assert job.ref_acid_sp_job is job.ref_sp_jobs[0]
        assert job.ref_conjugate_base_sp_job is job.ref_sp_jobs[1]
        assert job.ref_acid_sp_job.label == "test_gaussian_pka_ref_sp_HRef_sp"
        assert (
            job.ref_conjugate_base_sp_job.label
            == "test_gaussian_pka_ref_sp_Ref_sp"
        )
        assert list(job.ref_acid_sp_job.molecule.symbols) == list(
            href_optimized.symbols
        )
        assert list(job.ref_conjugate_base_sp_job.molecule.symbols) == list(
            ref_optimized.symbols
        )
        assert job.ref_acid_sp_job.settings.charge == settings.reference_charge
        assert (
            job.ref_conjugate_base_sp_job.settings.charge
            == settings.reference_charge - 1
        )
        assert run_spy.call_count == 2

    def test_gaussian_pka_parallel_stops_before_sp_on_opt_failure(
        self, tmp_path, pbs_server, gaussian_jobrunner_no_scratch, mocker
    ):
        from chemsmart.jobs.gaussian.pka import GaussianpKaJob

        reference_file = tmp_path / "reference.xyz"
        self._write_reference_xyz(reference_file)
        settings = self._build_gaussian_pka_settings(reference_file)
        target_molecule = Molecule(
            symbols=["C", "H"],
            positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.1]],
            charge=0,
            multiplicity=1,
        )

        gaussian_jobrunner_no_scratch.run_in_serial = False
        job = GaussianpKaJob(
            molecule=target_molecule,
            settings=settings,
            label="test_gaussian_pka_parallel_gate",
            jobrunner=gaussian_jobrunner_no_scratch,
            parallel=True,
        )

        mocker.patch.object(
            job, "_run_opt_worker", return_value={"success": True}
        )
        mocker.patch.object(
            job,
            "_collect_futures",
            return_value=([], ["opt phase failure"]),
        )
        create_sp_spy = mocker.patch.object(job, "_create_sp_jobs")

        with pytest.raises(RuntimeError, match="Optimization phase failed"):
            job._run_parallel()

        create_sp_spy.assert_not_called()

    def test_dias_job_uses_batch_parallel_per_stage(
        self, pbs_server, gaussian_jobrunner_no_scratch, mocker
    ):
        from chemsmart.jobs.gaussian.dias import GaussianDIASJob
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings

        mock_molecules = [MockMolecule(), MockMolecule(), MockMolecule()]
        settings = GaussianJobSettings()

        gaussian_jobrunner_no_scratch.run_in_serial = False
        job = GaussianDIASJob(
            molecules=mock_molecules,
            settings=settings,
            label="test_dias",
            jobrunner=gaussian_jobrunner_no_scratch,
            fragment_indices="1",
            every_n_points=1,
            mode="irc",
        )

        molecule_jobs = [Mock(label="mol_1"), Mock(label="mol_2")]
        fragment1_jobs = [Mock(label="f1_1")]
        fragment2_jobs = [Mock(label="f2_1")]

        mocker.patch.object(
            type(job),
            "all_molecules_jobs",
            new_callable=PropertyMock,
            return_value=molecule_jobs,
        )
        mocker.patch.object(
            type(job),
            "fragment1_jobs",
            new_callable=PropertyMock,
            return_value=fragment1_jobs,
        )
        mocker.patch.object(
            type(job),
            "fragment2_jobs",
            new_callable=PropertyMock,
            return_value=fragment2_jobs,
        )

        mock_batch_cls = mocker.patch(
            "chemsmart.jobs.gaussian.dias.GaussianBatchJob"
        )

        job._run()

        assert mock_batch_cls.call_count == 3
        first_call = mock_batch_cls.call_args_list[0][1]
        second_call = mock_batch_cls.call_args_list[1][1]
        third_call = mock_batch_cls.call_args_list[2][1]

        assert first_call["jobs"] == molecule_jobs
        assert first_call["label"] == "test_dias_molecules_batch"
        assert second_call["jobs"] == fragment1_jobs
        assert second_call["label"] == "test_dias_fragment1_batch"
        assert third_call["jobs"] == fragment2_jobs
        assert third_call["label"] == "test_dias_fragment2_batch"

        for call in mock_batch_cls.call_args_list:
            assert call[1]["run_in_serial"] is False
            assert call[1]["jobrunner"] == gaussian_jobrunner_no_scratch

        assert mock_batch_cls.return_value.run.call_count == 3

    def test_dias_job_serial_fragment1_stops_on_incomplete(
        self, pbs_server, gaussian_jobrunner_no_scratch, mocker
    ):
        from chemsmart.jobs.gaussian.dias import GaussianDIASJob
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings

        mock_molecules = [MockMolecule(), MockMolecule(), MockMolecule()]
        settings = GaussianJobSettings()

        gaussian_jobrunner_no_scratch.run_in_serial = True
        job = GaussianDIASJob(
            molecules=mock_molecules,
            settings=settings,
            label="test_dias_serial",
            jobrunner=gaussian_jobrunner_no_scratch,
            fragment_indices="1",
            every_n_points=1,
            mode="irc",
        )

        mock_job1 = Mock()
        mock_job1.label = "f1_1"
        mock_job1.is_complete.return_value = False

        mock_job2 = Mock()
        mock_job2.label = "f1_2"
        mock_job2.is_complete.return_value = True

        mocker.patch.object(
            type(job),
            "fragment1_jobs",
            new_callable=PropertyMock,
            return_value=[mock_job1, mock_job2],
        )

        job._run_fragment1_jobs()

        mock_job1.run.assert_called_once()
        mock_job2.run.assert_not_called()


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


class TestOrcaPkaBatchExecution:
    """Tests for ORCA pKa staged execution via ORCApKaBatchJob."""

    def test_orca_pka_job_parallel_mode_dispatches_to_run_parallel(
        self, pbs_server, orca_jobrunner_no_scratch, mocker
    ):
        from chemsmart.jobs.orca.pka import ORCApKaJob
        from chemsmart.jobs.orca.settings import ORCApKaJobSettings

        mock_molecule = MockMolecule()
        settings = ORCApKaJobSettings(
            proton_index=1,
            scheme="direct",
            charge=0,
            multiplicity=1,
            functional="B3LYP",
            basis="def2-SVP",
        )
        settings.scheme = "proton exchange"
        settings.reference_file = "reference_acid.xyz"

        orca_jobrunner_no_scratch.run_in_serial = False
        job = ORCApKaJob(
            molecule=mock_molecule,
            settings=settings,
            label="test_orca_pka",
            jobrunner=orca_jobrunner_no_scratch,
            parallel=True,
        )

        run_parallel_spy = mocker.patch.object(job, "_run_parallel")
        run_opt_spy = mocker.patch.object(job, "_run_opt_jobs")
        run_sp_spy = mocker.patch.object(job, "_run_sp_jobs")

        job._run()

        run_parallel_spy.assert_called_once()
        run_opt_spy.assert_not_called()
        run_sp_spy.assert_not_called()

    def test_orca_pka_job_serial_stops_before_sp_when_opt_incomplete(
        self, pbs_server, orca_jobrunner_no_scratch, mocker
    ):
        from chemsmart.jobs.orca.pka import ORCApKaJob
        from chemsmart.jobs.orca.settings import ORCApKaJobSettings

        mock_molecule = MockMolecule()
        settings = ORCApKaJobSettings(
            proton_index=1,
            scheme="direct",
            charge=0,
            multiplicity=1,
            functional="B3LYP",
            basis="def2-SVP",
        )

        orca_jobrunner_no_scratch.run_in_serial = True
        job = ORCApKaJob(
            molecule=mock_molecule,
            settings=settings,
            label="test_orca_pka_serial",
            jobrunner=orca_jobrunner_no_scratch,
            parallel=True,
        )

        opt_job1 = Mock(label="opt_ha")
        opt_job1.is_complete.return_value = False
        opt_job2 = Mock(label="opt_a")
        opt_job2.is_complete.return_value = True

        mocker.patch.object(
            type(job),
            "opt_jobs",
            new_callable=PropertyMock,
            return_value=[opt_job1, opt_job2],
        )
        run_opt_spy = mocker.patch.object(job, "_run_opt_jobs")
        run_sp_spy = mocker.patch.object(job, "_run_sp_jobs")
        run_parallel_spy = mocker.patch.object(job, "_run_parallel")

        job._run()

        # Parallel mode is disabled when run_in_serial=True.
        assert job.parallel is False
        run_parallel_spy.assert_not_called()
        run_opt_spy.assert_called_once()
        run_sp_spy.assert_not_called()

    def test_orca_pka_parallel_reference_sp_jobs_use_reference_payloads(
        self, pbs_server, orca_jobrunner_no_scratch, mocker
    ):
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.orca.pka import ORCApKaJob
        from chemsmart.jobs.orca.settings import (
            ORCAJobSettings,
            ORCApKaJobSettings,
        )

        target_ha_fallback = Molecule(
            symbols=["H"],
            positions=[[0.0, 0.0, 0.0]],
            charge=0,
            multiplicity=1,
        )
        target_a_fallback = Molecule(
            symbols=["He"],
            positions=[[1.0, 0.0, 0.0]],
            charge=-1,
            multiplicity=1,
        )
        ref_hb_fallback = Molecule(
            symbols=["Li"],
            positions=[[2.0, 0.0, 0.0]],
            charge=1,
            multiplicity=1,
        )
        ref_b_fallback = Molecule(
            symbols=["Be"],
            positions=[[3.0, 0.0, 0.0]],
            charge=0,
            multiplicity=1,
        )

        target_ha_optimized = Molecule(
            symbols=["B"],
            positions=[[4.0, 0.0, 0.0]],
            charge=0,
            multiplicity=1,
        )
        target_a_optimized = Molecule(
            symbols=["C"],
            positions=[[5.0, 0.0, 0.0]],
            charge=-1,
            multiplicity=1,
        )
        ref_hb_optimized = Molecule(
            symbols=["N"],
            positions=[[6.0, 0.0, 0.0]],
            charge=1,
            multiplicity=1,
        )
        ref_b_optimized = Molecule(
            symbols=["O"],
            positions=[[7.0, 0.0, 0.0]],
            charge=0,
            multiplicity=1,
        )

        settings = ORCApKaJobSettings(
            proton_index=1,
            scheme="direct",
            charge=0,
            multiplicity=1,
            functional="B3LYP",
            basis="def2-SVP",
        )

        orca_jobrunner_no_scratch.run_in_serial = False
        job = ORCApKaJob(
            molecule=target_ha_fallback,
            settings=settings,
            label="test_orca_pka_parallel_refs",
            jobrunner=orca_jobrunner_no_scratch,
            parallel=True,
        )

        # Mock opt jobs with vibrational_frequencies (all positive → no imaginary)
        ha_opt_job = Mock(label="ha_opt")
        ha_opt_job._output.return_value = Mock(
            normal_termination=True,
            molecule=target_ha_optimized,
            vibrational_frequencies=[100.0, 200.0],
        )
        a_opt_job = Mock(label="a_opt")
        a_opt_job._output.return_value = Mock(
            normal_termination=True,
            molecule=target_a_optimized,
            vibrational_frequencies=[100.0, 200.0],
        )
        hb_opt_job = Mock(label="hb_opt")
        hb_opt_job._output.return_value = Mock(
            normal_termination=True,
            molecule=ref_hb_optimized,
            vibrational_frequencies=[100.0, 200.0],
        )
        b_opt_job = Mock(label="b_opt")
        b_opt_job._output.return_value = Mock(
            normal_termination=True,
            molecule=ref_b_optimized,
            vibrational_frequencies=[100.0, 200.0],
        )

        ha_sp_settings = ORCAJobSettings(functional="B3LYP", basis="def2-SVP")
        a_sp_settings = ORCAJobSettings(functional="PBE0", basis="def2-SVP")
        hb_sp_settings = ORCAJobSettings(
            functional="M06-2X", basis="def2-TZVP"
        )
        b_sp_settings = ORCAJobSettings(
            functional="wB97X-D", basis="ma-def2-SVP"
        )

        mocker.patch.object(
            type(job),
            "has_reference_jobs",
            new_callable=PropertyMock,
            return_value=True,
        )
        mocker.patch.object(
            type(job),
            "opt_jobs",
            new_callable=PropertyMock,
            return_value=[ha_opt_job, a_opt_job],
        )
        mocker.patch.object(
            type(job),
            "ref_opt_jobs",
            new_callable=PropertyMock,
            return_value=[hb_opt_job, b_opt_job],
        )
        mocker.patch.object(
            type(job),
            "protonated_molecule",
            new_callable=PropertyMock,
            return_value=target_ha_fallback,
        )
        mocker.patch.object(
            type(job),
            "conjugate_base_molecule",
            new_callable=PropertyMock,
            return_value=target_a_fallback,
        )
        mocker.patch.object(
            type(job),
            "reference_molecule",
            new_callable=PropertyMock,
            return_value=ref_hb_fallback,
        )
        mocker.patch.object(
            type(job),
            "reference_conjugate_base_molecule",
            new_callable=PropertyMock,
            return_value=ref_b_fallback,
        )
        mocker.patch.object(
            job.settings,
            "_create_solution_phase_sp_settings",
            return_value=(ha_sp_settings, a_sp_settings),
        )
        mocker.patch.object(
            job.settings,
            "reference_pair_sp_job_settings",
            return_value=(hb_sp_settings, b_sp_settings),
        )
        mocker.patch("chemsmart.jobs.orca.pka.ORCASinglePointJob.run")

        job._run_parallel()

        assert [sp_job.label for sp_job in job.sp_jobs] == [
            f"{job._acid_basename}_sp",
            f"{job._conjugate_base_label}_sp",
        ]
        assert [sp_job.label for sp_job in job.ref_sp_jobs] == [
            f"{job._ref_basename}_sp",
            f"{job._ref_conjugate_base_label}_sp",
        ]

        assert job.sp_jobs[0].molecule.symbols == target_ha_optimized.symbols
        assert job.sp_jobs[1].molecule.symbols == target_a_optimized.symbols
        assert job.ref_sp_jobs[0].molecule.symbols == ref_hb_optimized.symbols
        assert job.ref_sp_jobs[1].molecule.symbols == ref_b_optimized.symbols

        assert job.sp_jobs[0].settings.functional == ha_sp_settings.functional
        assert job.sp_jobs[1].settings.functional == a_sp_settings.functional
        assert (
            job.ref_sp_jobs[0].settings.functional == hb_sp_settings.functional
        )
        assert (
            job.ref_sp_jobs[1].settings.functional == b_sp_settings.functional
        )
        assert job.sp_jobs[0].settings.basis == ha_sp_settings.basis
        assert job.sp_jobs[1].settings.basis == a_sp_settings.basis
        assert job.ref_sp_jobs[0].settings.basis == hb_sp_settings.basis
        assert job.ref_sp_jobs[1].settings.basis == b_sp_settings.basis

        assert job.ref_sp_jobs[0].label != job.sp_jobs[0].label
        assert job.ref_sp_jobs[1].label != job.sp_jobs[1].label

    def test_orca_pka_parallel_stops_before_sp_on_opt_failure(
        self, pbs_server, orca_jobrunner_no_scratch, mocker
    ):
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.orca.pka import ORCApKaJob
        from chemsmart.jobs.orca.settings import ORCApKaJobSettings

        mock_molecule = Molecule(
            symbols=["H", "C"],
            positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.1]],
            charge=0,
            multiplicity=1,
        )
        settings = ORCApKaJobSettings(
            proton_index=1,
            scheme="direct",
            charge=0,
            multiplicity=1,
            functional="B3LYP",
            basis="def2-SVP",
        )

        orca_jobrunner_no_scratch.run_in_serial = False
        job = ORCApKaJob(
            molecule=mock_molecule,
            settings=settings,
            label="test_orca_pka_parallel_gate",
            jobrunner=orca_jobrunner_no_scratch,
            parallel=True,
        )

        mocker.patch.object(
            job, "_run_opt_worker", return_value={"success": True}
        )
        mocker.patch.object(
            job,
            "_collect_futures",
            return_value=([], ["opt phase failure"]),
        )
        prepare_sp_spy = mocker.patch.object(
            job,
            "_prepare_sp_jobs",
            return_value=(),
        )

        with pytest.raises(RuntimeError, match="Optimization phase failed"):
            job._run_parallel()

        prepare_sp_spy.assert_not_called()
