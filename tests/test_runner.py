"""Tests for JobRunner class and no_run_in_parallel flag."""

from unittest.mock import Mock

import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.runner import (
    JobRunner,
    run_phase_jobs,
)


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

    def test_run_phase_jobs_passes_stop_on_incomplete(
        self, pbs_server, mocker
    ):
        from chemsmart.jobs.job import Job

        runner = JobRunner(server=pbs_server, fake=True)
        mock_child = Mock()
        mock_child.run.return_value = None
        mock_child.is_complete.return_value = True
        mock_execute = mocker.patch.object(Job, "_execute_phase_jobs")

        run_phase_jobs(
            parent_runner=runner,
            jobs=[mock_child],
            stop_on_incomplete=True,
            phase_label="test phase",
        )

        mock_execute.assert_called_once()
        assert mock_execute.call_args.kwargs["stop_on_incomplete"] is True


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
    """Tests for shared BatchJob orchestration."""

    @staticmethod
    def _dummy_batch_cls():
        from chemsmart.jobs.batch import BatchJob

        class DummyBatchJob(BatchJob):
            PROGRAM = "dummy"

        return DummyBatchJob

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

    def test_parallel_request_falls_back_to_serial(self, pbs_server):
        """In-process parallel is disabled; children run serially."""
        from chemsmart.jobs.batch import BatchExecutionError

        dummy_batch_cls = self._dummy_batch_cls()
        runner = JobRunner(
            server=pbs_server,
            fake=True,
            num_cores=16,
            mem_gb=32,
        )

        fail_job = Mock(label="fail_job")
        fail_job.run.side_effect = RuntimeError("fail")
        fail_job.is_complete.return_value = False

        later_job = Mock(label="later_job")
        later_job.run.return_value = None
        later_job.is_complete.return_value = True

        batch = dummy_batch_cls(
            jobs=[fail_job, later_job],
            fail_fast=True,
            jobrunner=runner,
        )

        with pytest.raises(
            BatchExecutionError, match="1 attempted, 1 failed, 1 not started"
        ):
            batch.run()

        fail_job.run.assert_called_once()
        later_job.run.assert_not_called()
        assert fail_job.jobrunner.num_cores == 16
        assert fail_job.jobrunner.mem_gb == 32

    def test_run_child_jobs_as_batch_always_serial(self, pbs_server):
        """Nested batches always run children serially with full resources."""
        from chemsmart.jobs.batch import run_child_jobs_as_batch

        dummy_batch_cls = self._dummy_batch_cls()
        runner = JobRunner(
            server=pbs_server,
            fake=True,
            num_cores=16,
            mem_gb=32,
        )
        parent = Mock()
        parent.label = "parent"
        parent.jobrunner = runner

        child_a = Mock(label="child_a")
        child_a.run.return_value = None
        child_a.is_complete.return_value = True
        child_b = Mock(label="child_b")
        child_b.run.return_value = None
        child_b.is_complete.return_value = True

        batch = run_child_jobs_as_batch(
            batch_cls=dummy_batch_cls,
            jobs=[child_a, child_b],
            parent=parent,
            label_suffix="_batch",
            fail_fast=False,
        )

        assert batch.fail_fast is False
        assert batch.label == "parent_batch"
        child_a.run.assert_called_once()
        child_b.run.assert_called_once()
        assert child_a.jobrunner.num_cores == 16
        assert child_a.jobrunner.mem_gb == 32
        assert child_b.jobrunner.num_cores == 16
        assert child_b.jobrunner.mem_gb == 32


class TestGaussianBatchDelegation:
    """Tests for Gaussian multi-subjob workflows using GaussianBatchJob."""

    def test_crest_job_nested_batch_always_serial(
        self, pbs_server, gaussian_jobrunner_no_scratch, mocker
    ):
        from chemsmart.jobs.gaussian.crest import GaussianCrestJob
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings

        mock_molecules = [MockMolecule(), MockMolecule(), MockMolecule()]
        settings = GaussianJobSettings()

        # Nested crest children run serially regardless of runner policy.
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
        assert call_kwargs["fail_fast"] is False
        assert call_kwargs["label"] == "test_crest_batch"
        assert call_kwargs["jobrunner"] == gaussian_jobrunner_no_scratch
        mock_batch.run.assert_called_once()

    def test_traj_job_nested_batch_always_serial(
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
        assert call_kwargs["fail_fast"] is False
        assert call_kwargs["label"] == "test_traj_batch"
        assert call_kwargs["jobrunner"] == gaussian_jobrunner_no_scratch
        mock_batch.run.assert_called_once()

    def test_traj_array_task_maps_stable_end_slice(
        self, pbs_server, gaussian_jobrunner_no_scratch, mocker, monkeypatch
    ):
        """SLURM task id indexes the stable last-N set, not incompletes."""
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings
        from chemsmart.jobs.gaussian.traj import GaussianTrajJob

        monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "2")
        molecules = [MockMolecule() for _ in range(5)]
        for index, mol in enumerate(molecules):
            mol.energy = float(index)

        settings = GaussianJobSettings()
        job = GaussianTrajJob(
            molecules=molecules,
            settings=settings,
            label="traj_array",
            jobrunner=gaussian_jobrunner_no_scratch,
            proportion_structures_to_use=1.0,
            num_structures_to_run=3,
            skip_completed=False,
        )

        prepared = []
        for index in range(5):
            child = Mock(label=f"traj_array_c{index + 1}")
            child.run.return_value = None
            child.is_complete.return_value = index == 3  # c4 complete
            prepared.append(child)
        mocker.patch.object(job, "_prepare_all_jobs", return_value=prepared)
        mock_batch_cls = mocker.patch(
            "chemsmart.jobs.gaussian.traj.GaussianBatchJob"
        )

        job._run_all_jobs()

        # Stable last-3 is c3,c4,c5; task 2 runs c4 only.
        prepared[0].run.assert_not_called()
        prepared[1].run.assert_not_called()
        prepared[2].run.assert_not_called()
        prepared[3].run.assert_called_once()
        prepared[4].run.assert_not_called()
        mock_batch_cls.return_value.run.assert_not_called()

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
        assert call_kwargs["fail_fast"] is False
        mock_batch.run.assert_called_once()

    def test_qrc_job_nested_batch_always_serial(
        self, pbs_server, gaussian_jobrunner_no_scratch, mocker
    ):
        from chemsmart.jobs.gaussian.qrc import GaussianQRCJob
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings

        settings = GaussianJobSettings()
        gaussian_jobrunner_no_scratch.no_run_in_parallel = False
        job = GaussianQRCJob(
            molecule=MockMolecule(),
            settings=settings,
            label="test_qrc",
            jobrunner=gaussian_jobrunner_no_scratch,
        )

        mock_jobs = [Mock(label="qrc_f"), Mock(label="qrc_r")]
        mocker.patch.object(
            type(job),
            "both_qrc_jobs",
            new_callable=mocker.PropertyMock,
            return_value=mock_jobs,
        )

        mock_batch_cls = mocker.patch(
            "chemsmart.jobs.gaussian.qrc.GaussianBatchJob"
        )
        mock_batch = mock_batch_cls.return_value

        job._run_both_jobs()

        mock_batch_cls.assert_called_once()
        call_kwargs = mock_batch_cls.call_args.kwargs
        assert call_kwargs["jobs"] == mock_jobs
        assert call_kwargs["fail_fast"] is False
        assert call_kwargs["label"] == "test_qrc_batch"
        mock_batch.run.assert_called_once()

    def test_qrc_array_task_runs_selected_child_only(
        self, pbs_server, gaussian_jobrunner_no_scratch, mocker, monkeypatch
    ):
        """SLURM_ARRAY_TASK_ID selects one QRC child via nestable helper."""
        from chemsmart.jobs.gaussian.qrc import GaussianQRCJob
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings

        monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "2")
        gaussian_jobrunner_no_scratch.num_cores = 16
        gaussian_jobrunner_no_scratch.mem_gb = 32

        settings = GaussianJobSettings()
        job = GaussianQRCJob(
            molecule=MockMolecule(),
            settings=settings,
            label="test_qrc_array",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        child_f = Mock(label="qrc_f")
        child_f.run.return_value = None
        child_f.is_complete.return_value = True
        child_r = Mock(label="qrc_r")
        child_r.run.return_value = None
        child_r.is_complete.return_value = True
        mocker.patch.object(
            type(job),
            "both_qrc_jobs",
            new_callable=mocker.PropertyMock,
            return_value=[child_f, child_r],
        )
        mock_batch_cls = mocker.patch(
            "chemsmart.jobs.gaussian.qrc.GaussianBatchJob"
        )

        job._run_both_jobs()

        child_f.run.assert_not_called()
        child_r.run.assert_called_once()
        assert child_r.jobrunner is not gaussian_jobrunner_no_scratch
        assert child_r.jobrunner.num_cores == 16
        assert child_r.jobrunner.mem_gb == 32
        mock_batch_cls.return_value.run.assert_not_called()

    def test_qrc_array_task_raises_when_selected_child_incomplete(
        self, pbs_server, gaussian_jobrunner_no_scratch, mocker, monkeypatch
    ):
        """Nestable array path fails if the selected child did not complete."""
        from chemsmart.jobs.batch import BatchExecutionError
        from chemsmart.jobs.gaussian.qrc import GaussianQRCJob
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings

        monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "1")
        settings = GaussianJobSettings()
        job = GaussianQRCJob(
            molecule=MockMolecule(),
            settings=settings,
            label="test_qrc_incomplete",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        child_f = Mock(label="qrc_f")
        child_f.run.return_value = None
        child_f.is_complete.return_value = False
        child_r = Mock(label="qrc_r")
        child_r.run.return_value = None
        child_r.is_complete.return_value = True
        mocker.patch.object(
            type(job),
            "both_qrc_jobs",
            new_callable=mocker.PropertyMock,
            return_value=[child_f, child_r],
        )

        with pytest.raises(
            BatchExecutionError, match="incomplete after execution"
        ):
            job._run_both_jobs()

        child_f.run.assert_called_once()
        child_r.run.assert_not_called()

    def test_multi_molecule_orca_qrc_outer_array_runs_both_directions(
        self, pbs_server, orca_jobrunner_no_scratch, mocker, monkeypatch
    ):
        """Outer SLURM array must not steal nestable QRC child selection.

        ``chemsmart sub ... orca -i 1,2,3 qrc`` yields ORCABatchJob of three
        ORCAQRCJobs. Each array task runs one molecule; that QRC must still
        run both forward and reverse. Task id 3 previously raised ValueError
        against the two nested QRC children.
        """
        import os

        from chemsmart.jobs.orca.batch import ORCABatchJob
        from chemsmart.jobs.orca.qrc import ORCAQRCJob
        from chemsmart.jobs.orca.settings import ORCAJobSettings

        monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "3")
        settings = ORCAJobSettings()
        qrc_parents = []
        nested_pairs = []
        for index in range(3):
            parent = ORCAQRCJob(
                molecule=MockMolecule(),
                settings=settings,
                label=f"mol{index}_qrc",
                jobrunner=orca_jobrunner_no_scratch,
                skip_completed=False,
            )
            forward = Mock(label=f"mol{index}_f")
            forward.run.return_value = None
            forward.is_complete.return_value = True
            reverse = Mock(label=f"mol{index}_r")
            reverse.run.return_value = None
            reverse.is_complete.return_value = True
            nested = [forward, reverse]
            nested_pairs.append(nested)
            mocker.patch.object(
                parent, "_prepare_both_qrc_jobs", return_value=nested
            )
            qrc_parents.append(parent)

        batch = ORCABatchJob(
            jobs=qrc_parents,
            label="mols_qrc_batch",
            jobrunner=orca_jobrunner_no_scratch,
        )
        batch.run()

        # Outer array selected molecule 3 only.
        for index, nested in enumerate(nested_pairs):
            if index == 2:
                nested[0].run.assert_called_once()
                nested[1].run.assert_called_once()
            else:
                nested[0].run.assert_not_called()
                nested[1].run.assert_not_called()

        # Outer array env restored after the selected child finishes.
        assert os.environ.get("SLURM_ARRAY_TASK_ID") == "3"


class TestRunListFailureAggregation:
    """List execution aggregates failures like BatchJob."""

    def test_serial_list_execution_raises_aggregated_failures(
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

        from chemsmart.jobs.job import Job

        ok_job = Mock(spec=Job, label="ok_job", TYPE="g16opt")
        ok_job.run.return_value = None
        fail_job = Mock(spec=Job, label="fail_job", TYPE="g16opt")
        fail_job.run.side_effect = RuntimeError("boom")

        mocker.patch.object(JobRunner, "from_job", return_value=jobrunner)

        with pytest.warns(DeprecationWarning, match="bare list of Job"):
            with pytest.raises(
                BatchExecutionError, match="fail_job"
            ) as exc_info:
                process_pipeline.__wrapped__(ctx, [ok_job, fail_job])

        assert "1 of 2 list job(s) failed" in str(exc_info.value)
        ok_job.run.assert_called_once()
        fail_job.run.assert_called_once()

    def test_sub_list_emits_deprecation_and_submits_individually(
        self, pbs_server, monkeypatch
    ):
        import click

        from chemsmart.cli.sub import process_pipeline, sub
        from chemsmart.jobs.job import Job
        from chemsmart.settings.server import Server

        jobrunner = JobRunner(
            server=pbs_server, fake=True, no_run_in_parallel=True
        )
        ctx = click.Context(sub)
        ctx.ensure_object(dict)
        ctx.obj["jobrunner"] = jobrunner
        ctx.obj["subcommand"] = [
            {
                "name": "sub",
                "kwargs": {
                    "server": "dummy",
                    "test": True,
                    "print_command": False,
                    "time_hours": None,
                    "queue": None,
                    "verbose": False,
                },
                "args": (),
                "params": {},
            }
        ]

        captured = {"submissions": []}
        fake_server = Server(name="dummy")

        def _fake_submit(job, test=False, cli_args=None, **kw):
            captured["submissions"].append((job, test, cli_args))

        fake_server.submit = _fake_submit
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )
        monkeypatch.setattr(
            "chemsmart.cli.sub.CtxObjArguments.reconstruct_command_line",
            lambda self: ["sub", "gaussian", "opt"],
        )

        job_a = Mock(spec=Job, label="a", TYPE="g16opt")
        job_b = Mock(spec=Job, label="b", TYPE="g16opt")

        with pytest.warns(DeprecationWarning, match="bare list of Job"):
            process_pipeline.__wrapped__(
                ctx,
                [job_a, job_b],
                server="dummy",
                test=True,
                print_command=False,
            )

        assert len(captured["submissions"]) == 2
        assert captured["submissions"][0][0] is job_a
        assert captured["submissions"][1][0] is job_b
        assert job_a.jobrunner is jobrunner
        assert job_b.jobrunner is jobrunner

    def test_run_forces_batchjob_serial_despite_parallel_flag(
        self, pbs_server, mocker
    ):
        """chemsmart run serializes top-level BatchJob children."""
        import click

        from chemsmart.cli.run import process_pipeline, run
        from chemsmart.jobs.batch import BatchJob

        jobrunner = JobRunner(
            server=pbs_server,
            fake=True,
            num_cores=16,
            mem_gb=32,
        )
        ctx = click.Context(run)
        ctx.ensure_object(dict)
        ctx.obj["jobrunner"] = jobrunner

        from chemsmart.jobs.job import Job

        child = Mock(spec=Job, label="child", TYPE="g16opt")
        child.run.return_value = None
        child.is_complete.return_value = True

        class _DummyBatch(BatchJob):
            PROGRAM = "test"

        batch = _DummyBatch(
            jobs=[child],
            jobrunner=jobrunner,
            label="top_batch",
        )
        mocker.patch.object(JobRunner, "from_job", return_value=jobrunner)

        process_pipeline.__wrapped__(ctx, batch)

        child.run.assert_called_once()
        assert child.jobrunner.num_cores == 16
        assert child.jobrunner.mem_gb == 32

    def test_serial_batch_propagates_full_cores_and_memory(self, pbs_server):
        """Serial BatchJob children receive full parent num_cores / mem_gb."""
        from chemsmart.jobs.batch import BatchJob

        runner = JobRunner(
            server=pbs_server,
            fake=True,
            num_cores=16,
            mem_gb=32,
        )

        child_a = Mock(label="child_a")
        child_a.run.return_value = None
        child_a.is_complete.return_value = True
        child_b = Mock(label="child_b")
        child_b.run.return_value = None
        child_b.is_complete.return_value = True

        class _DummyBatch(BatchJob):
            PROGRAM = "test"

        batch = _DummyBatch(
            jobs=[child_a, child_b],
            jobrunner=runner,
            label="resource_batch",
        )
        batch.run()

        assert child_a.jobrunner.num_cores == 16
        assert child_a.jobrunner.mem_gb == 32
        assert child_b.jobrunner.num_cores == 16
        assert child_b.jobrunner.mem_gb == 32
        assert child_a.jobrunner is not child_b.jobrunner
        assert child_a.jobrunner is not runner


class TestBatchSerialExecutionPolicy:
    """Serial local/nested BatchJob policy with full per-child resources."""

    @staticmethod
    def _dummy_batch_cls():
        from chemsmart.jobs.batch import BatchJob

        class DummyBatchJob(BatchJob):
            PROGRAM = "dummy"

        return DummyBatchJob

    def test_toplevel_run_batch_is_serial_with_full_resources(
        self, pbs_server, mocker
    ):
        """Top-level chemsmart run BatchJob: serial children, full cores/mem."""
        import click

        from chemsmart.cli.run import process_pipeline, run
        from chemsmart.jobs.job import Job

        jobrunner = JobRunner(
            server=pbs_server,
            fake=True,
            num_cores=16,
            mem_gb=32,
        )
        ctx = click.Context(run)
        ctx.ensure_object(dict)
        ctx.obj["jobrunner"] = jobrunner

        call_order = []
        cores_seen = []
        mem_seen = []
        children = []
        for index in range(4):
            child = Mock(spec=Job, label=f"mol_{index}", TYPE="g16opt")
            child.is_complete.return_value = True

            def _run(child_ref=child):
                call_order.append(child_ref.label)
                cores_seen.append(child_ref.jobrunner.num_cores)
                mem_seen.append(child_ref.jobrunner.mem_gb)

            child.run.side_effect = _run
            children.append(child)

        batch_cls = self._dummy_batch_cls()
        batch = batch_cls(
            jobs=children,
            jobrunner=jobrunner,
            label="mols_batch",
        )
        mocker.patch.object(JobRunner, "from_job", return_value=jobrunner)

        process_pipeline.__wrapped__(ctx, batch)

        assert call_order == ["mol_0", "mol_1", "mol_2", "mol_3"]
        assert cores_seen == [16, 16, 16, 16]
        assert mem_seen == [32, 32, 32, 32]

    def test_run_in_parallel_does_not_oversubscribe_cores(
        self, pbs_server, mocker
    ):
        """--run-in-parallel on run must not assign N× full cores or split cores."""
        import click

        from chemsmart.cli.run import process_pipeline, run
        from chemsmart.jobs.job import Job

        parent_cores = 16
        parent_mem = 32
        num_children = 4
        jobrunner = JobRunner(
            server=pbs_server,
            fake=True,
            num_cores=parent_cores,
            mem_gb=parent_mem,
        )
        ctx = click.Context(run)
        ctx.ensure_object(dict)
        ctx.obj["jobrunner"] = jobrunner

        children = []
        for index in range(num_children):
            child = Mock(spec=Job, label=f"job_{index}", TYPE="g16opt")
            child.run.return_value = None
            child.is_complete.return_value = True
            children.append(child)

        batch_cls = self._dummy_batch_cls()
        batch = batch_cls(
            jobs=children,
            jobrunner=jobrunner,
            label="parallel_flag_batch",
        )
        mocker.patch.object(JobRunner, "from_job", return_value=jobrunner)

        process_pipeline.__wrapped__(ctx, batch)

        split_cores = parent_cores // num_children
        for child in children:
            assert child.jobrunner.num_cores == parent_cores
            assert child.jobrunner.mem_gb == parent_mem
            assert child.jobrunner.num_cores != split_cores

    def test_nested_crest_qrc_serial_despite_parallel_runner(
        self, pbs_server, gaussian_jobrunner_no_scratch, mocker
    ):
        """Nested crest/QRC batches stay serial when runner allows parallel."""
        from chemsmart.jobs.gaussian.crest import GaussianCrestJob
        from chemsmart.jobs.gaussian.qrc import GaussianQRCJob
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings
        from chemsmart.jobs.orca.qrc import ORCAQRCJob
        from chemsmart.jobs.orca.settings import ORCAJobSettings

        gaussian_jobrunner_no_scratch.no_run_in_parallel = False
        settings = GaussianJobSettings()

        crest = GaussianCrestJob(
            molecules=[MockMolecule(), MockMolecule()],
            settings=settings,
            label="crest_nested",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        crest_jobs = [Mock(label="c1"), Mock(label="c2")]
        mocker.patch.object(
            crest, "_prepare_all_jobs", return_value=crest_jobs
        )
        crest._run_all_jobs()

        qrc = GaussianQRCJob(
            molecule=MockMolecule(),
            settings=settings,
            label="qrc_nested",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        qrc_jobs = [Mock(label="qf"), Mock(label="qr")]
        mocker.patch.object(
            type(qrc),
            "both_qrc_jobs",
            new_callable=mocker.PropertyMock,
            return_value=qrc_jobs,
        )
        qrc._run_both_jobs()

        orca_settings = ORCAJobSettings()
        orca_qrc = ORCAQRCJob(
            molecule=MockMolecule(),
            settings=orca_settings,
            label="orca_qrc_nested",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        orca_jobs = [Mock(label="of"), Mock(label="or")]
        mocker.patch.object(
            type(orca_qrc),
            "both_qrc_jobs",
            new_callable=mocker.PropertyMock,
            return_value=orca_jobs,
        )
        orca_qrc._run_both_jobs()

    def test_single_child_keeps_full_engine_resources(self, pbs_server):
        """A single child still receives the parent full core/memory allocation."""
        runner = JobRunner(
            server=pbs_server,
            fake=True,
            num_cores=64,
            mem_gb=128,
        )
        child = Mock(label="only_child")
        child.run.return_value = None
        child.is_complete.return_value = True

        batch = self._dummy_batch_cls()(
            jobs=[child],
            jobrunner=runner,
            label="single_child_batch",
        )
        batch.run()

        assert child.jobrunner.num_cores == 64
        assert child.jobrunner.mem_gb == 128

    def test_batchjob_has_no_inprocess_parallel_runner(self):
        """Unsplit in-process parallel path is removed from BatchJob."""
        from chemsmart.jobs.batch import BatchJob

        assert "_run_jobs_in_parallel" not in vars(BatchJob)


class TestBatchExecutionModes:
    """local_batch vs array_task mode resolution for BatchJob.run()."""

    @staticmethod
    def _dummy_batch_cls():
        from chemsmart.jobs.batch import BatchJob

        class DummyBatchJob(BatchJob):
            PROGRAM = "dummy"

        return DummyBatchJob

    def test_batch_execution_mode_enum_values(self):
        from chemsmart.jobs.batch import BatchExecutionMode

        assert BatchExecutionMode.LOCAL_BATCH.value == "local_batch"
        assert BatchExecutionMode.ARRAY_TASK.value == "array_task"
        assert not hasattr(BatchExecutionMode, "MULTI_NODE")

    def test_resolve_batch_execution_mode_defaults_to_local_batch(
        self, monkeypatch
    ):
        from chemsmart.jobs.batch import (
            BatchExecutionMode,
            resolve_batch_execution_mode,
        )

        for key in ("SLURM_ARRAY_TASK_ID", "PBS_ARRAYID", "LSB_JOBINDEX"):
            monkeypatch.delenv(key, raising=False)

        assert resolve_batch_execution_mode() is BatchExecutionMode.LOCAL_BATCH

    def test_resolve_batch_execution_mode_array_task_from_slurm(
        self, monkeypatch
    ):
        from chemsmart.jobs.batch import (
            BatchExecutionMode,
            resolve_array_task_id,
            resolve_batch_execution_mode,
        )

        monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "2")
        assert resolve_array_task_id() == 2
        assert resolve_batch_execution_mode() is BatchExecutionMode.ARRAY_TASK

    def test_array_task_logs_execution_mode(
        self, pbs_server, monkeypatch, caplog
    ):
        import logging

        from chemsmart.jobs.batch import BatchExecutionMode

        monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "2")
        runner = JobRunner(
            server=pbs_server,
            fake=True,
            num_cores=16,
            mem_gb=32,
        )
        children = []
        for index in range(4):
            child = Mock(label=f"child_{index}")
            child.run.return_value = None
            child.is_complete.return_value = True
            children.append(child)

        batch = self._dummy_batch_cls()(
            jobs=children,
            jobrunner=runner,
            label="mols_batch",
        )
        with caplog.at_level(logging.INFO):
            batch.run()

        assert any(
            f"execution={BatchExecutionMode.ARRAY_TASK.value}"
            in record.message
            and "task=2/4" in record.message
            and "cores=16" in record.message
            for record in caplog.records
        )

    def test_array_task_runs_only_selected_child(
        self, pbs_server, monkeypatch
    ):
        """SLURM_ARRAY_TASK_ID=2 runs child index 1 with full resources."""
        monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "2")

        runner = JobRunner(
            server=pbs_server,
            fake=True,
            num_cores=16,
            mem_gb=32,
        )
        children = []
        for index in range(4):
            child = Mock(label=f"child_{index}")
            child.run.return_value = None
            child.is_complete.return_value = True
            children.append(child)

        batch = self._dummy_batch_cls()(
            jobs=children,
            jobrunner=runner,
            label="array_batch",
        )
        batch.run()

        children[0].run.assert_not_called()
        children[1].run.assert_called_once()
        children[2].run.assert_not_called()
        children[3].run.assert_not_called()
        assert children[1].jobrunner.num_cores == 16
        assert children[1].jobrunner.mem_gb == 32

    def test_array_task_clears_scheduler_env_during_child_run(
        self, pbs_server, monkeypatch
    ):
        """Outer array task id must not be visible inside child.run()."""
        import os

        from chemsmart.jobs.batch import resolve_array_task_id

        monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "2")
        seen = {}

        def _record_env(**kwargs):
            seen["task_id"] = resolve_array_task_id()
            seen["slurm"] = os.environ.get("SLURM_ARRAY_TASK_ID")

        runner = JobRunner(
            server=pbs_server, fake=True, no_run_in_parallel=True
        )
        children = []
        for index in range(3):
            child = Mock(label=f"child_{index}")
            child.run.side_effect = _record_env if index == 1 else None
            child.is_complete.return_value = True
            children.append(child)

        batch = self._dummy_batch_cls()(
            jobs=children,
            jobrunner=runner,
            label="array_batch",
        )
        batch.run()

        assert seen["task_id"] is None
        assert seen["slurm"] is None
        assert os.environ.get("SLURM_ARRAY_TASK_ID") == "2"

    def test_array_task_rejects_out_of_range_task_id(
        self, pbs_server, monkeypatch
    ):
        monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "5")
        runner = JobRunner(server=pbs_server, fake=True)
        children = [
            Mock(
                label="child_0",
                run=Mock(),
                is_complete=Mock(return_value=True),
            ),
            Mock(
                label="child_1",
                run=Mock(),
                is_complete=Mock(return_value=True),
            ),
        ]
        batch = self._dummy_batch_cls()(
            jobs=children,
            jobrunner=runner,
            label="array_batch",
        )
        with pytest.raises(ValueError, match="out of range"):
            batch.run()

    def test_local_batch_without_array_env_runs_all_children(
        self, pbs_server, monkeypatch
    ):
        for key in ("SLURM_ARRAY_TASK_ID", "PBS_ARRAYID", "LSB_JOBINDEX"):
            monkeypatch.delenv(key, raising=False)

        runner = JobRunner(
            server=pbs_server, fake=True, no_run_in_parallel=True
        )
        children = []
        for index in range(3):
            child = Mock(label=f"child_{index}")
            child.run.return_value = None
            child.is_complete.return_value = True
            children.append(child)

        batch = self._dummy_batch_cls()(
            jobs=children,
            jobrunner=runner,
            label="local_batch",
        )
        batch.run()

        for child in children:
            child.run.assert_called_once()

    def test_local_batch_logs_serial_policy(
        self, pbs_server, monkeypatch, caplog
    ):
        import logging

        from chemsmart.jobs.batch import BatchExecutionMode

        for key in ("SLURM_ARRAY_TASK_ID", "PBS_ARRAYID", "LSB_JOBINDEX"):
            monkeypatch.delenv(key, raising=False)

        runner = JobRunner(
            server=pbs_server,
            fake=True,
            num_cores=8,
            mem_gb=16,
        )
        children = [
            Mock(label="a", run=Mock(), is_complete=Mock(return_value=True)),
            Mock(label="b", run=Mock(), is_complete=Mock(return_value=True)),
            Mock(label="c", run=Mock(), is_complete=Mock(return_value=True)),
            Mock(label="d", run=Mock(), is_complete=Mock(return_value=True)),
        ]
        for child in children:
            child.run.return_value = None

        batch = self._dummy_batch_cls()(
            jobs=children,
            jobrunner=runner,
            label="mols_batch",
        )
        with caplog.at_level(logging.INFO):
            batch.run()

        assert any(
            f"execution={BatchExecutionMode.LOCAL_BATCH.value}"
            in record.message
            and "children=4" in record.message
            and "policy=serial" in record.message
            for record in caplog.records
        )

    def test_nested_serial_logs_serial_nested_policy(
        self, pbs_server, monkeypatch, caplog
    ):
        import logging

        from chemsmart.jobs.batch import BatchExecutionMode

        for key in ("SLURM_ARRAY_TASK_ID", "PBS_ARRAYID", "LSB_JOBINDEX"):
            monkeypatch.delenv(key, raising=False)

        runner = JobRunner(
            server=pbs_server, fake=True, no_run_in_parallel=True
        )
        children = []
        for index in range(8):
            child = Mock(label=f"conf_{index}")
            child.run.return_value = None
            child.is_complete.return_value = True
            children.append(child)

        batch = self._dummy_batch_cls()(
            jobs=children,
            jobrunner=runner,
            label="crest_children",
            nested_serial=True,
        )
        with caplog.at_level(logging.INFO):
            batch.run()

        assert any(
            f"execution={BatchExecutionMode.LOCAL_BATCH.value}"
            in record.message
            and "children=8" in record.message
            and "policy=serial_nested" in record.message
            for record in caplog.records
        )

    def test_run_child_jobs_as_batch_sets_nested_serial(
        self, pbs_server, monkeypatch
    ):
        from chemsmart.jobs.batch import run_child_jobs_as_batch

        for key in ("SLURM_ARRAY_TASK_ID", "PBS_ARRAYID", "LSB_JOBINDEX"):
            monkeypatch.delenv(key, raising=False)

        parent = Mock()
        parent.label = "crest_parent"
        parent.jobrunner = JobRunner(
            server=pbs_server, fake=True, no_run_in_parallel=True
        )
        children = []
        for index in range(2):
            child = Mock(label=f"c{index}")
            child.run.return_value = None
            child.is_complete.return_value = True
            children.append(child)

        batch = run_child_jobs_as_batch(
            batch_cls=self._dummy_batch_cls(),
            jobs=children,
            parent=parent,
            label_suffix="_children",
        )
        assert batch.nested_serial is True
        assert batch.label == "crest_parent_children"

    def test_slurm_nodelist_still_runs_serial_local_batch(
        self, pbs_server, monkeypatch, caplog
    ):
        """Multi-node SLURM allocation does not enable in-process fan-out."""
        import logging

        from chemsmart.jobs.batch import BatchExecutionMode

        for key in ("SLURM_ARRAY_TASK_ID", "PBS_ARRAYID", "LSB_JOBINDEX"):
            monkeypatch.delenv(key, raising=False)
        monkeypatch.setenv("SLURM_JOB_NODELIST", "nodeA,nodeB")

        runner = JobRunner(
            server=pbs_server, fake=True, no_run_in_parallel=True
        )
        children = []
        for index in range(2):
            child = Mock(label=f"child_{index}")
            child.run.return_value = None
            child.is_complete.return_value = True
            children.append(child)

        batch = self._dummy_batch_cls()(
            jobs=children,
            jobrunner=runner,
            label="multi_node_batch",
        )
        with caplog.at_level(logging.INFO):
            batch.run()

        assert any(
            f"execution={BatchExecutionMode.LOCAL_BATCH.value}"
            in record.message
            and "children=2" in record.message
            for record in caplog.records
        )
        assert not any(
            "execution=multi_node" in record.message
            for record in caplog.records
        )
        children[0].run.assert_called_once()
        children[1].run.assert_called_once()

    def test_pbs_arrayid_selects_array_task_mode(
        self, pbs_server, monkeypatch
    ):
        monkeypatch.delenv("SLURM_ARRAY_TASK_ID", raising=False)
        monkeypatch.setenv("PBS_ARRAYID", "1")

        runner = JobRunner(server=pbs_server, fake=True)
        child_a = Mock(label="child_a")
        child_a.run.return_value = None
        child_a.is_complete.return_value = True
        child_b = Mock(label="child_b")
        child_b.run.return_value = None
        child_b.is_complete.return_value = True

        batch = self._dummy_batch_cls()(
            jobs=[child_a, child_b],
            jobrunner=runner,
            label="pbs_array_batch",
        )
        batch.run()

        child_a.run.assert_called_once()
        child_b.run.assert_not_called()
