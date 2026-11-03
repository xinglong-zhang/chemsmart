"""Tests for PySCF job objects and the fake job runner.

The fake runner writes the real driver script and a synthetic results file
without importing pyscf, so these tests exercise the
settings -> writer -> output -> parse path in an environment where PySCF is
not installed.
"""

import glob
import os
import sys
from unittest.mock import patch

import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.io.pyscf.output import PySCFOutput
from chemsmart.jobs.pyscf.hess import PySCFHessJob
from chemsmart.jobs.pyscf.job import PySCFJob
from chemsmart.jobs.pyscf.opt import PySCFOptJob
from chemsmart.jobs.pyscf.runner import PySCFJobRunner
from chemsmart.jobs.pyscf.settings import PySCFJobSettings
from chemsmart.jobs.pyscf.singlepoint import PySCFSinglePointJob
from chemsmart.jobs.pyscf.writer import PySCFScriptWriter
from chemsmart.settings.executable import PySCFExecutable
from tests.helpers import make_pyscf_settings


class TestPySCFJobConstruction:
    """Job identity, file paths, and constructor guards."""

    def test_job_file_paths_follow_the_label(
        self, pyscf_water_molecule, temporary_working_dir
    ):
        job = PySCFSinglePointJob(
            molecule=pyscf_water_molecule,
            settings=make_pyscf_settings("sp"),
            label="water_sp",
        )
        folder = str(temporary_working_dir)
        assert job.PROGRAM == "PySCF"
        assert job.inputfile == os.path.join(folder, "water_sp.py")
        assert job.outputfile == os.path.join(folder, "water_sp.out")
        assert job.resultsfile == os.path.join(folder, "water_sp.h5")

    def test_job_without_results_is_not_complete(
        self, pyscf_water_molecule, temporary_working_dir
    ):
        job = PySCFSinglePointJob(
            molecule=pyscf_water_molecule,
            settings=make_pyscf_settings("sp"),
            label="water_sp",
        )
        assert not os.path.exists(job.resultsfile)
        assert job.is_complete() is False

    def test_settings_must_be_pyscf_settings(
        self, pyscf_water_molecule, temporary_working_dir
    ):
        with pytest.raises(ValueError, match="Settings must be instance of"):
            PySCFSinglePointJob(
                molecule=pyscf_water_molecule,
                settings={"functional": "b3lyp"},
                label="water_sp",
            )

    def test_molecule_must_be_a_molecule(self, temporary_working_dir):
        with pytest.raises(ValueError):
            PySCFSinglePointJob(
                molecule="water",
                settings=make_pyscf_settings("sp"),
                label="water_sp",
            )


class TestFakePySCFJobRunner:
    """The preview runner must produce a parseable, clearly fake result."""

    @staticmethod
    def _run(job_class, jobtype, molecule, jobrunner, **extra):
        job = job_class(
            molecule=molecule,
            settings=make_pyscf_settings(jobtype, **extra),
            label=f"water_{jobtype}",
            jobrunner=jobrunner,
        )
        assert jobrunner.run(job) == 0
        return job

    def test_preview_writes_script_log_and_results(
        self,
        pyscf_water_molecule,
        pyscf_jobrunner_no_scratch,
        temporary_working_dir,
    ):
        job = self._run(
            PySCFSinglePointJob,
            "sp",
            pyscf_water_molecule,
            pyscf_jobrunner_no_scratch,
        )
        assert os.path.exists(job.inputfile)
        assert os.path.exists(job.outputfile)
        assert os.path.exists(job.resultsfile)

    def test_preview_result_is_marked_fake(
        self,
        pyscf_water_molecule,
        pyscf_jobrunner_no_scratch,
        temporary_working_dir,
    ):
        job = self._run(
            PySCFSinglePointJob,
            "sp",
            pyscf_water_molecule,
            pyscf_jobrunner_no_scratch,
        )
        output = PySCFOutput(job.outputfile)
        # A preview must never be mistaken for a calculation: the committed
        # real artifacts report the producing PySCF version instead.
        assert output.version == "fake"
        assert output.engine == "cpu"

    @pytest.mark.parametrize(
        ("job_class", "jobtype", "extra"),
        [
            (PySCFOptJob, "opt", {}),
            (PySCFHessJob, "hess", {"freq": True}),
        ],
    )
    def test_preview_generates_a_driver_script_per_jobtype(
        self,
        job_class,
        jobtype,
        extra,
        pyscf_water_molecule,
        pyscf_jobrunner_no_scratch,
        temporary_working_dir,
    ):
        job = self._run(
            job_class,
            jobtype,
            pyscf_water_molecule,
            pyscf_jobrunner_no_scratch,
            **extra,
        )
        with open(job.inputfile) as handle:
            script = handle.read()
        assert "CHEMSMART PySCF driver" in script
        assert "def2-svp" in script

    def test_preview_in_scratch_copies_results_back(
        self,
        pyscf_water_molecule,
        pyscf_jobrunner_scratch,
        temporary_working_dir,
    ):
        job = self._run(
            PySCFSinglePointJob,
            "sp",
            pyscf_water_molecule,
            pyscf_jobrunner_scratch,
        )
        # The calculation runs in the scratch directory, but the artifacts
        # must end up next to the job.
        assert os.path.exists(job.resultsfile)
        assert os.path.exists(job.outputfile)
        assert PySCFOutput(job.outputfile).version == "fake"

    def test_preview_records_the_requested_solvent(
        self,
        pyscf_water_molecule,
        pyscf_jobrunner_no_scratch,
        temporary_working_dir,
    ):
        job = self._run(
            PySCFSinglePointJob,
            "sp",
            pyscf_water_molecule,
            pyscf_jobrunner_no_scratch,
            solvent_model="cpcm",
            solvent_id="water",
        )
        output = PySCFOutput(job.outputfile)
        assert output.solvent_on
        assert output.solvent_model == "cpcm"
        assert output.solvent_id == "water"


class TestPySCFRerun:
    """A rerun must never report the results of an earlier run."""

    def test_rerun_that_dies_does_not_report_the_old_result(
        self,
        pyscf_water_molecule,
        pbs_server,
        pyscf_jobrunner_no_scratch,
        temporary_working_dir,
    ):
        first = PySCFSinglePointJob(
            molecule=pyscf_water_molecule,
            settings=make_pyscf_settings("sp"),
            label="rerun",
            jobrunner=pyscf_jobrunner_no_scratch,
            skip_completed=False,
        )
        first.run()
        assert first.is_complete()

        changed = make_pyscf_settings("sp")
        changed.functional = "pbe"
        rerun = PySCFSinglePointJob(
            molecule=pyscf_water_molecule,
            settings=changed,
            label="rerun",
            jobrunner=PySCFJobRunner(server=pbs_server, scratch=False),
            skip_completed=False,
        )
        dies = f"{sys.executable} -c 'raise SystemExit(1)'"
        with patch.object(PySCFJobRunner, "_get_command", return_value=dies):
            rerun.run()
        assert not rerun.is_complete()
        assert not os.path.exists(rerun.resultsfile)


class TestPySCFForcesRequest:
    """A gradient is computed only when its cost is justified."""

    @pytest.mark.parametrize(
        ("job_class", "jobtype", "extra", "expected"),
        [
            (PySCFSinglePointJob, "sp", {}, False),
            (PySCFSinglePointJob, "sp", {"forces": True}, True),
            (PySCFHessJob, "hess", {}, False),
            (PySCFOptJob, "opt", {}, True),
        ],
    )
    def test_forces_follow_the_request(
        self,
        job_class,
        jobtype,
        extra,
        expected,
        pyscf_water_molecule,
        temporary_working_dir,
    ):
        job = job_class(
            molecule=pyscf_water_molecule,
            settings=make_pyscf_settings(jobtype, **extra),
            label="forces",
        )
        assert PySCFScriptWriter(job=job).build_config()["forces"] is expected


def test_optimized_structure_of_a_completed_opt_job(
    pyscf_water_molecule, pyscf_jobrunner_no_scratch, temporary_working_dir
):
    job = PySCFOptJob(
        molecule=pyscf_water_molecule,
        settings=make_pyscf_settings("opt", freq=False),
        label="opt",
        jobrunner=pyscf_jobrunner_no_scratch,
    )
    job.run()
    assert job.optimized_structure().is_optimized_structure


class TestPySCFJobFactories:
    def test_from_filename_reads_a_real_structure(
        self,
        pyscf_water_xyz,
        pyscf_jobrunner_no_scratch,
        temporary_working_dir,
    ):
        job = PySCFSinglePointJob.from_filename(
            pyscf_water_xyz,
            settings=make_pyscf_settings("sp"),
            label="from_file",
            jobrunner=pyscf_jobrunner_no_scratch,
        )
        assert list(job.molecule.chemical_symbols) == ["O", "H", "H"]
        assert job.label == "from_file"

    def test_from_pubchem_uses_the_fetched_structure(
        self,
        pyscf_water_molecule,
        pyscf_jobrunner_no_scratch,
        temporary_working_dir,
    ):
        with patch.object(
            Molecule, "from_pubchem", return_value=pyscf_water_molecule
        ) as fetch:
            job = PySCFSinglePointJob.from_pubchem(
                "962",
                settings=make_pyscf_settings("sp"),
                label="fetched",
                jobrunner=pyscf_jobrunner_no_scratch,
            )
        fetch.assert_called_once_with(identifier="962")
        assert job.molecule.num_atoms == 3


class TestPySCFFactoryRunnerSelection:
    def test_job_from_a_file_can_run_on_its_selected_runner(
        self, pyscf_water_xyz, pbs_server, temporary_working_dir
    ):
        job = PySCFSinglePointJob.from_filename(
            pyscf_water_xyz,
            settings=make_pyscf_settings("sp"),
            label="selected",
            server=pbs_server,
            fake=True,
        )
        job.run()
        assert job.is_complete()

    def test_given_runner_is_kept(
        self,
        pyscf_water_xyz,
        pbs_server,
        pyscf_jobrunner_no_scratch,
        temporary_working_dir,
    ):
        job = PySCFSinglePointJob.from_filename(
            pyscf_water_xyz,
            settings=make_pyscf_settings("sp"),
            label="given",
            jobrunner=pyscf_jobrunner_no_scratch,
            server=pbs_server,
        )
        assert job.jobrunner is pyscf_jobrunner_no_scratch


class TestPySCFJobContract:
    def test_label_defaults_to_the_empirical_formula(
        self, pyscf_water_molecule, temporary_working_dir
    ):
        job = PySCFSinglePointJob(
            molecule=pyscf_water_molecule, settings=make_pyscf_settings("sp")
        )
        assert job.label == "H2O"
        assert job.errfile.endswith("H2O.err")

    def test_settings_class(self):
        assert PySCFSinglePointJob.settings_class() is PySCFJobSettings

    def test_base_job_has_no_stages_of_its_own(
        self, pyscf_water_molecule, temporary_working_dir
    ):
        job = PySCFJob(
            molecule=pyscf_water_molecule,
            settings=make_pyscf_settings("sp"),
            label="base",
        )
        with pytest.raises(NotImplementedError, match="PySCFJob"):
            job.stages

    def test_opt_with_freq_appends_the_hessian_stage(
        self, pyscf_water_molecule, temporary_working_dir
    ):
        job = PySCFOptJob(
            molecule=pyscf_water_molecule,
            settings=make_pyscf_settings("opt", freq=True),
            label="opt",
        )
        assert job.stages == ["scf", "opt", "hess"]

    def test_backup_keeps_script_log_and_results_together(
        self,
        pyscf_water_molecule,
        pyscf_jobrunner_no_scratch,
        temporary_working_dir,
    ):
        job = PySCFSinglePointJob(
            molecule=pyscf_water_molecule,
            settings=make_pyscf_settings("sp"),
            label="kept",
            jobrunner=pyscf_jobrunner_no_scratch,
        )
        job.run()
        job.backup()
        backed_up = {
            os.path.basename(path)
            for path in glob.glob(
                os.path.join(job.folder, "**", "kept.*"), recursive=True
            )
            if os.path.dirname(path) != job.folder
        }
        assert backed_up == {"kept.py", "kept.out", "kept.h5"}


class TestPySCFJobRunner:
    def test_real_job_runs_end_to_end_in_a_subprocess(
        self, pbs_server, pyscf_h2_molecule, temporary_working_dir
    ):
        pytest.importorskip("pyscf")
        job = PySCFSinglePointJob(
            molecule=pyscf_h2_molecule,
            settings=PySCFJobSettings(
                ab_initio="hf",
                basis="sto-3g",
                density_fit=False,
                charge=0,
                multiplicity=1,
                jobtype="sp",
            ),
            label="h2",
            jobrunner=PySCFJobRunner(server=pbs_server, scratch=False),
        )
        job.run()
        assert job.is_complete()
        output = PySCFOutput(job.outputfile)
        assert output.energies == pytest.approx([-1.116759], abs=1e-5)
        assert output.provenance["pyscf_version"] != "fake"
        # A successful run leaves no error log behind.
        assert not os.path.exists(job.errfile)

    def test_scratch_defaults_to_off(self, pbs_server):
        assert PySCFJobRunner(server=pbs_server).scratch is False

    def test_command_uses_the_configured_interpreter(
        self, pbs_server, pyscf_water_molecule, temporary_working_dir
    ):
        runner = PySCFJobRunner(server=pbs_server, scratch=False)
        job = PySCFSinglePointJob(
            molecule=pyscf_water_molecule,
            settings=make_pyscf_settings("sp"),
            label="command",
            jobrunner=runner,
        )
        runner._prerun(job)
        assert runner._get_command(job) == (
            f"{runner.executable.get_executable()} {runner.job_inputfile}"
        )

    @pytest.mark.parametrize(
        ("engine", "expected"), [("cpu", ""), ("gpu", "0,1")]
    )
    def test_device_visibility_follows_the_engine(
        self,
        engine,
        expected,
        pbs_server,
        pyscf_water_molecule,
        temporary_working_dir,
        monkeypatch,
    ):
        monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
        runner = PySCFJobRunner(server=pbs_server, scratch=False, num_gpus=2)
        job = PySCFSinglePointJob(
            molecule=pyscf_water_molecule,
            settings=make_pyscf_settings("sp", engine=engine),
            label="devices",
            jobrunner=runner,
        )
        env = runner._update_os_environ(job)
        assert env["CUDA_VISIBLE_DEVICES"] == expected

    def test_missing_server_file_is_reported(self, pbs_server, caplog):
        runner = PySCFJobRunner(server=pbs_server, scratch=False)
        with patch.object(
            PySCFExecutable,
            "from_servername",
            side_effect=FileNotFoundError("no such server"),
        ):
            with pytest.raises(FileNotFoundError, match="no such server"):
                runner.executable

    def test_failed_copy_from_scratch_is_logged_not_raised(
        self,
        pyscf_water_molecule,
        pyscf_jobrunner_scratch,
        temporary_working_dir,
        caplog,
    ):
        job = PySCFSinglePointJob(
            molecule=pyscf_water_molecule,
            settings=make_pyscf_settings("sp"),
            label="copyfail",
            jobrunner=pyscf_jobrunner_scratch,
        )
        with patch(
            "chemsmart.jobs.pyscf.runner.copy",
            side_effect=PermissionError("read-only"),
        ):
            job.run()
        assert "read-only" in caplog.text


def test_charge_and_multiplicity_can_come_from_the_molecule(
    pyscf_water_molecule, temporary_working_dir
):
    settings = PySCFJobSettings(
        functional="b3lyp", basis="def2-svp", jobtype="sp", defgrid="defgrid3"
    )
    job = PySCFSinglePointJob(
        molecule=pyscf_water_molecule, settings=settings, label="fallback"
    )
    config = PySCFScriptWriter(job=job).build_config()
    assert (config["charge"], config["multiplicity"]) == (0, 1)
    assert config["atom_grid"] == [99, 590]
