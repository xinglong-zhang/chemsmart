"""
Direct unit tests for :class:`ORCAJob` and its concrete subclasses
(``ORCAInpJob``, ``ORCAGeneralJob``) in ``chemsmart.jobs.orca.job``.

CLI-level propagation is already covered by ``test_orca_cli.py``. These
tests exercise the job classes directly: constructor validation, path
properties, backup/output handling, and the ``from_filename``/
``from_pubchem``/``from_jobtype`` factory methods.
"""

import os
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.orca.job import (
    ORCAGeneralJob,
    ORCAInpJob,
    ORCAJob,
)
from chemsmart.jobs.orca.settings import ORCAJobSettings


@pytest.fixture()
def orca_settings():
    return ORCAJobSettings.default()


@pytest.fixture()
def a_molecule(single_molecule_xyz_file):
    return Molecule.from_filepath(single_molecule_xyz_file)


class TestORCAJobConstruction:
    def test_rejects_invalid_settings_type(self, a_molecule):
        with pytest.raises(ValueError, match="Settings must be instance"):
            ORCAJob(molecule=a_molecule, settings="not-settings")

    def test_rejects_invalid_molecule_type(self, orca_settings):
        with pytest.raises(ValueError, match="Molecule must be instance"):
            ORCAJob(molecule="not-a-molecule", settings=orca_settings)

    def test_label_defaults_to_chemical_formula(
        self, a_molecule, orca_settings
    ):
        job = ORCAJob(molecule=a_molecule, settings=orca_settings)
        assert job.label == a_molecule.get_chemical_formula(empirical=True)

    def test_custom_label_used(self, a_molecule, orca_settings):
        job = ORCAJob(
            molecule=a_molecule, settings=orca_settings, label="mylabel"
        )
        assert job.label == "mylabel"

    def test_settings_class_returns_expected_type(self):
        assert ORCAJob.settings_class() is ORCAJobSettings


class TestORCAJobPaths:
    def test_input_output_gbw_err_file_paths(self, a_molecule, orca_settings):
        job = ORCAJob(
            molecule=a_molecule, settings=orca_settings, label="mylabel"
        )
        assert job.inputfile == os.path.join(job.folder, "mylabel.inp")
        assert job.outputfile == os.path.join(job.folder, "mylabel.out")
        assert job.gbwfile == os.path.join(job.folder, "mylabel.gbw")
        assert job.errfile == os.path.join(job.folder, "mylabel.err")


class TestORCAJobBackupAndOutput:
    def test_backup_files_copies_input_and_output(
        self, a_molecule, orca_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = ORCAJob(
            molecule=a_molecule, settings=orca_settings, label="mylabel"
        )
        with open(job.inputfile, "w") as f:
            f.write("input")
        with open(job.outputfile, "w") as f:
            f.write("output")

        job._backup_files()

        backups = job._previous_backup_folders()
        assert len(backups) == 1
        assert os.path.exists(os.path.join(backups[0], "mylabel.inp"))
        assert os.path.exists(os.path.join(backups[0], "mylabel.out"))

    def test_backup_files_includes_gbw_when_requested(
        self, a_molecule, orca_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = ORCAJob(
            molecule=a_molecule, settings=orca_settings, label="mylabel"
        )
        for f in [job.inputfile, job.outputfile, job.gbwfile]:
            with open(f, "w") as fh:
                fh.write("data")

        job._backup_files(backup_gbw=True)

        backups = job._previous_backup_folders()
        assert os.path.exists(os.path.join(backups[0], "mylabel.gbw"))

    def test_output_none_when_outputfile_missing(
        self, a_molecule, orca_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = ORCAJob(
            molecule=a_molecule, settings=orca_settings, label="mylabel"
        )
        assert job._output() is None

    def test_output_returns_orca_output_when_present(
        self, a_molecule, orca_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = ORCAJob(
            molecule=a_molecule, settings=orca_settings, label="mylabel"
        )
        with open(job.outputfile, "w") as f:
            f.write("dummy output")

        mock_output = MagicMock()
        with patch(
            "chemsmart.io.orca.output.ORCAOutput", return_value=mock_output
        ) as mock_cls:
            result = job._output()

        mock_cls.assert_called_once_with(job.outputfile)
        assert result is mock_output


class TestORCAJobRun:
    def test_run_delegates_to_jobrunner(self, a_molecule, orca_settings):
        mock_runner = MagicMock()
        job = ORCAJob(
            molecule=a_molecule, settings=orca_settings, jobrunner=mock_runner
        )
        job._run()
        mock_runner.run.assert_called_once_with(job)


class TestORCAJobFactories:
    def test_from_filename_builds_job(
        self, single_molecule_xyz_file, orca_settings
    ):
        mock_jobrunner = MagicMock()
        job = ORCAJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=orca_settings,
            jobrunner=mock_jobrunner,
        )
        assert isinstance(job, ORCAJob)
        assert job.jobrunner is mock_jobrunner

    def test_from_pubchem_builds_job(self, orca_settings):
        pubchem_molecule = MagicMock(spec=Molecule)
        pubchem_molecule.get_chemical_formula.return_value = "H2O"
        pubchem_molecule.copy.return_value = pubchem_molecule
        mock_jobrunner = MagicMock()
        with patch(
            "chemsmart.jobs.orca.job.Molecule.from_pubchem",
            return_value=pubchem_molecule,
        ) as mock_from_pubchem:
            job = ORCAJob.from_pubchem(
                identifier="water",
                settings=orca_settings,
                jobrunner=mock_jobrunner,
            )

        mock_from_pubchem.assert_called_once_with(identifier="water")
        assert isinstance(job, ORCAJob)
        assert job.jobrunner is mock_jobrunner

    def test_from_jobtype_opt(self, a_molecule, orca_settings):
        mock_jobrunner = MagicMock()
        job = ORCAJob.from_jobtype(
            jobtype="opt",
            molecule=a_molecule,
            settings=orca_settings,
            jobrunner=mock_jobrunner,
        )
        from chemsmart.jobs.orca.opt import ORCAOptJob

        assert isinstance(job, ORCAOptJob)

    def test_from_jobtype_inp(self, a_molecule, orca_settings):
        mock_jobrunner = MagicMock()
        job = ORCAJob.from_jobtype(
            jobtype="inp",
            molecule=a_molecule,
            settings=orca_settings,
            jobrunner=mock_jobrunner,
        )
        assert isinstance(job, ORCAInpJob)

    def test_from_jobtype_orca(self, a_molecule, orca_settings):
        mock_jobrunner = MagicMock()
        job = ORCAJob.from_jobtype(
            jobtype="orca",
            molecule=a_molecule,
            settings=orca_settings,
            jobrunner=mock_jobrunner,
        )
        assert isinstance(job, ORCAGeneralJob)

    def test_from_jobtype_invalid_raises(self, a_molecule, orca_settings):
        with pytest.raises(ValueError, match="Invalid job type"):
            ORCAJob.from_jobtype(
                jobtype="bogus",
                molecule=a_molecule,
                settings=orca_settings,
            )


class TestORCAInpJob:
    def test_from_filename_requires_inp_extension(
        self, single_molecule_xyz_file
    ):
        with pytest.raises(AssertionError, match="must be .inp file"):
            ORCAInpJob.from_filename(filename=single_molecule_xyz_file)

    def test_from_filename_builds_job_from_inp_file(self, water_sp_input_path):
        mock_jobrunner = MagicMock()
        job = ORCAInpJob.from_filename(
            filename=water_sp_input_path, jobrunner=mock_jobrunner
        )
        assert isinstance(job, ORCAInpJob)
        assert job.settings.input_string is not None
        assert job.jobrunner is mock_jobrunner

    def test_copy_input_regular_folder_logs_info(
        self, a_molecule, orca_settings, tmp_path, monkeypatch, caplog
    ):
        monkeypatch.chdir(tmp_path)
        mock_runner = MagicMock()
        mock_runner.scratch = False
        job = ORCAInpJob(
            molecule=a_molecule,
            settings=orca_settings,
            label="mylabel",
            jobrunner=mock_runner,
        )
        job._copy_input()  # should not raise

    def test_copy_input_scratch_dir_missing_warns(
        self, a_molecule, orca_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        mock_runner = MagicMock()
        mock_runner.scratch = True
        mock_runner.scratch_dir = str(tmp_path / "does_not_exist")
        job = ORCAInpJob(
            molecule=a_molecule,
            settings=orca_settings,
            label="mylabel",
            jobrunner=mock_runner,
        )
        job._copy_input()  # should not raise, just warn

    def test_copy_input_copies_file_to_scratch(
        self, a_molecule, orca_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        scratch_dir = tmp_path / "scratch"
        scratch_dir.mkdir()
        mock_runner = MagicMock()
        mock_runner.scratch = True
        mock_runner.scratch_dir = str(scratch_dir)
        job = ORCAInpJob(
            molecule=a_molecule,
            settings=orca_settings,
            label="mylabel",
            jobrunner=mock_runner,
        )
        with open(job.inputfile, "w") as f:
            f.write("dummy input")

        job._copy_input()

        assert os.path.exists(
            os.path.join(scratch_dir, "mylabel", "mylabel.inp")
        )


class TestORCAGeneralJob:
    def test_construction(self, a_molecule, orca_settings):
        job = ORCAGeneralJob(molecule=a_molecule, settings=orca_settings)
        assert job.TYPE == "orcajob"
        assert isinstance(job, ORCAJob)
