"""
Direct unit tests for :class:`XTBJob` in ``chemsmart.jobs.xtb.job``.

CLI-level propagation is already covered by ``test_xtb_cli.py``, and full
run-integration tests live in ``test_XTBJobs.py``. These tests exercise
the base job class directly: constructor validation, path properties,
backup/output handling, and the ``from_filename`` factory.
"""

import os
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.xtb.job import XTBJob
from chemsmart.jobs.xtb.settings import XTBJobSettings


@pytest.fixture()
def xtb_settings():
    return XTBJobSettings.default()


@pytest.fixture()
def a_molecule(single_molecule_xyz_file):
    return Molecule.from_filepath(single_molecule_xyz_file)


class TestXTBJobConstruction:
    def test_rejects_invalid_settings_type(self, a_molecule):
        with pytest.raises(ValueError, match="Settings must be instance"):
            XTBJob(molecule=a_molecule, settings="not-settings")

    def test_rejects_invalid_molecule_type(self, xtb_settings):
        with pytest.raises(ValueError, match="Molecule must be instance"):
            XTBJob(molecule="not-a-molecule", settings=xtb_settings)

    def test_label_defaults_to_chemical_formula(
        self, a_molecule, xtb_settings
    ):
        job = XTBJob(molecule=a_molecule, settings=xtb_settings)
        assert job.label == a_molecule.get_chemical_formula(empirical=True)

    def test_custom_label_used(self, a_molecule, xtb_settings):
        job = XTBJob(
            molecule=a_molecule, settings=xtb_settings, label="mylabel"
        )
        assert job.label == "mylabel"

    def test_settings_class_returns_expected_type(self):
        assert XTBJob.settings_class() is XTBJobSettings


class TestXTBJobPaths:
    def test_all_path_properties(self, a_molecule, xtb_settings):
        job = XTBJob(
            molecule=a_molecule, settings=xtb_settings, label="mylabel"
        )
        assert job.xyzfile == os.path.join(job.folder, "mylabel.xyz")
        assert job.inputfile == os.path.join(job.folder, "mylabel.inp")
        assert job.outputfile == os.path.join(job.folder, "mylabel.out")
        assert job.errfile == os.path.join(job.folder, "mylabel.err")
        assert job.restartfile == os.path.join(job.folder, "xtbrestart")

    def test_base_folder_returns_determine_folder(
        self, a_molecule, xtb_settings
    ):
        job = XTBJob(
            molecule=a_molecule, settings=xtb_settings, label="mylabel"
        )
        assert job.base_folder() == job._determine_folder()

    def test_determine_folder_skips_label_subfolder_when_label_falsy(
        self, a_molecule, xtb_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = XTBJob(
            molecule=a_molecule, settings=xtb_settings, label="mylabel"
        )
        job.label = ""
        assert job._determine_folder() == os.path.abspath(str(tmp_path))


class TestXTBJobBackupAndOutput:
    def test_backup_files_copies_all_tracked_files(
        self, a_molecule, xtb_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = XTBJob(
            molecule=a_molecule, settings=xtb_settings, label="mylabel"
        )
        for f in [
            job.xyzfile,
            job.inputfile,
            job.outputfile,
            job.restartfile,
        ]:
            with open(f, "w") as fh:
                fh.write("data")

        job._backup_files()

        backups = job._previous_backup_folders()
        assert len(backups) == 1
        assert os.path.exists(os.path.join(backups[0], "mylabel.xyz"))
        assert os.path.exists(os.path.join(backups[0], "mylabel.inp"))
        assert os.path.exists(os.path.join(backups[0], "mylabel.out"))
        assert os.path.exists(os.path.join(backups[0], "xtbrestart"))

    def test_output_none_when_outputfile_missing(
        self, a_molecule, xtb_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = XTBJob(
            molecule=a_molecule, settings=xtb_settings, label="mylabel"
        )
        assert job._output() is None

    def test_output_returns_xtb_output_when_present(
        self, a_molecule, xtb_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = XTBJob(
            molecule=a_molecule, settings=xtb_settings, label="mylabel"
        )
        with open(job.outputfile, "w") as f:
            f.write("dummy output")

        mock_output = MagicMock()
        with patch(
            "chemsmart.io.xtb.output.XTBOutput",
            return_value=mock_output,
        ) as mock_cls:
            result = job._output()

        mock_cls.assert_called_once_with(folder=job.folder)
        assert result is mock_output

    def test_output_none_on_exception(
        self, a_molecule, xtb_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = XTBJob(
            molecule=a_molecule, settings=xtb_settings, label="mylabel"
        )
        with open(job.outputfile, "w") as f:
            f.write("dummy output")

        with patch(
            "chemsmart.io.xtb.output.XTBOutput",
            side_effect=RuntimeError("boom"),
        ):
            assert job._output() is None


class TestXTBJobRun:
    def test_run_delegates_to_jobrunner(self, a_molecule, xtb_settings):
        mock_runner = MagicMock()
        job = XTBJob(
            molecule=a_molecule,
            settings=xtb_settings,
            jobrunner=mock_runner,
        )
        job._run()
        mock_runner.run.assert_called_once_with(job)


class TestXTBJobFactories:
    def test_from_filename_builds_job(
        self, single_molecule_xyz_file, xtb_settings
    ):
        mock_jobrunner = MagicMock()
        job = XTBJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=xtb_settings,
            jobrunner=mock_jobrunner,
        )
        assert isinstance(job, XTBJob)
        assert job.jobrunner is mock_jobrunner

    def test_from_filename_without_jobrunner_leaves_it_none(
        self, single_molecule_xyz_file, xtb_settings
    ):
        job = XTBJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=xtb_settings,
        )
        assert isinstance(job, XTBJob)
        assert job.jobrunner is None
