"""
Direct unit tests for :class:`CRESTJob` in ``chemsmart.jobs.crest.job``.

CLI-level propagation is already covered by ``test_crest_cli.py``, and a
full run-integration test lives in ``test_CRESTJobs.py``. These tests
exercise the base job class directly: constructor validation, path
properties, backup/output handling, and the ``from_filename`` factory.
"""

import os
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.crest.job import CRESTJob
from chemsmart.jobs.crest.settings import CRESTJobSettings


@pytest.fixture()
def crest_settings():
    return CRESTJobSettings.default()


@pytest.fixture()
def a_molecule(single_molecule_xyz_file):
    return Molecule.from_filepath(single_molecule_xyz_file)


class TestCRESTJobConstruction:
    def test_rejects_invalid_settings_type(self, a_molecule):
        with pytest.raises(ValueError, match="Settings must be instance"):
            CRESTJob(molecule=a_molecule, settings="not-settings")

    def test_rejects_invalid_molecule_type(self, crest_settings):
        with pytest.raises(ValueError, match="Molecule must be instance"):
            CRESTJob(molecule="not-a-molecule", settings=crest_settings)

    def test_label_defaults_to_chemical_formula(
        self, a_molecule, crest_settings
    ):
        job = CRESTJob(molecule=a_molecule, settings=crest_settings)
        assert job.label == a_molecule.get_chemical_formula(empirical=True)

    def test_custom_label_used(self, a_molecule, crest_settings):
        job = CRESTJob(
            molecule=a_molecule, settings=crest_settings, label="mylabel"
        )
        assert job.label == "mylabel"

    def test_settings_class_returns_expected_type(self):
        assert CRESTJob.settings_class() is CRESTJobSettings


class TestCRESTJobPaths:
    def test_all_path_properties(self, a_molecule, crest_settings):
        job = CRESTJob(
            molecule=a_molecule, settings=crest_settings, label="mylabel"
        )
        assert job.xyzfile == os.path.join(job.folder, "mylabel.xyz")
        assert job.outputfile == os.path.join(job.folder, "mylabel.out")
        assert job.errfile == os.path.join(job.folder, "mylabel.err")
        assert job.constraints_file == os.path.join(
            job.folder, "constraints.inp"
        )
        assert job.toml_file == os.path.join(job.folder, "mylabel.toml")
        assert job.conformers_file == os.path.join(
            job.folder, "crest_conformers.xyz"
        )
        assert job.best_file == os.path.join(job.folder, "crest_best.xyz")
        assert job.rotamers_file == os.path.join(
            job.folder, "crest_rotamers.xyz"
        )

    def test_base_folder_returns_determine_folder(
        self, a_molecule, crest_settings
    ):
        job = CRESTJob(
            molecule=a_molecule, settings=crest_settings, label="mylabel"
        )
        assert job.base_folder() == job._determine_folder()

    def test_determine_folder_skips_label_subfolder_when_label_falsy(
        self, a_molecule, crest_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = CRESTJob(
            molecule=a_molecule, settings=crest_settings, label="mylabel"
        )
        job.label = ""
        assert job._determine_folder() == os.path.abspath(str(tmp_path))


class TestCRESTJobBackupAndOutput:
    def test_backup_files_copies_xyz_and_output(
        self, a_molecule, crest_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = CRESTJob(
            molecule=a_molecule, settings=crest_settings, label="mylabel"
        )
        with open(job.xyzfile, "w") as f:
            f.write("xyz")
        with open(job.outputfile, "w") as f:
            f.write("output")

        job._backup_files()

        backups = job._previous_backup_folders()
        assert len(backups) == 1
        assert os.path.exists(os.path.join(backups[0], "mylabel.xyz"))
        assert os.path.exists(os.path.join(backups[0], "mylabel.out"))

    def test_output_none_when_outputfile_missing(
        self, a_molecule, crest_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = CRESTJob(
            molecule=a_molecule, settings=crest_settings, label="mylabel"
        )
        assert job._output() is None

    def test_output_returns_crest_output_when_present(
        self, a_molecule, crest_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = CRESTJob(
            molecule=a_molecule, settings=crest_settings, label="mylabel"
        )
        with open(job.outputfile, "w") as f:
            f.write("dummy output")

        mock_output = MagicMock()
        with patch(
            "chemsmart.io.crest.output.CRESTOutput",
            return_value=mock_output,
        ) as mock_cls:
            result = job._output()

        mock_cls.assert_called_once_with(folder=job.folder)
        assert result is mock_output

    def test_output_none_on_exception(
        self, a_molecule, crest_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = CRESTJob(
            molecule=a_molecule, settings=crest_settings, label="mylabel"
        )
        with open(job.outputfile, "w") as f:
            f.write("dummy output")

        with patch(
            "chemsmart.io.crest.output.CRESTOutput",
            side_effect=RuntimeError("boom"),
        ):
            assert job._output() is None


class TestCRESTJobRun:
    def test_run_delegates_to_jobrunner(self, a_molecule, crest_settings):
        mock_runner = MagicMock()
        job = CRESTJob(
            molecule=a_molecule,
            settings=crest_settings,
            jobrunner=mock_runner,
        )
        job._run()
        mock_runner.run.assert_called_once_with(job)


class TestCRESTJobFactories:
    def test_from_filename_builds_job(
        self, single_molecule_xyz_file, crest_settings
    ):
        mock_jobrunner = MagicMock()
        job = CRESTJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=crest_settings,
            jobrunner=mock_jobrunner,
        )
        assert isinstance(job, CRESTJob)
        assert job.jobrunner is mock_jobrunner

    def test_from_filename_without_jobrunner_leaves_it_none(
        self, single_molecule_xyz_file, crest_settings
    ):
        """Unlike GaussianJob.from_filename, CRESTJob.from_filename has
        no jobrunner-auto-creation logic -- omitting jobrunner simply
        leaves it None on the constructed job."""
        job = CRESTJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=crest_settings,
        )
        assert isinstance(job, CRESTJob)
        assert job.jobrunner is None
