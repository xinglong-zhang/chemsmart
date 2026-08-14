"""
Direct unit tests for :class:`GaussianJob` and its concrete subclasses
(``GaussianComJob``, ``GaussianGeneralJob``) in
``chemsmart.jobs.gaussian.job``.

CLI-level propagation is already covered by ``test_gaussian_cli.py``.
These tests exercise the job classes directly: constructor validation,
path properties, backup/output handling, and the ``from_filename``/
``from_pubchem``/``from_jobtype`` factory methods.
"""

import os
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.job import (
    GaussianComJob,
    GaussianGeneralJob,
    GaussianJob,
)
from chemsmart.jobs.gaussian.settings import GaussianJobSettings


@pytest.fixture()
def gaussian_settings():
    return GaussianJobSettings.default()


@pytest.fixture()
def a_molecule(single_molecule_xyz_file):
    return Molecule.from_filepath(single_molecule_xyz_file)


class TestGaussianJobConstruction:
    def test_rejects_invalid_settings_type(self, a_molecule):
        with pytest.raises(ValueError, match="Settings must be instance"):
            GaussianJob(molecule=a_molecule, settings="not-settings")

    def test_rejects_invalid_molecule_type(self, gaussian_settings):
        with pytest.raises(ValueError, match="Molecule must be instance"):
            GaussianJob(molecule="not-a-molecule", settings=gaussian_settings)

    def test_label_defaults_to_chemical_formula(
        self, a_molecule, gaussian_settings
    ):
        job = GaussianJob(molecule=a_molecule, settings=gaussian_settings)
        assert job.label == a_molecule.get_chemical_formula(empirical=True)

    def test_custom_label_used(self, a_molecule, gaussian_settings):
        job = GaussianJob(
            molecule=a_molecule, settings=gaussian_settings, label="mylabel"
        )
        assert job.label == "mylabel"

    def test_settings_class_returns_expected_type(self):
        assert GaussianJob.settings_class() is GaussianJobSettings


class TestGaussianJobPaths:
    def test_input_output_chk_err_file_paths(
        self, a_molecule, gaussian_settings
    ):
        job = GaussianJob(
            molecule=a_molecule, settings=gaussian_settings, label="mylabel"
        )
        assert job.inputfile == os.path.join(job.folder, "mylabel.com")
        assert job.outputfile == os.path.join(job.folder, "mylabel.log")
        assert job.chkfile == os.path.join(job.folder, "mylabel.chk")
        assert job.errfile == os.path.join(job.folder, "mylabel.err")


class TestGaussianJobBackupAndOutput:
    def test_backup_files_copies_input_and_output(
        self, a_molecule, gaussian_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = GaussianJob(
            molecule=a_molecule, settings=gaussian_settings, label="mylabel"
        )
        with open(job.inputfile, "w") as f:
            f.write("input")
        with open(job.outputfile, "w") as f:
            f.write("output")

        job._backup_files()

        backups = job._previous_backup_folders()
        assert len(backups) == 1
        assert os.path.exists(os.path.join(backups[0], "mylabel.com"))
        assert os.path.exists(os.path.join(backups[0], "mylabel.log"))

    def test_backup_files_includes_chk_when_requested(
        self, a_molecule, gaussian_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = GaussianJob(
            molecule=a_molecule, settings=gaussian_settings, label="mylabel"
        )
        for f in [job.inputfile, job.outputfile, job.chkfile]:
            with open(f, "w") as fh:
                fh.write("data")

        job._backup_files(backup_chk=True)

        backups = job._previous_backup_folders()
        assert os.path.exists(os.path.join(backups[0], "mylabel.chk"))

    def test_output_none_when_outputfile_missing(
        self, a_molecule, gaussian_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = GaussianJob(
            molecule=a_molecule, settings=gaussian_settings, label="mylabel"
        )
        assert job._output() is None

    def test_output_returns_gaussian_output_when_present(
        self, a_molecule, gaussian_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = GaussianJob(
            molecule=a_molecule, settings=gaussian_settings, label="mylabel"
        )
        with open(job.outputfile, "w") as f:
            f.write("dummy output")

        mock_output = MagicMock()
        with patch(
            "chemsmart.io.gaussian.output.Gaussian16Output",
            return_value=mock_output,
        ) as mock_cls:
            result = job._output()

        mock_cls.assert_called_once_with(filename=job.outputfile)
        assert result is mock_output

    def test_output_falls_back_to_pbc_output_on_value_error(
        self, a_molecule, gaussian_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = GaussianJob(
            molecule=a_molecule, settings=gaussian_settings, label="mylabel"
        )
        with open(job.outputfile, "w") as f:
            f.write("dummy output")

        mock_pbc_output = MagicMock()
        with (
            patch(
                "chemsmart.io.gaussian.output.Gaussian16Output",
                side_effect=ValueError("not a standard output"),
            ),
            patch(
                "chemsmart.io.gaussian.output.Gaussian16OutputWithPBC",
                return_value=mock_pbc_output,
            ) as mock_pbc_cls,
        ):
            result = job._output()

        mock_pbc_cls.assert_called_once_with(filename=job.outputfile)
        assert result is mock_pbc_output


class TestGaussianJobRun:
    def test_run_delegates_to_jobrunner(self, a_molecule, gaussian_settings):
        mock_runner = MagicMock()
        job = GaussianJob(
            molecule=a_molecule,
            settings=gaussian_settings,
            jobrunner=mock_runner,
        )
        job._run()
        mock_runner.run.assert_called_once_with(job)


class TestGaussianJobFactories:
    def test_from_filename_builds_job(
        self, single_molecule_xyz_file, gaussian_settings
    ):
        mock_jobrunner = MagicMock()
        job = GaussianJob.from_filename(
            filename=single_molecule_xyz_file,
            settings=gaussian_settings,
            jobrunner=mock_jobrunner,
        )
        assert isinstance(job, GaussianJob)
        assert job.jobrunner is mock_jobrunner

    def test_from_pubchem_builds_job(self, gaussian_settings):
        pubchem_molecule = MagicMock(spec=Molecule)
        pubchem_molecule.get_chemical_formula.return_value = "H2O"
        pubchem_molecule.copy.return_value = pubchem_molecule
        mock_jobrunner = MagicMock()
        with patch(
            "chemsmart.jobs.gaussian.job.Molecule.from_pubchem",
            return_value=pubchem_molecule,
        ) as mock_from_pubchem:
            job = GaussianJob.from_pubchem(
                identifier="water",
                settings=gaussian_settings,
                jobrunner=mock_jobrunner,
            )

        mock_from_pubchem.assert_called_once_with(identifier="water")
        assert isinstance(job, GaussianJob)
        assert job.jobrunner is mock_jobrunner

    def test_from_jobtype_opt(self, a_molecule, gaussian_settings):
        mock_jobrunner = MagicMock()
        job = GaussianJob.from_jobtype(
            jobtype="opt",
            molecule=a_molecule,
            settings=gaussian_settings,
            jobrunner=mock_jobrunner,
        )
        from chemsmart.jobs.gaussian.opt import GaussianOptJob

        assert isinstance(job, GaussianOptJob)

    def test_from_jobtype_com(self, a_molecule, gaussian_settings):
        mock_jobrunner = MagicMock()
        job = GaussianJob.from_jobtype(
            jobtype="com",
            molecule=a_molecule,
            settings=gaussian_settings,
            jobrunner=mock_jobrunner,
        )
        assert isinstance(job, GaussianComJob)

    def test_from_jobtype_g16_without_explicit_jobrunner_works(
        self, a_molecule, gaussian_settings
    ):
        """When ``jobrunner`` is *not* supplied, the "g16" branch creates
        one and correctly returns a ``GaussianGeneralJob``."""
        job = GaussianJob.from_jobtype(
            jobtype="g16",
            molecule=a_molecule,
            settings=gaussian_settings,
        )
        assert isinstance(job, GaussianGeneralJob)

    def test_from_jobtype_g16_with_explicit_jobrunner_is_buggy(
        self, a_molecule, gaussian_settings
    ):
        """Documents a real bug in ``GaussianJob.from_jobtype``: the
        ``return GaussianGeneralJob(...)`` for the "g16" branch is nested
        inside ``if jobrunner is None:`` instead of the ``elif
        jobtype.lower() == "g16":`` block. When a caller *already* has a
        jobrunner (the common case — see e.g. the ``qmmm`` subcommands
        that build parent jobs with a jobrunner in hand), this falls
        through to the ``else: raise ValueError(f"Invalid job type:
        {jobtype}")`` clause even though "g16" is a valid job type."""
        mock_jobrunner = MagicMock()
        with pytest.raises(ValueError, match="Invalid job type: g16"):
            GaussianJob.from_jobtype(
                jobtype="g16",
                molecule=a_molecule,
                settings=gaussian_settings,
                jobrunner=mock_jobrunner,
            )

    def test_from_jobtype_invalid_silently_returns_none(
        self, a_molecule, gaussian_settings
    ):
        """Documents a second bug in ``from_jobtype``: the "invalid job
        type" ``raise ValueError`` is only reachable from inside the
        "g16" branch (see the "buggy" test above). There is no top-level
        ``else`` for jobtypes that are not "opt"/"com"/"g16" at all, so a
        genuinely unrecognized jobtype silently falls through and returns
        ``None`` instead of raising."""
        result = GaussianJob.from_jobtype(
            jobtype="bogus",
            molecule=a_molecule,
            settings=gaussian_settings,
        )
        assert result is None


class TestGaussianComJob:
    def test_from_filename_crashes_on_none_molecule(
        self, gaussian_opt_inputfile
    ):
        """Documents a third bug: ``GaussianComJob.from_filename`` always
        constructs the job with ``molecule=None`` (it only reads route
        info, not coordinates, from the ``.com`` file), but the parent
        ``GaussianJob.__init__`` unconditionally requires
        ``isinstance(molecule, Molecule)``. So this factory method
        currently cannot succeed at all."""
        mock_jobrunner = MagicMock()
        with pytest.raises(ValueError, match="Molecule must be instance"):
            GaussianComJob.from_filename(
                filename=gaussian_opt_inputfile, jobrunner=mock_jobrunner
            )


class TestGaussianGeneralJob:
    def test_construction(self, a_molecule, gaussian_settings):
        job = GaussianGeneralJob(
            molecule=a_molecule, settings=gaussian_settings
        )
        assert job.TYPE == "g16"
        assert isinstance(job, GaussianJob)
