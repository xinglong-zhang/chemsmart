import os
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.jobs.crest.conformers import CRESTConformerSearchJob
from chemsmart.jobs.crest.job import CRESTJob
from chemsmart.jobs.crest.settings import CRESTJobSettings


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


class TestCRESTJobConstructionValidation:
    def test_rejects_non_crest_job_settings(self, water_molecule):
        with pytest.raises(ValueError, match="Settings must be instance"):
            CRESTConformerSearchJob(
                molecule=water_molecule,
                settings="not-a-settings-object",
                label="water_conformers",
            )

    def test_rejects_non_molecule(self):
        settings = CRESTJobSettings(
            charge=0, multiplicity=1, jobtype="conformers"
        )
        with pytest.raises(ValueError, match="Molecule must be instance"):
            CRESTConformerSearchJob(
                molecule="not-a-molecule",
                settings=settings,
                label="water_conformers",
            )

    def test_label_defaults_to_empirical_chemical_formula(
        self, water_molecule
    ):
        settings = CRESTJobSettings(
            charge=0, multiplicity=1, jobtype="conformers"
        )
        job = CRESTConformerSearchJob(
            molecule=water_molecule, settings=settings, label=None
        )
        assert job.label == water_molecule.get_chemical_formula(empirical=True)


class TestCRESTJobs:
    def test_conformers_job_with_constraints(
        self, temporary_working_dir, water_molecule, crest_jobrunner_no_scratch
    ):
        settings = CRESTJobSettings(
            charge=0,
            multiplicity=1,
            jobtype="conformers",
            solvent_model="alpb",
            solvent_id="water",
            constraints=[[1, 2]],
            force_constant=0.5,
            gfn_version="gfn2",
        )
        job = CRESTConformerSearchJob(
            molecule=water_molecule,
            settings=settings,
            label="water_conformers",
            jobrunner=crest_jobrunner_no_scratch,
        )
        assert isinstance(job, CRESTConformerSearchJob)
        assert crest_jobrunner_no_scratch.run(job) == 0

        folder = Path(job.folder)
        assert (folder / f"{job.label}.xyz").is_file()
        assert (folder / "constraints.inp").is_file()
        assert (folder / f"{job.label}.out").is_file()

        inp = (folder / "constraints.inp").read_text()
        assert "$constrain" in inp
        assert "distance: 1, 2" in inp
        assert "force constant=0.5" in inp

        out = (folder / f"{job.label}.out").read_text()
        assert "--cinp constraints.inp" in out
        assert "--gfn2" in out
        assert "--alpb water" in out


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


class TestCRESTConformerSearchJobRun:
    def test_conformers_job_with_constraints(
        self, temporary_working_dir, water_molecule, crest_jobrunner_no_scratch
    ):
        settings = CRESTJobSettings(
            charge=0,
            multiplicity=1,
            jobtype="conformers",
            solvent_model="alpb",
            solvent_id="water",
            constraints=[[1, 2]],
            force_constant=0.5,
            gfn_version="gfn2",
        )
        job = CRESTConformerSearchJob(
            molecule=water_molecule,
            settings=settings,
            label="water_conformers",
            jobrunner=crest_jobrunner_no_scratch,
        )
        assert isinstance(job, CRESTConformerSearchJob)
        assert crest_jobrunner_no_scratch.run(job) == 0

        folder = Path(job.folder)
        assert (folder / f"{job.label}.xyz").is_file()
        assert (folder / "constraints.inp").is_file()
        assert (folder / f"{job.label}.out").is_file()

        inp = (folder / "constraints.inp").read_text()
        assert "$constrain" in inp
        assert "distance: 1, 2" in inp
        assert "force constant=0.5" in inp

        out = (folder / f"{job.label}.out").read_text()
        assert "--cinp constraints.inp" in out
        assert "--gfn2" in out
        assert "--alpb water" in out
