"""
Direct unit tests for :class:`NCIPLOTJob` in
``chemsmart.jobs.nciplot.job``.

Covers the filenames/molecule XOR validation, label derivation, file
path properties, backup, output detection, and completion checking.
Also documents a crash bug when ``settings=None`` (the documented
default) is used (see ``BUGS_FOUND.md``).
"""

import os
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.jobs.nciplot.job import NCIPLOTJob
from chemsmart.jobs.nciplot.settings import NCIPLOTJobSettings


@pytest.fixture()
def nciplot_settings():
    return NCIPLOTJobSettings()


class TestNCIPLOTJobSettingsDefaultCrashes:
    def test_settings_none_always_crashes(self):
        """Documents a bug (see BUGS_FOUND.md): ``settings`` defaults to
        ``None`` and is only *conditionally* type-checked
        (``if settings is not None and not isinstance(...)``), but
        ``self.settings = settings.copy()`` runs unconditionally, so the
        documented default value crashes with an unrelated
        ``AttributeError`` instead of the class's own ``ValueError``
        validation path."""
        with pytest.raises(AttributeError, match="'NoneType' object has no"):
            NCIPLOTJob(filenames=["a.xyz"], settings=None)


class TestNCIPLOTJobValidation:
    def test_wrong_settings_type_raises(self):
        with pytest.raises(ValueError, match="Settings must be instance"):
            NCIPLOTJob(filenames=["a.xyz"], settings="not-settings")

    def test_both_filenames_and_molecule_raises(
        self, nciplot_settings, ethanol_molecule
    ):
        with pytest.raises(ValueError, match="not both"):
            NCIPLOTJob(
                filenames=["a.xyz"],
                molecule=ethanol_molecule,
                settings=nciplot_settings,
            )

    def test_neither_filenames_nor_molecule_raises(self, nciplot_settings):
        with pytest.raises(ValueError, match="both are None"):
            NCIPLOTJob(settings=nciplot_settings)

    def test_empty_filenames_list_raises(self, nciplot_settings):
        with pytest.raises(ValueError, match="No filenames provided"):
            NCIPLOTJob(filenames=[], settings=nciplot_settings)

    def test_molecule_wrong_type_raises(self, nciplot_settings):
        with pytest.raises(ValueError, match="Molecule must be instance"):
            NCIPLOTJob(molecule="not-a-molecule", settings=nciplot_settings)


class TestNCIPLOTJobLabelDerivation:
    def test_single_filename_label_is_stem(self, nciplot_settings):
        job = NCIPLOTJob(filenames=["complex.xyz"], settings=nciplot_settings)
        assert job.label == "complex"

    def test_single_filename_explicit_label_wins(self, nciplot_settings):
        job = NCIPLOTJob(
            filenames=["complex.xyz"],
            settings=nciplot_settings,
            label="custom",
        )
        assert job.label == "custom"

    def test_multiple_filenames_label_joins_stems(self, nciplot_settings):
        job = NCIPLOTJob(
            filenames=["frag1.xyz", "frag2.xyz"], settings=nciplot_settings
        )
        assert job.label == "frag1_and_frag2"

    def test_molecule_label_is_none_by_default(
        self, nciplot_settings, ethanol_molecule
    ):
        job = NCIPLOTJob(molecule=ethanol_molecule, settings=nciplot_settings)
        assert job.label is None

    def test_molecule_explicit_label_used(
        self, nciplot_settings, ethanol_molecule
    ):
        job = NCIPLOTJob(
            molecule=ethanol_molecule,
            settings=nciplot_settings,
            label="mymol",
        )
        assert job.label == "mymol"


class TestNCIPLOTJobFilePaths:
    def test_inputfile_outputfile_errfile(self, nciplot_settings):
        job = NCIPLOTJob(
            filenames=["complex.xyz"],
            settings=nciplot_settings,
            label="mynci",
        )
        assert job.inputfile == f"{job.folder}/mynci.nci"
        assert job.outputfile == f"{job.folder}/mynci.nciout"
        assert job.errfile == f"{job.folder}/mynci.ncierr"


class TestNCIPLOTJobBackup:
    def test_backup_files_calls_backup_file_with_created_folder(
        self, nciplot_settings
    ):
        job = NCIPLOTJob(
            filenames=["complex.xyz"],
            settings=nciplot_settings,
            label="mynci",
        )
        with (
            patch.object(
                job, "_create_backup_folder_name", return_value="/bk/folder"
            ),
            patch.object(job, "backup_file") as mock_backup_file,
        ):
            job._backup_files()

        mock_backup_file.assert_called_once_with(
            job.outputfile, folder="/bk/folder"
        )


class TestNCIPLOTJobOutput:
    def test_output_returns_none_when_missing(self, nciplot_settings):
        job = NCIPLOTJob(
            filenames=["complex.xyz"],
            settings=nciplot_settings,
            label="mynci_missing",
        )
        assert job._output() is None

    def test_output_returns_abspath_when_present(
        self, nciplot_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = NCIPLOTJob(
            filenames=["complex.xyz"],
            settings=nciplot_settings,
            label="mynci_present",
        )
        with open(job.outputfile, "w") as f:
            f.write("some output\n")
        assert job._output() == os.path.abspath(job.outputfile)


class TestNCIPLOTJobIsComplete:
    def test_incomplete_when_file_missing(
        self, nciplot_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = NCIPLOTJob(
            filenames=["complex.xyz"],
            settings=nciplot_settings,
            label="mynci",
        )
        assert job._job_is_complete() is False

    def test_incomplete_when_file_present_but_no_end_marker(
        self, nciplot_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = NCIPLOTJob(
            filenames=["complex.xyz"],
            settings=nciplot_settings,
            label="mynci",
        )
        with open(job.outputfile, "w") as f:
            f.write("Some line\nAnother line\n")
        assert job._job_is_complete() is False

    def test_complete_when_last_line_starts_with_end(
        self, nciplot_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = NCIPLOTJob(
            filenames=["complex.xyz"],
            settings=nciplot_settings,
            label="mynci",
        )
        with open(job.outputfile, "w") as f:
            f.write("Some line\nEnd of NCIPLOT run\n")
        assert job._job_is_complete() is True

    def test_incomplete_when_file_is_empty(
        self, nciplot_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = NCIPLOTJob(
            filenames=["complex.xyz"],
            settings=nciplot_settings,
            label="mynci",
        )
        open(job.outputfile, "w").close()
        assert job._job_is_complete() is False


class TestNCIPLOTJobRun:
    def test_run_delegates_to_jobrunner(self, nciplot_settings):
        mock_runner = MagicMock()
        job = NCIPLOTJob(
            filenames=["complex.xyz"],
            settings=nciplot_settings,
            label="mynci",
            jobrunner=mock_runner,
        )
        job._run(extra="kwarg")
        mock_runner.run.assert_called_once_with(job, extra="kwarg")

    def test_settings_class_returns_nciplot_settings(self):
        assert NCIPLOTJob.settings_class() is NCIPLOTJobSettings
