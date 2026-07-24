"""
Direct unit tests for :class:`NCIPLOTJobRunner`.

Covers scratch vs. job-directory path setup, the input-file preparation
branches (PubChem molecule, supported-format copy, format conversion,
missing/invalid filenames), command generation, and postrun scratch
file copying. ``FakeNCIPLOTJobRunner``/``FakeNCIPLOT`` are already
exercised indirectly by many existing fixtures across the test suite.
"""

import os
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.jobs.nciplot.job import NCIPLOTJob
from chemsmart.jobs.nciplot.runner import NCIPLOTJobRunner
from chemsmart.jobs.nciplot.settings import NCIPLOTJobSettings


@pytest.fixture()
def nciplot_settings():
    return NCIPLOTJobSettings()


@pytest.fixture()
def nciplot_job_from_files(
    single_molecule_xyz_file, nciplot_settings, tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    return NCIPLOTJob(
        filenames=[single_molecule_xyz_file],
        settings=nciplot_settings,
        label="ncitest",
    )


class TestNCIPLOTJobRunnerSetup:
    def test_defaults_scratch_to_true(self, pbs_server):
        # Constructing with scratch left unset (None) would trigger the
        # real ``_set_scratch`` executable/scratch-dir resolution, which
        # the minimal test server fixture doesn't support. The relevant
        # behavior — falling back to the class default — is exercised via
        # the class attribute directly.
        assert NCIPLOTJobRunner.SCRATCH is True
        runner = NCIPLOTJobRunner(server=pbs_server, scratch=False)
        assert runner.scratch is False

    def test_prerun_job_directory_paths(
        self, pbs_server, nciplot_job_from_files
    ):
        runner = NCIPLOTJobRunner(server=pbs_server, scratch=False)
        runner._assign_variables(nciplot_job_from_files)

        assert runner.running_directory == nciplot_job_from_files.folder
        assert runner.job_inputfile == os.path.abspath(
            nciplot_job_from_files.inputfile
        )
        assert runner.job_outputfile == os.path.abspath(
            nciplot_job_from_files.outputfile
        )
        assert runner.job_errfile == os.path.abspath(
            nciplot_job_from_files.errfile
        )

    def test_prerun_scratch_directory_paths(
        self, pbs_server, nciplot_job_from_files, tmp_path
    ):
        scratch_dir = tmp_path / "scratch"
        runner = NCIPLOTJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        runner._assign_variables(nciplot_job_from_files)

        expected_dir = os.path.join(
            str(scratch_dir), nciplot_job_from_files.label
        )
        assert runner.running_directory == expected_dir
        assert os.path.isdir(expected_dir)
        assert runner.job_inputfile.endswith("ncitest.nci")
        assert runner.job_outputfile.endswith("ncitest.nciout")
        assert runner.job_errfile.endswith("ncitest.ncierr")


class TestNCIPLOTJobRunnerPrepareFiles:
    def test_prepare_files_requires_filenames_when_no_molecule(
        self, pbs_server, nciplot_settings, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = NCIPLOTJob.__new__(NCIPLOTJob)
        job.molecule = None
        job.filenames = None
        runner = NCIPLOTJobRunner(server=pbs_server, scratch=False)
        with pytest.raises(AssertionError, match="No molecule provided"):
            runner._prepare_files(job)

    def test_prepare_files_rejects_non_list_filenames(
        self, pbs_server, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = NCIPLOTJob.__new__(NCIPLOTJob)
        job.molecule = None
        job.filenames = "not-a-list"
        runner = NCIPLOTJobRunner(server=pbs_server, scratch=False)
        with pytest.raises(TypeError, match="Expected filenames"):
            runner._prepare_files(job)

    def test_prepare_files_rejects_empty_filenames(
        self, pbs_server, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = NCIPLOTJob.__new__(NCIPLOTJob)
        job.molecule = None
        job.filenames = []
        runner = NCIPLOTJobRunner(server=pbs_server, scratch=False)
        with pytest.raises(ValueError, match="No filenames provided"):
            runner._prepare_files(job)

    def test_prepare_files_copies_supported_formats(
        self, pbs_server, nciplot_job_from_files, single_molecule_xyz_file
    ):
        runner = NCIPLOTJobRunner(server=pbs_server, scratch=False)
        runner._assign_variables(nciplot_job_from_files)
        os.makedirs(runner.running_directory, exist_ok=True)

        with patch.object(runner, "_copy_input_files") as mock_copy:
            runner._prepare_files(nciplot_job_from_files)

        mock_copy.assert_called_once_with(nciplot_job_from_files)

    def test_prepare_files_converts_unsupported_formats(
        self,
        pbs_server,
        gaussian_opt_inputfile,
        nciplot_settings,
        tmp_path,
        monkeypatch,
    ):
        monkeypatch.chdir(tmp_path)
        job = NCIPLOTJob(
            filenames=[gaussian_opt_inputfile],
            settings=nciplot_settings,
            label="conv_test",
        )
        runner = NCIPLOTJobRunner(server=pbs_server, scratch=False)
        runner._assign_variables(job)
        os.makedirs(runner.running_directory, exist_ok=True)

        with patch.object(runner, "_write_xyz_from_input_files") as mock_conv:
            runner._prepare_files(job)

        mock_conv.assert_called_once_with(job)

    def test_prepare_files_writes_xyz_from_pubchem_molecule(
        self,
        pbs_server,
        single_molecule_xyz_file,
        nciplot_settings,
        tmp_path,
        monkeypatch,
    ):
        monkeypatch.chdir(tmp_path)
        from chemsmart.io.molecules.structure import Molecule

        molecule = Molecule.from_filepath(single_molecule_xyz_file)
        job = NCIPLOTJob(
            molecule=molecule, settings=nciplot_settings, label="pubchem_test"
        )
        runner = NCIPLOTJobRunner(server=pbs_server, scratch=False)
        runner._assign_variables(job)
        os.makedirs(runner.running_directory, exist_ok=True)

        with patch.object(runner, "_write_xyz_from_pubchem") as mock_write:
            runner._prepare_files(job)

        mock_write.assert_called_once_with(job)


class TestNCIPLOTJobRunnerCopyInputFiles:
    def test_copy_input_files_raises_for_missing_file(
        self, pbs_server, nciplot_job_from_files, tmp_path
    ):
        runner = NCIPLOTJobRunner(server=pbs_server, scratch=False)
        runner._assign_variables(nciplot_job_from_files)
        os.makedirs(runner.running_directory, exist_ok=True)
        nciplot_job_from_files.filenames = [
            str(tmp_path / "does_not_exist.xyz")
        ]

        with pytest.raises(FileNotFoundError):
            runner._copy_input_files(nciplot_job_from_files)

    def test_copy_input_files_copies_existing_file(
        self, pbs_server, nciplot_job_from_files, single_molecule_xyz_file
    ):
        runner = NCIPLOTJobRunner(server=pbs_server, scratch=False)
        runner._assign_variables(nciplot_job_from_files)
        os.makedirs(runner.running_directory, exist_ok=True)
        nciplot_job_from_files.filenames = [single_molecule_xyz_file]

        runner._copy_input_files(nciplot_job_from_files)

        copied = os.path.join(
            runner.running_directory,
            os.path.basename(single_molecule_xyz_file),
        )
        assert os.path.exists(copied)


class TestNCIPLOTJobRunnerWriteXyzFromPubchem:
    def test_write_xyz_from_pubchem_calls_molecule_write(
        self, pbs_server, nciplot_job_from_files
    ):
        runner = NCIPLOTJobRunner(server=pbs_server, scratch=False)
        runner._assign_variables(nciplot_job_from_files)
        os.makedirs(runner.running_directory, exist_ok=True)

        mock_molecule = MagicMock()
        nciplot_job_from_files.molecule = mock_molecule

        runner._write_xyz_from_pubchem(nciplot_job_from_files)

        expected_path = os.path.join(
            runner.running_directory, f"{nciplot_job_from_files.label}.xyz"
        )
        mock_molecule.write_xyz.assert_called_once_with(
            filename=expected_path, mode="w"
        )


class TestNCIPLOTJobRunnerCommand:
    def test_get_command_uses_executable_and_inputfile(
        self, pbs_server, nciplot_job_from_files
    ):
        runner = NCIPLOTJobRunner(server=pbs_server, scratch=False)
        runner._assign_variables(nciplot_job_from_files)

        with patch.object(
            runner, "_get_executable", return_value="/path/to/nciplot"
        ):
            command = runner._get_command(nciplot_job_from_files)

        assert command == f"/path/to/nciplot {runner.job_inputfile}"


class TestNCIPLOTJobRunnerPostrun:
    def test_postrun_no_op_without_scratch(
        self, pbs_server, nciplot_job_from_files
    ):
        runner = NCIPLOTJobRunner(server=pbs_server, scratch=False)
        runner._assign_variables(nciplot_job_from_files)
        # Should not raise.
        runner._postrun(nciplot_job_from_files)

    def test_postrun_copies_non_tmp_files_from_scratch(
        self, pbs_server, nciplot_job_from_files, tmp_path
    ):
        scratch_dir = tmp_path / "scratch"
        runner = NCIPLOTJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        runner._assign_variables(nciplot_job_from_files)
        os.makedirs(runner.running_directory, exist_ok=True)
        os.makedirs(nciplot_job_from_files.folder, exist_ok=True)

        keep_file = os.path.join(runner.running_directory, "result.nciout")
        with open(keep_file, "w") as f:
            f.write("data")
        tmp_file = os.path.join(runner.running_directory, "scratch.tmp")
        with open(tmp_file, "w") as f:
            f.write("temp")

        runner._postrun(nciplot_job_from_files)

        assert os.path.exists(
            os.path.join(nciplot_job_from_files.folder, "result.nciout")
        )
        assert not os.path.exists(
            os.path.join(nciplot_job_from_files.folder, "scratch.tmp")
        )
