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
from chemsmart.jobs.nciplot.runner import (
    FakeNCIPLOT,
    FakeNCIPLOTJobRunner,
    NCIPLOTJobRunner,
)
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

    def test_scratch_none_falls_back_to_class_default(self, pbs_server):
        """scratch=None (the constructor default) should resolve to the
        class-level SCRATCH default (True) before being forwarded to the
        base JobRunner constructor. The base __init__ is mocked out since
        a real scratch=True setup would trigger executable resolution
        that this minimal test server config doesn't support."""

        def _fake_init(self, **kwargs):
            self.server = kwargs.get("server")
            self.scratch = kwargs.get("scratch")
            self.fake = kwargs.get("fake")
            self.delete_scratch = None

        with patch(
            "chemsmart.jobs.nciplot.runner.JobRunner.__init__",
            _fake_init,
        ):
            runner = NCIPLOTJobRunner(server=pbs_server, scratch=None)
        assert runner.scratch is True

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

    def test_prerun_scratch_directory_already_exists(
        self, pbs_server, nciplot_job_from_files, tmp_path
    ):
        """The scratch job dir may already exist from a prior run; the
        os.makedirs call should be skipped rather than raising."""
        scratch_dir = tmp_path / "scratch"
        expected_dir = os.path.join(
            str(scratch_dir), nciplot_job_from_files.label
        )
        os.makedirs(expected_dir)

        runner = NCIPLOTJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(scratch_dir)
        )
        runner._assign_variables(nciplot_job_from_files)

        assert runner.running_directory == expected_dir


class TestNCIPLOTJobRunnerExecutable:
    def test_executable_success(self, pbs_server):
        runner = NCIPLOTJobRunner(server=pbs_server, scratch=False)
        mock_exe = MagicMock()
        with patch(
            "chemsmart.jobs.nciplot.runner.NCIPLOTExecutable.from_servername",
            return_value=mock_exe,
        ) as mock_from_servername:
            assert runner.executable is mock_exe
        mock_from_servername.assert_called_once_with(
            servername=pbs_server.name
        )

    def test_executable_raises_file_not_found(self, pbs_server):
        runner = NCIPLOTJobRunner(server=pbs_server, scratch=False)
        with patch(
            "chemsmart.jobs.nciplot.runner.NCIPLOTExecutable.from_servername",
            side_effect=FileNotFoundError("no server file"),
        ):
            with pytest.raises(FileNotFoundError):
                runner.executable

    def test_get_executable_returns_path(self, pbs_server):
        runner = NCIPLOTJobRunner(server=pbs_server, scratch=False)
        mock_exe = MagicMock()
        mock_exe.get_executable.return_value = "/path/to/nciplot"
        with patch.object(
            NCIPLOTJobRunner,
            "executable",
            new_callable=lambda: property(lambda self: mock_exe),
        ):
            assert runner._get_executable() == "/path/to/nciplot"


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


class TestNCIPLOTJobRunnerWriteXyzFromInputFiles:
    def test_raises_for_missing_file(
        self, pbs_server, nciplot_job_from_files, tmp_path
    ):
        runner = NCIPLOTJobRunner(server=pbs_server, scratch=False)
        runner._assign_variables(nciplot_job_from_files)
        os.makedirs(runner.running_directory, exist_ok=True)
        nciplot_job_from_files.filenames = [
            str(tmp_path / "does_not_exist.com")
        ]

        with pytest.raises(FileNotFoundError):
            runner._write_xyz_from_input_files(nciplot_job_from_files)

    def test_converts_and_copies_promolecular_xyz(
        self, pbs_server, nciplot_job_from_files, tmp_path
    ):
        runner = NCIPLOTJobRunner(server=pbs_server, scratch=False)
        runner._assign_variables(nciplot_job_from_files)
        os.makedirs(runner.running_directory, exist_ok=True)

        input_file = tmp_path / "molecule.com"
        input_file.write_text("dummy input")
        nciplot_job_from_files.filenames = [str(input_file)]

        expected_xyz = tmp_path / "molecule.xyz"

        def _fake_convert_files(self):
            expected_xyz.write_text("3\nmol\nC 0 0 0\nH 1 0 0\nH 0 1 0\n")

        with patch(
            "chemsmart.jobs.nciplot.runner.FileConverter.convert_files",
            _fake_convert_files,
        ):
            runner._write_xyz_from_input_files(nciplot_job_from_files)

        copied = os.path.join(
            runner.running_directory, "molecule_promolecular.xyz"
        )
        assert os.path.exists(copied)

    def test_wraps_conversion_failure_as_value_error(
        self, pbs_server, nciplot_job_from_files, tmp_path
    ):
        runner = NCIPLOTJobRunner(server=pbs_server, scratch=False)
        runner._assign_variables(nciplot_job_from_files)
        os.makedirs(runner.running_directory, exist_ok=True)

        input_file = tmp_path / "molecule.com"
        input_file.write_text("dummy input")
        nciplot_job_from_files.filenames = [str(input_file)]

        with patch(
            "chemsmart.jobs.nciplot.runner.FileConverter.convert_files",
            side_effect=RuntimeError("bad format"),
        ):
            with pytest.raises(ValueError, match="Could not convert file"):
                runner._write_xyz_from_input_files(nciplot_job_from_files)


class TestNCIPLOTJobRunnerWriteInput:
    def test_write_input_delegates_to_writer(
        self, pbs_server, nciplot_job_from_files
    ):
        runner = NCIPLOTJobRunner(server=pbs_server, scratch=False)
        runner._assign_variables(nciplot_job_from_files)
        os.makedirs(runner.running_directory, exist_ok=True)

        with patch(
            "chemsmart.jobs.nciplot.writer.NCIPLOTInputWriter"
        ) as mock_writer_cls:
            mock_writer = mock_writer_cls.return_value
            runner._write_input(nciplot_job_from_files)

        mock_writer_cls.assert_called_once_with(job=nciplot_job_from_files)
        mock_writer.write.assert_called_once_with(
            target_directory=runner.running_directory
        )


class TestNCIPLOTJobRunnerCreateProcess:
    def test_create_process_invokes_popen(
        self, pbs_server, nciplot_job_from_files
    ):
        runner = NCIPLOTJobRunner(server=pbs_server, scratch=False)
        runner._assign_variables(nciplot_job_from_files)
        os.makedirs(runner.running_directory, exist_ok=True)

        mock_exe = MagicMock()
        mock_exe.env = {"PATH": "/usr/bin"}
        mock_process = MagicMock()
        with (
            patch.object(
                NCIPLOTJobRunner,
                "executable",
                new_callable=lambda: property(lambda self: mock_exe),
            ),
            patch(
                "chemsmart.jobs.nciplot.runner.subprocess.Popen",
                return_value=mock_process,
            ) as mock_popen,
        ):
            result = runner._create_process(
                nciplot_job_from_files, "echo hello", env={"PATH": "/usr/bin"}
            )

        assert result is mock_process
        args, kwargs = mock_popen.call_args
        assert args[0] == ["echo", "hello"]
        assert kwargs["cwd"] == runner.running_directory
        assert os.path.exists(runner.job_outputfile)
        assert os.path.exists(runner.job_errfile)


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


class TestFakeNCIPLOTJobRunner:
    def test_run_writes_output_and_appends_fake_suffix_to_label(
        self,
        pbs_server,
        single_molecule_xyz_file,
        nciplot_settings,
        tmp_path,
        monkeypatch,
    ):
        monkeypatch.chdir(tmp_path)
        runner = FakeNCIPLOTJobRunner(
            server=pbs_server, scratch=False, fake=True
        )
        job = NCIPLOTJob(
            filenames=[single_molecule_xyz_file],
            settings=nciplot_settings,
            label="ncitest",
            jobrunner=runner,
        )
        original_label = job.label

        returncode = runner.run(job)

        assert returncode is None or returncode == 0
        assert job.label == original_label + "_fake"
        assert os.path.exists(runner.job_outputfile)
        with open(runner.job_outputfile) as f:
            content = f.read()
        assert "O   R   C   A" not in content  # sanity: not ORCA's fake output
        assert "NCIPLOT" in content


class TestFakeNCIPLOT:
    def test_raises_for_missing_file(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="not found"):
            FakeNCIPLOT(str(tmp_path / "missing.nci"))

    def test_properties(self, tmp_path):
        input_file = tmp_path / "job.nci"
        input_file.write_text("&nciplot\n&end\n")
        fake = FakeNCIPLOT(str(input_file))

        assert fake.file_folder == str(tmp_path)
        assert fake.filename == "job.nci"
        assert fake.input_filepath == str(input_file)
        assert fake.output_filepath == os.path.join(
            str(tmp_path), "job.nciout"
        )

    def test_run_writes_fake_output_with_header_and_input_contents(
        self, tmp_path
    ):
        input_file = tmp_path / "job.nci"
        input_file.write_text("&nciplot\ndensity job.wfn\n&end\n")
        fake = FakeNCIPLOT(str(input_file))

        returncode = fake.run()

        assert returncode is None
        assert os.path.exists(fake.output_filepath)
        with open(fake.output_filepath) as f:
            content = f.read()
        assert "NCIPLOT" in content
        assert "INPUT INFORMATION" in content
        assert "density job.wfn" in content
        assert "Start --" in content
        assert "End --" in content
