"""Direct unit coverage for XTBJobRunner/FakeXTBJobRunner's file-management
and execution-plumbing helpers, which are otherwise only exercised
indirectly (and incompletely) through full job-run integration tests
elsewhere."""

import os
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.jobs.xtb.runner import FakeXTB, FakeXTBJobRunner, XTBJobRunner
from chemsmart.jobs.xtb.settings import XTBJobSettings
from chemsmart.settings.executable import XTBExecutable


class DummyXTBJob:
    def __init__(self, folder, label, settings=None, molecule=None):
        self.folder = str(folder)
        self.label = label
        self.settings = settings or XTBJobSettings(charge=0, multiplicity=1)
        self.molecule = molecule if molecule is not None else MagicMock()

    @property
    def execution_name(self):
        return f"xtb_{self.label}"

    @property
    def xyzfile(self):
        return str(Path(self.folder) / f"{self.label}.xyz")

    @property
    def inputfile(self):
        return str(Path(self.folder) / f"{self.label}.inp")

    @property
    def errfile(self):
        return str(Path(self.folder) / f"{self.label}.err")

    @property
    def outputfile(self):
        return str(Path(self.folder) / f"{self.label}.out")

    def is_complete(self):
        return False


class TestXTBJobRunnerInit:
    def test_init_defaults_scratch_to_class_attribute(self, pbs_server):
        runner = XTBJobRunner(server=pbs_server, scratch=None, fake=True)
        assert runner.scratch is XTBJobRunner.SCRATCH

    def test_init_scratch_true_without_resolvable_dir_falls_back_to_false(
        self, pbs_server
    ):
        # The base JobRunner._set_scratch() can't determine a scratch dir
        # for this server fixture, so it downgrades scratch back to False.
        runner = XTBJobRunner(server=pbs_server, scratch=True, fake=True)
        assert runner.scratch is False


class TestXTBJobRunnerExecutable:
    def test_executable_property_returns_executable(
        self, pbs_server, monkeypatch
    ):
        runner = XTBJobRunner(server=pbs_server, scratch=False, fake=True)
        fake_executable = XTBExecutable(executable_folder=None)
        monkeypatch.setattr(
            XTBExecutable,
            "from_servername",
            staticmethod(lambda servername: fake_executable),
        )
        assert runner.executable is fake_executable

    def test_executable_property_wraps_filenotfounderror(
        self, pbs_server, monkeypatch
    ):
        runner = XTBJobRunner(server=pbs_server, scratch=False, fake=True)
        runner.server = SimpleNamespace(name="no-such-server")

        def _raise(servername):
            raise FileNotFoundError(f"no config for {servername}")

        monkeypatch.setattr(
            XTBExecutable, "from_servername", staticmethod(_raise)
        )
        with pytest.raises(FileNotFoundError):
            runner.executable


class TestAssignVariables:
    def test_assign_variables_dispatches_to_scratch_setup(
        self, pbs_server, tmp_path, monkeypatch
    ):
        runner = XTBJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(tmp_path),
            fake=True,
        )
        calls = []

        def _fake_scratch(job):
            calls.append("scratch")
            runner.running_directory = "dummy"
            runner.job_xyzfile = "dummy.xyz"
            runner.job_inputfile = "dummy.inp"
            runner.job_errfile = "dummy.err"

        def _fake_job_directory(job):
            calls.append("job_directory")
            runner.running_directory = "dummy"
            runner.job_xyzfile = "dummy.xyz"
            runner.job_inputfile = "dummy.inp"
            runner.job_errfile = "dummy.err"

        monkeypatch.setattr(
            runner, "_set_up_variables_in_scratch", _fake_scratch
        )
        monkeypatch.setattr(
            runner, "_set_up_variables_in_job_directory", _fake_job_directory
        )
        monkeypatch.setattr(
            type(runner),
            "executable",
            property(lambda self: SimpleNamespace(local_run=None)),
        )
        job = DummyXTBJob(folder=tmp_path / "jobfolder", label="assign_test")
        runner._assign_variables(job)
        assert calls == ["scratch"]
        assert not hasattr(job, "local")
        assert os.path.isdir(job.folder)

    def test_assign_variables_dispatches_to_job_directory_setup(
        self, pbs_server, tmp_path, monkeypatch
    ):
        runner = XTBJobRunner(server=pbs_server, scratch=False, fake=True)
        calls = []

        def _fake_scratch(job):
            calls.append("scratch")
            runner.running_directory = "dummy"
            runner.job_xyzfile = "dummy.xyz"
            runner.job_inputfile = "dummy.inp"
            runner.job_errfile = "dummy.err"

        def _fake_job_directory(job):
            calls.append("job_directory")
            runner.running_directory = "dummy"
            runner.job_xyzfile = "dummy.xyz"
            runner.job_inputfile = "dummy.inp"
            runner.job_errfile = "dummy.err"

        monkeypatch.setattr(
            runner, "_set_up_variables_in_scratch", _fake_scratch
        )
        monkeypatch.setattr(
            runner, "_set_up_variables_in_job_directory", _fake_job_directory
        )
        monkeypatch.setattr(
            type(runner),
            "executable",
            property(lambda self: SimpleNamespace(local_run=True)),
        )
        job = DummyXTBJob(folder=tmp_path / "jobfolder2", label="assign_test")
        runner._assign_variables(job)
        assert calls == ["job_directory"]
        assert job.local is True

    def test_assign_variables_scratch_true_but_no_scratch_dir_uses_job_directory(
        self, pbs_server, tmp_path, monkeypatch
    ):
        # scratch=True alone isn't enough -- scratch_dir must also be set.
        runner = XTBJobRunner(server=pbs_server, scratch=True, fake=True)
        calls = []

        def _fake_scratch(job):
            calls.append("scratch")
            runner.running_directory = "dummy"
            runner.job_xyzfile = "dummy.xyz"
            runner.job_inputfile = "dummy.inp"
            runner.job_errfile = "dummy.err"

        def _fake_job_directory(job):
            calls.append("job_directory")
            runner.running_directory = "dummy"
            runner.job_xyzfile = "dummy.xyz"
            runner.job_inputfile = "dummy.inp"
            runner.job_errfile = "dummy.err"

        monkeypatch.setattr(
            runner, "_set_up_variables_in_scratch", _fake_scratch
        )
        monkeypatch.setattr(
            runner, "_set_up_variables_in_job_directory", _fake_job_directory
        )
        monkeypatch.setattr(
            type(runner),
            "executable",
            property(lambda self: SimpleNamespace(local_run=None)),
        )
        job = DummyXTBJob(folder=tmp_path / "jobfolder3", label="assign_test")
        runner._assign_variables(job)
        assert calls == ["job_directory"]


class TestSetUpVariablesInScratch:
    def test_creates_fresh_scratch_directory_and_sets_paths(
        self, pbs_server, tmp_path
    ):
        scratch_root = tmp_path / "scratch_root"
        runner = XTBJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(scratch_root),
            fake=True,
        )
        job = DummyXTBJob(folder=tmp_path / "jobfolder", label="myjob")
        runner._set_up_variables_in_scratch(job)

        expected_dir = os.path.join(str(scratch_root), job.execution_name)
        assert runner.running_directory == expected_dir
        assert os.path.isdir(expected_dir)
        assert runner.job_xyzfile == os.path.abspath(
            os.path.join(expected_dir, f"{job.label}.xyz")
        )
        assert runner.job_inputfile == os.path.abspath(
            os.path.join(expected_dir, f"{job.label}.inp")
        )
        assert runner.job_errfile == os.path.abspath(job.errfile)

    def test_reuses_existing_scratch_directory(self, pbs_server, tmp_path):
        scratch_root = tmp_path / "scratch_root2"
        job = DummyXTBJob(folder=tmp_path / "jobfolder2", label="existing")
        expected_dir = os.path.join(str(scratch_root), job.execution_name)
        os.makedirs(expected_dir)

        runner = XTBJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(scratch_root),
            fake=True,
        )
        runner._set_up_variables_in_scratch(job)
        assert runner.running_directory == expected_dir


class TestSetUpVariablesInJobDirectory:
    def test_sets_paths_from_job_folder(self, pbs_server, tmp_path):
        runner = XTBJobRunner(server=pbs_server, scratch=False, fake=True)
        job = DummyXTBJob(folder=tmp_path, label="myjob")
        runner._set_up_variables_in_job_directory(job)
        assert runner.running_directory == job.folder
        assert runner.job_xyzfile == os.path.abspath(job.xyzfile)
        assert runner.job_inputfile == os.path.abspath(job.inputfile)
        assert runner.job_errfile == os.path.abspath(job.errfile)


class TestWriteInput:
    def test_writes_xyz_via_molecule(self, pbs_server, tmp_path):
        runner = XTBJobRunner(server=pbs_server, scratch=False, fake=True)
        job = DummyXTBJob(folder=tmp_path, label="myjob")
        runner._set_up_variables_in_job_directory(job)
        runner._write_input(job)
        job.molecule.write_xyz.assert_called_once_with(
            runner.job_xyzfile, mode="w"
        )

    def test_writes_constraints_input_when_present(self, pbs_server, tmp_path):
        settings = XTBJobSettings(
            charge=0, multiplicity=1, constraints=[[1, 2]], force_constant=0.5
        )
        job = DummyXTBJob(folder=tmp_path, label="myjob", settings=settings)
        runner = XTBJobRunner(server=pbs_server, scratch=False, fake=True)
        runner._set_up_variables_in_job_directory(job)
        with patch(
            "chemsmart.jobs.xtb.writer.XTBInputWriter"
        ) as mock_writer_cls:
            mock_writer = MagicMock()
            mock_writer_cls.return_value = mock_writer
            runner._write_input(job)
        mock_writer_cls.assert_called_once_with(job=job)
        mock_writer.write.assert_called_once_with(
            target_directory=runner.running_directory
        )

    def test_no_constraints_skips_input_writer(self, pbs_server, tmp_path):
        job = DummyXTBJob(folder=tmp_path, label="myjob")
        runner = XTBJobRunner(server=pbs_server, scratch=False, fake=True)
        runner._set_up_variables_in_job_directory(job)
        with patch(
            "chemsmart.jobs.xtb.writer.XTBInputWriter"
        ) as mock_writer_cls:
            runner._write_input(job)
        mock_writer_cls.assert_not_called()


class TestGetCommand:
    def test_get_command_uses_executable_and_settings_args(
        self, pbs_server, tmp_path, monkeypatch
    ):
        settings = XTBJobSettings(charge=0, multiplicity=1, jobtype="sp")
        job = DummyXTBJob(folder=tmp_path, label="myjob", settings=settings)
        runner = XTBJobRunner(server=pbs_server, scratch=False, fake=True)
        runner._set_up_variables_in_job_directory(job)
        monkeypatch.setattr(
            type(runner),
            "executable",
            property(
                lambda self: SimpleNamespace(get_executable=lambda: "xtb")
            ),
        )
        command = runner._get_command(job)
        assert command[0] == "xtb"
        assert command[1] == runner.job_xyzfile
        assert "--chrg" in command
        assert "0" in command


class TestGetSettingsArgs:
    def _runner(self, pbs_server):
        runner = XTBJobRunner(server=pbs_server, scratch=False, fake=True)
        runner.job_inputfile = "/tmp/fake_job.inp"
        return runner

    def test_opt_jobtype_without_freq_uses_opt_flag(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = XTBJobSettings(
            charge=0,
            multiplicity=1,
            jobtype="opt",
            freq=False,
            optimization_level="tight",
        )
        args = runner.get_settings_args(settings)
        assert "--opt" in args
        assert "--ohess" not in args
        assert "tight" in args

    def test_opt_jobtype_with_freq_uses_ohess_flag(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = XTBJobSettings(
            charge=0, multiplicity=1, jobtype="opt", freq=True
        )
        args = runner.get_settings_args(settings)
        assert "--ohess" in args
        assert "--opt" not in args

    def test_opt_jobtype_without_optimization_level_omits_level_arg(
        self, pbs_server
    ):
        runner = self._runner(pbs_server)
        settings = XTBJobSettings(
            charge=0, multiplicity=1, jobtype="opt", optimization_level=None
        )
        args = runner.get_settings_args(settings)
        assert args.count("--opt") == 1
        # no bare optimization-level token follows --opt
        idx = args.index("--opt")
        assert args[idx + 1] == "--chrg"

    def test_hess_jobtype_uses_hess_flag(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = XTBJobSettings(charge=0, multiplicity=1, jobtype="hess")
        args = runner.get_settings_args(settings)
        assert "--hess" in args

    def test_sp_jobtype_adds_no_jobtype_flag(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = XTBJobSettings(charge=0, multiplicity=1, jobtype="sp")
        args = runner.get_settings_args(settings)
        assert "--opt" not in args
        assert "--hess" not in args

    def test_unsupported_jobtype_raises(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = XTBJobSettings(charge=0, multiplicity=1, jobtype="path")
        with pytest.raises(ValueError, match="Unsupported xTB jobtype"):
            runner.get_settings_args(settings)

    def test_charge_and_uhf_from_multiplicity(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = XTBJobSettings(charge=-1, multiplicity=3, jobtype="sp")
        args = runner.get_settings_args(settings)
        chrg_idx = args.index("--chrg")
        assert args[chrg_idx + 1] == "-1"
        uhf_idx = args.index("--uhf")
        assert args[uhf_idx + 1] == "2"

    def test_solvent_model_and_id_both_set_adds_solvent_flag(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = XTBJobSettings(
            charge=0,
            multiplicity=1,
            jobtype="sp",
            solvent_model="gbsa",
            solvent_id="water",
        )
        args = runner.get_settings_args(settings)
        assert "--gbsa" in args
        idx = args.index("--gbsa")
        assert args[idx + 1] == "water"

    def test_solvent_model_without_id_omits_solvent_flag(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = XTBJobSettings(
            charge=0, multiplicity=1, jobtype="sp", solvent_model="gbsa"
        )
        args = runner.get_settings_args(settings)
        assert "--gbsa" not in args

    def test_grad_flag_added_when_true(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = XTBJobSettings(
            charge=0, multiplicity=1, jobtype="sp", grad=True
        )
        args = runner.get_settings_args(settings)
        assert "--grad" in args

    def test_grad_flag_omitted_when_false(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = XTBJobSettings(
            charge=0, multiplicity=1, jobtype="sp", grad=False
        )
        args = runner.get_settings_args(settings)
        assert "--grad" not in args

    def test_constraints_add_input_flag_with_basename(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = XTBJobSettings(
            charge=0, multiplicity=1, jobtype="sp", constraints=[[1, 2]]
        )
        args = runner.get_settings_args(settings)
        assert "--input" in args
        idx = args.index("--input")
        assert args[idx + 1] == "fake_job.inp"

    def test_no_constraints_omits_input_flag(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = XTBJobSettings(charge=0, multiplicity=1, jobtype="sp")
        args = runner.get_settings_args(settings)
        assert "--input" not in args

    def test_additional_flags_appended_split_on_whitespace(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = XTBJobSettings(
            charge=0,
            multiplicity=1,
            jobtype="sp",
            additional_flags="--squick --foo bar",
        )
        args = runner.get_settings_args(settings)
        assert args[-3:] == ["--squick", "--foo", "bar"]

    def test_no_additional_flags_adds_nothing_extra(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = XTBJobSettings(charge=0, multiplicity=1, jobtype="sp")
        args = runner.get_settings_args(settings)
        assert args[-1] != "None"


class TestGfnArgs:
    def test_none_returns_empty_list(self, pbs_server):
        assert self._runner(pbs_server)._gfn_args(None) == []

    def test_gfn_digit_version_uses_gfn_flag_with_digit(self, pbs_server):
        assert self._runner(pbs_server)._gfn_args("gfn2") == ["--gfn", "2"]

    def test_gfnff_uses_dedicated_flag(self, pbs_server):
        assert self._runner(pbs_server)._gfn_args("gfnff") == ["--gfnff"]

    def test_other_value_uses_generic_dashdash_flag(self, pbs_server):
        assert self._runner(pbs_server)._gfn_args("gfn2//gfnff") == [
            "--gfn2//gfnff"
        ]

    def _runner(self, pbs_server):
        return XTBJobRunner(server=pbs_server, scratch=False, fake=True)


class TestCreateProcess:
    def test_create_process_opens_files_and_calls_popen(
        self, pbs_server, tmp_path
    ):
        runner = XTBJobRunner(server=pbs_server, scratch=False, fake=True)
        runner.running_directory = str(tmp_path)
        runner.job_outputfile = str(tmp_path / "job.out")
        runner.job_errfile = str(tmp_path / "job.err")
        fake_executable = SimpleNamespace(env={"FOO": "bar"})
        with (
            patch.object(
                type(runner),
                "executable",
                property(lambda self: fake_executable),
            ),
            patch("chemsmart.jobs.xtb.runner.subprocess.Popen") as mock_popen,
        ):
            mock_process = MagicMock()
            mock_popen.return_value = mock_process
            job = DummyXTBJob(folder=tmp_path, label="myjob")
            result = runner._create_process(job, ["xtb", "mol.xyz"], {})

        assert result is mock_process
        mock_popen.assert_called_once()
        call_args, call_kwargs = mock_popen.call_args
        assert call_args[0] == ["xtb", "mol.xyz"]
        assert call_kwargs["cwd"] == str(tmp_path)
        assert os.path.exists(runner.job_outputfile)
        assert os.path.exists(runner.job_errfile)


class TestPostrun:
    def test_postrun_no_op_when_not_scratch(self, pbs_server, tmp_path):
        runner = XTBJobRunner(server=pbs_server, scratch=False, fake=True)
        job = DummyXTBJob(folder=tmp_path, label="myjob")
        # Should return immediately without touching running_directory.
        runner._postrun(job)

    def test_postrun_copies_scratch_files_to_job_folder(
        self, pbs_server, tmp_path
    ):
        scratch_dir = tmp_path / "scratch"
        scratch_dir.mkdir()
        (scratch_dir / "myjob.xyz").write_text("3\n\nO 0 0 0\n")
        job_folder = tmp_path / "jobfolder"
        job_folder.mkdir()

        runner = XTBJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(tmp_path),
            fake=True,
        )
        runner.running_directory = str(scratch_dir)
        job = DummyXTBJob(folder=job_folder, label="myjob")

        runner._postrun(job)

        assert (job_folder / "myjob.xyz").exists()

    def test_postrun_skips_file_already_at_destination(
        self, pbs_server, tmp_path
    ):
        # running_directory == job.folder, so source and destination
        # resolve to the same path and the copy is skipped.
        job_folder = tmp_path / "jobfolder"
        job_folder.mkdir()
        (job_folder / "myjob.xyz").write_text("3\n\nO 0 0 0\n")

        runner = XTBJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(tmp_path),
            fake=True,
        )
        runner.running_directory = str(job_folder)
        job = DummyXTBJob(folder=job_folder, label="myjob")

        runner._postrun(job)
        assert (job_folder / "myjob.xyz").exists()

    def test_postrun_suppresses_isadirectoryerror_for_subdirectories(
        self, pbs_server, tmp_path
    ):
        scratch_dir = tmp_path / "scratch2"
        (scratch_dir / "subdir").mkdir(parents=True)
        job_folder = tmp_path / "jobfolder2"
        job_folder.mkdir()

        runner = XTBJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(tmp_path),
            fake=True,
        )
        runner.running_directory = str(scratch_dir)
        job = DummyXTBJob(folder=job_folder, label="myjob")

        # Should not raise even though "subdir" is a directory, not a file.
        runner._postrun(job)


class TestFakeXTBJobRunner:
    def test_executable_property_is_local_fake(self, pbs_server):
        runner = FakeXTBJobRunner(server=pbs_server, scratch=False, fake=True)
        executable = runner.executable
        assert isinstance(executable, XTBExecutable)
        assert executable.local_run is True
        assert executable.executable_folder is None

    def test_set_up_variables_in_scratch_appends_fake_suffix(
        self, pbs_server, tmp_path
    ):
        runner = FakeXTBJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(tmp_path),
            fake=True,
        )
        job = DummyXTBJob(folder=tmp_path / "jobfolder", label="myjob")
        # The scratch directory name is derived from the program-prefixed
        # label as it was *before* the "_fake" suffix rename, since
        # running_directory is computed first and the rename happens only
        # afterward.
        expected_dir = os.path.join(str(tmp_path), job.execution_name)
        runner._set_up_variables_in_scratch(job)
        assert job.label == "myjob_fake"
        assert runner.running_directory == expected_dir
        assert runner.job_outputfile == os.path.abspath(job.outputfile)
        assert runner.job_errfile == os.path.abspath(job.errfile)

    def test_set_up_variables_in_scratch_reuses_existing_directory(
        self, pbs_server, tmp_path
    ):
        runner = FakeXTBJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(tmp_path),
            fake=True,
        )
        job = DummyXTBJob(folder=tmp_path / "jobfolder", label="existing")
        expected_dir = os.path.join(str(tmp_path), job.execution_name)
        os.makedirs(expected_dir)

        runner._set_up_variables_in_scratch(job)
        assert runner.running_directory == expected_dir

    def test_set_up_variables_in_job_directory_appends_fake_suffix(
        self, pbs_server, tmp_path
    ):
        runner = FakeXTBJobRunner(server=pbs_server, scratch=False, fake=True)
        job = DummyXTBJob(folder=tmp_path, label="myjob")
        runner._set_up_variables_in_job_directory(job)
        assert job.label == "myjob_fake"
        assert runner.running_directory == job.folder
        assert runner.job_xyzfile == os.path.abspath(job.xyzfile)

    def test_run_executes_full_fake_pipeline(self, pbs_server, tmp_path):
        runner = FakeXTBJobRunner(server=pbs_server, scratch=False, fake=True)
        settings = XTBJobSettings(charge=0, multiplicity=1, jobtype="sp")

        def _write_xyz(path, mode="w"):
            with open(path, "w") as f:
                f.write("1\n\nH 0 0 0\n")

        molecule = MagicMock()
        molecule.write_xyz.side_effect = _write_xyz
        job = DummyXTBJob(
            folder=tmp_path,
            label="myjob",
            settings=settings,
            molecule=molecule,
        )

        exit_code = runner.run(job)

        assert exit_code == 0
        assert job.label == "myjob_fake"
        assert os.path.exists(job.outputfile)
        assert os.path.exists(job.errfile)
        assert "fake xtb" in Path(job.outputfile).read_text()


class TestFakeXTB:
    def test_raises_when_xyzfile_missing(self, tmp_path):
        fake = FakeXTB(
            xyzfile=str(tmp_path / "missing.xyz"),
            outputfile=str(tmp_path / "out.out"),
            errfile=str(tmp_path / "err.err"),
            command=["xtb", "missing.xyz"],
        )
        with pytest.raises(FileNotFoundError, match="not found"):
            fake.run()

    def test_writes_output_and_empty_err_file(self, tmp_path):
        xyzfile = tmp_path / "mol.xyz"
        xyzfile.write_text("1\n\nH 0 0 0\n")
        outputfile = tmp_path / "job.out"
        errfile = tmp_path / "job.err"
        fake = FakeXTB(
            xyzfile=str(xyzfile),
            outputfile=str(outputfile),
            errfile=str(errfile),
            command=["xtb", "mol.xyz", "--gfn", "2"],
        )
        result = fake.run()
        assert result == 0
        content = outputfile.read_text()
        assert "xtb version 0.0.0 (Fake)" in content
        assert "xtb mol.xyz --gfn 2" in content
        assert "finished run (fake xtb)" in content
        assert errfile.read_text() == ""
