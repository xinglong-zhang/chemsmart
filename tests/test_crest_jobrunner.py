"""Direct unit coverage for CRESTJobRunner/FakeCRESTJobRunner's
file-management and execution-plumbing helpers, which are otherwise only
exercised indirectly (and incompletely) through full job-run integration
tests elsewhere."""

import os
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.jobs.crest.runner import (
    CRESTJobRunner,
    FakeCREST,
    FakeCRESTJobRunner,
)
from chemsmart.jobs.crest.settings import CRESTJobSettings
from chemsmart.settings.executable import CRESTExecutable


class DummyCRESTJob:
    def __init__(self, folder, label, settings=None, molecule=None):
        self.folder = str(folder)
        self.label = label
        self.settings = settings or CRESTJobSettings(charge=0, multiplicity=1)
        self.molecule = molecule if molecule is not None else MagicMock()

    @property
    def execution_name(self):
        return f"crest_{self.label}"

    @property
    def xyzfile(self):
        return str(Path(self.folder) / f"{self.label}.xyz")

    @property
    def errfile(self):
        return str(Path(self.folder) / f"{self.label}.err")

    @property
    def outputfile(self):
        return str(Path(self.folder) / f"{self.label}.out")

    def is_complete(self):
        return False


class TestCRESTJobRunnerInit:
    def test_init_defaults_scratch_to_class_attribute(self, pbs_server):
        runner = CRESTJobRunner(server=pbs_server, scratch=None, fake=True)
        assert runner.scratch is CRESTJobRunner.SCRATCH

    def test_init_scratch_true_without_resolvable_dir_falls_back_to_false(
        self, pbs_server
    ):
        runner = CRESTJobRunner(server=pbs_server, scratch=True, fake=True)
        assert runner.scratch is False


class TestCRESTJobRunnerExecutable:
    def test_executable_property_returns_executable(
        self, pbs_server, monkeypatch
    ):
        runner = CRESTJobRunner(server=pbs_server, scratch=False, fake=True)
        fake_executable = CRESTExecutable(executable_folder=None)
        monkeypatch.setattr(
            CRESTExecutable,
            "from_servername",
            staticmethod(lambda servername: fake_executable),
        )
        assert runner.executable is fake_executable

    def test_executable_property_wraps_filenotfounderror(
        self, pbs_server, monkeypatch
    ):
        runner = CRESTJobRunner(server=pbs_server, scratch=False, fake=True)
        runner.server = SimpleNamespace(name="no-such-server")

        def _raise(servername):
            raise FileNotFoundError(f"no config for {servername}")

        monkeypatch.setattr(
            CRESTExecutable, "from_servername", staticmethod(_raise)
        )
        with pytest.raises(FileNotFoundError):
            runner.executable


class TestAssignVariables:
    def _patched(self, runner, monkeypatch, local_run):
        calls = []

        def _fake_scratch(job):
            calls.append("scratch")
            runner.running_directory = "dummy"
            runner.job_xyzfile = "dummy.xyz"
            runner.job_errfile = "dummy.err"

        def _fake_job_directory(job):
            calls.append("job_directory")
            runner.running_directory = "dummy"
            runner.job_xyzfile = "dummy.xyz"
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
            property(lambda self: SimpleNamespace(local_run=local_run)),
        )
        return calls

    def test_assign_variables_dispatches_to_scratch_setup(
        self, pbs_server, tmp_path, monkeypatch
    ):
        runner = CRESTJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(tmp_path),
            fake=True,
        )
        calls = self._patched(runner, monkeypatch, local_run=None)
        job = DummyCRESTJob(folder=tmp_path / "jobfolder", label="assign_test")
        runner._assign_variables(job)
        assert calls == ["scratch"]
        assert not hasattr(job, "local")
        assert os.path.isdir(job.folder)

    def test_assign_variables_dispatches_to_job_directory_setup(
        self, pbs_server, tmp_path, monkeypatch
    ):
        runner = CRESTJobRunner(server=pbs_server, scratch=False, fake=True)
        calls = self._patched(runner, monkeypatch, local_run=True)
        job = DummyCRESTJob(
            folder=tmp_path / "jobfolder2", label="assign_test"
        )
        runner._assign_variables(job)
        assert calls == ["job_directory"]
        assert job.local is True


class TestSetUpVariablesInScratch:
    def test_creates_fresh_scratch_directory_and_sets_paths(
        self, pbs_server, tmp_path
    ):
        scratch_root = tmp_path / "scratch_root"
        runner = CRESTJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(scratch_root),
            fake=True,
        )
        job = DummyCRESTJob(folder=tmp_path / "jobfolder", label="myjob")
        runner._set_up_variables_in_scratch(job)

        expected_dir = os.path.join(str(scratch_root), job.execution_name)
        assert runner.running_directory == expected_dir
        assert os.path.isdir(expected_dir)
        assert runner.job_xyzfile == os.path.abspath(
            os.path.join(expected_dir, f"{job.label}.xyz")
        )
        assert runner.job_errfile == os.path.abspath(job.errfile)

    def test_reuses_existing_scratch_directory(self, pbs_server, tmp_path):
        scratch_root = tmp_path / "scratch_root2"
        job = DummyCRESTJob(folder=tmp_path / "jobfolder2", label="existing")
        expected_dir = os.path.join(str(scratch_root), job.execution_name)
        os.makedirs(expected_dir)

        runner = CRESTJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(scratch_root),
            fake=True,
        )
        runner._set_up_variables_in_scratch(job)
        assert runner.running_directory == expected_dir


class TestSetUpVariablesInJobDirectory:
    def test_sets_paths_from_job_folder(self, pbs_server, tmp_path):
        runner = CRESTJobRunner(server=pbs_server, scratch=False, fake=True)
        job = DummyCRESTJob(folder=tmp_path, label="myjob")
        runner._set_up_variables_in_job_directory(job)
        assert runner.running_directory == job.folder
        assert runner.job_xyzfile == os.path.abspath(job.xyzfile)
        assert runner.job_errfile == os.path.abspath(job.errfile)


class TestWriteInput:
    def test_writes_xyz_and_always_invokes_input_writer(
        self, pbs_server, tmp_path
    ):
        job = DummyCRESTJob(folder=tmp_path, label="myjob")
        runner = CRESTJobRunner(server=pbs_server, scratch=False, fake=True)
        runner._set_up_variables_in_job_directory(job)
        with patch(
            "chemsmart.jobs.crest.writer.CRESTInputWriter"
        ) as mock_writer_cls:
            mock_writer = MagicMock()
            mock_writer_cls.return_value = mock_writer
            runner._write_input(job)
        job.molecule.write_xyz.assert_called_once_with(
            runner.job_xyzfile, mode="w"
        )
        mock_writer_cls.assert_called_once_with(job=job)
        mock_writer.write.assert_called_once_with(
            target_directory=runner.running_directory
        )


class TestGetCommand:
    def test_get_command_uses_executable_and_settings_args(
        self, pbs_server, tmp_path, monkeypatch
    ):
        settings = CRESTJobSettings(charge=0, multiplicity=1)
        job = DummyCRESTJob(folder=tmp_path, label="myjob", settings=settings)
        runner = CRESTJobRunner(server=pbs_server, scratch=False, fake=True)
        runner._set_up_variables_in_job_directory(job)
        monkeypatch.setattr(
            type(runner),
            "executable",
            property(
                lambda self: SimpleNamespace(get_executable=lambda: "crest")
            ),
        )
        command = runner._get_command(job)
        assert command[0] == "crest"
        assert command[1] == runner.job_xyzfile
        assert "--chrg" in command
        assert "0" in command


class TestGetSettingsArgs:
    def _runner(self, pbs_server):
        return CRESTJobRunner(server=pbs_server, scratch=False, fake=True)

    def test_charge_and_uhf_from_multiplicity(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = CRESTJobSettings(charge=-1, multiplicity=3)
        args = runner.get_settings_args(settings)
        chrg_idx = args.index("--chrg")
        assert args[chrg_idx + 1] == "-1"
        uhf_idx = args.index("--uhf")
        assert args[uhf_idx + 1] == "2"

    def test_no_charge_or_multiplicity_omits_flags(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = CRESTJobSettings()
        args = runner.get_settings_args(settings)
        assert "--chrg" not in args
        assert "--uhf" not in args

    def test_solvent_model_and_id_both_set_adds_solvent_flag(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = CRESTJobSettings(
            charge=0,
            multiplicity=1,
            solvent_model="gbsa",
            solvent_id="water",
        )
        args = runner.get_settings_args(settings)
        assert "--gbsa" in args
        idx = args.index("--gbsa")
        assert args[idx + 1] == "water"

    def test_solvent_model_without_id_omits_solvent_flag(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = CRESTJobSettings(
            charge=0, multiplicity=1, solvent_model="gbsa"
        )
        args = runner.get_settings_args(settings)
        assert "--gbsa" not in args

    @pytest.mark.parametrize(
        "field,flag,value",
        [
            ("energy_window", "--ewin", 6.0),
            ("rmsd_threshold", "--rthr", 0.5),
            ("energy_threshold", "--ethr", 0.1),
            ("bconst_threshold", "--bthr", 0.02),
            ("population_threshold", "--pthr", 0.05),
            ("temperature", "--temp", 298.0),
            ("md_timestep", "--tstep", 4),
            ("md_length", "--mdlen", 10),
            ("md_dump_step", "--mddump", 500),
            ("vbias_dump_interval", "--vbdump", 100),
            ("additional_md_temperature", "--tnmd", 400),
        ],
    )
    def test_optional_numeric_flag_added_when_set(
        self, pbs_server, field, flag, value
    ):
        runner = self._runner(pbs_server)
        settings = CRESTJobSettings(**{field: value})
        args = runner.get_settings_args(settings)
        assert flag in args
        idx = args.index(flag)
        assert args[idx + 1] == str(value)

    @pytest.mark.parametrize(
        "field,flag",
        [
            ("energy_window", "--ewin"),
            ("rmsd_threshold", "--rthr"),
            ("energy_threshold", "--ethr"),
            ("bconst_threshold", "--bthr"),
            ("population_threshold", "--pthr"),
            ("temperature", "--temp"),
            ("md_timestep", "--tstep"),
            ("md_length", "--mdlen"),
            ("md_dump_step", "--mddump"),
            ("vbias_dump_interval", "--vbdump"),
            ("additional_md_temperature", "--tnmd"),
        ],
    )
    def test_optional_numeric_flag_omitted_when_unset(
        self, pbs_server, field, flag
    ):
        runner = self._runner(pbs_server)
        settings = CRESTJobSettings()
        args = runner.get_settings_args(settings)
        assert flag not in args

    def test_no_topology_check_flag(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = CRESTJobSettings(no_topology_check=True)
        args = runner.get_settings_args(settings)
        assert "--notopo" in args

    def test_no_reference_topology_check_flag(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = CRESTJobSettings(no_reference_topology_check=True)
        args = runner.get_settings_args(settings)
        assert "--noreftopo" in args

    def test_topology_check_flags_omitted_by_default(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = CRESTJobSettings()
        args = runner.get_settings_args(settings)
        assert "--notopo" not in args
        assert "--noreftopo" not in args

    def test_optimization_level_flag(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = CRESTJobSettings(optimization_level="tight")
        args = runner.get_settings_args(settings)
        idx = args.index("--optlev")
        assert args[idx + 1] == "tight"

    def test_quick_mode_flag(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = CRESTJobSettings(quick_mode="squick")
        args = runner.get_settings_args(settings)
        assert "--squick" in args

    def test_nci_flag(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = CRESTJobSettings(nci=True)
        args = runner.get_settings_args(settings)
        assert "--nci" in args

    def test_nci_flag_omitted_when_false(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = CRESTJobSettings(nci=False)
        args = runner.get_settings_args(settings)
        assert "--nci" not in args

    def test_nprocs_flag(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = CRESTJobSettings(nprocs=8)
        args = runner.get_settings_args(settings)
        idx = args.index("--T")
        assert args[idx + 1] == "8"

    def test_constraints_add_cinp_flag(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = CRESTJobSettings(constraints=[[1, 2]])
        args = runner.get_settings_args(settings)
        idx = args.index("--cinp")
        assert args[idx + 1] == "constraints.inp"

    def test_no_constraints_omits_cinp_flag(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = CRESTJobSettings()
        args = runner.get_settings_args(settings)
        assert "--cinp" not in args

    def test_additional_flags_appended_split_on_whitespace(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = CRESTJobSettings(additional_flags="--foo --bar baz")
        args = runner.get_settings_args(settings)
        assert args[-3:] == ["--foo", "--bar", "baz"]

    def test_no_additional_flags_adds_nothing_extra(self, pbs_server):
        runner = self._runner(pbs_server)
        settings = CRESTJobSettings()
        args = runner.get_settings_args(settings)
        assert "None" not in args


class TestGfnArgs:
    def test_none_returns_empty_list(self, pbs_server):
        runner = CRESTJobRunner(server=pbs_server, scratch=False, fake=True)
        assert runner._gfn_args(None) == []

    def test_value_uses_generic_dashdash_flag(self, pbs_server):
        runner = CRESTJobRunner(server=pbs_server, scratch=False, fake=True)
        assert runner._gfn_args("gfn2") == ["--gfn2"]
        assert runner._gfn_args("gfnff") == ["--gfnff"]


class TestCreateProcess:
    def test_create_process_opens_files_and_calls_popen(
        self, pbs_server, tmp_path
    ):
        runner = CRESTJobRunner(server=pbs_server, scratch=False, fake=True)
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
            patch(
                "chemsmart.jobs.crest.runner.subprocess.Popen"
            ) as mock_popen,
        ):
            mock_process = MagicMock()
            mock_popen.return_value = mock_process
            job = DummyCRESTJob(folder=tmp_path, label="myjob")
            result = runner._create_process(job, ["crest", "mol.xyz"], {})

        assert result is mock_process
        call_args, call_kwargs = mock_popen.call_args
        assert call_args[0] == ["crest", "mol.xyz"]
        assert call_kwargs["cwd"] == str(tmp_path)
        assert os.path.exists(runner.job_outputfile)
        assert os.path.exists(runner.job_errfile)


class TestPostrun:
    def test_postrun_no_op_when_not_scratch(self, pbs_server, tmp_path):
        runner = CRESTJobRunner(server=pbs_server, scratch=False, fake=True)
        job = DummyCRESTJob(folder=tmp_path, label="myjob")
        runner._postrun(job)

    def test_postrun_copies_scratch_files_to_job_folder(
        self, pbs_server, tmp_path
    ):
        scratch_dir = tmp_path / "scratch"
        scratch_dir.mkdir()
        (scratch_dir / "myjob.xyz").write_text("3\n\nO 0 0 0\n")
        job_folder = tmp_path / "jobfolder"
        job_folder.mkdir()

        runner = CRESTJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(tmp_path),
            fake=True,
        )
        runner.running_directory = str(scratch_dir)
        job = DummyCRESTJob(folder=job_folder, label="myjob")

        runner._postrun(job)
        assert (job_folder / "myjob.xyz").exists()

    def test_postrun_skips_file_already_at_destination(
        self, pbs_server, tmp_path
    ):
        job_folder = tmp_path / "jobfolder"
        job_folder.mkdir()
        (job_folder / "myjob.xyz").write_text("3\n\nO 0 0 0\n")

        runner = CRESTJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(tmp_path),
            fake=True,
        )
        runner.running_directory = str(job_folder)
        job = DummyCRESTJob(folder=job_folder, label="myjob")

        runner._postrun(job)
        assert (job_folder / "myjob.xyz").exists()

    def test_postrun_suppresses_isadirectoryerror_for_subdirectories(
        self, pbs_server, tmp_path
    ):
        scratch_dir = tmp_path / "scratch2"
        (scratch_dir / "subdir").mkdir(parents=True)
        job_folder = tmp_path / "jobfolder2"
        job_folder.mkdir()

        runner = CRESTJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(tmp_path),
            fake=True,
        )
        runner.running_directory = str(scratch_dir)
        job = DummyCRESTJob(folder=job_folder, label="myjob")

        runner._postrun(job)


class TestFakeCRESTJobRunner:
    def test_executable_property_is_local_fake(self, pbs_server):
        runner = FakeCRESTJobRunner(
            server=pbs_server, scratch=False, fake=True
        )
        executable = runner.executable
        assert isinstance(executable, CRESTExecutable)
        assert executable.local_run is True
        assert executable.executable_folder is None

    def test_set_up_variables_in_scratch_appends_fake_suffix(
        self, pbs_server, tmp_path
    ):
        runner = FakeCRESTJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(tmp_path),
            fake=True,
        )
        job = DummyCRESTJob(folder=tmp_path / "jobfolder", label="myjob")
        expected_dir = os.path.join(str(tmp_path), job.execution_name)
        runner._set_up_variables_in_scratch(job)
        assert job.label == "myjob_fake"
        assert runner.running_directory == expected_dir
        assert runner.job_outputfile == os.path.abspath(job.outputfile)
        assert runner.job_errfile == os.path.abspath(job.errfile)

    def test_set_up_variables_in_scratch_reuses_existing_directory(
        self, pbs_server, tmp_path
    ):
        runner = FakeCRESTJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(tmp_path),
            fake=True,
        )
        job = DummyCRESTJob(folder=tmp_path / "jobfolder", label="existing")
        expected_dir = os.path.join(str(tmp_path), job.execution_name)
        os.makedirs(expected_dir)

        runner._set_up_variables_in_scratch(job)
        assert runner.running_directory == expected_dir

    def test_set_up_variables_in_job_directory_appends_fake_suffix(
        self, pbs_server, tmp_path
    ):
        runner = FakeCRESTJobRunner(
            server=pbs_server, scratch=False, fake=True
        )
        job = DummyCRESTJob(folder=tmp_path, label="myjob")
        runner._set_up_variables_in_job_directory(job)
        assert job.label == "myjob_fake"
        assert runner.running_directory == job.folder
        assert runner.job_xyzfile == os.path.abspath(job.xyzfile)

    def test_run_executes_full_fake_pipeline(self, pbs_server, tmp_path):
        runner = FakeCRESTJobRunner(
            server=pbs_server, scratch=False, fake=True
        )
        settings = CRESTJobSettings(charge=0, multiplicity=1)

        def _write_xyz(path, mode="w"):
            with open(path, "w") as f:
                f.write("1\n\nH 0 0 0\n")

        molecule = MagicMock()
        molecule.write_xyz.side_effect = _write_xyz
        job = DummyCRESTJob(
            folder=tmp_path,
            label="myjob",
            settings=settings,
            molecule=molecule,
        )

        with patch("chemsmart.jobs.crest.writer.CRESTInputWriter"):
            exit_code = runner.run(job)

        assert exit_code == 0
        assert job.label == "myjob_fake"
        assert os.path.exists(job.outputfile)
        assert os.path.exists(job.errfile)
        assert "CREST terminated normally" in Path(job.outputfile).read_text()


class TestFakeCREST:
    def test_raises_when_xyzfile_missing(self, tmp_path):
        fake = FakeCREST(
            xyzfile=str(tmp_path / "missing.xyz"),
            outputfile=str(tmp_path / "out.out"),
            errfile=str(tmp_path / "err.err"),
            command=["crest", "missing.xyz"],
        )
        with pytest.raises(FileNotFoundError, match="not found"):
            fake.run()

    def test_writes_output_and_empty_err_file(self, tmp_path):
        xyzfile = tmp_path / "mol.xyz"
        xyzfile.write_text("1\n\nH 0 0 0\n")
        outputfile = tmp_path / "job.out"
        errfile = tmp_path / "job.err"
        fake = FakeCREST(
            xyzfile=str(xyzfile),
            outputfile=str(outputfile),
            errfile=str(errfile),
            command=["crest", "mol.xyz", "--gfn2"],
        )
        result = fake.run()
        assert result == 0
        content = outputfile.read_text()
        assert "command: crest mol.xyz --gfn2" in content
        assert "CREST terminated normally" in content
        assert errfile.read_text() == ""
