import os
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from chemsmart.cli.run import run
from chemsmart.cli.sub import sub
from chemsmart.jobs.crest.runner import CRESTJobRunner, FakeCRESTJobRunner
from chemsmart.jobs.gaussian.runner import (
    FakeGaussian,
    FakeGaussianJobRunner,
    GaussianJobRunner,
)
from chemsmart.jobs.iterate.runner import IterateJobRunner
from chemsmart.jobs.job import Job
from chemsmart.jobs.orca.runner import FakeORCAJobRunner, ORCAJobRunner
from chemsmart.jobs.runner import (
    JobRunner,
    _executable_class_for_program,
    decide_phase_transition,
    get_configured_max_submitters,
    get_submitter_worker_count,
    run_phase_jobs,
)
from chemsmart.jobs.thermochemistry.runner import ThermochemistryJobRunner
from chemsmart.jobs.xtb.runner import FakeXTBJobRunner, XTBJobRunner
from chemsmart.settings.server import Server


class ConcreteJobRunner(JobRunner):
    JOBTYPES = ["concrete_runner_type"]
    PROGRAM = "concrete"
    FAKE = False
    SCRATCH = False

    @property
    def executable(self):
        return None

    def _get_command(self, job):
        return "echo hi"

    def _create_process(self, job, command, env):
        return MagicMock()


@pytest.fixture()
def runner(pbs_server):
    return ConcreteJobRunner(server=pbs_server, scratch=False)


class DummyORCAJob:
    PROGRAM = "ORCA"
    execution_name = Job.execution_name

    def __init__(self, folder, label):
        self.folder = str(folder)
        self.label = label

    @property
    def inputfile(self):
        return str(Path(self.folder) / f"{self.label}.inp")

    @property
    def gbwfile(self):
        return str(Path(self.folder) / f"{self.label}.gbw")

    @property
    def errfile(self):
        return str(Path(self.folder) / f"{self.label}.err")

    @property
    def outputfile(self):
        return str(Path(self.folder) / f"{self.label}.out")


class DummyGaussianJob:
    PROGRAM = "Gaussian"
    execution_name = Job.execution_name

    def __init__(self, folder, label):
        self.folder = str(folder)
        self.label = label

    @property
    def inputfile(self):
        return str(Path(self.folder) / f"{self.label}.com")

    @property
    def chkfile(self):
        return str(Path(self.folder) / f"{self.label}.chk")

    @property
    def errfile(self):
        return str(Path(self.folder) / f"{self.label}.err")

    @property
    def outputfile(self):
        return str(Path(self.folder) / f"{self.label}.log")


class TestProgramAwareNames:
    @pytest.mark.parametrize(
        "program,label,expected",
        [
            ("ORCA", "molecule_opt", "orca_molecule_opt"),
            ("Gaussian", "molecule_opt", "gaussian_molecule_opt"),
            (None, "molecule_opt", "molecule_opt"),
            ("ORCA", None, None),
        ],
    )
    def test_execution_name(self, program, label, expected):
        """Derive a program prefix without mutating the calculation label.

        Missing program or label values retain the original label.
        """
        job = Job(molecule=None, label=label, jobrunner=None)
        job.PROGRAM = program

        assert job.execution_name == expected
        assert job.label == label

    @pytest.mark.parametrize(
        "runner_cls,program,input_attribute,expected_filename",
        [
            (
                GaussianJobRunner,
                "Gaussian",
                "job_inputfile",
                "molecule_opt.com",
            ),
            (
                FakeGaussianJobRunner,
                "Gaussian",
                "job_inputfile",
                "molecule_opt_fake.com",
            ),
            (ORCAJobRunner, "ORCA", "job_inputfile", "molecule_opt.inp"),
            (
                FakeORCAJobRunner,
                "ORCA",
                "job_inputfile",
                "molecule_opt_fake.inp",
            ),
            (XTBJobRunner, "xtb", "job_inputfile", "molecule_opt.inp"),
            (
                FakeXTBJobRunner,
                "xtb",
                "job_inputfile",
                "molecule_opt_fake.inp",
            ),
            (CRESTJobRunner, "crest", "job_xyzfile", "molecule_opt.xyz"),
            (
                FakeCRESTJobRunner,
                "crest",
                "job_xyzfile",
                "molecule_opt_fake.xyz",
            ),
            (
                ThermochemistryJobRunner,
                "Thermochemistry",
                "job_inputfile",
                "molecule_opt.log",
            ),
        ],
        ids=lambda value: value.__name__ if isinstance(value, type) else None,
    )
    def test_scratch_directory_uses_execution_name(
        self,
        runner_cls,
        program,
        input_attribute,
        expected_filename,
        pbs_server,
        tmp_path,
    ):
        """Prefix scratch directories while preserving calculation filenames.

        Only configure paths; no external program is executed. Fake runners
        keep their existing ``_fake`` label suffix inside the prefixed directory.
        """
        runner = runner_cls(
            server=pbs_server, scratch=True, scratch_dir=str(tmp_path)
        )
        job = Job(molecule=None, label="molecule_opt", jobrunner=runner)
        job.PROGRAM = program
        job.inputfile = str(tmp_path / "molecule_opt.log")
        job.outputfile = str(tmp_path / "molecule_opt.out")
        job.errfile = str(tmp_path / "molecule_opt.err")

        runner._set_up_variables_in_scratch(job)

        expected_dir = tmp_path / f"{program.lower()}_molecule_opt"
        assert expected_dir.is_dir()
        assert Path(runner.running_directory) == expected_dir
        assert (
            Path(getattr(runner, input_attribute))
            == expected_dir / expected_filename
        )
        assert job.label == Path(expected_filename).stem


class TestJobRunnerSelection:
    def test_fake_gaussian_runner_selected_when_fake_enabled(self, pbs_server):
        job = SimpleNamespace(TYPE="g16opt")
        runner = JobRunner.from_job(
            job=job, server=pbs_server, scratch=False, fake=True
        )
        assert isinstance(runner, FakeGaussianJobRunner)
        assert runner.fake is True

    def test_fake_orca_runner_selected_when_fake_enabled(self, pbs_server):
        job = SimpleNamespace(TYPE="orcasp")
        runner = JobRunner.from_job(
            job=job, server=pbs_server, scratch=False, fake=True
        )
        assert isinstance(runner, FakeORCAJobRunner)
        assert runner.fake is True

    def test_fake_xtb_runner_selected_when_fake_enabled(self, pbs_server):
        job = SimpleNamespace(TYPE="xtbhess")
        runner = JobRunner.from_job(
            job=job, server=pbs_server, scratch=False, fake=True
        )
        assert isinstance(runner, FakeXTBJobRunner)
        assert runner.fake is True

    def test_real_runner_selected_when_fake_disabled(self, pbs_server):
        gaussian_job = SimpleNamespace(TYPE="g16opt")
        gaussian_runner = JobRunner.from_job(
            job=gaussian_job, server=pbs_server, scratch=False, fake=False
        )
        assert isinstance(gaussian_runner, GaussianJobRunner)

        orca_job = SimpleNamespace(TYPE="orcasp")
        orca_runner = JobRunner.from_job(
            job=orca_job, server=pbs_server, scratch=False, fake=False
        )
        assert isinstance(orca_runner, ORCAJobRunner)

        xtb_job = SimpleNamespace(TYPE="xtbsp")
        xtb_runner = JobRunner.from_job(
            job=xtb_job, server=pbs_server, scratch=False, fake=False
        )
        assert isinstance(xtb_runner, XTBJobRunner)

    def test_fake_flag_propagates_when_no_fake_runner_exists(self, pbs_server):
        job = SimpleNamespace(TYPE="iterate")
        runner = JobRunner.from_job(
            job=job, server=pbs_server, scratch=False, fake=True
        )
        assert isinstance(runner, IterateJobRunner)
        assert runner.fake is True

    def test_fake_orca_appends_fake_suffix_in_job_directory(
        self, pbs_server, tmp_path
    ):
        runner = FakeORCAJobRunner(server=pbs_server, scratch=False, fake=True)
        job = DummyORCAJob(folder=tmp_path, label="orca_opt")

        runner._set_up_variables_in_job_directory(job)

        assert job.label == "orca_opt_fake"
        assert Path(runner.job_inputfile).name == "orca_opt_fake.inp"
        assert Path(runner.job_gbwfile).name == "orca_opt_fake.gbw"
        assert Path(runner.job_errfile).name == "orca_opt_fake.err"
        assert Path(runner.job_outputfile).name == "orca_opt_fake.out"

    def test_fake_orca_appends_fake_suffix_in_scratch(
        self, pbs_server, tmp_path
    ):
        runner = FakeORCAJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(tmp_path),
            fake=True,
        )
        job = DummyORCAJob(folder=tmp_path, label="orca_opt")

        runner._set_up_variables_in_scratch(job)

        assert Path(runner.running_directory) == tmp_path / "orca_orca_opt"
        assert job.label == "orca_opt_fake"
        assert Path(runner.job_inputfile).name == "orca_opt_fake.inp"
        assert Path(runner.job_gbwfile).name == "orca_opt_fake.gbw"
        assert Path(runner.job_errfile).name == "orca_opt_fake.err"
        assert Path(runner.job_outputfile).name == "orca_opt_fake.out"

    def test_fake_gaussian_does_not_duplicate_fake_suffix(
        self, pbs_server, tmp_path
    ):
        runner = FakeGaussianJobRunner(
            server=pbs_server, scratch=False, fake=True
        )
        job = DummyGaussianJob(folder=tmp_path, label="gaussian_opt_fake")

        runner._set_up_variables_in_job_directory(job)

        assert job.label == "gaussian_opt_fake"
        assert Path(runner.job_inputfile).name == "gaussian_opt_fake.com"
        assert Path(runner.job_chkfile).name == "gaussian_opt_fake.chk"
        assert Path(runner.job_errfile).name == "gaussian_opt_fake.err"

    def test_fake_orca_does_not_duplicate_fake_suffix(
        self, pbs_server, tmp_path
    ):
        runner = FakeORCAJobRunner(server=pbs_server, scratch=False, fake=True)
        job = DummyORCAJob(folder=tmp_path, label="orca_opt_fake")

        runner._set_up_variables_in_job_directory(job)

        assert job.label == "orca_opt_fake"
        assert Path(runner.job_inputfile).name == "orca_opt_fake.inp"
        assert Path(runner.job_gbwfile).name == "orca_opt_fake.gbw"
        assert Path(runner.job_errfile).name == "orca_opt_fake.err"
        assert Path(runner.job_outputfile).name == "orca_opt_fake.out"


class TestScratchCLI:
    """Scratch CLI wiring: omit vs explicit flags."""

    def test_run_omitted_scratch_leaves_none_for_from_job(
        self,
        monkeypatch,
        single_molecule_xyz_file,
        gaussian_project_config_dir,
    ):
        monkeypatch.setenv(
            "CHEMSMART_CONFIG_DIR", str(gaussian_project_config_dir)
        )
        observed = {"scratch_arg": "unset"}

        def _from_job(cls, job, server, scratch=None, fake=False, **kwargs):
            observed["scratch_arg"] = scratch
            return type("R", (), {"scratch": scratch})()

        monkeypatch.setattr(
            "chemsmart.jobs.runner.JobRunner.from_job",
            classmethod(_from_job),
        )
        monkeypatch.setattr("chemsmart.jobs.job.Job.run", lambda self: None)

        result = CliRunner().invoke(
            run,
            [
                "--fake",
                "gaussian",
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "opt",
            ],
        )
        assert result.exit_code == 0, result.output
        assert observed["scratch_arg"] is None

    def test_run_explicit_scratch_reaches_from_job(
        self,
        monkeypatch,
        single_molecule_xyz_file,
        gaussian_project_config_dir,
    ):
        monkeypatch.setenv(
            "CHEMSMART_CONFIG_DIR", str(gaussian_project_config_dir)
        )
        observed = {"scratch_arg": "unset"}

        def _from_job(cls, job, server, scratch=None, fake=False, **kwargs):
            observed["scratch_arg"] = scratch
            return type("R", (), {"scratch": scratch})()

        monkeypatch.setattr(
            "chemsmart.jobs.runner.JobRunner.from_job",
            classmethod(_from_job),
        )
        monkeypatch.setattr("chemsmart.jobs.job.Job.run", lambda self: None)

        result = CliRunner().invoke(
            run,
            [
                "--fake",
                "--scratch",
                "gaussian",
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "opt",
            ],
        )
        assert result.exit_code == 0, result.output
        assert observed["scratch_arg"] is True

    def test_sub_omitted_scratch_does_not_reconstruct_no_scratch(
        self,
        monkeypatch,
        single_molecule_xyz_file,
        gaussian_project_config_dir,
        captured,
    ):
        monkeypatch.setenv(
            "CHEMSMART_CONFIG_DIR", str(gaussian_project_config_dir)
        )
        fake_server = Server(name="dummy")
        captured["cli_args"] = None
        fake_server.submit = (
            lambda job, test=False, cli_args=None, **kw: captured.update(
                cli_args=cli_args
            )
        )
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )

        result = CliRunner().invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "gaussian",
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "opt",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "--no-scratch" not in captured["cli_args"]
        assert "--scratch" not in captured["cli_args"]

    def test_sub_time_hours_and_queue_override_server_settings(
        self,
        monkeypatch,
        single_molecule_xyz_file,
        gaussian_project_config_dir,
    ):
        monkeypatch.setenv(
            "CHEMSMART_CONFIG_DIR", str(gaussian_project_config_dir)
        )
        fake_server = Server(name="dummy")
        fake_server.submit = lambda job, test=False, cli_args=None, **kw: None
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )

        result = CliRunner().invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--time-hours",
                "48.0",
                "--queue",
                "gpu",
                "gaussian",
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "opt",
            ],
        )
        assert result.exit_code == 0, result.output
        assert fake_server.num_hours == 48.0
        assert fake_server.queue_name == "gpu"

    def test_sub_verbose_flag_enables_stream_and_debug_logging(
        self,
        monkeypatch,
        single_molecule_xyz_file,
        gaussian_project_config_dir,
    ):
        monkeypatch.setenv(
            "CHEMSMART_CONFIG_DIR", str(gaussian_project_config_dir)
        )
        fake_server = Server(name="dummy")
        fake_server.submit = lambda job, test=False, cli_args=None, **kw: None
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )
        captured = {}
        monkeypatch.setattr(
            "chemsmart.cli.sub.create_logger",
            lambda **kw: captured.update(kw),
        )

        result = CliRunner().invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--verbose",
                "gaussian",
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "opt",
            ],
        )
        assert result.exit_code == 0, result.output
        assert captured == {"stream": True, "debug": True}


def _write_server_yaml(
    path, *, gaussian_scratch, orca_scratch, xtb_scratch=None
):
    """Write a minimal server YAML with optional program SCRATCH keys."""
    gaussian_line = (
        f"    SCRATCH: {gaussian_scratch}\n"
        if gaussian_scratch is not None
        else ""
    )
    orca_line = (
        f"    SCRATCH: {orca_scratch}\n" if orca_scratch is not None else ""
    )
    xtb_line = (
        f"    SCRATCH: {xtb_scratch}\n" if xtb_scratch is not None else ""
    )
    path.write_text(
        "SERVER:\n"
        "    SCHEDULER: PBS\n"
        "    MEM_GB: 8\n"
        "    NUM_CORES: 2\n"
        "    NUM_GPUS: 0\n"
        "    NUM_THREADS: 2\n"
        "    SUBMIT_COMMAND: qsub\n"
        "    SCRATCH_DIR: null\n"
        "GAUSSIAN:\n"
        "    EXEFOLDER: ~/programs/g16\n"
        "    LOCAL_RUN: True\n"
        f"{gaussian_line}"
        "    ENVARS: |\n"
        "        export SCRATCH=~/scratch\n"
        "ORCA:\n"
        "    EXEFOLDER: ~/programs/orca\n"
        "    LOCAL_RUN: False\n"
        f"{orca_line}"
        "    ENVARS: |\n"
        "        export SCRATCH=~/scratch\n"
        "XTB:\n"
        "    EXEFOLDER: null\n"
        "    LOCAL_RUN: True\n"
        f"{xtb_line}"
        "    ENVARS: |\n"
        "        export SCRATCH=~/scratch\n"
    )
    return path


class TestScratchYamlOverride:
    """CLI → YAML program SCRATCH → runner class SCRATCH."""

    def test_program_scratch_from_servername_reads_yaml(self, tmp_path):
        from chemsmart.settings.executable import (
            GaussianExecutable,
            ORCAExecutable,
            XTBExecutable,
        )

        yaml_path = _write_server_yaml(
            tmp_path / "mixed.yaml",
            gaussian_scratch=True,
            orca_scratch=False,
            xtb_scratch=True,
        )
        assert (
            GaussianExecutable.program_scratch_from_servername(str(yaml_path))
            is True
        )
        assert (
            ORCAExecutable.program_scratch_from_servername(str(yaml_path))
            is False
        )
        assert (
            XTBExecutable.program_scratch_from_servername(str(yaml_path))
            is True
        )

    def test_program_scratch_from_servername_missing_key_is_none(
        self, tmp_path
    ):
        from chemsmart.settings.executable import ORCAExecutable

        yaml_path = _write_server_yaml(
            tmp_path / "no_orca_scratch.yaml",
            gaussian_scratch=True,
            orca_scratch=None,
        )
        assert (
            ORCAExecutable.program_scratch_from_servername(str(yaml_path))
            is None
        )

    def test_omit_uses_yaml_false_over_orca_class_true(self, tmp_path):
        yaml_path = _write_server_yaml(
            tmp_path / "orca_off.yaml",
            gaussian_scratch=True,
            orca_scratch=False,
        )
        server = Server.from_yaml(str(yaml_path))
        runner = JobRunner.from_job(
            job=SimpleNamespace(TYPE="orcasp"),
            server=server,
            scratch=None,
            fake=True,
        )
        assert isinstance(runner, FakeORCAJobRunner)
        assert ORCAJobRunner.SCRATCH is True
        assert runner.scratch is False

    def test_omit_uses_yaml_true_for_gaussian(self, tmp_path):
        yaml_path = _write_server_yaml(
            tmp_path / "g16_on.yaml",
            gaussian_scratch=True,
            orca_scratch=False,
        )
        server = Server.from_yaml(str(yaml_path))
        scratch_dir = tmp_path / "scratch"
        scratch_dir.mkdir()
        runner = JobRunner.from_job(
            job=SimpleNamespace(TYPE="g16opt"),
            server=server,
            scratch=None,
            fake=True,
            scratch_dir=str(scratch_dir),
        )
        assert isinstance(runner, FakeGaussianJobRunner)
        assert runner.scratch is True

    def test_omit_falls_back_to_class_when_yaml_key_absent(self, tmp_path):
        yaml_path = _write_server_yaml(
            tmp_path / "orca_absent.yaml",
            gaussian_scratch=True,
            orca_scratch=None,
        )
        server = Server.from_yaml(str(yaml_path))
        scratch_dir = tmp_path / "scratch"
        scratch_dir.mkdir()
        runner = JobRunner.from_job(
            job=SimpleNamespace(TYPE="orcasp"),
            server=server,
            scratch=None,
            fake=True,
            scratch_dir=str(scratch_dir),
        )
        assert runner.scratch is ORCAJobRunner.SCRATCH

    def test_cli_scratch_true_overrides_yaml_false(self, tmp_path):
        yaml_path = _write_server_yaml(
            tmp_path / "orca_off.yaml",
            gaussian_scratch=True,
            orca_scratch=False,
        )
        server = Server.from_yaml(str(yaml_path))
        scratch_dir = tmp_path / "scratch"
        scratch_dir.mkdir()
        runner = JobRunner.from_job(
            job=SimpleNamespace(TYPE="orcasp"),
            server=server,
            scratch=True,
            fake=True,
            scratch_dir=str(scratch_dir),
        )
        assert runner.scratch is True

    def test_cli_scratch_false_overrides_yaml_true(self, tmp_path):
        yaml_path = _write_server_yaml(
            tmp_path / "g16_on.yaml",
            gaussian_scratch=True,
            orca_scratch=False,
        )
        server = Server.from_yaml(str(yaml_path))
        runner = JobRunner.from_job(
            job=SimpleNamespace(TYPE="g16opt"),
            server=server,
            scratch=False,
            fake=True,
        )
        assert runner.scratch is False

    def test_missing_scratch_dir_falls_back_to_job_folder(
        self, pbs_server, tmp_path, monkeypatch
    ):
        runner = FakeGaussianJobRunner(
            server=pbs_server, scratch=False, fake=True
        )
        runner.scratch = True
        runner._scratch_dir = None
        runner._set_scratch.cache_clear()
        fake_exe = SimpleNamespace(
            scratch_dir=str(tmp_path / "does_not_exist"),
            local_run=None,
        )
        monkeypatch.setattr(
            type(runner),
            "executable",
            property(lambda self: fake_exe),
        )
        runner._set_scratch()
        job = DummyGaussianJob(folder=tmp_path, label="gaussian_opt")
        runner._assign_variables(job)

        assert runner.scratch is False
        assert runner.running_directory == job.folder


class TestGaussianJobRunnerDirectHelpers:
    """Direct unit coverage for GaussianJobRunner's file-management and
    execution-plumbing helpers, which are otherwise only exercised
    indirectly (and incompletely) through full job-run integration
    tests elsewhere."""

    def test_init_defaults_scratch_to_class_attribute(self, pbs_server):
        runner = GaussianJobRunner(server=pbs_server, scratch=None, fake=True)
        assert runner.scratch is GaussianJobRunner.SCRATCH

    def test_executable_property_wraps_filenotfounderror(
        self, pbs_server, monkeypatch
    ):
        runner = GaussianJobRunner(server=pbs_server, scratch=False, fake=True)
        runner.server = SimpleNamespace(name="no-such-server")

        from chemsmart.jobs.gaussian import runner as runner_module

        def _raise(servername):
            raise FileNotFoundError(f"no config for {servername}")

        monkeypatch.setattr(
            runner_module.GaussianExecutable,
            "from_servername",
            staticmethod(_raise),
        )
        with pytest.raises(FileNotFoundError):
            runner.executable

    def test_assign_variables_dispatches_to_scratch_setup(
        self, pbs_server, tmp_path, monkeypatch
    ):
        runner = GaussianJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(tmp_path),
            fake=True,
        )
        calls = []
        monkeypatch.setattr(
            runner,
            "_set_up_variables_in_scratch",
            lambda job: calls.append("scratch"),
        )
        monkeypatch.setattr(
            runner,
            "_set_up_variables_in_job_directory",
            lambda job: calls.append("job_directory"),
        )
        monkeypatch.setattr(
            type(runner),
            "executable",
            property(lambda self: SimpleNamespace(local_run=None)),
        )
        job = DummyGaussianJob(folder=tmp_path, label="assign_test")
        runner._assign_variables(job)
        assert calls == ["scratch"]
        assert not hasattr(job, "local")

    def test_assign_variables_dispatches_to_job_directory_setup(
        self, pbs_server, tmp_path, monkeypatch
    ):
        runner = GaussianJobRunner(server=pbs_server, scratch=False, fake=True)
        calls = []
        monkeypatch.setattr(
            runner,
            "_set_up_variables_in_scratch",
            lambda job: calls.append("scratch"),
        )
        monkeypatch.setattr(
            runner,
            "_set_up_variables_in_job_directory",
            lambda job: calls.append("job_directory"),
        )
        monkeypatch.setattr(
            type(runner),
            "executable",
            property(lambda self: SimpleNamespace(local_run=True)),
        )
        job = DummyGaussianJob(folder=tmp_path, label="assign_test")
        runner._assign_variables(job)
        assert calls == ["job_directory"]
        assert job.local is True

    def test_set_up_variables_in_scratch_creates_fresh_directory(
        self, pbs_server, tmp_path
    ):
        runner = GaussianJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(tmp_path / "scratch_root"),
            fake=True,
        )
        job = DummyGaussianJob(
            folder=tmp_path / "jobfolder", label="mytest_job"
        )
        runner._set_up_variables_in_scratch(job)

        import os

        expected_dir = os.path.join(
            str(tmp_path / "scratch_root"), job.execution_name
        )
        assert os.path.isdir(expected_dir)
        assert runner.running_directory == expected_dir
        assert runner.job_inputfile == os.path.abspath(
            os.path.join(expected_dir, f"{job.label}.com")
        )
        assert runner.job_chkfile == os.path.abspath(
            os.path.join(expected_dir, f"{job.label}.chk")
        )
        assert runner.job_errfile == os.path.abspath(
            os.path.join(expected_dir, f"{job.label}.err")
        )
        assert runner.job_outputfile == os.path.abspath(
            os.path.join(expected_dir, f"{job.label}.log")
        )

    def test_set_up_variables_in_scratch_reuses_existing_directory(
        self, pbs_server, tmp_path
    ):
        import os

        scratch_root = tmp_path / "scratch_root2"
        job = DummyGaussianJob(
            folder=tmp_path / "jobfolder2", label="existing_job"
        )
        os.makedirs(os.path.join(str(scratch_root), job.execution_name))

        runner = GaussianJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(scratch_root),
            fake=True,
        )
        runner._set_up_variables_in_scratch(job)
        assert runner.running_directory == os.path.join(
            str(scratch_root), job.execution_name
        )

    def test_set_up_variables_in_job_directory(self, pbs_server, tmp_path):
        import os

        runner = GaussianJobRunner(server=pbs_server, scratch=False, fake=True)
        job = DummyGaussianJob(folder=tmp_path / "jobdir", label="direct_job")
        runner._set_up_variables_in_job_directory(job)

        assert runner.running_directory == job.folder
        assert runner.job_inputfile == os.path.abspath(job.inputfile)
        assert runner.job_chkfile == os.path.abspath(job.chkfile)
        assert runner.job_errfile == os.path.abspath(job.errfile)
        assert runner.job_outputfile == os.path.abspath(job.outputfile)

    def test_prerun_delegates_to_assign_variables(
        self, pbs_server, tmp_path, monkeypatch
    ):
        runner = GaussianJobRunner(server=pbs_server, scratch=False, fake=True)
        calls = []
        monkeypatch.setattr(
            runner, "_assign_variables", lambda job: calls.append(job)
        )
        job = DummyGaussianJob(folder=tmp_path, label="prerun_job")
        runner._prerun(job)
        assert calls == [job]

    def test_write_input_uses_gaussian_input_writer(
        self, pbs_server, tmp_path, monkeypatch
    ):
        runner = GaussianJobRunner(server=pbs_server, scratch=False, fake=True)
        runner.running_directory = str(tmp_path)
        job = DummyGaussianJob(folder=tmp_path, label="write_input_job")

        captured = {}

        class FakeInputWriter:
            def __init__(self, job):
                captured["job"] = job

            def write(self, target_directory):
                captured["target_directory"] = target_directory

        monkeypatch.setattr(
            "chemsmart.jobs.gaussian.writer.GaussianInputWriter",
            FakeInputWriter,
        )
        runner._write_input(job)
        assert captured["job"] is job
        assert captured["target_directory"] == str(tmp_path)

    def test_get_command_combines_executable_and_inputfile(
        self, pbs_server, monkeypatch
    ):
        runner = GaussianJobRunner(server=pbs_server, scratch=False, fake=True)
        monkeypatch.setattr(runner, "_get_executable", lambda: "/opt/g16")
        runner.job_inputfile = "/scratch/job/mytest.com"
        assert (
            runner._get_command(job=None) == "/opt/g16 /scratch/job/mytest.com"
        )

    def test_get_executable_delegates_to_executable_object(
        self, pbs_server, monkeypatch
    ):
        runner = GaussianJobRunner(server=pbs_server, scratch=False, fake=True)
        monkeypatch.setattr(
            type(runner),
            "executable",
            property(
                lambda self: SimpleNamespace(
                    get_executable=lambda: "/usr/bin/g16"
                )
            ),
        )
        assert runner._get_executable() == "/usr/bin/g16"

    def test_create_process_opens_files_and_spawns_popen(
        self, pbs_server, tmp_path, monkeypatch
    ):
        runner = GaussianJobRunner(server=pbs_server, scratch=False, fake=True)
        runner.running_directory = str(tmp_path)
        runner.job_outputfile = str(tmp_path / "job.log")
        runner.job_errfile = str(tmp_path / "job.err")
        monkeypatch.setattr(
            type(runner),
            "executable",
            property(lambda self: SimpleNamespace(env={"PATH": "/bin"})),
        )

        captured = {}

        class FakeProcess:
            pass

        def fake_popen(args, stdout, stderr, env, cwd):
            captured["args"] = args
            captured["env"] = env
            captured["cwd"] = cwd
            captured["stdout_closed"] = stdout.closed
            captured["stderr_closed"] = stderr.closed
            return FakeProcess()

        from chemsmart.jobs.gaussian import runner as runner_module

        monkeypatch.setattr(runner_module.subprocess, "Popen", fake_popen)

        result = runner._create_process(
            job=None, command="g16 job.com", env={"PATH": "/bin"}
        )

        assert isinstance(result, FakeProcess)
        assert captured["args"] == ["g16", "job.com"]
        assert captured["cwd"] == str(tmp_path)
        assert (tmp_path / "job.log").exists()
        assert (tmp_path / "job.err").exists()

    def test_postrun_noop_when_scratch_disabled(self, pbs_server, tmp_path):
        runner = GaussianJobRunner(server=pbs_server, scratch=False, fake=True)
        job = DummyGaussianJob(folder=tmp_path, label="noop_job")
        runner._postrun(job)  # should not raise, nothing to copy

    def test_postrun_copies_matching_files_to_job_folder(
        self, pbs_server, tmp_path
    ):
        import os

        running_dir = tmp_path / "scratch_run"
        running_dir.mkdir()
        job_folder = tmp_path / "job_folder"
        job_folder.mkdir()

        (running_dir / "copy_job.log").write_text("fake log contents")

        runner = GaussianJobRunner(server=pbs_server, scratch=True, fake=True)
        runner.running_directory = str(running_dir)
        job = DummyGaussianJob(folder=job_folder, label="copy_job")

        runner._postrun(job)

        assert (job_folder / "copy_job.log").exists()
        assert (job_folder / "copy_job.log").read_text() == (
            "fake log contents"
        )
        assert os.path.exists(str(running_dir / "copy_job.log"))

    def test_postrun_skips_files_whose_full_path_starts_with_gau_prefix(
        self, pbs_server, tmp_path, monkeypatch
    ):
        """Defensive branch: the "except files starting with Gau-" check
        tests the full glob'd path (not the basename), so it only ever
        skips a file when the running_directory string itself begins
        with "Gau-" (an edge case, not Gaussian's own scratch files,
        which never match the job.label* glob pattern in the first
        place)."""
        monkeypatch.chdir(tmp_path)
        import os

        os.makedirs("Gau-scratch")
        with open(os.path.join("Gau-scratch", "prefixed_job.log"), "w") as f:
            f.write("should not be copied")
        job_folder = tmp_path / "job_folder"
        job_folder.mkdir()

        runner = GaussianJobRunner(server=pbs_server, scratch=True, fake=True)
        runner.running_directory = "Gau-scratch"
        job = DummyGaussianJob(folder=job_folder, label="prefixed_job")

        runner._postrun(job)

        assert not (job_folder / "prefixed_job.log").exists()

    def test_postrun_logs_and_continues_on_copy_failure(
        self, pbs_server, tmp_path, monkeypatch
    ):
        running_dir = tmp_path / "scratch_run2"
        running_dir.mkdir()
        job_folder = tmp_path / "job_folder2"
        job_folder.mkdir()
        (running_dir / "fail_job.log").write_text("contents")

        from chemsmart.jobs.gaussian import runner as runner_module

        def _raise_copy(src, dst):
            raise OSError("disk full")

        monkeypatch.setattr(runner_module, "copy", _raise_copy)

        runner = GaussianJobRunner(server=pbs_server, scratch=True, fake=True)
        runner.running_directory = str(running_dir)
        job = DummyGaussianJob(folder=job_folder, label="fail_job")

        runner._postrun(job)  # should not raise; failure is logged

        assert not (job_folder / "fail_job.log").exists()


class TestFakeGaussianJobRunnerDirectHelpers:
    def test_set_up_variables_in_scratch_reuses_existing_directory(
        self, pbs_server, tmp_path
    ):
        import os

        scratch_root = tmp_path / "fake_scratch_root"
        job = DummyGaussianJob(
            folder=tmp_path / "jobfolder", label="fake_existing_job"
        )
        os.makedirs(os.path.join(str(scratch_root), job.execution_name))

        expected_dir = os.path.join(str(scratch_root), job.execution_name)
        runner = FakeGaussianJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(scratch_root),
            fake=True,
        )
        # _set_up_variables_in_scratch appends "_fake" to job.label as a
        # side effect, so the expected directory must be captured first.
        runner._set_up_variables_in_scratch(job)
        assert runner.running_directory == expected_dir

    def test_set_up_variables_in_scratch_creates_fresh_directory(
        self, pbs_server, tmp_path
    ):
        import os

        scratch_root = tmp_path / "fake_scratch_root_fresh"
        job = DummyGaussianJob(
            folder=tmp_path / "jobfolder_fresh", label="fake_fresh_job"
        )
        expected_dir = os.path.join(str(scratch_root), job.execution_name)

        runner = FakeGaussianJobRunner(
            server=pbs_server,
            scratch=True,
            scratch_dir=str(scratch_root),
            fake=True,
        )
        runner._set_up_variables_in_scratch(job)
        assert os.path.isdir(expected_dir)
        assert runner.running_directory == expected_dir

    def test_run_executes_full_fake_pipeline_in_order(
        self, pbs_server, tmp_path, monkeypatch
    ):
        runner = FakeGaussianJobRunner(
            server=pbs_server, scratch=False, fake=True
        )
        calls = []
        monkeypatch.setattr(
            runner, "_prerun", lambda job: calls.append("prerun")
        )
        monkeypatch.setattr(
            runner, "_write_input", lambda job: calls.append("write_input")
        )
        monkeypatch.setattr(
            runner, "_postrun", lambda job: calls.append("postrun")
        )
        monkeypatch.setattr(
            runner,
            "_postrun_cleanup",
            lambda job: calls.append("postrun_cleanup"),
        )
        runner.job_inputfile = "unused.com"

        class FakeFakeGaussian:
            def __init__(self, file_to_run):
                calls.append(f"init:{file_to_run}")

            def run(self):
                calls.append("fake_run")
                return 0

        from chemsmart.jobs.gaussian import runner as runner_module

        monkeypatch.setattr(runner_module, "FakeGaussian", FakeFakeGaussian)

        job = DummyGaussianJob(folder=tmp_path, label="pipeline_job")
        result = runner.run(job)

        assert result == 0
        assert calls == [
            "prerun",
            "write_input",
            "init:unused.com",
            "fake_run",
            "postrun",
            "postrun_cleanup",
        ]


def _write_gaussian_com(
    tmp_path, filename="test.com", route="# opt freq b3lyp/6-31g(d)", mult=1
):
    lines = [
        "%chk=test.chk",
        "%mem=4GB",
        route,
        "",
        "title",
        "",
        f"0 {mult}",
        "C 0.0 0.0 0.0",
        "H 0.0 0.0 1.0",
        "",
    ]
    path = tmp_path / filename
    path.write_text("\n".join(lines))
    return str(path)


class TestFakeGaussianDirect:
    """Direct unit coverage for the FakeGaussian simulator (used by
    FakeGaussianJobRunner.run) that isn't reached by the higher-level
    fake-job-run tests elsewhere: the missing-file error, odd-multiplicity
    spin, extra input blocks, and the semiempirical/ab-initio/unknown-method
    branches of the fake SCF-energy line."""

    def test_missing_input_file_raises_filenotfounderror(self):
        with pytest.raises(FileNotFoundError):
            FakeGaussian("/no/such/path/does_not_exist.com")

    def test_spin_is_u_for_non_singlet_multiplicity(self, tmp_path):
        path = _write_gaussian_com(tmp_path, mult=2)
        fake = FakeGaussian(path)
        assert fake.spin == "U"

    def test_run_writes_additional_input_blocks(self, tmp_path):
        path = _write_gaussian_com(tmp_path, route="# opt freq b3lyp/gen")
        # append a gen/genecp-style basis-set block (content_groups[3:]),
        # separated from the coordinate block by a blank line.
        with open(path, "a") as f:
            f.write("\nC 0\n6-31g(d)\n****\n\n")

        fake = FakeGaussian(path)
        assert len(fake.input_blocks) > 3
        fake.run()
        output = Path(fake.output_filepath).read_text()
        assert "6-31g(d)" in output

    def test_run_uses_semiempirical_energy_line(self, tmp_path):
        path = _write_gaussian_com(tmp_path, route="# opt freq pm6")
        fake = FakeGaussian(path)
        fake.run()
        output = Path(fake.output_filepath).read_text()
        assert "E(PM6)" in output

    def test_run_uses_ab_initio_energy_line(self, tmp_path):
        # A bare ab initio method with no "/basis" suffix is required so
        # that .functional stays None (any "method/basis" token is always
        # parsed as a functional first, regardless of the method).
        path = _write_gaussian_com(tmp_path, route="# opt freq hf")
        fake = FakeGaussian(path)
        assert fake.input_object.functional is None
        assert fake.input_object.ab_initio == "hf"
        fake.run()
        output = Path(fake.output_filepath).read_text()
        assert "E(hf)" in output

    def test_run_uses_unknown_method_energy_line(self, tmp_path):
        path = _write_gaussian_com(tmp_path, route="# opt freq uff")
        fake = FakeGaussian(path)
        assert fake.input_object.functional is None
        assert fake.input_object.semiempirical is None
        assert fake.input_object.ab_initio is None
        fake.run()
        output = Path(fake.output_filepath).read_text()
        assert "E(Unknow Method)" in output


class TestSubResourceOverrides:
    def test_cli_resources_reach_submission_server(
        self,
        monkeypatch,
        single_molecule_xyz_file,
        gaussian_project_config_dir,
        captured,
    ):
        monkeypatch.setenv(
            "CHEMSMART_CONFIG_DIR", str(gaussian_project_config_dir)
        )
        fake_server = Server(
            name="dummy",
            NUM_CORES=2,
            NUM_GPUS=0,
            MEM_GB=8,
            NUM_HOURS=1,
            QUEUE_NAME="normal",
        )

        def _capture_script(self, job, cli_args, **kwargs):
            captured.update(
                server=self,
                num_cores=self.num_cores,
                num_gpus=self.num_gpus,
                mem_gb=self.mem_gb,
                num_hours=self.num_hours,
                queue_name=self.queue_name,
            )

        monkeypatch.setattr(
            Server,
            "_check_running_jobs",
            lambda self, job: None,
        )
        monkeypatch.setattr(
            Server,
            "_write_submission_script",
            _capture_script,
        )
        monkeypatch.setattr(
            Server,
            "current",
            classmethod(lambda cls: fake_server),
        )

        result = CliRunner().invoke(
            sub,
            [
                "--test",
                "--num-cores",
                "7",
                "--num-gpus",
                "2",
                "--mem-gb",
                "23",
                "--time-hours",
                "12.5",
                "--queue",
                "debug",
                "gaussian",
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "opt",
            ],
        )

        assert result.exit_code == 0, result.output
        assert captured["server"] is fake_server
        assert captured["num_cores"] == 7
        assert captured["num_gpus"] == 2
        assert captured["mem_gb"] == 23
        assert captured["num_hours"] == 12.5
        assert captured["queue_name"] == "debug"


class TestJobRunnerConstruction:
    def test_rejects_invalid_server_type(self):
        with pytest.raises(ValueError, match="server must be instance"):
            ConcreteJobRunner(server=12345, scratch=False)

    def test_accepts_server_instance(self, pbs_server):
        runner = ConcreteJobRunner(server=pbs_server, scratch=False)
        assert runner.server is pbs_server

    def test_num_cores_defaults_from_server(self, pbs_server):
        runner = ConcreteJobRunner(server=pbs_server, scratch=False)
        assert runner.num_cores == pbs_server.num_cores

    def test_num_cores_explicit_override(self, pbs_server):
        runner = ConcreteJobRunner(
            server=pbs_server, scratch=False, num_cores=99
        )
        assert runner.num_cores == 99

    def test_repr_contains_class_and_server(self, runner):
        text = repr(runner)
        assert "ConcreteJobRunner" in text

    def test_server_as_string_resolves_via_from_servername(self):
        fake_server = MagicMock(spec=Server)
        fake_server.num_cores = 8
        fake_server.num_gpus = 0
        fake_server.mem_gb = 32
        with patch.object(
            Server, "from_servername", return_value=fake_server
        ) as mock_from_servername:
            runner = ConcreteJobRunner(server="myserver", scratch=False)
        mock_from_servername.assert_called_once_with("myserver")
        assert runner.server is fake_server

    def test_num_gpus_explicit_override(self, pbs_server):
        runner = ConcreteJobRunner(
            server=pbs_server, scratch=False, num_gpus=3
        )
        assert runner.num_gpus == 3

    def test_mem_gb_explicit_override(self, pbs_server):
        runner = ConcreteJobRunner(
            server=pbs_server, scratch=False, mem_gb=256
        )
        assert runner.mem_gb == 256

    def test_servername_property_delegates_to_server(self, runner):
        assert runner.servername == runner.server.name


class TestJobRunnerScratchDirSetter:
    def test_setter_rejects_nonexistent_path(self, runner):
        with pytest.raises(FileNotFoundError):
            runner.scratch_dir = "/path/does/not/exist/anywhere"

    def test_setter_accepts_existing_path(self, runner, tmp_path):
        runner.scratch_dir = str(tmp_path)
        assert runner._scratch_dir == str(tmp_path)

    def test_setter_expands_user_home(self, runner):
        # Should not raise for a bare "~" (always exists) and should be
        # expanded to an absolute path.
        runner.scratch_dir = "~"
        assert runner._scratch_dir == os.path.expanduser("~")

    def test_setter_accepts_none_without_validation(self, runner):
        """Covers the setter's `if value is not None:` False arm: a
        None value skips expanduser/exists-check entirely."""
        runner._scratch_dir = "/some/previous/value"
        runner.scratch_dir = None
        assert runner._scratch_dir is None

    def test_getter_computes_and_caches_when_unset(self, pbs_server, tmp_path):
        """Covers the getter's `if self._scratch_dir is None:` branch:
        the first access to .scratch_dir must resolve and cache it via
        _set_scratch(), rather than staying None forever."""
        with patch.object(type(pbs_server), "scratch_dir", new=str(tmp_path)):
            runner = ConcreteJobRunner(server=pbs_server, scratch=True)
            assert runner._scratch_dir is None
            resolved = runner.scratch_dir
        assert resolved == str(tmp_path)
        assert runner._scratch_dir == str(tmp_path)


class TestJobRunnerSetScratch:
    def test_uses_explicit_scratch_dir_when_set(self, pbs_server, tmp_path):
        runner = ConcreteJobRunner(
            server=pbs_server, scratch=True, scratch_dir=str(tmp_path)
        )
        assert runner._scratch_dir == str(tmp_path)

    def test_disables_scratch_when_no_source_available(self, pbs_server):
        with (
            patch.object(type(pbs_server), "scratch_dir", new=None),
            patch("chemsmart.jobs.runner.user_settings") as mock_user_settings,
        ):
            mock_user_settings.scratch = None
            runner = ConcreteJobRunner(server=pbs_server, scratch=True)
        assert runner.scratch is False

    def test_disables_scratch_when_resolved_directory_missing(
        self, pbs_server
    ):
        with patch.object(
            type(pbs_server), "scratch_dir", new="/does/not/exist/at/all"
        ):
            runner = ConcreteJobRunner(server=pbs_server, scratch=True)
        assert runner.scratch is False


class TestJobRunnerRunOrchestration:
    def test_run_calls_hooks_in_order(self, runner):
        job = MagicMock()
        calls = []
        runner._prerun = MagicMock(
            side_effect=lambda j: calls.append("prerun")
        )
        runner._write_input = MagicMock(
            side_effect=lambda j: calls.append("write_input")
        )
        runner._get_command = MagicMock(
            side_effect=lambda j: calls.append("get_command") or "cmd"
        )
        runner._update_os_environ = MagicMock(
            side_effect=lambda j: calls.append("update_env") or {}
        )
        runner._create_process = MagicMock(
            side_effect=lambda j, command, env: calls.append("create_process")
            or MagicMock()
        )
        runner._run = MagicMock(
            side_effect=lambda p, **kw: calls.append("run")
        )
        runner._postrun = MagicMock(
            side_effect=lambda j: calls.append("postrun")
        )
        runner._postrun_cleanup = MagicMock(
            side_effect=lambda j: calls.append("postrun_cleanup")
        )

        runner.run(job)

        assert calls == [
            "prerun",
            "write_input",
            "get_command",
            "update_env",
            "create_process",
            "run",
            "postrun",
            "postrun_cleanup",
        ]

    def test_default_run_communicates_and_polls_process(self, runner):
        process = MagicMock()
        process.poll.return_value = 0
        result = runner._run(process)
        process.communicate.assert_called_once()
        process.poll.assert_called_once()
        assert result == 0

    def test_prerun_write_input_postrun_are_no_ops_by_default(self, runner):
        """ConcreteJobRunner doesn't override these optional hooks, so
        the base class's own pass-through bodies run for real."""
        job = MagicMock()
        assert runner._prerun(job) is None
        assert runner._write_input(job) is None
        assert runner._postrun(job) is None


class TestUpdateOsEnviron:
    def test_no_executable_returns_plain_environ_copy(self, runner):
        job = MagicMock()
        env = runner._update_os_environ(job)
        assert env == dict(os.environ)

    def test_empty_executable_env_returns_plain_environ_copy(self, runner):
        job = MagicMock()
        with patch.object(
            type(runner),
            "executable",
            new=property(lambda self: MagicMock(env={})),
        ):
            env = runner._update_os_environ(job)
        assert env == dict(os.environ)

    def test_executable_env_vars_are_applied_with_user_expansion(self, runner):
        job = MagicMock()
        with patch.object(
            type(runner),
            "executable",
            new=property(
                lambda self: MagicMock(
                    env={"MY_STR_VAR": "~/mydir", "MY_INT_VAR": 5}
                )
            ),
        ):
            env = runner._update_os_environ(job)
        assert env["MY_STR_VAR"] == os.path.expanduser("~/mydir")
        assert env["MY_INT_VAR"] == "5"


class TestJobRunnerCopy:
    def test_copy_returns_shallow_copy(self, runner):
        copied = runner.copy()
        assert copied is not runner
        assert copied.server is runner.server


class TestJobRunnerFromJob:
    def test_from_job_dispatches_to_registered_runner(self, pbs_server):
        from chemsmart.jobs.thermochemistry.job import ThermochemistryJob
        from chemsmart.jobs.thermochemistry.runner import (
            ThermochemistryJobRunner,
        )

        job = MagicMock(spec=ThermochemistryJob)
        job.TYPE = "thermochemistry"

        result = JobRunner.from_job(job, server=pbs_server, fake=True)
        assert isinstance(result, ThermochemistryJobRunner)

    def test_from_job_raises_for_unknown_jobtype(self, pbs_server):
        job = MagicMock()
        job.TYPE = "totally_unregistered_jobtype_xyz"
        with pytest.raises(ValueError, match="Could not find any runners"):
            JobRunner.from_job(job, server=pbs_server)

    def test_from_job_treats_notimplemented_jobtypes_as_empty_list(
        self, pbs_server
    ):
        """A registered runner subclass that never overrides JOBTYPES
        (still NotImplemented) must be treated as supporting no
        jobtypes, not crash the `jobtype in runner_jobtypes` check."""

        class _RunnerWithoutJobtypes(JobRunner):
            PROGRAM = "concrete_no_jobtypes"
            FAKE = False
            SCRATCH = False

            @property
            def executable(self):
                return None

            def _get_command(self, job):
                return "echo hi"

            def _create_process(self, job, command, env):
                return MagicMock()

        job = MagicMock()
        job.TYPE = "concrete_runner_type"
        result = JobRunner.from_job(job, server=pbs_server, fake=True)
        assert isinstance(result, ConcreteJobRunner)


class TestJobRunnerErrFileCleanup:
    def test_remove_err_files_removes_existing_suffixes(
        self, runner, tmp_path
    ):
        job = MagicMock()
        job.folder = str(tmp_path)
        job.label = "myjob"
        for suffix in [".err", ".pbserr", ".slurmerr"]:
            (tmp_path / f"myjob{suffix}").write_text("x")

        runner._remove_err_files(job)

        for suffix in [".err", ".pbserr", ".slurmerr"]:
            assert not (tmp_path / f"myjob{suffix}").exists()

    def test_remove_err_files_no_op_when_missing(self, runner, tmp_path):
        job = MagicMock()
        job.folder = str(tmp_path)
        job.label = "myjob"
        # Should not raise even though no err files exist.
        runner._remove_err_files(job)

    def test_append_suffix_to_job_label_appends_once(self, runner):
        job = MagicMock()
        job.label = "myjob"
        runner._append_suffix_to_job_label(job, "_fake")
        assert job.label == "myjob_fake"
        runner._append_suffix_to_job_label(job, "_fake")
        assert job.label == "myjob_fake"

    def test_append_suffix_no_op_for_falsy_suffix(self, runner):
        job = MagicMock()
        job.label = "myjob"
        runner._append_suffix_to_job_label(job, "")
        assert job.label == "myjob"


class TestJobRunnerPostrunCleanup:
    def test_cleanup_removes_err_files_when_complete(self, runner):
        job = MagicMock()
        job.is_complete.return_value = True
        runner._remove_err_files = MagicMock()
        runner._delete_scratch_directory = MagicMock()
        runner.scratch = False

        runner._postrun_cleanup(job)

        runner._remove_err_files.assert_called_once_with(job)
        runner._delete_scratch_directory.assert_not_called()

    def test_cleanup_no_op_when_incomplete(self, runner):
        job = MagicMock()
        job.is_complete.return_value = False
        runner._remove_err_files = MagicMock()

        runner._postrun_cleanup(job)

        runner._remove_err_files.assert_not_called()

    def test_cleanup_deletes_scratch_when_enabled(self, runner):
        job = MagicMock()
        job.is_complete.return_value = True
        runner._remove_err_files = MagicMock()
        runner._delete_scratch_directory = MagicMock()
        runner.scratch = True
        runner.delete_scratch = True

        runner._postrun_cleanup(job)

        runner._delete_scratch_directory.assert_called_once()


class TestJobRunnerDeleteScratchDirectory:
    def test_no_op_without_running_directory_attr(self, runner):
        # runner has no `running_directory` set at all.
        runner._delete_scratch_directory()  # should not raise

    def test_refuses_when_running_directory_equals_scratch_root(
        self, runner, tmp_path
    ):
        runner.running_directory = str(tmp_path)
        runner._scratch_dir = str(tmp_path)
        with patch("chemsmart.jobs.runner.rmtree") as mock_rmtree:
            runner._delete_scratch_directory()
        mock_rmtree.assert_not_called()

    def test_refuses_when_running_directory_outside_scratch(
        self, runner, tmp_path
    ):
        scratch_root = tmp_path / "scratch"
        scratch_root.mkdir()
        outside_dir = tmp_path / "outside"
        outside_dir.mkdir()

        runner.running_directory = str(outside_dir)
        runner._scratch_dir = str(scratch_root)
        with patch("chemsmart.jobs.runner.rmtree") as mock_rmtree:
            runner._delete_scratch_directory()
        mock_rmtree.assert_not_called()
        assert outside_dir.exists()

    def test_deletes_when_running_directory_inside_scratch(
        self, runner, tmp_path
    ):
        scratch_root = tmp_path / "scratch"
        job_dir = scratch_root / "myjob"
        job_dir.mkdir(parents=True)

        runner.running_directory = str(job_dir)
        runner._scratch_dir = str(scratch_root)

        runner._delete_scratch_directory()

        assert not job_dir.exists()
        assert scratch_root.exists()

    def test_logs_error_when_scratch_dir_missing(
        self, runner, tmp_path, caplog
    ):
        """Covers the `if not sd.exists() or not sd.is_dir():` branch:
        running_directory exists, but the resolved scratch_dir itself
        doesn't -- refuse to proceed instead of raising."""
        running_dir = tmp_path / "running"
        running_dir.mkdir()

        runner.running_directory = str(running_dir)
        runner._scratch_dir = str(tmp_path / "does_not_exist_scratch")

        with (
            patch("chemsmart.jobs.runner.rmtree") as mock_rmtree,
            caplog.at_level("ERROR"),
        ):
            runner._delete_scratch_directory()

        mock_rmtree.assert_not_called()
        assert "doesn't exist or is not a directory" in caplog.text

    def test_rmtree_failure_is_logged_not_raised(self, runner, tmp_path):
        scratch_root = tmp_path / "scratch"
        job_dir = scratch_root / "myjob"
        job_dir.mkdir(parents=True)

        runner.running_directory = str(job_dir)
        runner._scratch_dir = str(scratch_root)

        with patch(
            "chemsmart.jobs.runner.rmtree",
            side_effect=OSError("simulated rmtree failure"),
        ):
            runner._delete_scratch_directory()  # should not raise


class TestModuleLevelHelpers:
    def test_get_configured_max_submitters_from_env(self, monkeypatch):
        monkeypatch.setenv("CHEMSMART_MAX_SUBMITTERS", "7")
        assert get_configured_max_submitters(jobrunner=None) == 7

    def test_get_configured_max_submitters_from_runner_num_cores(
        self, monkeypatch
    ):
        monkeypatch.delenv("CHEMSMART_MAX_SUBMITTERS", raising=False)
        runner = MagicMock(spec=[])
        runner.num_cores = 4
        assert get_configured_max_submitters(jobrunner=runner) == 4

    def test_get_configured_max_submitters_zero_env_falls_through(
        self, monkeypatch
    ):
        """A non-positive env value is treated as unset (covers
        _positive_int_or_none's `parsed <= 0` branch), falling through
        to the next resolution step instead of returning 0."""
        monkeypatch.setenv("CHEMSMART_MAX_SUBMITTERS", "0")
        runner = MagicMock(spec=[])
        runner.num_cores = 5
        assert get_configured_max_submitters(jobrunner=runner) == 5

    def test_get_configured_max_submitters_from_runner_max_submitters(
        self, monkeypatch
    ):
        monkeypatch.delenv("CHEMSMART_MAX_SUBMITTERS", raising=False)
        runner = MagicMock(spec=["max_submitters"])
        runner.max_submitters = 9
        assert get_configured_max_submitters(jobrunner=runner) == 9

    def test_get_configured_max_submitters_from_server_max_submitters(
        self, monkeypatch
    ):
        monkeypatch.delenv("CHEMSMART_MAX_SUBMITTERS", raising=False)
        runner = MagicMock(spec=["server"])
        runner.server = MagicMock(spec=["max_submitters"])
        runner.server.max_submitters = 6
        assert get_configured_max_submitters(jobrunner=runner) == 6

    def test_get_configured_max_submitters_from_server_num_cores(
        self, monkeypatch
    ):
        monkeypatch.delenv("CHEMSMART_MAX_SUBMITTERS", raising=False)
        runner = MagicMock(spec=["server"])
        runner.server = MagicMock(spec=["num_cores"])
        runner.server.num_cores = 12
        assert get_configured_max_submitters(jobrunner=runner) == 12

    def test_get_configured_max_submitters_falls_back_to_cpu_count(
        self, monkeypatch
    ):
        monkeypatch.delenv("CHEMSMART_MAX_SUBMITTERS", raising=False)
        with patch("os.cpu_count", return_value=16):
            assert get_configured_max_submitters(jobrunner=None) == 16

    def test_get_submitter_worker_count_bounds_by_num_jobs(self, monkeypatch):
        monkeypatch.setenv("CHEMSMART_MAX_SUBMITTERS", "10")
        assert get_submitter_worker_count(None, num_jobs=3) == 3

    def test_get_submitter_worker_count_bounds_by_max_submitters(
        self, monkeypatch
    ):
        monkeypatch.setenv("CHEMSMART_MAX_SUBMITTERS", "2")
        assert get_submitter_worker_count(None, num_jobs=10) == 2

    def test_get_submitter_worker_count_at_least_one_for_zero_jobs(self):
        assert get_submitter_worker_count(None, num_jobs=0) == 1


class TestDecidePhaseTransition:
    def test_proceeds_with_no_failures(self):
        decision = decide_phase_transition(phase_name="opt")
        assert decision.proceed is True
        assert decision.should_raise is False

    def test_raises_on_failures(self):
        decision = decide_phase_transition(
            phase_name="opt", failures=["worker1 failed"]
        )
        assert decision.proceed is False
        assert decision.should_raise is True
        assert "worker1 failed" in decision.message

    def test_halts_without_raising_when_incomplete_required(self):
        decision = decide_phase_transition(
            phase_name="opt", require_complete=True, is_complete=False
        )
        assert decision.proceed is False
        assert decision.should_raise is False
        assert "incomplete" in decision.message


class TestExecutableClassForProgram:
    def test_falsy_program_returns_none(self):
        assert _executable_class_for_program(None) is None
        assert _executable_class_for_program("") is None

    def test_notimplemented_program_returns_none(self):
        assert _executable_class_for_program(NotImplemented) is None


class TestRunPhaseJobs:
    def test_delegates_to_job_execute_phase_jobs(self):
        with patch(
            "chemsmart.jobs.runner.Job._execute_phase_jobs"
        ) as mock_execute:
            run_phase_jobs(
                parent_runner="runner-sentinel",
                jobs=["job1"],
                stop_on_incomplete=True,
                phase_label="opt",
            )
        mock_execute.assert_called_once_with(
            parent_runner="runner-sentinel",
            jobs=["job1"],
            jobs_factory=None,
            stop_on_incomplete=True,
            before_run=None,
            logger_obj=None,
            phase_label="opt",
        )
