from pathlib import Path
from types import SimpleNamespace

from click.testing import CliRunner

from chemsmart.cli.run import run
from chemsmart.cli.sub import sub
from chemsmart.jobs.gaussian.runner import (
    FakeGaussianJobRunner,
    GaussianJobRunner,
)
from chemsmart.jobs.iterate.runner import IterateJobRunner
from chemsmart.jobs.orca.runner import FakeORCAJobRunner, ORCAJobRunner
from chemsmart.jobs.runner import JobRunner
from chemsmart.jobs.xtb.runner import FakeXTBJobRunner, XTBJobRunner
from chemsmart.settings.server import Server


class DummyORCAJob:
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
    ):
        monkeypatch.setenv(
            "CHEMSMART_CONFIG_DIR", str(gaussian_project_config_dir)
        )
        fake_server = Server(name="dummy")
        captured = {"cli_args": None}
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
