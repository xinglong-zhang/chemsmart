import os
from io import StringIO

from click.testing import CliRunner

from chemsmart.cli.run import run
from chemsmart.cli.sub import sub
from chemsmart.settings.executable import GaussianExecutable, ORCAExecutable
from chemsmart.settings.server import Server
from chemsmart.settings.submitters import PBSSubmitter, SLURMSubmitter


class TestServer:
    def test_server_yaml(self, server_yaml_file):
        assert os.path.exists(server_yaml_file)
        assert os.path.isfile(server_yaml_file)
        server = Server.from_yaml(name=server_yaml_file)
        assert server.scheduler.lower() == "pbs"
        assert server.queue_name == "normal"
        assert server.num_hours == 24
        assert server.mem_gb == 375
        assert server.num_cores == 64
        assert server.num_gpus == 0
        assert server.num_threads == 64
        assert server.submit_command == "qsub"
        assert server.scratch_dir is None
        assert server.use_hosts is True
        assert (
            server.extra_commands == """export PATH=$HOME/bin/chemsmart:$PATH
export PATH=$HOME/bin/chemsmart/chemsmart/cli:$PATH
export PATH=$HOME/bin/chemsmart/chemsmart/scripts:$PATH
export PYTHONPATH=$HOME/bin/chemsmart:$PYTHONPATH
"""
        )
        assert server.extra_scheduler_directives == "#PBS -m abe\n"

    def test_gaussian_executable(self, server_yaml_file):
        gaussian_executable = GaussianExecutable.from_servername(
            server_yaml_file
        )
        assert gaussian_executable.executable_folder == os.path.expanduser(
            "~/programs/g16"
        )
        assert gaussian_executable.local_run is True

        gaussian_conda_env = """source ~/anaconda3/etc/profile.d/conda.sh
conda activate ~/anaconda3/envs/chemsmart
"""
        assert gaussian_executable.conda_env == gaussian_conda_env

        gaussian_modules = """module purge
module load craype-x86-rome
module load libfabric/1.11.0.4.125
"""
        assert gaussian_executable.modules == gaussian_modules

        assert (
            gaussian_executable.scripts
            == 'tcsh -c "source ~/programs/g16/bsd/g16.login"\n'
        )

        gassian_envars = """export SCRATCH=~/scratch
export GAUSS_EXEDIR=~/programs/g16
export g16root=~/programs/g16

"""
        assert gaussian_executable.envars == gassian_envars

    def test_orca_executable(self, server_yaml_file):
        orca_executable = ORCAExecutable.from_servername(server_yaml_file)
        assert orca_executable.executable_folder == os.path.expanduser(
            "~/programs/orca_6_0_0"
        )
        assert orca_executable.local_run is False

        assert orca_executable.conda_env is None

        assert orca_executable.modules is None

        assert orca_executable.scripts is None

        orca_envars = """export PATH=~/programs/openmpi-4.1.6/build/bin:$PATH
export LD_LIBRARY_PATH=~/programs/openmpi-4.1.6/build/lib:$LD_LIBRARY_PATH
"""
        assert orca_executable.envars == orca_envars

    def test_slurm_submitter_writes_extra_scheduler_directives(self):
        server = Server(
            "custom-slurm",
            SCHEDULER="SLURM",
            NUM_CORES=8,
            MEM_GB=24,
            NUM_GPUS=0,
            EXTRA_SCHEDULER_DIRECTIVES="#SBATCH --reservation=xlzhang_1\n",
        )
        job = type("DummyJob", (), {"label": "job1"})()
        submitter = SLURMSubmitter(job=job, server=server)

        buffer = StringIO()
        submitter._write_scheduler_options(buffer)
        assert "#SBATCH --reservation=xlzhang_1\n" in buffer.getvalue()

    def test_pbs_submitter_writes_extra_scheduler_directives(self):
        server = Server(
            "custom-pbs",
            SCHEDULER="PBS",
            NUM_CORES=8,
            MEM_GB=24,
            NUM_GPUS=0,
            EXTRA_SCHEDULER_DIRECTIVES="#PBS -m abe\n",
        )
        job = type("DummyJob", (), {"label": "job1"})()
        submitter = PBSSubmitter(job=job, server=server)

        buffer = StringIO()
        submitter._write_scheduler_options(buffer)
        assert "#PBS -m abe\n" in buffer.getvalue()


def _write_gaussian_project(tmp_path):
    config_root = tmp_path / "chemsmart_cfg"
    gaussian_cfg = config_root / "gaussian"
    gaussian_cfg.mkdir(parents=True)
    (gaussian_cfg / "test.yaml").write_text(
        "gas:\n  functional: B3LYP\n  basis: def2-SVP\n"
        "solv:\n  functional: B3LYP\n  basis: def2-SVP\n"
        "  solvent_model: smd\n  solvent_id: water\n"
    )
    return config_root


class TestScratchCLIDefaults:
    """Omitted --scratch should leave program defaults via from_job."""

    def test_run_omitted_scratch_leaves_none_for_program_default(
        self, tmp_path, monkeypatch, single_molecule_xyz_file
    ):
        monkeypatch.setenv(
            "CHEMSMART_CONFIG_DIR", str(_write_gaussian_project(tmp_path))
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

    def test_sub_omitted_scratch_does_not_force_no_scratch(
        self, tmp_path, monkeypatch, single_molecule_xyz_file
    ):
        monkeypatch.setenv(
            "CHEMSMART_CONFIG_DIR", str(_write_gaussian_project(tmp_path))
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

    def test_run_explicit_scratch_is_not_silently_disabled(
        self, tmp_path, monkeypatch, single_molecule_xyz_file
    ):
        """Placeholder JobRunner must not clear --scratch before from_job."""
        monkeypatch.setenv(
            "CHEMSMART_CONFIG_DIR", str(_write_gaussian_project(tmp_path))
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

    def test_sub_explicit_scratch_keeps_jobrunner_scratch_true(
        self, tmp_path, monkeypatch, single_molecule_xyz_file
    ):
        monkeypatch.setenv(
            "CHEMSMART_CONFIG_DIR", str(_write_gaussian_project(tmp_path))
        )
        fake_server = Server(name="dummy")
        captured = {"job": None, "cli_args": None}

        def _submit(job, test=False, cli_args=None, **kw):
            captured["job"] = job
            captured["cli_args"] = cli_args

        fake_server.submit = _submit
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )

        result = CliRunner().invoke(
            sub,
            [
                "--test",
                "--scratch",
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
        assert "--scratch" in captured["cli_args"]
        assert captured["job"].jobrunner.scratch is True
