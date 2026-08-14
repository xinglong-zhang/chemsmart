"""
Direct unit tests for ``chemsmart.settings.submitters``: ``RunScript``,
the abstract ``Submitter`` base class, and its scheduler-specific
subclasses (``PBSSubmitter``, ``SLURMSubmitter``, ``SLFSubmitter``,
``FUGAKUSubmitter``).

These write plain text scripts to the current working directory with no
subprocess execution, so most tests just run inside ``tmp_path`` (via
``monkeypatch.chdir``) and inspect the written content.
"""

import os
from io import StringIO
from unittest.mock import MagicMock, patch

import pytest

from chemsmart.settings.server import Server
from chemsmart.settings.submitters import (
    FUGAKUSubmitter,
    PBSSubmitter,
    RunScript,
    SLFSubmitter,
    SLURMSubmitter,
    Submitter,
    user_settings,
)


def make_job(label="job1", program="Gaussian", complete=False, folder="."):
    job = MagicMock()
    job.label = label
    job.PROGRAM = program
    job.folder = folder
    job.is_complete.return_value = complete
    return job


@pytest.fixture()
def pbs_test_server():
    return Server(
        "custom-pbs",
        SCHEDULER="PBS",
        NUM_CORES=8,
        MEM_GB=24,
        NUM_GPUS=0,
        NUM_HOURS=24,
        QUEUE_NAME="normal",
    )


@pytest.fixture()
def slurm_test_server():
    return Server(
        "custom-slurm",
        SCHEDULER="SLURM",
        NUM_CORES=8,
        MEM_GB=24,
        NUM_GPUS=0,
        NUM_HOURS=24,
        QUEUE_NAME="normal",
    )


class TestRunScript:
    def test_write_creates_file_with_expected_content(
        self, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        script = RunScript(filename="run.py", cli_args=["gaussian", "opt"])
        script.write()

        assert os.path.exists("run.py")
        with open("run.py") as f:
            content = f.read()
        assert "from chemsmart.cli.run import run" in content
        assert "run(['gaussian', 'opt'])" in content
        assert "if __name__ == '__main__':" in content


class TestSubmitterDunderMethods:
    def test_eq_compares_by_name(self, pbs_test_server):
        job = make_job()
        s1 = PBSSubmitter(job=job, server=pbs_test_server)
        s2 = PBSSubmitter(name="PBS", job=job, server=pbs_test_server)
        assert s1 == s2

    def test_hash_based_on_name(self, pbs_test_server):
        job = make_job()
        s1 = PBSSubmitter(job=job, server=pbs_test_server)
        assert hash(s1) == hash("PBS")

    def test_str_and_repr(self, pbs_test_server):
        job = make_job()
        s1 = PBSSubmitter(job=job, server=pbs_test_server)
        assert str(s1) == "Submitter: PBS"
        assert repr(s1) == "Submitter(name=PBS)"


class TestSubmitterPaths:
    def test_submit_folder_is_job_folder(self, pbs_test_server, tmp_path):
        job = make_job(folder=str(tmp_path))
        submitter = PBSSubmitter(job=job, server=pbs_test_server)
        assert submitter.submit_folder == str(tmp_path)

    def test_submit_script_uses_label(self, pbs_test_server):
        job = make_job(label="mylabel")
        submitter = PBSSubmitter(job=job, server=pbs_test_server)
        assert submitter.submit_script == "chemsmart_sub_mylabel.sh"

    def test_submit_script_default_without_label(self, pbs_test_server):
        job = make_job(label=None)
        submitter = PBSSubmitter(job=job, server=pbs_test_server)
        assert submitter.submit_script == "chemsmart_sub.sh"

    def test_array_submit_script_uses_label(self, pbs_test_server):
        job = make_job(label="mylabel")
        submitter = PBSSubmitter(job=job, server=pbs_test_server)
        assert (
            submitter.array_submit_script == "chemsmart_sub_array_mylabel.sh"
        )

    def test_run_script_uses_label(self, pbs_test_server):
        job = make_job(label="mylabel")
        submitter = PBSSubmitter(job=job, server=pbs_test_server)
        assert submitter.run_script == "chemsmart_run_mylabel.py"


class TestSubmitterExecutableDispatch:
    def test_unsupported_program_raises(self, pbs_test_server):
        job = make_job(program="SomeUnsupportedProgram")
        submitter = PBSSubmitter(job=job, server=pbs_test_server)
        with pytest.raises(ValueError, match="not supported"):
            _ = submitter.executable

    def test_gaussian_program_resolves_executable(
        self, pbs_test_server, server_yaml_file
    ):
        job = make_job(program="Gaussian")
        server = Server.from_yaml(name=server_yaml_file)
        submitter = PBSSubmitter(job=job, server=server)
        executable = submitter.executable
        assert executable is not None


class TestSubmitterWrite:
    def test_write_warns_when_job_already_complete(
        self, pbs_test_server, tmp_path, monkeypatch, caplog
    ):
        monkeypatch.chdir(tmp_path)
        job = make_job(complete=True)
        submitter = PBSSubmitter(job=job, server=pbs_test_server)
        with (
            patch.object(submitter, "_write_runscript") as mock_run,
            patch.object(submitter, "_write_submitscript") as mock_submit,
        ):
            submitter.write(cli_args=["gaussian", "opt"])

        mock_run.assert_called_once_with(["gaussian", "opt"])
        mock_submit.assert_called_once()

    def test_write_submitscript_produces_full_script(
        self, pbs_test_server, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = make_job(label="mylabel")
        submitter = PBSSubmitter(job=job, server=pbs_test_server)
        with patch.object(
            type(submitter),
            "executable",
            new=MagicMock(
                conda_env=None, modules=None, scripts=None, envars=None
            ),
        ):
            submitter._write_submitscript()

        assert os.path.exists("chemsmart_sub_mylabel.sh")
        with open("chemsmart_sub_mylabel.sh") as f:
            content = f.read()
        assert content.startswith("#!/bin/bash")
        assert "#PBS -o mylabel.pbsout" in content
        assert "cd $PBS_O_WORKDIR" in content
        assert "./chemsmart_run_mylabel.py &" in content


class TestSubmitterArrayJob:
    def test_write_array_job_no_jobs_warns_and_returns(
        self, pbs_test_server, caplog
    ):
        job = make_job()
        submitter = PBSSubmitter(job=job, server=pbs_test_server)
        with caplog.at_level("WARNING"):
            submitter.write_array_job(jobs=[])
        assert "No jobs provided" in caplog.text

    def test_write_array_job_shared_cli_args(
        self, pbs_test_server, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = make_job(label="mylabel")
        submitter = PBSSubmitter(job=job, server=pbs_test_server)
        jobs = [make_job(label=f"j{i}") for i in range(3)]

        with patch.object(
            submitter, "_write_array_submitscript"
        ) as mock_submit:
            submitter.write_array_job(jobs=jobs, cli_args=["gaussian", "opt"])

        assert submitter.jobs == jobs
        mock_submit.assert_called_once()
        for i in range(3):
            assert os.path.exists(f"chemsmart_run_array_{i}.py")
            with open(f"chemsmart_run_array_{i}.py") as f:
                assert "run(['gaussian', 'opt'])" in f.read()

    def test_write_array_job_per_job_cli_args(
        self, pbs_test_server, tmp_path, monkeypatch
    ):
        monkeypatch.chdir(tmp_path)
        job = make_job(label="mylabel")
        submitter = PBSSubmitter(job=job, server=pbs_test_server)
        jobs = [make_job(label=f"j{i}") for i in range(2)]
        per_job_args = [["gaussian", "opt", "0"], ["gaussian", "opt", "1"]]

        with patch.object(submitter, "_write_array_submitscript"):
            submitter.write_array_job(jobs=jobs, cli_args=per_job_args)

        with open("chemsmart_run_array_0.py") as f:
            assert "run(['gaussian', 'opt', '0'])" in f.read()
        with open("chemsmart_run_array_1.py") as f:
            assert "run(['gaussian', 'opt', '1'])" in f.read()

    def test_write_array_job_command_format(self, pbs_test_server):
        job = make_job()
        submitter = PBSSubmitter(job=job, server=pbs_test_server)
        buf = StringIO()
        submitter._write_array_job_command(buf)
        content = buf.getvalue()
        assert "SLURM_ARRAY_TASK_ID" in content
        assert "PBS_ARRAYID" in content
        assert "LSB_JOBINDEX" in content
        assert "python chemsmart_run_array_${TASK_ID}.py" in content


class TestSubmitterFromSchedulerType:
    def test_dispatches_to_pbs(self, pbs_test_server):
        job = make_job()
        result = Submitter.from_scheduler_type(
            "PBS", job=job, server=pbs_test_server
        )
        assert isinstance(result, PBSSubmitter)

    def test_dispatches_to_slurm(self, slurm_test_server):
        job = make_job()
        result = Submitter.from_scheduler_type(
            "SLURM", job=job, server=slurm_test_server
        )
        assert isinstance(result, SLURMSubmitter)

    def test_unknown_scheduler_raises(self, pbs_test_server):
        job = make_job()
        with pytest.raises(ValueError, match="Could not find any submitters"):
            Submitter.from_scheduler_type(
                "BOGUS_SCHEDULER", job=job, server=pbs_test_server
            )


class TestPBSSubmitter:
    def test_scheduler_options_basic(self, pbs_test_server):
        job = make_job(label="mylabel")
        submitter = PBSSubmitter(job=job, server=pbs_test_server)
        buf = StringIO()
        with patch.object(user_settings, "data", {}):
            submitter._write_scheduler_options(buf)
        content = buf.getvalue()
        assert "#PBS -o mylabel.pbsout" in content
        assert "#PBS -e mylabel.pbserr" in content
        assert "#PBS -q normal" in content
        assert "#PBS -l walltime=24:00:00" in content
        assert "gpus=" not in content

    def test_scheduler_options_includes_gpus_when_requested(self, tmp_path):
        server = Server(
            "gpu-pbs", SCHEDULER="PBS", NUM_CORES=8, MEM_GB=24, NUM_GPUS=2
        )
        job = make_job(label="mylabel")
        submitter = PBSSubmitter(job=job, server=server)
        buf = StringIO()
        with patch.object(user_settings, "data", {}):
            submitter._write_scheduler_options(buf)
        assert "#PBS -l gpus=2" in buf.getvalue()

    def test_scheduler_options_includes_project_and_email(
        self, pbs_test_server
    ):
        job = make_job(label="mylabel")
        submitter = PBSSubmitter(job=job, server=pbs_test_server)
        buf = StringIO()
        with patch.object(
            user_settings,
            "data",
            {"PROJECT": "proj123", "EMAIL": "me@example.com"},
        ):
            submitter._write_scheduler_options(buf)
        content = buf.getvalue()
        assert "#PBS -P proj123" in content
        assert "#PBS -M me@example.com" in content
        assert "#PBS -m abe" in content

    def test_change_to_job_directory(self, pbs_test_server):
        job = make_job()
        submitter = PBSSubmitter(job=job, server=pbs_test_server)
        buf = StringIO()
        submitter._write_change_to_job_directory(buf)
        assert buf.getvalue() == "cd $PBS_O_WORKDIR\n\n"


class TestSLURMSubmitter:
    def test_scheduler_options_basic(self, slurm_test_server):
        job = make_job(label="mylabel")
        submitter = SLURMSubmitter(job=job, server=slurm_test_server)
        buf = StringIO()
        with patch.object(user_settings, "data", {}):
            submitter._write_scheduler_options(buf)
        content = buf.getvalue()
        assert "#SBATCH --job-name=mylabel" in content
        assert "#SBATCH --partition=normal" in content
        assert "#SBATCH --time=24:00:00" in content

    def test_array_scheduler_options_default_throttle(self, slurm_test_server):
        job = make_job(label="mylabel")
        submitter = SLURMSubmitter(job=job, server=slurm_test_server)
        submitter.jobs = [make_job() for _ in range(5)]
        buf = StringIO()
        with patch.object(user_settings, "data", {}):
            submitter._write_array_scheduler_options(buf, num_nodes=None)
        assert "#SBATCH --array=1-5\n" in buf.getvalue()

    def test_array_scheduler_options_with_node_throttle(
        self, slurm_test_server
    ):
        job = make_job(label="mylabel")
        submitter = SLURMSubmitter(job=job, server=slurm_test_server)
        submitter.jobs = [make_job() for _ in range(5)]
        buf = StringIO()
        with patch.object(user_settings, "data", {}):
            submitter._write_array_scheduler_options(buf, num_nodes=2)
        assert "#SBATCH --array=1-5%2\n" in buf.getvalue()

    def test_change_to_job_directory(self, slurm_test_server):
        job = make_job()
        submitter = SLURMSubmitter(job=job, server=slurm_test_server)
        buf = StringIO()
        submitter._write_change_to_job_directory(buf)
        assert buf.getvalue() == "cd $SLURM_SUBMIT_DIR\n\n"


class TestSLFSubmitterBug:
    """Documents a real bug: ``Server`` has no ``num_nodes`` attribute at
    all, but ``SLFSubmitter._write_scheduler_options`` unconditionally
    reads ``self.server.num_nodes``. Any attempt to write an LSF/SLF
    submission script crashes. See BUGS_FOUND.md."""

    def test_change_to_job_directory_works(self):
        server = Server("slf-server", SCHEDULER="SLF", NUM_CORES=8, MEM_GB=24)
        job = make_job()
        submitter = SLFSubmitter(job=job, server=server)
        buf = StringIO()
        submitter._write_change_to_job_directory(buf)
        assert buf.getvalue() == "cd $LS_SUBCWD\n\n"

    def test_scheduler_options_crashes_on_missing_num_nodes(self):
        server = Server(
            "slf-server",
            SCHEDULER="SLF",
            NUM_CORES=8,
            MEM_GB=24,
            NUM_GPUS=0,
            NUM_HOURS=24,
        )
        job = make_job(label="mylabel")
        submitter = SLFSubmitter(job=job, server=server)
        buf = StringIO()
        with patch.object(user_settings, "data", {}):
            with pytest.raises(AttributeError, match="num_nodes"):
                submitter._write_scheduler_options(buf)


class TestFUGAKUSubmitterBug:
    """Documents a real bug: ``FUGAKUSubmitter._write_scheduler_options``
    references ``self.project``, which is never assigned anywhere in the
    class or its base. Any attempt to write a FUGAKU submission script
    crashes (after also requiring ``RSCGRP`` in user settings). See
    BUGS_FOUND.md."""

    def test_change_to_job_directory_works(self):
        server = Server(
            "fugaku-server", SCHEDULER="FUGAKU", NUM_CORES=8, MEM_GB=24
        )
        job = make_job()
        submitter = FUGAKUSubmitter(job=job, server=server)
        buf = StringIO()
        submitter._write_change_to_job_directory(buf)
        assert buf.getvalue() == "cd $PJM_O_WORKDIR\n\n"

    def test_scheduler_options_crashes_on_missing_project_attr(self):
        server = Server(
            "fugaku-server",
            SCHEDULER="FUGAKU",
            NUM_CORES=8,
            MEM_GB=24,
            NUM_HOURS=24,
        )
        job = make_job(label="mylabel")
        submitter = FUGAKUSubmitter(job=job, server=server)
        buf = StringIO()
        with patch.object(user_settings, "data", {"RSCGRP": "small"}):
            with pytest.raises(AttributeError, match="project"):
                submitter._write_scheduler_options(buf)
