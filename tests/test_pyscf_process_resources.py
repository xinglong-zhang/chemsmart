"""Chemistry-free checks for PySCF process resource evidence."""

import subprocess
import sys
from unittest.mock import patch

from chemsmart.jobs.pyscf.runner import PySCFJobRunner


def test_pyscf_runner_observes_exact_default_limits():
    runner = object.__new__(PySCFJobRunner)
    runner.mem_gb = 4
    process = subprocess.Popen(
        [sys.executable, "-c", "print('complete')"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )

    returncode = runner._run(process)

    assert returncode == 0
    assert runner._process_observation["state"] == "exited"
    assert runner._process_observation["timeout_seconds"] == 600
    assert runner._process_observation["memory_limit_mb"] == 4096
    assert runner._process_observation["peak_rss_mb"] > 0


def test_pyscf_process_creation_owns_a_new_session(tmp_path):
    runner = object.__new__(PySCFJobRunner)
    runner.job_errfile = str(tmp_path / "water.err")
    runner.job_outputfile = str(tmp_path / "water.out")
    runner.job_resultsfile = str(tmp_path / "water.h5")
    runner.running_directory = str(tmp_path)
    command = (sys.executable, str(tmp_path / "water.py"))

    with patch("chemsmart.jobs.pyscf.runner.subprocess.Popen") as popen:
        runner._create_process(None, command, {"PATH": "bounded"})

    assert popen.call_args.args[0] == command
    assert popen.call_args.kwargs["start_new_session"] is True
    assert popen.call_args.kwargs.get("shell", False) is False


def test_pyscf_launch_failure_records_boundaries_without_retry():
    runner = object.__new__(PySCFJobRunner)
    runner.mem_gb = 4
    command = ("/missing/pyscf/python", "water.py")
    with (
        patch.object(runner, "_prerun"),
        patch.object(runner, "_write_input"),
        patch.object(runner, "_get_command", return_value=command),
        patch.object(runner, "_update_os_environ", return_value={}),
        patch.object(
            runner,
            "_create_process",
            side_effect=FileNotFoundError("missing interpreter"),
        ) as create_process,
        patch.object(runner, "_postrun"),
        patch.object(runner, "_postrun_cleanup"),
    ):
        returncode = runner.run(object())

    assert returncode is None
    assert create_process.call_count == 1
    assert runner._process_observation["state"] == "launch_failed"
    assert runner._process_observation["timeout_seconds"] == 600
    assert runner._process_observation["memory_limit_mb"] == 4096
    assert runner._process_observation["termination_requested"] is False
