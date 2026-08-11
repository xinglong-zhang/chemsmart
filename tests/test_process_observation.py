"""Deterministic, chemistry-free checks for bounded process execution."""

import os
import subprocess
import sys

from chemsmart.utils.process_observation import (
    launch_failure_observation,
    observe_process,
)


def _python_process(source):
    return subprocess.Popen(
        [sys.executable, "-c", source],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )


def test_process_observation_records_normal_exit_and_memory():
    process = _python_process(
        "import time; value = bytearray(2 * 1024 * 1024); "
        "print(len(value), flush=True); time.sleep(0.15)"
    )

    result = observe_process(
        process,
        timeout_seconds=2,
        memory_limit_mb=256,
        sample_interval_seconds=0.02,
    )

    assert result.stdout.strip() == "2097152"
    assert result.observation.state == "exited"
    assert result.observation.returncode == 0
    assert result.observation.process_group_owned is True
    assert result.observation.timed_out is False
    assert result.observation.memory_limit_exceeded is False
    assert result.observation.memory_observation_state == "observed"
    assert result.observation.peak_rss_mb > 0
    assert result.observation.termination_requested is False
    assert result.observation.findings == ()


def test_timeout_terminates_the_owned_process_group():
    process = _python_process(
        "import subprocess, sys, time; "
        "child = subprocess.Popen([sys.executable, '-c', "
        "'import time; time.sleep(30)']); "
        "print(child.pid, flush=True); time.sleep(30)"
    )

    result = observe_process(
        process,
        timeout_seconds=0.2,
        memory_limit_mb=256,
        sample_interval_seconds=0.02,
        termination_grace_seconds=0.2,
    )

    child_pid = int(result.stdout.strip())
    assert result.observation.state == "timed_out_terminated"
    assert result.observation.timed_out is True
    assert result.observation.termination_requested is True
    assert result.observation.termination_confirmed is True
    assert result.observation.termination_scope == "process_group_tree"
    assert result.observation.remaining_process_ids == ()
    assert "process.timeout" in result.observation.findings
    try:
        os.kill(child_pid, 0)
    except ProcessLookupError:
        pass
    else:  # pragma: no cover - assertion includes the observed child PID
        raise AssertionError("timed-out descendant remains alive")


def test_observed_memory_limit_is_not_silently_relaxed():
    process = _python_process(
        "import time; value = bytearray(8 * 1024 * 1024); time.sleep(30)"
    )

    result = observe_process(
        process,
        timeout_seconds=2,
        memory_limit_mb=1,
        sample_interval_seconds=0.02,
        termination_grace_seconds=0.2,
    )

    assert result.observation.state == "memory_limit_exceeded_terminated"
    assert result.observation.memory_limit_mb == 1
    assert result.observation.memory_limit_exceeded is True
    assert result.observation.peak_rss_mb > 1
    assert result.observation.termination_confirmed is True
    assert "process.memory_limit_exceeded" in result.observation.findings


def test_launch_failure_is_explicit_and_digest_bound():
    observation = launch_failure_observation(
        timeout_seconds=600,
        memory_limit_mb=4096,
        error_type="FileNotFoundError",
    )

    assert observation.state == "launch_failed"
    assert observation.pid is None
    assert observation.timeout_seconds == 600
    assert observation.memory_limit_mb == 4096
    assert observation.termination_confirmed is None
    assert observation.findings == (
        "process.launch_failed.filenotfounderror",
    )
    assert len(observation.receipt_sha256) == 64
