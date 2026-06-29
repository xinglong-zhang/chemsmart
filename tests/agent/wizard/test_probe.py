import shlex
import subprocess

import pytest

from chemsmart.agent.wizard import ALL_PROBE_SPECS, ProbeError, ProbeRunner
from chemsmart.agent.wizard.probe import MAX_OUTPUT_BYTES, TRUNCATION_MARKER


def test_run_local_printenv():
    result = ProbeRunner.run_local(ALL_PROBE_SPECS["common.printenv_all"])

    assert result.returncode == 0
    assert result.mode == "local"
    assert result.host is None
    assert result.command == "printenv"
    assert isinstance(result.duration_s, float)


@pytest.mark.parametrize("command", ["rm -rf /", ["hostname"], []])
def test_run_local_rejects_raw_commands(command):
    with pytest.raises(ProbeError):
        ProbeRunner.run_local(command)


def test_run_local_truncates_large_output(monkeypatch):
    def fake_run(*args, **kwargs):
        return subprocess.CompletedProcess(
            args[0],
            0,
            stdout="x" * (MAX_OUTPUT_BYTES + 32),
            stderr="",
        )

    monkeypatch.setattr(subprocess, "run", fake_run)

    result = ProbeRunner.run_local(ALL_PROBE_SPECS["common.printenv_all"])

    assert result.returncode == 0
    assert result.truncated is True
    assert len(result.stdout) == MAX_OUTPUT_BYTES
    assert result.stdout.endswith(TRUNCATION_MARKER)


def test_run_local_timeout(monkeypatch):
    def fake_run(*args, **kwargs):
        raise subprocess.TimeoutExpired(cmd=args[0], timeout=kwargs["timeout"])

    monkeypatch.setattr(subprocess, "run", fake_run)

    result = ProbeRunner.run_local(ALL_PROBE_SPECS["common.printenv_all"])

    assert result.returncode == 124
    assert "[probe-timeout]" in result.stderr
    assert result.mode == "local"
    assert result.host is None


def test_run_ssh_uses_hardened_flags(monkeypatch):
    calls = []

    def fake_run(*args, **kwargs):
        calls.append(args[0])
        return subprocess.CompletedProcess(args[0], 0, stdout="ok", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    host = "cluster.example.edu"
    result = ProbeRunner.run_ssh(
        host,
        ALL_PROBE_SPECS["software.command_v"],
        exe_name="hostname",
    )

    ssh_args = calls[0]
    assert result.returncode == 0
    assert result.mode == "ssh"
    assert result.host == host
    assert ssh_args[0] == "ssh"
    assert "-A" not in ssh_args
    assert "StrictHostKeyChecking=no" not in ssh_args
    assert "accept-new" not in " ".join(ssh_args)
    assert ssh_args[1:9] == [
        "-o",
        "BatchMode=yes",
        "-o",
        "ClearAllForwardings=yes",
        "-o",
        "ConnectTimeout=10",
        "-o",
        "StrictHostKeyChecking=yes",
    ]
    assert ssh_args[9] == host
    assert ssh_args[10] == shlex.join(
        [
            "bash",
            "-lc",
            shlex.join(["command", "-v", "hostname"]),
        ]
    )


def test_run_ssh_rejects_raw_commands():
    with pytest.raises(ProbeError):
        ProbeRunner.run_ssh("cluster.example.edu", "rm -rf /")
