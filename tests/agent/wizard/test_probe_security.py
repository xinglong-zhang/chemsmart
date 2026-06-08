import shlex
import subprocess

import pytest

from chemsmart.agent.wizard import ALL_PROBE_SPECS, ProbeError, ProbeRunner


def test_raw_run_local_string_is_rejected():
    with pytest.raises(ProbeError, match="raw command string not allowed"):
        ProbeRunner.run_local("rm -rf /")


def test_shell_metacharacters_in_path_slot_are_rejected():
    with pytest.raises(ProbeError):
        ProbeRunner.run_local(
            ALL_PROBE_SPECS["software.readlink"],
            path="/tmp;rm",
        )


def test_command_substitution_in_path_slot_is_rejected():
    with pytest.raises(ProbeError):
        ProbeRunner.run_local(
            ALL_PROBE_SPECS["software.readlink"],
            path="/tmp/$(whoami)",
        )


def test_valid_probe_spec_runs_successfully(monkeypatch):
    def fake_run(*args, **kwargs):
        return subprocess.CompletedProcess(
            args[0],
            0,
            stdout="/tmp/ok\n",
            stderr="",
        )

    monkeypatch.setattr(subprocess, "run", fake_run)

    result = ProbeRunner.run_local(
        ALL_PROBE_SPECS["software.readlink"],
        path="/tmp/ok",
    )

    assert result.returncode == 0
    assert result.command == "readlink -f /tmp/ok"


def test_ssh_final_command_has_no_shell_composition(monkeypatch):
    calls = []

    def fake_run(*args, **kwargs):
        calls.append(args[0])
        return subprocess.CompletedProcess(
            args[0], 0, stdout="ok\n", stderr=""
        )

    monkeypatch.setattr(subprocess, "run", fake_run)

    ProbeRunner.run_ssh(
        "cluster.example.edu",
        ALL_PROBE_SPECS["software.readlink"],
        path="/tmp/ok",
    )

    final_cmd = calls[0][-1]
    assert final_cmd == shlex.join(
        ["bash", "-lc", shlex.join(["readlink", "-f", "/tmp/ok"])]
    )
    assert "&&" not in final_cmd
    assert "|" not in final_cmd
    assert "$(" not in final_cmd
