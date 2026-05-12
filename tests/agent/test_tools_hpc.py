from __future__ import annotations

from chemsmart.agent.tools_hpc import ssh_probe
from chemsmart.agent.transport import ExecResult, MockExecTransport
from chemsmart.agent.wizard.probe import ALL_PROBE_SPECS


def test_ssh_probe_known_probe_returns_parsed_payload(monkeypatch):
    transport = MockExecTransport(
        scripted_results=[
            ExecResult(
                returncode=0,
                stdout='{"nodes": []}',
                stderr="",
                command="sinfo --json",
                server="cluster-a",
                duration_s=0.25,
            )
        ]
    )
    monkeypatch.setattr(
        "chemsmart.agent.tools_hpc._TRANSPORT_FACTORY",
        lambda: transport,
    )

    result = ssh_probe(
        server="cluster-a",
        probe_name="survey.slurm.sinfo_json",
    )

    assert result == {
        "server": "cluster-a",
        "probe": "survey.slurm.sinfo_json",
        "returncode": 0,
        "stdout_truncated": '{"nodes": []}',
        "stderr_truncated": "",
        "parsed": {"nodes": []},
        "duration_s": 0.25,
    }
    assert transport.calls == [
        {
            "command": "sinfo --json",
            "server": "cluster-a",
            "timeout_s": 15,
        }
    ]


def test_ssh_probe_unknown_probe_returns_catalog():
    result = ssh_probe(server="cluster-a", probe_name="not-real")

    assert result == {
        "error": "unknown_probe",
        "available": sorted(ALL_PROBE_SPECS),
    }


def test_ssh_probe_captures_non_zero_returncode_and_stderr(monkeypatch):
    transport = MockExecTransport(
        scripted_results=[
            ExecResult(
                returncode=1,
                stdout="",
                stderr="permission denied",
                command="printenv",
                server="cluster-a",
                duration_s=0.5,
            )
        ]
    )
    monkeypatch.setattr(
        "chemsmart.agent.tools_hpc._TRANSPORT_FACTORY",
        lambda: transport,
    )

    result = ssh_probe(server="cluster-a", probe_name="common.printenv_all")

    assert result["returncode"] == 1
    assert result["stderr_truncated"] == "permission denied"
    assert result["parsed"] is None


def test_ssh_probe_timeout_returns_timeout_error(monkeypatch):
    transport = MockExecTransport(scripted_results=[TimeoutError("too slow")])
    monkeypatch.setattr(
        "chemsmart.agent.tools_hpc._TRANSPORT_FACTORY",
        lambda: transport,
    )

    result = ssh_probe(server="cluster-a", probe_name="common.printenv_all")

    assert result == {
        "error": "timeout",
        "server": "cluster-a",
        "probe": "common.printenv_all",
    }
