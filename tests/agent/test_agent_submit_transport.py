from __future__ import annotations

from pathlib import Path

from chemsmart.agent.transport import (
    LocalDryRunTransport,
    MockTransport,
    SshQsubTransport,
)
from chemsmart.settings.server import Server


def _server(name: str) -> Server:
    return Server.from_dict(
        {
            "name": name,
            "SCHEDULER": "PBS",
            "SUBMIT_COMMAND": "qsub",
        }
    )


def test_local_dry_run_transport_never_calls_subprocess_run(
    monkeypatch,
    tmp_path: Path,
):
    script_path = tmp_path / "chemsmart_sub_test.sh"
    script_path.write_text("#!/bin/bash\n")

    def fail(*args, **kwargs):
        raise AssertionError("subprocess.run should not be called")

    monkeypatch.setattr("subprocess.run", fail)

    result = LocalDryRunTransport().submit(
        script_path=str(script_path),
        working_dir=str(tmp_path),
        server=_server("local"),
    )

    assert result["returncode"] == 0
    assert result["job_id"] is None
    assert result["command_executed"] == f"qsub {script_path}"


def test_ssh_qsub_transport_uses_exact_remote_command(
    monkeypatch,
    tmp_path: Path,
):
    script_path = tmp_path / "chemsmart_sub_test.sh"
    script_path.write_text("#!/bin/bash\n")
    calls = []

    class FakeCompletedProcess:
        returncode = 0
        stdout = "12345.remote\n"
        stderr = ""

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        return FakeCompletedProcess()

    monkeypatch.setattr("subprocess.run", fake_run)

    result = SshQsubTransport().submit(
        script_path=str(script_path),
        working_dir=str(tmp_path),
        server=_server("remote-hpc"),
    )

    assert calls == [
        (
            [
                "ssh",
                "remote-hpc",
                f"cd {tmp_path} && qsub {script_path}",
            ],
            {
                "capture_output": True,
                "text": True,
                "check": False,
                "cwd": None,
            },
        )
    ]
    assert result["returncode"] == 0
    assert result["job_id"] == "12345.remote"
    assert (
        result["command_executed"]
        == f'ssh remote-hpc "cd {tmp_path} && qsub {script_path}"'
    )


def test_mock_transport_records_every_call(tmp_path: Path):
    script_path = tmp_path / "chemsmart_sub_test.sh"
    script_path.write_text("#!/bin/bash\n")
    transport = MockTransport()

    first = transport.submit(
        script_path=str(script_path),
        working_dir=str(tmp_path),
        server=_server("remote-hpc"),
    )
    second = transport.submit(
        script_path=str(script_path),
        working_dir=str(tmp_path),
        server=_server("remote-hpc"),
    )

    assert [call["server_name"] for call in transport.calls] == [
        "remote-hpc",
        "remote-hpc",
    ]
    assert first["job_id"] == "mock-job-0001"
    assert second["job_id"] == "mock-job-0002"
