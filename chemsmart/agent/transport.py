from __future__ import annotations

import os
import re
import shlex
import subprocess
from abc import ABC, abstractmethod
from typing import Any, TypedDict


class SubmitResult(TypedDict):
    submitted_at_path: str
    command_executed: str
    returncode: int
    job_id: str | None


class SubmitTransport(ABC):
    @abstractmethod
    def submit(self, script_path, working_dir, server) -> SubmitResult:
        raise NotImplementedError


class LocalDryRunTransport(SubmitTransport):
    def submit(self, script_path, working_dir, server) -> SubmitResult:
        return {
            "submitted_at_path": os.path.abspath(script_path),
            "command_executed": build_submit_command(
                script_path=script_path,
                working_dir=working_dir,
                server=server,
            ),
            "returncode": 0,
            "job_id": None,
        }


class SshQsubTransport(SubmitTransport):
    def submit(self, script_path, working_dir, server) -> SubmitResult:
        command, executed = build_submit_invocation(
            script_path=script_path,
            working_dir=working_dir,
            server=server,
        )
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False,
            cwd=None if command[0] == "ssh" else os.path.abspath(working_dir),
        )
        output = "\n".join(
            part.strip()
            for part in (result.stdout, result.stderr)
            if isinstance(part, str) and part.strip()
        )
        return {
            "submitted_at_path": os.path.abspath(script_path),
            "command_executed": executed,
            "returncode": int(result.returncode),
            "job_id": _parse_job_id(output),
        }


class MockTransport(SubmitTransport):
    def __init__(self):
        self.calls: list[dict[str, Any]] = []

    def submit(self, script_path, working_dir, server) -> SubmitResult:
        call_number = len(self.calls) + 1
        call = {
            "script_path": os.path.abspath(script_path),
            "working_dir": os.path.abspath(working_dir),
            "server_name": _server_host(server),
        }
        self.calls.append(call)
        return {
            "submitted_at_path": call["script_path"],
            "command_executed": build_submit_command(
                script_path=script_path,
                working_dir=working_dir,
                server=server,
            ),
            "returncode": 0,
            "job_id": f"mock-job-{call_number:04d}",
        }


def build_submit_command(script_path, working_dir, server) -> str:
    _, executed = build_submit_invocation(
        script_path=script_path,
        working_dir=working_dir,
        server=server,
    )
    return executed


def build_submit_invocation(script_path, working_dir, server):
    submit_command = getattr(server, "submit_command", None) or "qsub"
    script_path = os.path.abspath(script_path)
    working_dir = os.path.abspath(working_dir)
    host = _server_host(server)

    if _is_local_host(host):
        return [submit_command, script_path], f"{submit_command} {script_path}"

    remote_command = (
        f"cd {shlex.quote(working_dir)} && "
        f"{submit_command} {shlex.quote(script_path)}"
    )
    executed = f'ssh {host} "{remote_command}"'
    return ["ssh", host, remote_command], executed


def _server_host(server) -> str:
    name = getattr(server, "name", None)
    if isinstance(name, str) and name.strip():
        basename = os.path.basename(name.strip())
        return os.path.splitext(basename)[0]
    return "localhost"


def _is_local_host(host: str) -> bool:
    return host.lower() in {"localhost", "local"}


def _parse_job_id(output: str) -> str | None:
    if not output:
        return None

    match = re.search(r"\b(\d+(?:[._-][A-Za-z0-9._-]+)?)\b", output)
    if match is not None:
        return match.group(1)

    first_token = output.strip().split()
    if first_token:
        return first_token[0]
    return None
