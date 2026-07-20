from __future__ import annotations

import os
import re
import shlex
import subprocess
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, TypedDict


class SubmitResult(TypedDict):
    submitted_at_path: str
    command_executed: str
    returncode: int
    job_id: str | None


@dataclass(frozen=True)
class ExecResult:
    returncode: int
    stdout: str
    stderr: str
    command: str
    server: str | None
    duration_s: float


class SubmitTransport(ABC):
    @abstractmethod
    def submit(self, script_path, working_dir, server) -> SubmitResult:
        raise NotImplementedError


class ExecTransport(ABC):
    @abstractmethod
    def exec(
        self,
        command: str,
        server: str | None,
        timeout_s: int = 15,
    ) -> ExecResult:
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
        env = _sge_env(command)
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False,
            env=env,
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


class SshExecTransport(ExecTransport):
    def exec(
        self,
        command: str,
        server: str | None,
        timeout_s: int = 15,
    ) -> ExecResult:
        if not isinstance(server, str) or not server.strip():
            raise ValueError("SshExecTransport requires a non-empty server")

        host = server.strip()
        args = ["ssh", *_ssh_options(), host, command]
        return _run_exec_command(
            args=args,
            command=command,
            server=host,
            timeout_s=timeout_s,
        )


class LocalExecTransport(ExecTransport):
    def exec(
        self,
        command: str,
        server: str | None,
        timeout_s: int = 15,
    ) -> ExecResult:
        if server not in (None, "localhost"):
            raise ValueError(
                "LocalExecTransport only supports server=None or "
                "server='localhost'"
            )

        return _run_exec_command(
            args=shlex.split(command),
            command=command,
            server=server,
            timeout_s=timeout_s,
        )


class MockExecTransport(ExecTransport):
    def __init__(
        self,
        scripted_results: list[ExecResult | Exception] | None = None,
    ):
        self.calls: list[dict[str, Any]] = []
        self.scripted_results = list(scripted_results or [])

    def exec(
        self,
        command: str,
        server: str | None,
        timeout_s: int = 15,
    ) -> ExecResult:
        call = {
            "command": command,
            "server": server,
            "timeout_s": timeout_s,
        }
        self.calls.append(call)

        if self.scripted_results:
            result = self.scripted_results.pop(0)
            if isinstance(result, Exception):
                raise result
            return result

        return ExecResult(
            returncode=0,
            stdout="",
            stderr="",
            command=command,
            server=server,
            duration_s=0.0,
        )


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


def _target_absolute_path(value: str) -> str:
    """Absolutise a scheduler path without breaking POSIX target paths.

    Submit paths belong to the host that will run the job: a remote
    cluster over ssh, or a local POSIX scheduler. ``os.path.abspath``
    rewrote an already-absolute ``/tmp/job.sh`` into ``D:\\tmp\\job.sh``
    when the agent ran on Windows, so the emitted ``qsub``/``ssh`` command
    named a path the target could never resolve. Leave POSIX-absolute
    values untouched and absolutise only genuinely relative ones.
    """

    if str(value).startswith("/"):
        return str(value)
    return os.path.abspath(value)


def build_submit_invocation(script_path, working_dir, server):
    submit_command = getattr(server, "submit_command", None) or "qsub"
    script_path = _target_absolute_path(script_path)
    working_dir = _target_absolute_path(working_dir)
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
    kwargs = getattr(server, "kwargs", None)
    if isinstance(kwargs, dict):
        host = kwargs.get("HOST")
        if isinstance(host, str) and host.strip():
            return host.strip()

    name = getattr(server, "name", None)
    if isinstance(name, str) and name.strip():
        basename = os.path.basename(name.strip())
        return os.path.splitext(basename)[0]
    return "localhost"


def _is_local_host(host: str) -> bool:
    return host.lower() in {"localhost", "local"}


def _ssh_options() -> list[str]:
    return [
        "-o",
        "BatchMode=yes",
        "-o",
        "ClearAllForwardings=yes",
        "-o",
        "ConnectTimeout=10",
        "-o",
        "StrictHostKeyChecking=yes",
    ]


def _run_exec_command(
    *,
    args: list[str],
    command: str,
    server: str | None,
    timeout_s: int,
) -> ExecResult:
    start = time.monotonic()
    try:
        result = subprocess.run(
            args,
            capture_output=True,
            text=True,
            timeout=timeout_s,
            check=False,
        )
        stdout = _coerce_text(result.stdout)
        stderr = _coerce_text(result.stderr)
        returncode = int(result.returncode)
    except subprocess.TimeoutExpired as exc:
        raise TimeoutError(
            f"Command timed out after {timeout_s}s: {command}"
        ) from exc
    except FileNotFoundError as exc:
        stdout = ""
        stderr = str(exc)
        returncode = 127
    except OSError as exc:
        stdout = ""
        stderr = str(exc)
        returncode = 1

    return ExecResult(
        returncode=returncode,
        stdout=stdout,
        stderr=stderr,
        command=command,
        server=server,
        duration_s=time.monotonic() - start,
    )


def _coerce_text(value: str | bytes | None) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value


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


def _sge_env(command: list[str]) -> dict[str, str] | None:
    """Return an env dict with SGE_ROOT set when running qsub locally.

    qsub requires SGE_ROOT to be set. When the submit_command is an absolute
    path like /opt/sge/bin/lx-amd64/qsub and SGE_ROOT is missing from the
    current environment, infer it from the path (3 levels up from the binary).
    """
    if command and command[0] == "ssh":
        return None
    submit_cmd = command[0] if command else ""
    if "SGE_ROOT" in os.environ or "qsub" not in os.path.basename(submit_cmd):
        return None
    if not os.path.isabs(submit_cmd):
        return None
    from pathlib import Path

    parents = Path(submit_cmd).resolve().parents
    if len(parents) >= 3:
        inferred_root = str(parents[2])
        env = os.environ.copy()
        env["SGE_ROOT"] = inferred_root
        env.setdefault("SGE_CELL", "default")
        return env
    return None
