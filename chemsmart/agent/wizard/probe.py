"""Probe primitives for wizard transport-aware discovery."""

from __future__ import annotations

import shlex
import subprocess
import time
from dataclasses import dataclass
from typing import Sequence

TRUNCATION_MARKER = "…[truncated]"
MAX_OUTPUT_BYTES = 1024 * 1024
ALLOWED_COMMANDS = {
    "hostname",
    "env",
    "command",
    "type",
    "which",
    "sinfo",
    "scontrol",
    "squeue",
    "qstat",
    "qconf",
    "qhost",
    "bqueues",
    "bhosts",
    "lsid",
    "pbsnodes",
    "module",
    "conda",
    "df",
    "getent",
    "test",
    "readlink",
    "find",
    "lscpu",
    "printf",
    "echo",
    "uname",
    "id",
    "groups",
    "sacctmgr",
    "sacct",
}


class ProbeError(Exception):
    """Raised when a probe request is invalid."""


@dataclass(frozen=True)
class ProbeResult:
    command: str
    mode: str
    host: str | None
    returncode: int
    stdout: str
    stderr: str
    duration_s: float
    truncated: bool


class ProbeRunner:
    """Run allowlisted read-only probes locally or via SSH."""

    @classmethod
    def run_local(
        cls, command: list[str], timeout_s: int = 15
    ) -> ProbeResult:
        cls._validate_local_command(command)
        display_command = shlex.join(command)
        return cls._run(
            args=list(command),
            command_text=display_command,
            mode="local",
            host=None,
            timeout_s=timeout_s,
        )

    @classmethod
    def run_ssh(
        cls, host: str, command: str, timeout_s: int = 15
    ) -> ProbeResult:
        cls._validate_ssh_command(command)
        ssh_command = [
            "ssh",
            "-o",
            "BatchMode=yes",
            "-o",
            "ClearAllForwardings=yes",
            "-o",
            "ConnectTimeout=10",
            "-o",
            "StrictHostKeyChecking=yes",
            host,
            f"bash -lc {shlex.quote(command)}",
        ]
        return cls._run(
            args=ssh_command,
            command_text=command,
            mode="ssh",
            host=host,
            timeout_s=timeout_s,
        )

    @staticmethod
    def _validate_local_command(command: Sequence[str]) -> None:
        if not command:
            raise ProbeError("Probe command cannot be empty.")
        token = command[0].strip()
        ProbeRunner._validate_token(token)

    @staticmethod
    def _validate_ssh_command(command: str) -> None:
        token = command.strip().split(maxsplit=1)[0] if command.strip() else ""
        ProbeRunner._validate_token(token)

    @staticmethod
    def _validate_token(token: str) -> None:
        if not token:
            raise ProbeError("Probe command cannot be empty.")
        if token not in ALLOWED_COMMANDS:
            raise ProbeError(f"Probe command is not allowed: {token}")

    @classmethod
    def _run(
        cls,
        args: list[str],
        command_text: str,
        mode: str,
        host: str | None,
        timeout_s: int,
    ) -> ProbeResult:
        start = time.monotonic()
        try:
            result = subprocess.run(
                args,
                capture_output=True,
                text=True,
                timeout=timeout_s,
                check=False,
            )
            stdout = cls._coerce_text(result.stdout)
            stderr = cls._coerce_text(result.stderr)
            returncode = result.returncode
        except subprocess.TimeoutExpired as exc:
            stdout = cls._coerce_text(getattr(exc, "stdout", None))
            stderr = cls._coerce_text(getattr(exc, "stderr", None))
            timeout_marker = "[probe-timeout]"
            stderr = (
                f"{stderr.rstrip()}\n{timeout_marker}"
                if stderr.strip()
                else timeout_marker
            )
            returncode = 124
        except FileNotFoundError as exc:
            stdout = ""
            stderr = str(exc)
            returncode = 127
        except OSError as exc:
            stdout = ""
            stderr = str(exc)
            returncode = 1

        stdout, stdout_truncated = cls._truncate_output(stdout)
        stderr, stderr_truncated = cls._truncate_output(stderr)
        duration_s = time.monotonic() - start
        return ProbeResult(
            command=command_text,
            mode=mode,
            host=host,
            returncode=returncode,
            stdout=stdout,
            stderr=stderr,
            duration_s=duration_s,
            truncated=stdout_truncated or stderr_truncated,
        )

    @staticmethod
    def _coerce_text(value: str | bytes | None) -> str:
        if value is None:
            return ""
        if isinstance(value, bytes):
            return value.decode("utf-8", errors="replace")
        return value

    @staticmethod
    def _truncate_output(text: str) -> tuple[str, bool]:
        if len(text) <= MAX_OUTPUT_BYTES:
            return text, False
        limit = MAX_OUTPUT_BYTES - len(TRUNCATION_MARKER)
        return f"{text[:limit]}{TRUNCATION_MARKER}", True
