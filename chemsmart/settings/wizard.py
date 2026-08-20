"""The pure core of ``chemsmart wizard``: derive, render, write, verify.

Interactivity lives in the CLI layer; everything here takes typed facts
and returns typed results, so the whole setup path is testable without a
terminal or a cluster. The scheduler is the resource authority: defaults
come from the chosen queue's advertised limits, never from the login
node's own core count.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import shutil

from chemsmart.settings.probe import (
    HostFactsV1,
    QueueFactsV1,
    SchedulerFactsV1,
)

#: Schedulers the submission layer supports end-to-end today.
SUBMITTABLE_SCHEDULERS = ("SLURM", "PBS")


@dataclass(frozen=True)
class ServerChoicesV1:
    """The confirmed SERVER facts that become the YAML block."""

    scheduler: str | None
    queue_name: str | None
    num_hours: int
    mem_gb: int
    num_cores: int
    num_gpus: int
    num_threads: int
    submit_command: str | None
    scratch_dir: str | None
    provenance: tuple[str, ...] = ()


def default_hours(queue: QueueFactsV1 | None) -> int:
    if queue is None or queue.max_time_seconds is None:
        return 24
    return max(1, min(queue.max_time_seconds // 3600 or 1, 24))


def default_mem_gb(
    queue: QueueFactsV1 | None, host: HostFactsV1
) -> int:
    """floor(0.9 x advertised kB / 1024^2); queue first, host fallback."""

    source_kb = None
    if queue is not None and queue.mem_kb_per_node:
        source_kb = queue.mem_kb_per_node
    elif host.mem_kb:
        source_kb = host.mem_kb
    if not source_kb:
        return 4
    return max(1, int(0.9 * source_kb / (1024 * 1024)))


def default_cores(queue: QueueFactsV1 | None, host: HostFactsV1) -> int:
    if queue is not None and queue.cores_per_node:
        return queue.cores_per_node
    return max(1, host.cpu_count or 1)


def derive_choices(
    *,
    scheduler: SchedulerFactsV1 | None,
    queue: QueueFactsV1 | None,
    host: HostFactsV1,
    scratch_dir: str | None,
) -> ServerChoicesV1:
    submittable = (
        scheduler is not None
        and scheduler.scheduler in SUBMITTABLE_SCHEDULERS
    )
    gpus = 0
    if queue is not None and queue.gres.startswith("gpu"):
        tail = queue.gres.rsplit(":", 1)[-1]
        gpus = int(tail) if tail.isdigit() else 0
    cores = default_cores(queue if submittable else None, host)
    provenance = []
    if scheduler is not None:
        provenance.append(
            f"detected {scheduler.scheduler}"
            + (f" ({scheduler.version})" if scheduler.version else "")
        )
        if not submittable:
            provenance.append(
                f"{scheduler.scheduler} submission is not release-qualified;"
                " this configuration runs locally"
            )
    else:
        provenance.append("no batch scheduler detected; local execution")
    if queue is not None and submittable:
        limit = (
            f"{queue.max_time_seconds // 3600} h"
            if queue.max_time_seconds
            else "no stated limit"
        )
        provenance.append(
            f"queue '{queue.name}': {queue.cores_per_node} cores, "
            f"{(queue.mem_kb_per_node or 0) // 1024} MB per node, "
            f"max time {limit}"
        )
    return ServerChoicesV1(
        scheduler=scheduler.scheduler if submittable else None,
        queue_name=queue.name if (queue is not None and submittable) else None,
        num_hours=default_hours(queue if submittable else None),
        mem_gb=default_mem_gb(queue if submittable else None, host),
        num_cores=cores,
        num_gpus=gpus,
        num_threads=cores,
        submit_command=(
            scheduler.submit_path or None if submittable else None
        )
        if scheduler is not None
        else None,
        scratch_dir=scratch_dir,
        provenance=tuple(provenance),
    )


def _yaml_scalar(value: object) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def render_server_block(choices: ServerChoicesV1, *, host: HostFactsV1) -> str:
    """The SERVER: block, with its provenance stated as comments."""

    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    lines = [
        f"# Written by `chemsmart wizard` on {host.hostname or 'this host'} "
        f"at {stamp}.",
    ]
    for note in choices.provenance:
        lines.append(f"# {note}")
    lines += [
        "SERVER:",
        f"    SCHEDULER: {_yaml_scalar(choices.scheduler)}",
        f"    QUEUE_NAME: {_yaml_scalar(choices.queue_name)}",
        f"    NUM_HOURS: {choices.num_hours}",
        f"    MEM_GB: {choices.mem_gb}",
        f"    NUM_CORES: {choices.num_cores}",
        f"    NUM_GPUS: {choices.num_gpus}",
        f"    NUM_THREADS: {choices.num_threads}",
        f"    SUBMIT_COMMAND: {_yaml_scalar(choices.submit_command)}",
        f"    SCRATCH_DIR: {_yaml_scalar(choices.scratch_dir)}",
        "    USE_HOSTS: false",
        "    EXTRA_COMMANDS: |",
        "        # Host commands required before execution.",
        "",
    ]
    return "\n".join(lines)


def extract_top_level_block(text: str, key: str) -> str:
    """The ``key:`` block of a YAML text, comments-with-it, or ""."""

    lines = text.splitlines(keepends=True)
    start = None
    for index, line in enumerate(lines):
        if line.startswith(f"{key}:"):
            start = index
            break
    if start is None:
        return ""
    end = len(lines)
    for index in range(start + 1, len(lines)):
        line = lines[index]
        if line[:1] not in ("", " ", "\t", "\n", "\r", "#"):
            end = index
            break
    return "".join(lines[start:end])


def splice_top_level_block(
    existing_text: str, key: str, block: str
) -> str:
    """Replace one top-level block, preserving every other byte.

    A record starts at the column-0 ``key:`` line and ends before the
    next column-0 key; comment lines immediately above the old block are
    treated as part of it (they described the old block).
    """

    lines = existing_text.splitlines(keepends=True)
    start = None
    for index, line in enumerate(lines):
        if line.startswith(f"{key}:"):
            start = index
            break
    if start is None:
        separator = "" if existing_text.startswith(("\n", "")) else "\n"
        return block + separator + existing_text
    # Absorb the contiguous comment lines directly above the old block.
    first = start
    while first > 0 and lines[first - 1].lstrip().startswith("#"):
        first -= 1
    end = len(lines)
    for index in range(start + 1, len(lines)):
        line = lines[index]
        if line[:1] not in ("", " ", "\t", "\n", "\r", "#"):
            end = index
            break
        if line.lstrip().startswith("#") and index + 1 < len(lines) and lines[
            index + 1
        ][:1] not in ("", " ", "\t", "\n", "\r", "#"):
            # Comments directly above the next block belong to it.
            end = index
            break
    return "".join(lines[:first]) + block + "".join(lines[end:])


def splice_server_block(existing_text: str, server_block: str) -> str:
    return splice_top_level_block(existing_text, "SERVER", server_block)


def backup_path_for(path: Path) -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S")
    return path.with_name(path.name + f".bak-{stamp}")


def write_server_yaml(
    path: Path,
    server_block: str,
    *,
    template_text: str,
) -> Path | None:
    """Write the server YAML; returns the backup path when one was made.

    A new file is the template with its SERVER block replaced (program
    skeletons intact); an existing file is backed up first and only its
    SERVER block is spliced -- hand-tuned program blocks survive
    byte-for-byte.
    """

    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    backup = None
    if path.exists():
        backup = backup_path_for(path)
        shutil.copy2(path, backup)
        merged = splice_server_block(
            path.read_text(encoding="utf-8"), server_block
        )
    else:
        merged = splice_server_block(template_text, server_block)
    path.write_text(merged, encoding="utf-8")
    path.chmod(0o600)
    return backup


@dataclass(frozen=True)
class CheckResultV1:
    name: str
    status: str  # ok | warn | fail | skipped
    detail: str = ""


def run_verification(
    choices: ServerChoicesV1,
    *,
    runner=None,
    which=None,
) -> tuple[CheckResultV1, ...]:
    """The aiida-style setup checks, tabulated and non-fatal.

    Adapted from `verdi computer test` (aiida-core, MIT): clean shell
    output, scheduler responsiveness, identity, a byte-compared scratch
    round trip, login-shell cost, plus submit-binary presence.
    """

    import getpass
    import os
    import subprocess
    import time as time_module

    runner = subprocess.run if runner is None else runner
    which = shutil.which if which is None else which

    def run(command):
        try:
            return runner(
                command,
                capture_output=True,
                text=True,
                timeout=30,
                check=False,
            )
        except Exception as exc:  # noqa: BLE001 - reported, never raised
            class _Failed:
                returncode = 127
                stdout = ""
                stderr = str(exc)

            return _Failed()

    checks: list[CheckResultV1] = []

    probe = run(["bash", "-lc", "echo -n"])
    if probe.returncode == 0 and not probe.stdout and not probe.stderr:
        checks.append(CheckResultV1("clean login shell", "ok"))
    else:
        spurious = (probe.stdout + probe.stderr).strip()[:120]
        checks.append(
            CheckResultV1(
                "clean login shell",
                "fail",
                "the login shell prints output that can corrupt parsed "
                f"scheduler replies: {spurious!r}",
            )
        )

    if choices.scheduler == "SLURM":
        probe = run(["sinfo", "--version"])
    elif choices.scheduler == "PBS":
        probe = run(["qstat", "-Q"])
    else:
        probe = None
    if probe is None:
        checks.append(
            CheckResultV1(
                "scheduler responds", "skipped", "local execution only"
            )
        )
    elif probe.returncode == 0:
        checks.append(CheckResultV1("scheduler responds", "ok"))
    else:
        checks.append(
            CheckResultV1(
                "scheduler responds",
                "fail",
                f"exit {probe.returncode}: {probe.stderr.strip()[:120]}",
            )
        )

    try:
        identity = getpass.getuser()
        checks.append(CheckResultV1("identity", "ok", identity))
    except Exception:  # noqa: BLE001
        checks.append(CheckResultV1("identity", "fail", "no user name"))

    if choices.scratch_dir:
        try:
            scratch = Path(choices.scratch_dir).expanduser()
            scratch.mkdir(parents=True, exist_ok=True)
            token = f"chemsmart-wizard-{os.getpid()}-{time_module.time()}"
            probe_file = scratch / f".{token}"
            payload = token.encode("utf-8")
            probe_file.write_bytes(payload)
            read_back = probe_file.read_bytes()
            probe_file.unlink()
            if read_back == payload:
                checks.append(
                    CheckResultV1("scratch round trip", "ok", str(scratch))
                )
            else:
                checks.append(
                    CheckResultV1(
                        "scratch round trip",
                        "fail",
                        "read bytes differ from written bytes",
                    )
                )
        except OSError as exc:
            checks.append(
                CheckResultV1("scratch round trip", "fail", str(exc)[:120])
            )
    else:
        checks.append(
            CheckResultV1(
                "scratch round trip", "skipped", "no scratch configured"
            )
        )

    def timed(command):
        started = time_module.monotonic()
        for _ in range(3):
            run(command)
        return (time_module.monotonic() - started) / 3

    login = timed(["bash", "-lc", "true"])
    plain = timed(["bash", "-c", "true"])
    if login > 2 * plain + 0.1:
        checks.append(
            CheckResultV1(
                "login-shell cost",
                "warn",
                f"login shell {login:.2f}s vs plain {plain:.2f}s -- heavy "
                "startup files slow every submitted job",
            )
        )
    else:
        checks.append(CheckResultV1("login-shell cost", "ok"))

    if choices.submit_command:
        first = choices.submit_command.split()[0]
        present = Path(first).exists() or bool(which(first))
        checks.append(
            CheckResultV1(
                "submit command present",
                "ok" if present else "fail",
                choices.submit_command,
            )
        )
    else:
        checks.append(
            CheckResultV1(
                "submit command present", "skipped", "local execution only"
            )
        )
    return tuple(checks)


__all__ = [
    "CheckResultV1",
    "SUBMITTABLE_SCHEDULERS",
    "ServerChoicesV1",
    "backup_path_for",
    "default_cores",
    "default_hours",
    "default_mem_gb",
    "derive_choices",
    "extract_top_level_block",
    "run_verification",
    "render_server_block",
    "splice_server_block",
    "splice_top_level_block",
    "write_server_yaml",
]
