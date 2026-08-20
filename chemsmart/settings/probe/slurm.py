"""SLURM discovery: partitions as the controller advertises them.

aiida-core never queries partitions; this half is chemsmart's own. The
format string is pinned (never default human output), so the parser is a
pure function over a stable shape.
"""

from __future__ import annotations

import shutil
import subprocess
from typing import Callable

from .facts import QueueFactsV1, SchedulerFactsV1
from .units import ProbeUnitError, parse_memory_kb, parse_walltime_seconds

#: partition|avail|timelimit|cpus-per-node|memory-MB|nodes|gres
_SINFO_FORMAT = "%P|%a|%l|%c|%m|%D|%G"


def sinfo_command() -> tuple[str, ...]:
    return ("sinfo", "-h", "-o", _SINFO_FORMAT)


def sinfo_version_command() -> tuple[str, ...]:
    return ("sinfo", "--version")


def _int_or_none(text: str) -> int | None:
    cleaned = text.strip().rstrip("+")
    try:
        return int(cleaned)
    except ValueError:
        return None


def parse_sinfo(
    returncode: int, stdout: str, stderr: str
) -> tuple[QueueFactsV1, ...]:
    if returncode != 0:
        raise ProbeUnitError(
            f"sinfo exited {returncode}: {stderr.strip() or 'no output'}"
        )
    queues: list[QueueFactsV1] = []
    for line in stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split("|")
        if len(parts) != 7:
            continue  # an unexpected row is skipped, never a crash
        name, avail, limit, cpus, memory, nodes, gres = (
            part.strip() for part in parts
        )
        is_default = name.endswith("*")
        try:
            max_time = parse_walltime_seconds(limit)
        except ProbeUnitError:
            max_time = None
        try:
            mem_kb = parse_memory_kb(memory, default_unit="m")
        except ProbeUnitError:
            mem_kb = None
        queues.append(
            QueueFactsV1(
                name=name.rstrip("*"),
                is_default=is_default,
                available=avail.lower() == "up",
                max_time_seconds=max_time,
                cores_per_node=_int_or_none(cpus),
                mem_kb_per_node=mem_kb,
                node_count=_int_or_none(nodes),
                gres="" if gres in {"(null)", ""} else gres,
            )
        )
    return tuple(queues)


def parse_sinfo_version(returncode: int, stdout: str, stderr: str) -> str:
    if returncode != 0:
        return ""
    return stdout.strip().splitlines()[0] if stdout.strip() else ""


def probe_slurm(
    *,
    runner: Callable[..., "subprocess.CompletedProcess[str]"] = subprocess.run,
    which: Callable[[str], str | None] = shutil.which,
) -> SchedulerFactsV1:
    version_result = runner(
        list(sinfo_version_command()),
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    queues_result = runner(
        list(sinfo_command()),
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    return SchedulerFactsV1(
        scheduler="SLURM",
        version=parse_sinfo_version(
            version_result.returncode,
            version_result.stdout,
            version_result.stderr,
        ),
        submit_path=which("sbatch") or "",
        queues=parse_sinfo(
            queues_result.returncode,
            queues_result.stdout,
            queues_result.stderr,
        ),
        evidence=(f"sinfo -o '{_SINFO_FORMAT}'",),
    )


__all__ = [
    "parse_sinfo",
    "parse_sinfo_version",
    "probe_slurm",
    "sinfo_command",
    "sinfo_version_command",
]
