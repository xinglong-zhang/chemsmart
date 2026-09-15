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


def _merge_partition_rows(
    rows: list[QueueFactsV1],
) -> QueueFactsV1:
    """One partition's rows folded into the partition a user can pick.

    ``sinfo`` prints one row per distinct node configuration, not one per
    partition, so a heterogeneous partition arrives here many times over.
    Treating each row as its own queue offers the same name repeatedly in
    one choice list and then binds whichever row happened to come first,
    which is an arbitrary node rather than a described one.

    The facts kept are those of the node configuration **most of the
    partition's nodes have**, because the wizard is choosing a default job
    shape and the shape that fits the most nodes is the one most likely to
    be scheduled; the largest node would default every job to a size only
    a handful of nodes can satisfy. ``node_count`` is the partition total,
    and a ``gres`` is kept only when every row agrees, so a partition that
    mixes accelerated and plain nodes never claims a GPU nobody asked for.
    """

    if len(rows) == 1:
        return rows[0]
    typical = max(
        rows, key=lambda row: (row.node_count or 0, row.cores_per_node or 0)
    )
    gres = {row.gres for row in rows}
    limits = [row.max_time_seconds for row in rows]
    return QueueFactsV1(
        name=typical.name,
        is_default=any(row.is_default for row in rows),
        available=any(row.available for row in rows),
        # A row with no limit makes the partition's ceiling unlimited.
        max_time_seconds=(
            None
            if any(limit is None for limit in limits)
            else max(limit for limit in limits)
        ),
        cores_per_node=typical.cores_per_node,
        mem_kb_per_node=typical.mem_kb_per_node,
        node_count=sum(row.node_count or 0 for row in rows) or None,
        gres=gres.pop() if len(gres) == 1 else "",
    )


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
    # One entry per partition, in the order the controller first named it.
    by_name: dict[str, list[QueueFactsV1]] = {}
    for queue in queues:
        by_name.setdefault(queue.name, []).append(queue)
    return tuple(_merge_partition_rows(rows) for rows in by_name.values())


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
