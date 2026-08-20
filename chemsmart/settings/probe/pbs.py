"""PBS discovery: queue limits from qstat -Qf, node capability from pbsnodes.

The stanza-folding approach is adapted from aiida-core (MIT,
src/aiida/schedulers/plugins/pbsbaseclasses.py): a new record starts at an
unindented header line, indented ``key = value`` lines belong to it, and
deeper-indented lines continue the previous value.
"""

from __future__ import annotations

import shutil
import subprocess
from typing import Callable

from .facts import QueueFactsV1, SchedulerFactsV1
from .units import ProbeUnitError, parse_memory_kb, parse_walltime_seconds


def qstat_queues_command() -> tuple[str, ...]:
    return ("qstat", "-Q", "-f")


def pbsnodes_command() -> tuple[str, ...]:
    return ("pbsnodes", "-a")


def _stanzas(stdout: str, header: str) -> list[tuple[str, dict[str, str]]]:
    """Fold ``Header: name`` / bare-name stanzas into (name, fields) records.

    With a non-empty ``header``, a record starts at a line carrying that
    prefix; with ``header == ""`` (pbsnodes), any unindented non-empty line
    is a record name. Indented ``key = value`` lines belong to the current
    record; an indented line without ``=`` continues the previous value.
    """

    records: list[tuple[str, dict[str, str]]] = []
    name = ""
    fields: dict[str, str] = {}
    last_key = ""

    def is_header(raw: str) -> bool:
        if header:
            return raw.startswith(header)
        return bool(raw.strip()) and not raw[:1].isspace()

    for raw in stdout.splitlines():
        if is_header(raw):
            if name:
                records.append((name, fields))
            name = raw[len(header):].strip() if header else raw.strip()
            fields = {}
            last_key = ""
            continue
        if not name:
            continue
        stripped = raw.strip()
        if not stripped:
            continue
        if "=" in stripped:
            key, _, value = stripped.partition("=")
            last_key = key.strip().lower()
            fields[last_key] = value.strip()
        elif last_key:
            fields[last_key] = fields[last_key] + stripped
    if name:
        records.append((name, fields))
    return records


def parse_qstat_queues(
    returncode: int, stdout: str, stderr: str
) -> tuple[QueueFactsV1, ...]:
    if returncode != 0:
        raise ProbeUnitError(
            f"qstat -Q -f exited {returncode}: {stderr.strip() or 'no output'}"
        )
    queues: list[QueueFactsV1] = []
    for name, fields in _stanzas(stdout, "Queue:"):
        try:
            max_time = parse_walltime_seconds(
                fields.get("resources_max.walltime", "")
            )
        except ProbeUnitError:
            max_time = None
        try:
            mem_kb = parse_memory_kb(fields.get("resources_max.mem", ""))
        except ProbeUnitError:
            mem_kb = None
        cores = None
        raw_cores = fields.get("resources_max.ncpus", "").strip()
        if raw_cores.isdigit():
            cores = int(raw_cores)
        enabled = fields.get("enabled", "true").strip().lower() == "true"
        started = fields.get("started", "true").strip().lower() == "true"
        queues.append(
            QueueFactsV1(
                name=name,
                is_default=fields.get("queue_type", "").strip().lower()
                == "execution"
                and not queues,
                available=enabled and started,
                max_time_seconds=max_time,
                cores_per_node=cores,
                mem_kb_per_node=mem_kb,
            )
        )
    return tuple(queues)


def parse_pbsnodes(
    returncode: int, stdout: str, stderr: str
) -> tuple[int | None, int | None, int]:
    """-> (max cores per node, max memory kB per node, node count)."""

    if returncode != 0:
        return (None, None, 0)
    cores: list[int] = []
    mem_kb: list[int] = []
    count = 0
    for name, fields in _stanzas(stdout, ""):
        if not name or "=" in name:
            continue
        count += 1
        raw_cores = fields.get(
            "resources_available.ncpus", fields.get("np", "")
        ).strip()
        if raw_cores.isdigit():
            cores.append(int(raw_cores))
        raw_mem = fields.get(
            "resources_available.mem", fields.get("physmem", "")
        ).strip()
        if raw_mem:
            try:
                parsed = parse_memory_kb(raw_mem)
            except ProbeUnitError:
                parsed = None
            if parsed:
                mem_kb.append(parsed)
    return (
        max(cores) if cores else None,
        max(mem_kb) if mem_kb else None,
        count,
    )


def probe_pbs(
    *,
    runner: Callable[..., "subprocess.CompletedProcess[str]"] = subprocess.run,
    which: Callable[[str], str | None] = shutil.which,
) -> SchedulerFactsV1:
    queues_result = runner(
        list(qstat_queues_command()),
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    queues = parse_qstat_queues(
        queues_result.returncode, queues_result.stdout, queues_result.stderr
    )
    nodes_result = runner(
        list(pbsnodes_command()),
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    node_cores, node_mem_kb, node_count = parse_pbsnodes(
        nodes_result.returncode, nodes_result.stdout, nodes_result.stderr
    )
    # A queue that states no limits inherits the observed node capability.
    filled = tuple(
        QueueFactsV1(
            name=queue.name,
            is_default=queue.is_default,
            available=queue.available,
            max_time_seconds=queue.max_time_seconds,
            cores_per_node=queue.cores_per_node or node_cores,
            mem_kb_per_node=queue.mem_kb_per_node or node_mem_kb,
            node_count=queue.node_count or (node_count or None),
        )
        for queue in queues
    )
    return SchedulerFactsV1(
        scheduler="PBS",
        version="",
        submit_path=which("qsub") or "",
        queues=filled,
        evidence=("qstat -Q -f", "pbsnodes -a"),
    )


__all__ = [
    "parse_pbsnodes",
    "parse_qstat_queues",
    "pbsnodes_command",
    "probe_pbs",
    "qstat_queues_command",
]
