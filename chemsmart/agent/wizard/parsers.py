"""Scheduler output parsers for wizard queue discovery."""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass


@dataclass(frozen=True)
class QueueFacts:
    name: str
    default: bool
    max_walltime_hours: int | None
    default_walltime_hours: int | None
    default_mem_gb: int | None
    default_cores: int | None
    gpus_per_node: int | None
    enabled: bool
    started: bool


def parse_slurm_sinfo_json(payload: str) -> list[QueueFacts]:
    """Parse partitions from ``sinfo --json`` output."""

    data = json.loads(payload)
    queues: list[QueueFacts] = []
    for partition in data.get("partitions", []):
        name = (
            partition.get("name")
            or partition.get("partition")
            or partition.get("partition_name")
        )
        if not name:
            continue
        flags = [str(flag).lower() for flag in partition.get("flags", [])]
        state = str(
            partition.get("state") or partition.get("state_up") or "UP"
        )
        queues.append(
            QueueFacts(
                name=name,
                default=_as_bool(partition.get("default"))
                or "default" in flags,
                max_walltime_hours=_parse_slurm_json_time(
                    partition.get("maximums", {}).get("time")
                )
                or _parse_hours(partition.get("max_time")),
                default_walltime_hours=_parse_slurm_json_time(
                    partition.get("defaults", {}).get("time")
                )
                or _parse_hours(partition.get("default_time")),
                default_mem_gb=_parse_mem_gb(partition.get("default_mem")),
                default_cores=_parse_int(partition.get("default_cpus")),
                gpus_per_node=_parse_gpu_count(
                    partition.get("gres") or partition.get("tres")
                ),
                enabled=_is_active_state(state),
                started=_is_active_state(state),
            )
        )
    return queues


def parse_slurm_scontrol_partition_oneliner(payload: str) -> list[QueueFacts]:
    """Parse partitions from ``scontrol show partition --oneliner``."""

    queues: list[QueueFacts] = []
    for line in payload.splitlines():
        line = line.strip()
        if not line:
            continue
        fields = dict(re.findall(r"(\w+)=([^\s]+)", line))
        name = fields.get("PartitionName")
        if not name:
            continue
        default_cores = _first_int(
            fields.get("DefaultCPUsPerNode"),
            fields.get("MaxCPUsPerNode"),
            fields.get("TotalCPUs"),
        )
        default_mem_gb = _first_mem_gb(
            fields.get("DefMemPerNode"),
            _mem_times_cores(fields.get("DefMemPerCPU"), default_cores),
            _mem_times_cores(fields.get("MaxMemPerCPU"), default_cores),
        )
        state = fields.get("State", "UP")
        queues.append(
            QueueFacts(
                name=name,
                default=_as_bool(fields.get("Default")),
                max_walltime_hours=_parse_hours(fields.get("MaxTime")),
                default_walltime_hours=_parse_hours(fields.get("DefaultTime")),
                default_mem_gb=default_mem_gb,
                default_cores=default_cores,
                gpus_per_node=_parse_gpu_count(
                    fields.get("TRES") or fields.get("Gres")
                ),
                enabled=_is_active_state(state),
                started=_is_active_state(state),
            )
        )
    return queues


def parse_pbs_qstat_qf_json(payload: str) -> list[QueueFacts]:
    """Parse queues from ``qstat -Q -f -F json`` output."""

    data = json.loads(payload)
    queue_map = data.get("Queue") or data.get("Queues") or {}
    queues: list[QueueFacts] = []
    for name, fields in queue_map.items():
        default_cores = _first_int(
            fields.get("resources_default.ncpus"),
            fields.get("default_chunk.ncpus"),
        )
        queues.append(
            QueueFacts(
                name=name,
                default=_as_bool(
                    fields.get("default_queue")
                    or fields.get("resources_default.place") == "default"
                ),
                max_walltime_hours=_parse_hours(
                    fields.get("resources_max.walltime")
                ),
                default_walltime_hours=_parse_hours(
                    fields.get("resources_default.walltime")
                ),
                default_mem_gb=_parse_mem_gb(
                    fields.get("resources_default.mem")
                ),
                default_cores=default_cores,
                gpus_per_node=_first_int(
                    fields.get("resources_default.ngpus"),
                    fields.get("default_chunk.ngpus"),
                    fields.get("resources_available.ngpus"),
                ),
                enabled=_as_bool(fields.get("enabled"), default=True),
                started=_as_bool(fields.get("started"), default=True),
            )
        )
    return queues


def parse_pbs_qstat_qf_text(payload: str) -> list[QueueFacts]:
    """Parse queues from ``qstat -Q -f`` output."""

    blocks = re.split(r"(?m)^Queue\s+", payload)
    queues: list[QueueFacts] = []
    for block in blocks:
        block = block.strip()
        if not block:
            continue
        lines = block.splitlines()
        name = lines[0].strip()
        fields: dict[str, str] = {}
        for line in lines[1:]:
            match = re.match(r"\s*([\w.]+)\s*=\s*(.+)", line)
            if match:
                fields[match.group(1)] = match.group(2).strip()
        default_cores = _first_int(
            fields.get("resources_default.ncpus"),
            fields.get("default_chunk.ncpus"),
        )
        queues.append(
            QueueFacts(
                name=name,
                default=_as_bool(fields.get("default_queue")),
                max_walltime_hours=_parse_hours(
                    fields.get("resources_max.walltime")
                ),
                default_walltime_hours=_parse_hours(
                    fields.get("resources_default.walltime")
                ),
                default_mem_gb=_parse_mem_gb(
                    fields.get("resources_default.mem")
                ),
                default_cores=default_cores,
                gpus_per_node=_first_int(
                    fields.get("resources_default.ngpus"),
                    fields.get("default_chunk.ngpus"),
                    fields.get("resources_available.ngpus"),
                ),
                enabled=_as_bool(fields.get("enabled"), default=True),
                started=_as_bool(fields.get("started"), default=True),
            )
        )
    return queues


def parse_lsf_bqueues_l(payload: str) -> list[QueueFacts]:
    """Parse queues from ``bqueues -l`` output."""

    blocks = re.split(r"(?m)^-+\n", payload)
    queues: list[QueueFacts] = []
    for block in blocks:
        name_match = re.search(r"(?m)^QUEUE:\s*(\S+)", block)
        if not name_match:
            continue
        block = block.strip()
        name = name_match.group(1)
        status = " ".join(
            re.findall(r"(?:Open|Closed|Active|Inactive)", block)
        )
        enabled = "Closed" not in status
        started = "Inactive" not in status
        walltime = _parse_lsf_hours(block)
        queues.append(
            QueueFacts(
                name=name,
                default=_as_bool(
                    _search_group(block, r"DEFAULT_QUEUE\s*=\s*(\S+)")
                ),
                max_walltime_hours=walltime,
                default_walltime_hours=walltime,
                default_mem_gb=_parse_mem_gb(
                    _search_group(block, r"MEMLIMIT\s*=\s*([^\n]+)")
                ),
                default_cores=_first_int(
                    _search_group(block, r"PROCLIMIT\s*=\s*(\d+)"),
                    _search_group(
                        block, r"DEFAULT HOST SPECIFICATION:\s*(\d+)"
                    ),
                ),
                gpus_per_node=_parse_gpu_count(block),
                enabled=enabled,
                started=started,
            )
        )
    return queues


def parse_sge_qhost(payload: str) -> tuple[int | None, int | None]:
    """Parse ``qhost`` output and return (min_mem_gb, min_ncpu).

    Column layout (1-based):
      1=HOSTNAME 2=ARCH 3=NCPU 4=NSOC 5=NCOR 6=NTHR 7=LOAD 8=MEMTOT ...

    Returns the minimum values across all nodes so that the rendered YAML is
    conservative and safe for every queue / node type.
    """
    mem_values: list[int] = []
    cpu_values: list[int] = []
    for line in payload.splitlines():
        parts = line.split()
        # data rows have at least 8 columns; skip header and global
        if len(parts) < 8:
            continue
        hostname = parts[0]
        if hostname in ("HOSTNAME", "global", "-"):
            continue
        ncpu_raw = parts[2]
        mem_raw = parts[7]
        if ncpu_raw.isdigit():
            cpu_values.append(int(ncpu_raw))
        parsed_mem = _parse_mem_gb(mem_raw)
        if parsed_mem is not None:
            mem_values.append(parsed_mem)
    min_mem = min(mem_values) if mem_values else None
    min_cpu = min(cpu_values) if cpu_values else None
    return min_mem, min_cpu


def parse_sge_qconf_sql(payload: str) -> list[str]:
    """Parse queue names from ``qconf -sql`` output."""

    return [
        line.strip()
        for line in payload.splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]


def parse_sge_qconf_sq(payload: str) -> QueueFacts:
    """Parse one queue from ``qconf -sq <queue>`` output."""

    fields: dict[str, str] = {}
    for line in payload.splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split(None, 1)
        key = parts[0]
        value = parts[1].strip() if len(parts) > 1 else ""
        fields[key] = value
    name = fields.get("qname")
    if not name:
        raise ValueError("SGE queue configuration is missing qname.")
    state = fields.get("state", "enabled")
    complex_values = fields.get("complex_values", "")
    return QueueFacts(
        name=name,
        default=False,
        max_walltime_hours=_parse_sge_seconds(fields.get("h_rt")),
        default_walltime_hours=_parse_sge_seconds(fields.get("s_rt")),
        default_mem_gb=_first_mem_gb(
            fields.get("h_vmem"),
            _search_group(complex_values, r"h_vmem=([^,\s]+)"),
        ),
        default_cores=_parse_int(fields.get("slots")),
        gpus_per_node=_parse_gpu_count(complex_values),
        enabled="disabled" not in state.lower(),
        started="disabled" not in state.lower(),
    )


def _as_bool(value, default: bool = False) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if not text:
        return default
    return text in {"1", "true", "yes", "y", "on", "open", "active"}


def _first_int(*values) -> int | None:
    for value in values:
        parsed = _parse_int(value)
        if parsed is not None:
            return parsed
    return None


def _parse_int(value) -> int | None:
    if value is None:
        return None
    match = re.search(r"-?\d+", str(value))
    if not match:
        return None
    return int(match.group())


def _first_mem_gb(*values) -> int | None:
    for value in values:
        parsed = _parse_mem_gb(value)
        if parsed is not None:
            return parsed
    return None


def _mem_times_cores(value, cores: int | None) -> str | None:
    if value is None or cores is None:
        return None
    parsed = _parse_mem_gb(value)
    if parsed is None:
        return None
    return str(parsed * cores) + "G"


def _parse_mem_gb(value) -> int | None:
    if value is None:
        return None
    text = str(value).strip().replace(" ", "")
    if not text or text.upper() in {"UNLIMITED", "INFINITE", "INFINITY"}:
        return None
    match = re.search(r"(\d+(?:\.\d+)?)([KMGTP]B?|[KMGTP])?", text, re.I)
    if not match:
        return None
    number = float(match.group(1))
    unit = (match.group(2) or "M").upper().rstrip("B")
    scale = {
        "K": 1 / (1024 * 1024),
        "M": 1 / 1024,
        "G": 1,
        "T": 1024,
        "P": 1024 * 1024,
    }
    return max(1, math.ceil(number * scale.get(unit, 1 / 1024)))


def _parse_hours(value) -> int | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    upper = text.upper()
    if upper in {"UNLIMITED", "INFINITE", "INFINITY", "NOT_SET", "NONE", "-"}:
        return None
    if re.fullmatch(r"\d+-\d{2}:\d{2}:\d{2}", text):
        days, remainder = text.split("-", 1)
        hours, minutes, seconds = map(int, remainder.split(":"))
        total_hours = int(days) * 24 + hours
        if minutes or seconds:
            total_hours += 1
        return total_hours
    if re.fullmatch(r"\d{1,3}:\d{2}:\d{2}", text):
        hours, minutes, seconds = map(int, text.split(":"))
        return max(1, hours + int(bool(minutes or seconds)))
    if re.fullmatch(r"\d+(?:\.\d+)?h", text, re.I):
        return max(1, math.ceil(float(text[:-1])))
    if re.fullmatch(r"\d+", text):
        return int(text)
    return None


def _parse_sge_seconds(value) -> int | None:
    if value is None:
        return None
    text = str(value).strip()
    if re.fullmatch(r"\d+", text):
        seconds = int(text)
        return max(1, math.ceil(seconds / 3600))
    return _parse_hours(text)


def _parse_lsf_hours(block: str) -> int | None:
    value = _search_group(block, r"RUNLIMIT\s*=\s*([^\n]+)")
    if value is None:
        return None
    value = value.strip()
    if re.fullmatch(r"\d+(?:\.\d+)?\s*h", value, re.I):
        return max(1, math.ceil(float(value.split()[0])))
    return _parse_hours(value)


def _search_group(text: str, pattern: str) -> str | None:
    match = re.search(pattern, text, re.I)
    if not match:
        return None
    return match.group(1)


def _parse_gpu_count(value) -> int | None:
    if value is None:
        return None
    text = str(value)
    matches = [
        int(match)
        for match in re.findall(
            r"(?:gres/gpu|gpu(?:s)?|ngpus|num)\s*[=:]\s*(\d+)",
            text,
            re.I,
        )
    ]
    if matches:
        return max(matches)
    return None


def _parse_slurm_json_time(value) -> int | None:
    if not isinstance(value, dict):
        return None
    if value.get("infinite"):
        return None
    number = value.get("number")
    if isinstance(number, int):
        return max(1, math.ceil(number / 60))
    return None


def _is_active_state(value: str) -> bool:
    upper = str(value).upper()
    return not any(flag in upper for flag in ["DOWN", "DRAIN", "INACTIVE"])
