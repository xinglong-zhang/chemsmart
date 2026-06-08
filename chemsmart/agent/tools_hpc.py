"""HPC read-only tool wrappers."""

from __future__ import annotations

import json
import re
import shlex
from dataclasses import asdict
from typing import Any

from chemsmart.agent.error_summary import summarize_log
from chemsmart.agent.transport import ExecTransport, SshExecTransport
from chemsmart.agent.wizard.normalize import choose_queue
from chemsmart.agent.wizard.parsers import (
    parse_duration_seconds,
    parse_pbs_qstat_qf_json,
    parse_sge_qstat_gc,
    parse_slurm_scontrol_partition_oneliner,
)
from chemsmart.agent.wizard.probe import (
    ALL_PROBE_SPECS,
    ProbeError,
    resolve_probe_argv,
)

_OUTPUT_LIMIT = 4096
_LOG_TAIL_OUTPUT_LIMIT = 8192
_TRANSPORT_FACTORY = SshExecTransport
_SUPPORTED_SCHEDULERS = ("slurm", "pbs", "sge", "lsf")
_JOB_QUERY_PROBES = {
    "slurm": "query.slurm.squeue_job",
    "pbs": "query.pbs.qstat_job_json",
    "sge": "query.sge.qstat_job",
    "lsf": "query.lsf.bjobs_job_json",
}
_QUEUE_QUERY_PROBES = {
    "slurm": "survey.slurm.scontrol_partition",
    "pbs": "survey.pbs.qstat_json",
    "sge": "survey.sge.qstat_gc",
    "lsf": "survey.lsf.bqueues_json",
}
_DETECTION_PROBES = {
    "slurm": ("survey.slurm.scontrol_partition",),
    "pbs": ("survey.pbs.qstat_json", "survey.pbs.qstat_text"),
    "sge": ("survey.sge.qconf_sql",),
    "lsf": ("survey.lsf.bqueues_json",),
}
_PBS_STATE_MAP = {
    "B": "ARRAY_STARTED",
    "E": "EXITING",
    "F": "FINISHED",
    "H": "HELD",
    "M": "MOVED",
    "Q": "QUEUED",
    "R": "RUNNING",
    "S": "SUSPENDED",
    "T": "TRANSIT",
    "U": "USER_SUSPENDED",
    "W": "WAITING",
    "X": "FINISHED",
}
_LSF_STATE_MAP = {
    "DONE": "COMPLETED",
    "EXIT": "FAILED",
    "PEND": "PENDING",
    "PSUSP": "SUSPENDED",
    "RUN": "RUNNING",
    "SSUSP": "SUSPENDED",
    "UNKWN": "UNKNOWN",
    "USUSP": "SUSPENDED",
    "WAIT": "WAITING",
    "ZOMBI": "FAILED",
}


def ssh_probe(
    server: str,
    probe_name: str,
    timeout_s: int = 15,
) -> dict[str, Any]:
    """Run a predefined read-only probe on a remote HPC server."""

    if probe_name not in ALL_PROBE_SPECS:
        return {
            "error": "unknown_probe",
            "available": sorted(ALL_PROBE_SPECS),
        }

    spec = ALL_PROBE_SPECS[probe_name]
    transport = _make_transport()
    try:
        result = transport.exec(spec.command, server, timeout_s)
    except TimeoutError:
        return {
            "error": "timeout",
            "server": server,
            "probe": probe_name,
        }
    except ProbeError as exc:
        return {
            "error": "probe_requires_slots",
            "message": str(exc),
            "server": server,
            "probe": probe_name,
        }

    parsed = None
    if spec.parser is not None and result.stdout.strip():
        try:
            parsed = spec.parser(result.stdout)
        except Exception:
            parsed = None
    return {
        "server": server,
        "probe": probe_name,
        "returncode": result.returncode,
        "stdout_truncated": _truncate(result.stdout),
        "stderr_truncated": _truncate(result.stderr),
        "parsed": parsed,
        "duration_s": result.duration_s,
    }


def _make_transport() -> ExecTransport:
    return _TRANSPORT_FACTORY()


def scheduler_query(
    server: str,
    scheduler: str | None = None,
    job_id: str | None = None,
    timeout_s: int = 15,
) -> dict[str, Any]:
    """Inspect HPC scheduler queue state or one job through probe specs."""

    transport = _make_transport()
    scheduler_name = _normalize_scheduler_name(scheduler)
    if scheduler is not None and scheduler_name is None:
        return {
            "error": "unknown_scheduler",
            "supported": list(_SUPPORTED_SCHEDULERS),
        }

    if scheduler_name is None:
        detected = _detect_scheduler_kind(transport, server, timeout_s)
        if "error" in detected:
            return detected
        scheduler_name = detected["scheduler"]

    if job_id:
        probe_name = _JOB_QUERY_PROBES.get(scheduler_name)
        if probe_name is None:
            return {
                "error": "unknown_probe",
                "scheduler": scheduler_name,
            }
        outcome = _exec_probe(
            transport,
            server,
            probe_name,
            timeout_s,
            job_id=job_id,
        )
        if "error" in outcome:
            return outcome
        return _normalize_job_result(
            scheduler_name,
            outcome["parsed"],
        )

    probe_name = _QUEUE_QUERY_PROBES.get(scheduler_name)
    if probe_name is None:
        return {
            "error": "unknown_probe",
            "scheduler": scheduler_name,
        }
    outcome = _exec_probe(transport, server, probe_name, timeout_s)
    if "error" in outcome:
        return outcome
    return _normalize_queue_result(
        scheduler_name,
        outcome["parsed"],
        outcome["stdout"],
    )


def log_tail(
    server: str,
    path: str,
    lines: int = 200,
    grep: str | None = None,
    timeout_s: int = 15,
) -> dict[str, Any]:
    """Tail a remote log file and summarize common HPC error signatures."""

    if lines < 1 or lines > 10000:
        return {"error": "invalid_lines"}

    grep_pattern = _compile_grep(grep)
    if grep is not None and grep_pattern is None:
        return {"error": "invalid_grep"}

    transport = _make_transport()
    outcome = _exec_probe(
        transport,
        server,
        "query.log.tail",
        timeout_s,
        lines=str(lines),
        path=path,
    )
    if "error" in outcome:
        return outcome

    content = outcome["stdout"]
    filtered_lines = content.splitlines()
    if grep_pattern is not None:
        filtered_lines = [
            line
            for line in filtered_lines
            if grep_pattern.search(line) is not None
        ]
    filtered_text = "\n".join(filtered_lines)

    return {
        "server": server,
        "path": path,
        "lines_requested": lines,
        "lines_returned": len(filtered_lines),
        "content_truncated": _truncate(filtered_text, _LOG_TAIL_OUTPUT_LIMIT),
        "errors": [
            asdict(signature) for signature in summarize_log(filtered_text)
        ],
        "duration_s": outcome["duration_s"],
    }


def _truncate(text: str, limit: int = _OUTPUT_LIMIT) -> str:
    return text[:limit]


def _normalize_scheduler_name(scheduler: str | None) -> str | None:
    if scheduler is None:
        return None
    normalized = scheduler.strip().lower()
    if normalized in _SUPPORTED_SCHEDULERS:
        return normalized
    return None


def _detect_scheduler_kind(
    transport: ExecTransport,
    server: str,
    timeout_s: int,
) -> dict[str, Any]:
    matches: list[str] = []
    for scheduler, probe_names in _DETECTION_PROBES.items():
        for probe_name in probe_names:
            outcome = _exec_probe(transport, server, probe_name, timeout_s)
            if "error" in outcome:
                continue
            if outcome["returncode"] == 0 and outcome["stdout"].strip():
                matches.append(scheduler)
                break
    matches = sorted(set(matches))
    if len(matches) == 1:
        return {"scheduler": matches[0]}
    if len(matches) > 1:
        return {
            "error": "ambiguous_scheduler",
            "detected": matches,
            "supported": list(_SUPPORTED_SCHEDULERS),
        }
    return {
        "error": "scheduler_required",
        "supported": list(_SUPPORTED_SCHEDULERS),
    }


def _exec_probe(
    transport: ExecTransport,
    server: str,
    probe_name: str,
    timeout_s: int,
    **slots: str,
) -> dict[str, Any]:
    spec = ALL_PROBE_SPECS.get(probe_name)
    if spec is None:
        return {"error": "unknown_probe", "probe": probe_name}
    try:
        argv = (
            resolve_probe_argv(spec, **slots)
            if slots
            else resolve_probe_argv(spec)
        )
    except ProbeError as exc:
        return {
            "error": _probe_slot_error_name(slots, str(exc)),
            "message": str(exc),
            "probe": probe_name,
            "server": server,
        }
    command = shlex.join(argv)
    try:
        result = transport.exec(command, server, timeout_s)
    except TimeoutError:
        return {
            "error": "timeout",
            "server": server,
            "probe": probe_name,
            "scheduler": _scheduler_from_probe_name(probe_name),
        }
    if result.returncode != 0:
        return {
            "error": "probe_failed",
            "server": server,
            "probe": probe_name,
            "scheduler": _scheduler_from_probe_name(probe_name),
            "returncode": result.returncode,
            "stderr_truncated": _truncate(result.stderr),
        }

    parsed = None
    if spec.parser is not None and result.stdout.strip():
        try:
            parsed = spec.parser(result.stdout)
        except Exception as exc:
            return {
                "error": "parse_error",
                "message": str(exc),
                "server": server,
                "probe": probe_name,
                "scheduler": _scheduler_from_probe_name(probe_name),
                "stdout_truncated": _truncate(result.stdout),
            }

    return {
        "probe": probe_name,
        "returncode": result.returncode,
        "stdout": result.stdout,
        "parsed": parsed,
        "duration_s": result.duration_s,
    }


def _compile_grep(grep: str | None) -> re.Pattern[str] | None:
    if grep is None:
        return None
    if len(grep) > 128 or "\n" in grep or "\r" in grep:
        return None
    try:
        return re.compile(grep, re.IGNORECASE)
    except re.error:
        return None


def _probe_slot_error_name(
    slots: dict[str, str],
    message: str,
) -> str:
    lowered = message.lower()
    for slot_name in slots:
        if f"{slot_name} slot" in lowered or f"{slot_name!r}" in lowered:
            return f"invalid_{slot_name}"
    if len(slots) == 1:
        return f"invalid_{next(iter(slots))}"
    if "job_id" in slots:
        return "invalid_job_id"
    return "probe_error"


def _normalize_job_result(
    scheduler: str,
    raw: Any,
) -> dict[str, Any]:
    if scheduler == "slurm":
        return _normalize_slurm_job(raw)
    if scheduler == "pbs":
        return _normalize_pbs_job(raw)
    if scheduler == "sge":
        return _normalize_sge_job(raw)
    return _normalize_lsf_job(raw)


def _normalize_queue_result(
    scheduler: str,
    raw: Any,
    stdout: str,
) -> dict[str, Any]:
    if scheduler == "slurm":
        return _normalize_slurm_queue(stdout)
    if scheduler == "pbs":
        return _normalize_pbs_queue(raw)
    if scheduler == "sge":
        return _normalize_sge_queue(stdout)
    return _normalize_lsf_queue(raw)


def _normalize_slurm_job(raw: Any) -> dict[str, Any]:
    payload = dict(raw or {})
    return {
        "job_id": payload.get("job_id"),
        "state": _normalize_state("slurm", payload.get("state")),
        "queue": payload.get("queue"),
        "user": payload.get("user"),
        "name": payload.get("name"),
        "runtime_s": parse_duration_seconds(payload.get("runtime")),
        "node": _clean_node_name(payload.get("node_or_reason")),
        "scheduler": "slurm",
        "raw": raw,
    }


def _normalize_pbs_job(raw: Any) -> dict[str, Any]:
    jobs = (raw or {}).get("Jobs") or (raw or {}).get("jobs") or {}
    job_id, payload = next(iter(jobs.items()), (None, {}))
    return {
        "job_id": job_id or payload.get("Job_Id"),
        "state": _normalize_state("pbs", payload.get("job_state")),
        "queue": payload.get("queue"),
        "user": _before_separator(payload.get("Job_Owner"), "@"),
        "name": payload.get("Job_Name"),
        "runtime_s": parse_duration_seconds(
            payload.get("resources_used.walltime")
            or payload.get("Resource_List.walltime")
        ),
        "node": _parse_pbs_exec_host(payload.get("exec_host")),
        "scheduler": "pbs",
        "raw": raw,
    }


def _normalize_sge_job(raw: Any) -> dict[str, Any]:
    payload = dict(raw or {})
    queue_value = (
        payload.get("master_queue")
        or payload.get("queue_name")
        or payload.get("hard_queue_list")
        or payload.get("sge_o_queue")
    )
    return {
        "job_id": payload.get("job_number"),
        "state": _normalize_state(
            "sge",
            payload.get("state") or payload.get("job_state"),
        ),
        "queue": _before_separator(queue_value, "@"),
        "user": payload.get("owner"),
        "name": payload.get("job_name"),
        "runtime_s": parse_duration_seconds(
            payload.get("running_time")
            or payload.get("runtime")
            or payload.get("wallclock")
        ),
        "node": _after_separator(queue_value, "@"),
        "scheduler": "sge",
        "raw": raw,
    }


def _normalize_lsf_job(raw: Any) -> dict[str, Any]:
    records = (raw or {}).get("RECORDS") or (raw or {}).get("records") or []
    payload = records[0] if records else {}
    return {
        "job_id": _string_or_none(payload.get("JOBID")),
        "state": _normalize_state("lsf", payload.get("STAT")),
        "queue": payload.get("QUEUE"),
        "user": payload.get("USER"),
        "name": payload.get("JOB_NAME"),
        "runtime_s": parse_duration_seconds(
            payload.get("RUN_TIME") or payload.get("RUNTIME")
        ),
        "node": _parse_lsf_exec_host(payload.get("EXEC_HOST")),
        "scheduler": "lsf",
        "raw": raw,
    }


def _normalize_slurm_queue(stdout: str) -> dict[str, Any]:
    queues = parse_slurm_scontrol_partition_oneliner(stdout)
    chosen = choose_queue(queues) if queues else None
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    fields_by_name = {
        match.get("PartitionName"): match
        for match in (
            dict(re.findall(r"(\w+)=([^\s]+)", line)) for line in lines
        )
        if match.get("PartitionName")
    }
    selected = fields_by_name.get(chosen) if chosen else None
    if selected is None and fields_by_name:
        selected = next(iter(fields_by_name.values()))
        chosen = selected.get("PartitionName")
    return {
        "scheduler": "slurm",
        "partition_or_queue": chosen,
        "nodes": _parse_int(selected.get("TotalNodes")) if selected else None,
        "total_cpus": (
            _parse_int(selected.get("TotalCPUs")) if selected else None
        ),
        "raw": {
            "stdout": stdout,
            "queues": [queue.__dict__ for queue in queues],
        },
    }


def _normalize_pbs_queue(raw: Any) -> dict[str, Any]:
    queues = parse_pbs_qstat_qf_json(_stringify_json(raw))
    chosen = choose_queue(queues) if queues else None
    queue_map = (raw or {}).get("Queue") or (raw or {}).get("Queues") or {}
    payload = queue_map.get(chosen) if chosen else None
    if payload is None and queue_map:
        chosen, payload = next(iter(queue_map.items()))
    selected_queue = next(
        (queue for queue in queues if queue.name == chosen),
        None,
    )
    return {
        "scheduler": "pbs",
        "partition_or_queue": chosen,
        "nodes": _parse_int(
            (payload or {}).get("resources_assigned.nodect")
            or (payload or {}).get("resources_available.nodect")
        ),
        "total_cpus": (
            selected_queue.default_cores
            if selected_queue is not None
            else None
        ),
        "raw": raw,
    }


def _normalize_sge_queue(stdout: str) -> dict[str, Any]:
    totals = parse_sge_qstat_gc(stdout)
    queue_name, total_cpus = next(iter(totals.items()), (None, None))
    return {
        "scheduler": "sge",
        "partition_or_queue": queue_name,
        "nodes": None,
        "total_cpus": total_cpus,
        "raw": totals,
    }


def _normalize_lsf_queue(raw: Any) -> dict[str, Any]:
    records = (raw or {}).get("RECORDS") or (raw or {}).get("records") or []
    payload = records[0] if records else {}
    return {
        "scheduler": "lsf",
        "partition_or_queue": payload.get("QUEUE_NAME")
        or payload.get("QUEUE"),
        "nodes": None,
        "total_cpus": _parse_int(
            payload.get("PROCLIMIT") or payload.get("NJOBS")
        ),
        "raw": raw,
    }


def _normalize_state(scheduler: str, value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip().upper()
    if not text:
        return None
    if scheduler == "pbs":
        return _PBS_STATE_MAP.get(text, text)
    if scheduler == "lsf":
        return _LSF_STATE_MAP.get(text, text)
    if scheduler == "sge":
        return {
            "R": "RUNNING",
            "Q": "QUEUED",
            "H": "HELD",
            "S": "SUSPENDED",
            "T": "TRANSFERRING",
            "E": "ERROR",
            "D": "DELETING",
        }.get(text, text)
    return text


def _scheduler_from_probe_name(probe_name: str) -> str | None:
    parts = probe_name.split(".")
    if len(parts) >= 3:
        return parts[1]
    return None


def _parse_int(value: Any) -> int | None:
    if value is None:
        return None
    match = re.search(r"-?\d+", str(value))
    if match is None:
        return None
    return int(match.group())


def _parse_pbs_exec_host(value: Any) -> str | None:
    first = _string_or_none(value)
    if first is None:
        return None
    return _before_separator(_before_separator(first, "+"), "/")


def _parse_lsf_exec_host(value: Any) -> str | None:
    first = _string_or_none(value)
    if first is None:
        return None
    return _before_separator(_before_separator(first, ":"), "*")


def _clean_node_name(value: Any) -> str | None:
    text = _string_or_none(value)
    if text is None:
        return None
    if text.startswith("(") and text.endswith(")"):
        return None
    return text


def _before_separator(value: Any, separator: str) -> str | None:
    text = _string_or_none(value)
    if text is None:
        return None
    return text.split(separator, 1)[0].strip() or None


def _after_separator(value: Any, separator: str) -> str | None:
    text = _string_or_none(value)
    if text is None or separator not in text:
        return None
    return text.split(separator, 1)[1].strip() or None


def _string_or_none(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _stringify_json(value: Any) -> str:
    return json.dumps(value)
