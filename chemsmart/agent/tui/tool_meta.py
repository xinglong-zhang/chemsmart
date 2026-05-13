"""Shared tool metadata for approval and transcript UI."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from rich.text import Text

from chemsmart.agent.registry import ToolRegistry

_TOOL_META = {
    "build_molecule": {
        "risk": "read-only",
        "read_only": True,
        "style": "warning",
        "summary": "reads a structure file",
    },
    "recommend_method": {
        "risk": "read-only",
        "read_only": True,
        "style": "warning",
        "summary": "computes an advisory method recommendation",
    },
    "build_gaussian_settings": {
        "risk": "read-only",
        "read_only": True,
        "style": "warning",
        "summary": "builds validated Gaussian settings",
    },
    "build_orca_settings": {
        "risk": "read-only",
        "read_only": True,
        "style": "warning",
        "summary": "builds validated ORCA settings",
    },
    "build_job": {
        "risk": "local-state",
        "read_only": False,
        "style": "warning",
        "summary": "creates a job object handle",
    },
    "dry_run_input": {
        "risk": "mutates-state",
        "read_only": False,
        "style": "warning",
        "summary": "writes a local input file",
    },
    "validate_runtime": {
        "risk": "inspection",
        "read_only": True,
        "style": "warning",
        "summary": "checks local/runtime state",
    },
    "run_local": {
        "risk": "risky",
        "read_only": False,
        "style": "error",
        "summary": "starts a local job",
    },
    "extract_optimized_geometry": {
        "risk": "read-only",
        "read_only": True,
        "style": "warning",
        "summary": "reads output logs to extract geometry",
    },
    "submit_hpc": {
        "risk": "risky",
        "read_only": False,
        "style": "error",
        "summary": "may submit to a remote queue",
    },
    "wizard_probe": {
        "risk": "inspection",
        "read_only": True,
        "style": "warning",
        "summary": "probes local or remote server state",
    },
    "wizard_refresh": {
        "risk": "inspection",
        "read_only": True,
        "style": "warning",
        "summary": "refreshes or reuses the wizard node cache",
    },
    "wizard_verify": {
        "risk": "inspection",
        "read_only": True,
        "style": "warning",
        "summary": "verifies wizard/server transport wiring",
    },
    "wizard_write": {
        "risk": "risky",
        "read_only": False,
        "style": "error",
        "summary": "writes a server YAML config",
    },
    "read": {
        "risk": "inspection",
        "read_only": True,
        "style": "warning",
        "summary": "reads local file lines",
    },
    "ssh_probe": {
        "risk": "inspection",
        "read_only": True,
        "style": "warning",
        "summary": "runs a remote inspection probe",
    },
    "scheduler_query": {
        "risk": "inspection",
        "read_only": True,
        "style": "warning",
        "summary": "inspects remote scheduler state",
    },
    "log_tail": {
        "risk": "inspection",
        "read_only": True,
        "style": "warning",
        "summary": "tails remote logs for diagnostics",
    },
}

_STATUS_STYLE = {
    "pending": "warning",
    "approved": "success",
    "ok": "success",
    "denied": "error",
    "error": "error",
    "skipped": "error",
    "interrupted": "error",
}

_ASSUMPTION_CONFIDENCE = {
    "ok": "high",
    "partial": "med",
    "error": "low",
    "denied": "low",
    "skipped": "low",
    "ask_user": "low",
}


def tool_description(
    tool_name: str,
    *,
    registry: ToolRegistry | None = None,
) -> str:
    if registry is not None:
        return registry.describe_tool(tool_name)
    return tool_name


def tool_side_effect_summary(tool_name: str) -> str:
    return str(
        _TOOL_META.get(tool_name, {}).get(
            "summary", "may inspect or mutate local state"
        )
    )


def tool_risk_badge(tool_name: str) -> tuple[str, str]:
    meta = _TOOL_META.get(tool_name, {})
    return (
        str(meta.get("risk") or "unknown"),
        str(meta.get("style") or "warning"),
    )


def tool_read_only(tool_name: str) -> bool:
    return bool(_TOOL_META.get(tool_name, {}).get("read_only"))


def tool_status_style(status: str) -> str:
    return _STATUS_STYLE.get(status, "warning")


def format_assumptions_banner(
    entities: dict[str, Any] | None,
    recent_tool_status: str | None,
) -> str | None:
    if not entities:
        return None

    scheduler = _string_or_none(entities.get("last_scheduler"))
    server = _string_or_none(entities.get("last_server"))
    job_id = _string_or_none(entities.get("last_job_id"))
    log_path = _string_or_none(entities.get("last_log_path"))
    if not any((scheduler, server, job_id, log_path)):
        return None

    parts: list[str] = []
    if scheduler is not None:
        parts.append(scheduler.upper())
    if server is not None:
        parts.append(server)
    if job_id is not None:
        parts.append(f"job {job_id}")
    if log_path is not None:
        parts.append(f"log {Path(log_path).name}")
    confidence = (
        _ASSUMPTION_CONFIDENCE.get(recent_tool_status, "unknown")
        if recent_tool_status is not None
        else "unknown"
    )
    parts.append(f"conf={confidence}")

    line = " · ".join(parts)
    if len(line) > 120:
        return f"{line[:119]}…"
    return line


def pretty_tool_args(arguments: dict[str, Any] | None) -> str:
    payload = arguments or {}
    return json.dumps(payload, indent=2, sort_keys=True)


def render_tool_result_summary(
    tool_name: str,
    result: dict[str, Any] | None,
) -> str | None:
    payload = result or {}
    if not isinstance(payload, dict) or not payload:
        return None

    error = _result_error_summary(payload)
    if error is not None:
        return error

    if tool_name == "read":
        return _render_read_summary(payload)
    if tool_name == "ssh_probe":
        return _render_ssh_probe_summary(payload)
    if tool_name == "scheduler_query":
        return _render_scheduler_query_summary(payload)
    if tool_name == "log_tail":
        return _render_log_tail_summary(payload)
    return None


def render_tool_result_detail(
    tool_name: str,
    result: dict[str, Any] | None,
) -> list[Text] | None:
    payload = result or {}
    if not isinstance(payload, dict) or not payload:
        return None

    if payload.get("error") is not None:
        lines = []
        message = _string_or_none(payload.get("message"))
        if message is not None:
            lines.append(Text(message, style="error"))
        elif _string_or_none(payload.get("path")) is not None:
            lines.append(
                Text(
                    f"path: {_string_or_none(payload.get('path'))}",
                    style="dim",
                )
            )
        return lines or None

    if tool_name == "read":
        return _render_read_detail(payload)
    if tool_name == "ssh_probe":
        return _render_ssh_probe_detail(payload)
    if tool_name == "scheduler_query":
        return _render_scheduler_query_detail(payload)
    if tool_name == "log_tail":
        return _render_log_tail_detail(payload)
    return None


def _render_read_summary(payload: dict[str, Any]) -> str | None:
    start_line = payload.get("start_line")
    end_line = payload.get("end_line")
    total_lines = payload.get("total_lines")
    if start_line is None or end_line is None or total_lines is None:
        return None
    summary = f"L{start_line}-{end_line} of {total_lines}"
    if payload.get("truncated"):
        summary += ", truncated"
    return summary


def _render_ssh_probe_summary(payload: dict[str, Any]) -> str | None:
    returncode = payload.get("returncode")
    duration_s = payload.get("duration_s")
    if returncode is None or duration_s is None:
        return None
    return f"rc={returncode} in {float(duration_s):.2f}s"


def _render_scheduler_query_summary(payload: dict[str, Any]) -> str | None:
    scheduler = _string_or_none(payload.get("scheduler")) or "scheduler"
    state = _string_or_none(payload.get("state"))
    if state is not None:
        return f"{scheduler}: {state}"

    partition = _string_or_none(payload.get("partition_or_queue")) or "?"
    nodes = _string_or_none(payload.get("nodes")) or "?"
    total_cpus = _string_or_none(payload.get("total_cpus")) or "?"
    return f"{scheduler}: {partition} {nodes}n/{total_cpus}cpu"


def _render_log_tail_summary(payload: dict[str, Any]) -> str | None:
    lines_returned = payload.get("lines_returned")
    if lines_returned is None:
        return None

    errors = payload.get("errors") or []
    kinds: list[str] = []
    for item in errors[:3]:
        if not isinstance(item, dict):
            continue
        kind = _string_or_none(item.get("kind"))
        if kind is not None:
            kinds.append(kind)
    top_kinds = ", ".join(kinds) if kinds else "none"
    return f"{lines_returned}L, {len(errors)} errors: {top_kinds}"


def _render_read_detail(payload: dict[str, Any]) -> list[Text] | None:
    content = str(payload.get("content") or "")
    if not content:
        return None
    lines = [Text(line, style="dim") for line in content.splitlines()[:5]]
    truncated = bool(payload.get("truncated")) or len(content.splitlines()) > 5
    return _finalize_detail_lines(lines, truncated=truncated)


def _render_ssh_probe_detail(payload: dict[str, Any]) -> list[Text] | None:
    lines: list[Text] = []
    server = _string_or_none(payload.get("server"))
    probe = _string_or_none(payload.get("probe") or payload.get("probe_name"))
    if server is not None:
        lines.append(Text(f"server: {server}", style="dim"))
    if probe is not None:
        lines.append(Text(f"probe: {probe}", style="dim"))

    output_lines = _nonempty_lines(payload.get("stdout_truncated"))
    lines.extend(Text(line, style="dim") for line in output_lines[:3])

    stderr_lines = _nonempty_lines(payload.get("stderr_truncated"))
    if stderr_lines:
        lines.append(Text(stderr_lines[0], style="error"))

    truncated = len(output_lines) > 3 or len(stderr_lines) > 1
    return _finalize_detail_lines(lines, truncated=truncated)


def _render_scheduler_query_detail(
    payload: dict[str, Any],
) -> list[Text] | None:
    if payload.get("state") is not None:
        field_pairs = [
            ("job", payload.get("job_id")),
            ("state", payload.get("state")),
            ("queue", payload.get("queue")),
            ("user", payload.get("user")),
            ("node", payload.get("node")),
        ]
    else:
        field_pairs = [
            ("scheduler", payload.get("scheduler")),
            ("queue", payload.get("partition_or_queue")),
            ("nodes", payload.get("nodes")),
            ("cpus", payload.get("total_cpus")),
        ]

    lines = [
        Text(f"{label}: {value}", style="dim")
        for label, value in field_pairs
        if value is not None
    ]
    return _finalize_detail_lines(lines)


def _render_log_tail_detail(payload: dict[str, Any]) -> list[Text] | None:
    errors = payload.get("errors") or []
    if not errors:
        return [Text("No error signatures found.", style="dim")]

    lines = []
    for item in errors[:3]:
        if not isinstance(item, dict):
            continue
        kind = _string_or_none(item.get("kind")) or "error"
        line = _string_or_none(item.get("line")) or "(no line)"
        line_no = _string_or_none(item.get("line_no")) or "?"
        lines.append(Text(f"L{line_no} {kind}: {line}", style="error"))
    return _finalize_detail_lines(lines, truncated=len(errors) > 3)


def _finalize_detail_lines(
    lines: list[Text],
    *,
    truncated: bool = False,
) -> list[Text] | None:
    if truncated and (not lines or lines[-1].plain != "… truncated"):
        lines = [*lines, Text("… truncated", style="dim")]
    if len(lines) > 6:
        return [*lines[:5], Text("… truncated", style="dim")]
    return lines or None


def _result_error_summary(payload: dict[str, Any]) -> str | None:
    error = _string_or_none(payload.get("error"))
    if error is None:
        return None
    return f"error: {error}"


def _nonempty_lines(value: Any) -> list[str]:
    text = _string_or_none(value)
    if text is None:
        return []
    return [line for line in text.splitlines() if line.strip()]


def _string_or_none(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text or None
