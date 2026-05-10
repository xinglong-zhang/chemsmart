"""Shared tool metadata for approval and transcript UI."""

from __future__ import annotations

import json
from typing import Any

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
    "wizard_write": {
        "risk": "risky",
        "read_only": False,
        "style": "error",
        "summary": "writes a server YAML config",
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


def tool_description(
    tool_name: str,
    *,
    registry: ToolRegistry | None = None,
) -> str:
    if registry is not None:
        return registry.describe_tool(tool_name)
    return tool_name


def tool_side_effect_summary(tool_name: str) -> str:
    return _TOOL_META.get(tool_name, {}).get(
        "summary", "may inspect or mutate local state"
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


def pretty_tool_args(arguments: dict[str, Any] | None) -> str:
    payload = arguments or {}
    return json.dumps(payload, indent=2, sort_keys=True)
