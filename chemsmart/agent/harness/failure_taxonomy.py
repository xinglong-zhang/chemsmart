"""Stable deterministic classification for safe runtime failures."""

from __future__ import annotations

import re
from dataclasses import dataclass


@dataclass(frozen=True)
class RuntimeFailure:
    rule_id: str
    category: str
    message: str
    missing_info: tuple[str, ...] = ()


def classify_runtime_failure(
    *,
    stdout: str = "",
    stderr: str = "",
    returncode: int | None = None,
) -> RuntimeFailure:
    text = f"{stderr}\n{stdout}".strip()
    lowered = text.lower()
    if (
        "no project settings implemented" in lowered
        or "currently available projects:" in lowered
    ):
        missing = ["valid chemsmart project configuration"]
        available = _available_line(text, "Currently available projects:")
        if available:
            missing.append(f"available projects: {available}")
        return RuntimeFailure(
            "cmd.runtime.project_not_found",
            "workspace",
            "selected ChemSmart project YAML could not be resolved",
            tuple(missing),
        )
    if (
        "no server implemented" in lowered
        or "currently available servers:" in lowered
    ):
        missing = ["valid chemsmart server configuration"]
        available = _available_line(text, "Currently available servers:")
        if available:
            missing.append(f"available servers: {available}")
        return RuntimeFailure(
            "cmd.runtime.server_invalid",
            "workspace",
            "selected ChemSmart server YAML could not be resolved",
            tuple(missing),
        )
    if any(
        marker in lowered
        for marker in (
            "no such file or directory",
            "filenotfounderror",
            "does not exist",
        )
    ):
        return RuntimeFailure(
            "cmd.runtime.input_not_found",
            "input",
            "requested local input or endpoint file does not exist",
            ("existing local input file path",),
        )
    if any(
        marker in lowered
        for marker in (
            "no such option",
            "invalid value for",
            "usageerror",
            "badparameter",
        )
    ):
        return RuntimeFailure(
            "cmd.runtime.cli_value_error",
            "command",
            "ChemSmart rejected an option name, position, or value",
        )
    if any(
        marker in lowered
        for marker in (
            "open babel",
            "openbabel",
            "modulenotfounderror",
            "importerror",
        )
    ):
        return RuntimeFailure(
            "cmd.runtime.dependency_missing",
            "environment",
            "runtime dependency required by this workflow is unavailable",
        )
    if any(
        marker in lowered
        for marker in ("could not find any runners", "no runner", "runner")
    ):
        return RuntimeFailure(
            "cmd.runtime.runner_unavailable",
            "environment",
            "no compatible ChemSmart runner was available",
        )
    if re.search(r"assertionerror|valueerror|typeerror|traceback", lowered):
        return RuntimeFailure(
            "cmd.runtime.builder_error",
            "runtime",
            "ChemSmart settings or input builder rejected the command",
        )
    return RuntimeFailure(
        "cmd.runtime.unclassified",
        "runtime",
        f"safe ChemSmart runtime validation failed with exit code {returncode}",
    )


def _available_line(text: str, marker: str) -> str:
    if marker not in text:
        return ""
    return text.split(marker, 1)[1].splitlines()[0].strip()


__all__ = ["RuntimeFailure", "classify_runtime_failure"]
