"""Provider-independent projection of an adapted compact SPEC."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class CompactSpecProjection:
    intent: str
    message: str
    commands: tuple[str, ...]


def project_compact_spec(
    adapted: dict[str, Any],
    source: dict[str, Any],
) -> CompactSpecProjection:
    return CompactSpecProjection(
        intent=str(adapted.get("intent") or source.get("intent") or ""),
        message=str(adapted.get("message") or source.get("message") or ""),
        commands=tuple(
            command
            for command in adapted.get("commands", [])
            if isinstance(command, str) and command.strip()
        ),
    )


def non_workflow_result(
    projection: CompactSpecProjection,
) -> dict[str, Any]:
    status = (
        "infeasible"
        if projection.intent in {"decline", "chitchat"}
        else "needs_clarification"
    )
    return {
        "status": status,
        "command": "",
        "explanation": (
            projection.message or "No executable workflow was requested."
        ),
        "confidence": "high",
        "missing_info": [],
        "alternatives": [],
    }


__all__ = [
    "CompactSpecProjection",
    "non_workflow_result",
    "project_compact_spec",
]
