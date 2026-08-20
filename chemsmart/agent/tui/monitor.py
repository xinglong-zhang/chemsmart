"""View-only tail of the append-only runtime event stream.

The durable ``events.jsonl`` is the single source of live progress for both
a planning session and an approved execution; the interface reads the same
evidence the audit reads, so there is no second state channel to drift. The
reader takes no lock (the flock is writer-side), only ever parses
newline-terminated bytes, and treats a malformed line as a display gap
rather than a crash.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path
from typing import Any, Mapping

from .humanize import (
    HumanizedRowV1,
    humanize_provider_turn,
    humanize_tool_settled,
    humanize_tool_started,
)


@dataclass
class EventTailerV1:
    """Byte-offset incremental reader over an append-only JSONL file."""

    path: Path
    offset: int = 0

    def poll(self) -> tuple[dict[str, Any], ...]:
        if not self.path.exists():
            return ()
        size = self.path.stat().st_size
        if size <= self.offset:
            return ()
        with self.path.open("rb") as handle:
            handle.seek(self.offset)
            chunk = handle.read(size - self.offset)
        cut = chunk.rfind(b"\n")
        if cut < 0:
            # A torn tail: wait for the writer to finish the line.
            return ()
        events: list[dict[str, Any]] = []
        consumed = 0
        for line in chunk[: cut + 1].splitlines(keepends=True):
            consumed += len(line)
            stripped = line.strip()
            if not stripped:
                continue
            try:
                parsed = json.loads(stripped)
            except json.JSONDecodeError:
                events.append({"kind": "__display_gap__", "payload": {}})
                continue
            if isinstance(parsed, dict):
                events.append(parsed)
        self.offset += consumed
        return tuple(events)


# -- planning feed ----------------------------------------------------------

#: Event kinds worth a transcript row while a provider session runs.
_PLANNING_SETTLED = {"tool_succeeded", "tool_failed"}


def planning_feed_update(event: Mapping[str, Any]) -> HumanizedRowV1 | None:
    """Map one runtime event to a humanized planning row, or None."""

    kind = str(event.get("kind") or "")
    payload = event.get("payload")
    payload = payload if isinstance(payload, Mapping) else {}
    if kind == "tool_started":
        return humanize_tool_started(payload)
    if kind in _PLANNING_SETTLED:
        return humanize_tool_settled(kind, payload)
    if kind == "provider_turn_observed":
        return humanize_provider_turn(payload)
    if kind == "turn_blocked":
        reason = str(payload.get("reason") or "turn blocked")
        return HumanizedRowV1(icon="△", text=reason, state="failed")
    if kind == "tool_waiting":
        timeout = payload.get("timeout_seconds")
        text = "waiting on the engine"
        if timeout:
            text += f" (up to {timeout:g} s)"
        return HumanizedRowV1(icon="…", text=text, state="running")
    if kind == "tool_woke":
        waited = payload.get("waited_seconds")
        reason = payload.get("wake_reason") or "woke"
        text = f"engine wait ended: {reason}"
        if waited is not None:
            text += f" after {float(waited):.1f} s"
        return HumanizedRowV1(icon="…", text=text, state="finished")
    if kind == "__display_gap__":
        return HumanizedRowV1(
            icon="·", text="(one unreadable event line skipped)", state="note"
        )
    return None


#: request_id -> row identity, so a settle can mutate its started row.
def planning_row_key(event: Mapping[str, Any]) -> str:
    payload = event.get("payload")
    payload = payload if isinstance(payload, Mapping) else {}
    return str(payload.get("request_id") or "")


# -- execution feed ---------------------------------------------------------


@dataclass(frozen=True)
class ExecutionSignalV1:
    """One typed change a jobs/DAG panel can apply."""

    kind: str  # node_launched | engine_done | node_state | analysis_settled
    #        | report_rendered | terminal
    node_id: str = ""
    state: str = ""
    detail: str = ""
    payload: Mapping[str, Any] = field(default_factory=dict)


def execution_signal(event: Mapping[str, Any]) -> ExecutionSignalV1 | None:
    kind = str(event.get("kind") or "")
    payload = event.get("payload")
    payload = payload if isinstance(payload, Mapping) else {}
    if kind == "workflow_node_launch_reserved":
        return ExecutionSignalV1(
            kind="node_launched",
            node_id=str(payload.get("node_id") or ""),
            state="running",
            payload=payload,
        )
    if kind == "program_execution_observed":
        validated = bool(payload.get("validated"))
        return ExecutionSignalV1(
            kind="engine_done",
            node_id=str(payload.get("node_id") or ""),
            state=str(payload.get("execution_state") or ""),
            detail="validated" if validated else "awaiting validation",
            payload=payload,
        )
    if kind == "workflow_node_state_changed":
        return ExecutionSignalV1(
            kind="node_state",
            node_id=str(payload.get("node_id") or ""),
            state=str(payload.get("node_state") or ""),
            detail=str(payload.get("workflow_state") or ""),
            payload=payload,
        )
    if kind == "workflow_analysis_node_settled":
        return ExecutionSignalV1(
            kind="analysis_settled",
            node_id=str(payload.get("node_id") or ""),
            state=str(payload.get("state") or ""),
            detail=str(payload.get("reason") or ""),
            payload=payload,
        )
    if kind == "workflow_analysis_report_rendered":
        return ExecutionSignalV1(kind="report_rendered", payload=payload)
    if kind == "runtime_terminated":
        return ExecutionSignalV1(
            kind="terminal",
            state=str(payload.get("terminal_state") or ""),
            detail=str(payload.get("reason") or ""),
            payload=payload,
        )
    return None


__all__ = [
    "EventTailerV1",
    "ExecutionSignalV1",
    "execution_signal",
    "planning_feed_update",
    "planning_row_key",
]
