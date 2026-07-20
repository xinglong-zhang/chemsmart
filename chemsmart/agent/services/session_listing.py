"""Read-only projections used by the agent session-listing command."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from chemsmart.agent.services.session_store import current_session_dirs

UTC = timezone.utc


def load_session_snapshots(session_root: Path) -> list[dict[str, object]]:
    """Return newest-first summaries for current-schema session directories."""
    if not session_root.exists():
        return []
    snapshots = [
        _session_snapshot(session_dir)
        for session_dir in current_session_dirs(session_root)
    ]
    return sorted(
        snapshots,
        key=lambda snapshot: snapshot["timestamp"],
        reverse=True,
    )


def format_session_age(timestamp: datetime) -> str:
    """Format a UTC timestamp as a compact relative age."""
    delta_seconds = max(
        0, int((datetime.now(UTC) - timestamp).total_seconds())
    )
    if delta_seconds < 60:
        return f"{delta_seconds}s ago"
    minutes = delta_seconds // 60
    if minutes < 60:
        return f"{minutes}m ago"
    hours = minutes // 60
    if hours < 24:
        return f"{hours}h ago"
    days = hours // 24
    return f"{days}d ago"


def truncate_session_request(request: object, limit: int = 60) -> str:
    """Collapse whitespace and bound a request for table display."""
    compact = " ".join(str(request or "").split())
    if len(compact) <= limit:
        return compact
    return f"{compact[: limit - 3]}..."


def _session_snapshot(session_dir: Path) -> dict[str, object]:
    metadata = _load_json(session_dir / "session_metadata.json")
    state = _load_json(session_dir / "session.json")
    entries = _load_decision_entries(session_dir / "decision_log.jsonl")

    request = ""
    if isinstance(metadata, dict):
        request = str(metadata.get("request") or "")
    if not request and isinstance(state, dict):
        request = str(state.get("request") or "")
    if not request:
        for entry in entries:
            if entry.get("kind") == "request":
                request = str(entry.get("payload", {}).get("request") or "")
                break

    timestamp = _coerce_timestamp(
        (metadata or {}).get("ended_at")
        or (metadata or {}).get("started_at")
        or (entries[-1].get("ts") if entries else None)
        or (state or {}).get("started_at")
    )
    if timestamp is None:
        timestamp = datetime.fromtimestamp(session_dir.stat().st_mtime, tz=UTC)

    return {
        "session_id": session_dir.name,
        "request": request,
        "status": _session_status(metadata, state, entries),
        "timestamp": timestamp,
    }


def _session_status(
    metadata: dict | None,
    state: dict | None,
    entries: list[dict[str, object]],
) -> str:
    if isinstance(metadata, dict):
        if metadata.get("blocked"):
            return "blocked"
        if metadata.get("critic_verdict"):
            return "ok"

    summary = next(
        (
            entry.get("payload")
            for entry in reversed(entries)
            if entry.get("kind") == "session_summary"
        ),
        None,
    )
    if isinstance(summary, dict):
        if summary.get("blocked"):
            return "blocked"
        return "ok"

    if any(entry.get("kind") in {"error", "tool_error"} for entry in entries):
        return "error"

    if isinstance(state, dict):
        planned = int(state.get("total_steps_planned") or 0)
        current = int(state.get("current_step_index") or 0)
        if planned == 0 or current < planned:
            return "in-progress"
        return "ok"

    return "in-progress"


def _load_json(path: Path) -> dict | None:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _load_decision_entries(path: Path) -> list[dict[str, object]]:
    if not path.exists():
        return []
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def _coerce_timestamp(value: object) -> datetime | None:
    if not value:
        return None
    text = str(value)
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    try:
        return datetime.fromisoformat(text).astimezone(UTC)
    except ValueError:
        return None
