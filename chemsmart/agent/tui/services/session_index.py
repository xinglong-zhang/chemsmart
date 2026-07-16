"""Discovery helpers that keep runtime artifacts out of session UX."""

from __future__ import annotations

from pathlib import Path

_SESSION_MARKERS = (
    "session.json",
    "state.json",
    "session_metadata.json",
    "decision_log.jsonl",
)


def agent_session_dirs(session_root: str | Path) -> list[Path]:
    root = Path(session_root)
    if not root.is_dir():
        return []
    sessions = [
        path
        for path in root.iterdir()
        if path.is_dir()
        and not path.name.startswith(".")
        and any((path / marker).is_file() for marker in _SESSION_MARKERS)
    ]
    return sorted(sessions, reverse=True)


__all__ = ["agent_session_dirs"]
