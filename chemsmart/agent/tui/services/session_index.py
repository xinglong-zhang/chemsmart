"""Discovery helpers that keep runtime artifacts out of session UX."""

from __future__ import annotations

from pathlib import Path

from chemsmart.agent.services.session_store import current_session_dirs


def agent_session_dirs(session_root: str | Path) -> list[Path]:
    return current_session_dirs(session_root)


__all__ = ["agent_session_dirs"]
