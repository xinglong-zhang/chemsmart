"""Small immutable contracts shared by chat-screen presenters and actions."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ReadyCommand:
    command: str
    action: str
    workspace: Path
    semantic_verdict: str
    intent_verdict: str
    project_path: Path | None = None
    project_sha256: str | None = None
    source: str = ""
