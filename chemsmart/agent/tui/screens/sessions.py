"""Read-only session browser for Phase 1."""

from __future__ import annotations

import json
from pathlib import Path

from textual.app import ComposeResult
from textual.binding import Binding
from textual.screen import ModalScreen
from textual.widgets import Static


class SessionsScreen(ModalScreen[None]):
    BINDINGS = [
        Binding("escape", "close", "Close", priority=True),
        Binding("q", "close", "Close"),
        Binding("enter", "close", "Close"),
    ]

    DEFAULT_CSS = """
    SessionsScreen {
        align: center middle;
    }

    #sessions-modal {
        width: 80;
        height: auto;
        max-height: 24;
        border: round $primary;
        padding: 1 2;
        background: $surface;
    }
    """

    def __init__(self, session_root: Path) -> None:
        super().__init__()
        self.session_root = session_root

    def compose(self) -> ComposeResult:
        yield Static(self._render_text(), id="sessions-modal")

    def action_close(self) -> None:
        self.dismiss(None)

    def _render_text(self) -> str:
        lines = ["Sessions", "", "Use /resume <session-id> to load one."]
        session_dirs = (
            sorted(
                [
                    path
                    for path in self.session_root.iterdir()
                    if path.is_dir()
                ],
                reverse=True,
            )
            if self.session_root.exists()
            else []
        )
        if not session_dirs:
            lines.extend(["", "No sessions found."])
            return "\n".join(lines)

        for session_dir in session_dirs[:10]:
            request = ""
            metadata_path = session_dir / "session_metadata.json"
            session_path = session_dir / "session.json"
            if metadata_path.exists():
                data = json.loads(metadata_path.read_text())
                request = str(data.get("request") or "")
            elif session_path.exists():
                data = json.loads(session_path.read_text())
                request = str(data.get("request") or "")
            request = request.strip().replace("\n", " ")
            if len(request) > 48:
                request = f"{request[:45]}..."
            lines.append(
                f"- {session_dir.name}: {request or '(no request saved)'}"
            )
        return "\n".join(lines)
