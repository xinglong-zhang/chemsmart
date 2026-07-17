"""Read-only session browser for Phase 1."""

from __future__ import annotations

import json
from pathlib import Path

from textual.app import ComposeResult
from textual.binding import Binding
from textual.screen import ModalScreen
from textual.widgets import Static

from chemsmart.agent.tui.services.session_index import agent_session_dirs


class SessionsScreen(ModalScreen[str | None]):
    BINDINGS = [
        Binding("escape", "close", "Close", priority=True),
        Binding("q", "close", "Close"),
        Binding("up", "cursor_up", show=False),
        Binding("down", "cursor_down", show=False),
        Binding("enter", "select_session", show=False),
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
        self._session_ids = self._load_session_ids()
        self._selected_index = 0 if self._session_ids else -1

    def compose(self) -> ComposeResult:
        yield Static(self._render_text(), id="sessions-modal")

    def action_close(self) -> None:
        self.dismiss(None)

    def action_cursor_up(self) -> None:
        if not self._session_ids:
            return
        self._selected_index = max(0, self._selected_index - 1)
        self._refresh_text()

    def action_cursor_down(self) -> None:
        if not self._session_ids:
            return
        self._selected_index = min(
            len(self._session_ids) - 1,
            self._selected_index + 1,
        )
        self._refresh_text()

    def action_select_session(self) -> None:
        if not self._session_ids or self._selected_index < 0:
            self.dismiss(None)
            return
        self._resume_session(self._session_ids[self._selected_index])

    def _refresh_text(self) -> None:
        self.query_one("#sessions-modal", Static).update(self._render_text())

    def _resume_session(self, session_id: str) -> None:
        owner = (
            self.app.screen_stack[-2]
            if len(self.app.screen_stack) > 1
            else None
        )
        resume = getattr(owner, "_resume_or_prompt", None)
        self.dismiss(session_id)
        if callable(resume):
            self.app.call_after_refresh(resume, session_id)

    def _load_session_ids(self) -> list[str]:
        if not self.session_root.exists():
            return []
        return [
            path.name for path in agent_session_dirs(self.session_root)[:10]
        ]

    def _render_text(self) -> str:
        lines = ["Sessions", "", "Use /resume <session-id> to load one."]
        if not self._session_ids:
            lines.extend(["", "No sessions found."])
            return "\n".join(lines)

        for index, session_id in enumerate(self._session_ids):
            prefix = "▶ " if index == self._selected_index else "  "
            request = self._load_request_summary(session_id)
            lines.append(
                f"{prefix}{session_id}: {request or '(no request saved)'}"
            )
        return "\n".join(lines)

    def _load_request_summary(self, session_id: str) -> str:
        session_dir = self.session_root / session_id
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
        return request
