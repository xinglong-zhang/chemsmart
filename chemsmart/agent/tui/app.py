"""Application entrypoint for the Phase 1 agent TUI."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from textual.app import App

from chemsmart.agent.core import _default_session_root
from chemsmart.agent.services.command_logging import _silence_console_logging
from chemsmart.agent.tui.bindings import BINDINGS
from chemsmart.agent.tui.config import load_tui_config
from chemsmart.agent.tui.screens.chat import ChatScreen


class ChemsmartTuiApp(App[None]):
    CSS_PATH = "styles.tcss"
    BINDINGS = BINDINGS

    def __init__(
        self,
        *,
        plain: bool = False,
        session_root: str | Path | None = None,
        job_poll_interval: float = 5.0,
        runtime_v2: str = "active",
    ) -> None:
        super().__init__()
        self.tui_config = load_tui_config()
        for action, key in self.tui_config.keybindings.items():
            self.bind(key, action, description=action.replace("_", " "))
        self.plain = plain
        self.session_root = Path(session_root or _default_session_root())
        self.session_root.mkdir(parents=True, exist_ok=True)
        self.chat_screen = ChatScreen(
            session_root=self.session_root,
            job_poll_interval=job_poll_interval,
            runtime_v2=runtime_v2,
        )

    def on_mount(self) -> None:
        self.theme = "textual-dark"
        self.push_screen(self.chat_screen)

    def _delegate_to_chat(self, name: str) -> None:
        screen = self.chat_screen
        action = getattr(screen, name, None)
        if callable(action):
            action()

    def action_soft_cancel(self) -> None:
        self._delegate_to_chat("action_soft_cancel")

    def action_quit_if_empty(self) -> None:
        self._delegate_to_chat("action_quit_if_empty")

    def action_refresh_screen(self) -> None:
        self._delegate_to_chat("action_refresh_screen")

    def action_show_project_yaml(self) -> None:
        self._delegate_to_chat("action_show_project_yaml")

    def action_show_shortcuts(self) -> None:
        self._delegate_to_chat("action_show_shortcuts")

    def action_toggle_transcript(self) -> None:
        self._delegate_to_chat("action_toggle_transcript")

    def action_show_activity(self) -> None:
        self._delegate_to_chat("action_show_activity")

    def action_show_calculations(self) -> None:
        self._delegate_to_chat("action_show_calculations")

    def action_search_history(self) -> None:
        self._delegate_to_chat("action_search_history")

    def action_dismiss_overlay(self) -> None:
        self._delegate_to_chat("action_dismiss_overlay")


def launch_tui(
    *,
    plain: bool = False,
    session_root: str | Path | None = None,
    job_poll_interval: float = 5.0,
) -> None:
    app = ChemsmartTuiApp(
        plain=plain,
        session_root=session_root,
        job_poll_interval=job_poll_interval,
    )
    log_path = _silence_console_logging(app.session_root.parent)
    run_kwargs = {}
    if plain:
        app.animation_level = "none"
        run_kwargs = {"inline": True, "inline_no_clear": True, "mouse": False}
    try:
        app.run(**run_kwargs)
    finally:
        for handler in logging.getLogger().handlers:
            if isinstance(handler, logging.FileHandler):
                handler.flush()
        if log_path.exists() and log_path.stat().st_size:
            print(f"agent tui exited; debug log: {log_path}", file=sys.stderr)
