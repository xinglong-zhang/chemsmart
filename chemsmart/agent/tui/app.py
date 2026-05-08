"""Application entrypoint for the Phase 1 agent TUI."""

from __future__ import annotations

from pathlib import Path

from textual.app import App

from chemsmart.agent.core import _default_session_root
from chemsmart.agent.tui.bindings import BINDINGS
from chemsmart.agent.tui.screens.chat import ChatScreen


class ChemsmartTuiApp(App[None]):
    CSS_PATH = "styles.tcss"
    BINDINGS = BINDINGS

    def __init__(
        self,
        *,
        plain: bool = False,
        session_root: str | Path | None = None,
    ) -> None:
        super().__init__()
        self.plain = plain
        self.session_root = Path(session_root or _default_session_root())
        self.session_root.mkdir(parents=True, exist_ok=True)
        self.chat_screen = ChatScreen(session_root=self.session_root)

    def on_mount(self) -> None:
        self.push_screen(self.chat_screen)


def launch_tui(
    *,
    plain: bool = False,
    session_root: str | Path | None = None,
) -> None:
    app = ChemsmartTuiApp(plain=plain, session_root=session_root)
    if plain:
        app.animation_level = "none"
    app.run(inline=True, inline_no_clear=True, mouse=not plain)
