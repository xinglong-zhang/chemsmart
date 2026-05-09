from __future__ import annotations

import asyncio
from pathlib import Path

from textual.app import App, ComposeResult

from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.screens.sessions import SessionsScreen

from ._helpers import write_session_fixture


class _SessionsHarness(App[None]):
    def __init__(self, session_root: Path) -> None:
        super().__init__()
        self.session_root = session_root
        self.selected_session_id: str | None = None

    def compose(self) -> ComposeResult:
        return
        yield

    def on_mount(self) -> None:
        self.push_screen(
            SessionsScreen(self.session_root),
            self._capture_selection,
        )

    def _capture_selection(self, session_id: str | None) -> None:
        self.selected_session_id = session_id


def test_sessions_screen_down_enter_dismisses_selected_session(
    tmp_path: Path,
):
    session_root = tmp_path / "sessions"
    write_session_fixture(session_root, "session-001")
    write_session_fixture(session_root, "session-002")

    async def scenario() -> None:
        app = _SessionsHarness(session_root)
        async with app.run_test() as pilot:
            await pilot.pause()
            await pilot.press("down", "enter")
            await pilot.pause()
            assert app.selected_session_id == "session-001"

    asyncio.run(scenario())


def test_sessions_screen_calls_chat_resume_on_selection(tmp_path: Path):
    session_root = tmp_path / "sessions"
    write_session_fixture(session_root, "session-001")
    write_session_fixture(session_root, "session-002")
    resumed: list[str] = []

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=session_root)
        app.chat_screen._resume_or_prompt = resumed.append
        async with app.run_test() as pilot:
            await pilot.pause()
            app.chat_screen._handle_slash_command("/sessions")
            await pilot.pause()
            await pilot.press("down", "enter")
            await pilot.pause()
            assert resumed == ["session-001"]
            assert app.screen is app.chat_screen

    asyncio.run(scenario())
