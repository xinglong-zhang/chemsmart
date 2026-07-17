from __future__ import annotations

import asyncio
from pathlib import Path

from textual.app import App, ComposeResult
from textual.widgets import Static

from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.screens.sessions import SessionsScreen
from chemsmart.agent.tui.services.session_index import agent_session_dirs

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


def test_runtime_calculation_store_is_not_an_agent_session(tmp_path: Path):
    session_root = tmp_path / "sessions"
    write_session_fixture(session_root, "session-001")
    runtime = session_root / ".runtime" / "tui-one"
    runtime.mkdir(parents=True)
    (runtime / "decision_log.jsonl").write_text("{}\n", encoding="utf-8")

    assert [path.name for path in agent_session_dirs(session_root)] == [
        "session-001"
    ]


def test_legacy_state_only_directory_is_not_a_resumable_session(
    tmp_path: Path,
):
    session_root = tmp_path / "sessions"
    legacy = session_root / "legacy-001"
    legacy.mkdir(parents=True)
    (legacy / "state.json").write_text("{}", encoding="utf-8")

    assert agent_session_dirs(session_root) == []


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


def test_sessions_screen_updates_visible_selection_marker(tmp_path: Path):
    session_root = tmp_path / "sessions"
    write_session_fixture(session_root, "session-001")
    write_session_fixture(session_root, "session-002")

    async def scenario() -> None:
        app = _SessionsHarness(session_root)
        async with app.run_test() as pilot:
            await pilot.pause()
            modal = app.screen.query_one("#sessions-modal", Static)
            before = str(modal.renderable)
            assert "▶ session-002" in before
            assert "  session-001" in before

            await pilot.press("down")
            await pilot.pause()

            after = str(modal.renderable)
            assert "  session-002" in after
            assert "▶ session-001" in after

    asyncio.run(scenario())
