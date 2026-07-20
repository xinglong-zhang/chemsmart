from __future__ import annotations

import asyncio
from pathlib import Path

from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.widgets.composer import Composer
from chemsmart.agent.tui.widgets.popups import CwdMismatchOverlay

from ._helpers import write_session_fixture


def test_resume_cwd_mismatch_choices(monkeypatch, tmp_path: Path):
    session_root = tmp_path / "sessions"
    write_session_fixture(
        session_root,
        cwd="/tmp/recorded-chemsmart-session",
    )

    async def run_choice(choice: str):
        app = ChemsmartTuiApp(
            session_root=session_root, job_poll_interval=60.0
        )
        start_calls: list[tuple[str, str | None]] = []
        chdir_calls: list[str] = []
        monkeypatch.setattr(
            app.chat_screen,
            "start_resume",
            lambda session_id, cwd_override=None: start_calls.append(
                (session_id, cwd_override)
            ),
        )
        monkeypatch.setattr(
            "chemsmart.agent.tui.mixins.job_interaction.os.chdir",
            lambda path: chdir_calls.append(path),
        )
        async with app.run_test() as pilot:
            await pilot.pause()
            composer = app.query_one(Composer)
            composer.load_text("/resume session-001")
            await pilot.press("enter")
            await pilot.pause()
            assert isinstance(app.screen, CwdMismatchOverlay)
            await pilot.press(choice)
            await pilot.pause()
        return start_calls, chdir_calls

    async def scenario() -> None:
        start_calls, chdir_calls = await run_choice("n")
        assert start_calls == []
        assert chdir_calls == []

        start_calls, chdir_calls = await run_choice("c")
        assert chdir_calls == ["/tmp/recorded-chemsmart-session"]
        assert start_calls == [("session-001", None)]

        start_calls, chdir_calls = await run_choice("i")
        assert chdir_calls == []
        assert start_calls == [("session-001", str(Path.cwd().resolve()))]

    asyncio.run(scenario())
