from __future__ import annotations

import asyncio
from pathlib import Path

from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.widgets.footer import FooterWidget


def test_footer_spinner_visibility_and_animation(tmp_path: Path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            footer = app.query_one(FooterWidget)

            assert not footer.spinner_visible
            idle_text = str(footer.renderable)
            assert "Working…" not in idle_text
            assert all(
                label not in idle_text for label in footer.SPINNER_LABELS
            )

            footer.set_phase(Phase.PLANNING)
            footer.set_hint("Agent is planning…")
            await pilot.pause()

            assert footer.spinner_visible
            planning_text = str(footer.renderable)
            assert any(
                label in planning_text for label in footer.SPINNER_LABELS
            )

            first_spinner = footer.spinner_text
            await pilot.pause(0.3)
            second_spinner = footer.spinner_text
            assert first_spinner != second_spinner

            footer.set_phase(Phase.IDLE)
            await pilot.pause()

            assert not footer.spinner_visible
            idle_again = str(footer.renderable)
            assert "Working…" not in idle_again
            assert all(
                label not in idle_again for label in footer.SPINNER_LABELS
            )

    asyncio.run(scenario())


def test_footer_spinner_plain_mode_uses_static_working_text(tmp_path: Path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(
            plain=True,
            session_root=tmp_path / "sessions",
        )
        async with app.run_test() as pilot:
            await pilot.pause()
            footer = app.query_one(FooterWidget)
            footer.set_phase(Phase.RUNNING)
            footer.set_hint("Agent is checking tools…")
            await pilot.pause()

            assert footer.spinner_visible
            assert footer.spinner_text == "Working…"

            before = str(footer.renderable)
            await pilot.pause(0.3)
            after = str(footer.renderable)
            assert before == after
            assert "Working…" in after

    asyncio.run(scenario())
