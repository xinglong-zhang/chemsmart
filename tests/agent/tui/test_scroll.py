from __future__ import annotations

import asyncio
from pathlib import Path

from textual import events

from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.widgets.cells import AgentMessageCell
from chemsmart.agent.tui.widgets.transcript import Transcript


def test_transcript_scrolls_with_keyboard_and_mouse(tmp_path: Path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test(size=(100, 24)) as pilot:
            await pilot.pause()
            transcript = app.query_one(Transcript)
            transcript.clear_cells()
            for index in range(40):
                transcript.add_cell(
                    AgentMessageCell(f"message {index} " + "x" * 80)
                )
            await pilot.pause()
            transcript.scroll_end(animate=False)
            await pilot.pause()

            assert transcript.max_scroll_y > 0
            assert transcript.scroll_y == transcript.max_scroll_y

            last_cell = list(transcript.query_one("#cells").children)[-1]
            last_cell.focus()
            await pilot.pause()

            await pilot.press("pageup")
            await pilot.pause()
            assert transcript.scroll_y < transcript.max_scroll_y

            after_page_up = transcript.scroll_y
            last_cell.post_message(
                events.MouseScrollUp(
                    last_cell,
                    x=1,
                    y=1,
                    delta_x=0,
                    delta_y=-1,
                    button=0,
                    shift=False,
                    meta=False,
                    ctrl=False,
                    screen_x=1,
                    screen_y=1,
                )
            )
            await pilot.pause()
            assert transcript.scroll_y < after_page_up

            await pilot.press("home")
            await pilot.pause()
            assert transcript.scroll_y == 0

            await pilot.press("end")
            await pilot.pause()
            assert transcript.scroll_y == transcript.max_scroll_y

            await pilot.press("pageup")
            await pilot.pause()
            assert transcript.scroll_y < transcript.max_scroll_y

            await pilot.press("pagedown")
            await pilot.pause()
            assert transcript.scroll_y == transcript.max_scroll_y

    asyncio.run(scenario())
