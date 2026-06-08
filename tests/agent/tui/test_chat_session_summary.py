from __future__ import annotations

import asyncio
from pathlib import Path

from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.widgets.cells import AgentMessageCell
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.transcript import Transcript


def _summary_cells(app: ChemsmartTuiApp) -> list[AgentMessageCell]:
    transcript = app.query_one(Transcript).query_one("#cells")
    return [
        child
        for child in transcript.children
        if isinstance(child, AgentMessageCell)
        and child.border_title == "Summary"
    ]


def test_advisory_empty_plan_does_not_render_session_finished(
    tmp_path: Path,
):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            app.chat_screen._apply_log_entry(
                {
                    "kind": "plan",
                    "payload": {
                        "steps": [],
                        "rationale": "I'm the chemsmart agent.",
                        "estimated_cost": "none",
                        "intent": "advisory",
                    },
                }
            )
            app.chat_screen._apply_log_entry(
                {
                    "kind": "session_summary",
                    "payload": {
                        "blocked": False,
                        "block_reason": None,
                        "total_steps_executed": 0,
                        "total_steps_planned": 0,
                    },
                }
            )
            await pilot.pause()

            assert _summary_cells(app) == []
            footer = app.query_one(FooterWidget)
            assert footer.phase == Phase.IDLE
            assert footer.hint == "Ready"

    asyncio.run(scenario())


def test_blocked_summary_omits_unknown_block_reason_suffix(tmp_path: Path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            app.chat_screen._apply_log_entry(
                {
                    "kind": "session_summary",
                    "payload": {
                        "blocked": True,
                        "block_reason": "",
                        "total_steps_executed": 0,
                        "total_steps_planned": 1,
                    },
                }
            )
            await pilot.pause()

            summary_cells = _summary_cells(app)
            assert len(summary_cells) == 1
            assert (
                summary_cells[0].source_text == "Session finished (0/1 steps)."
            )
            assert "Block reason: unknown." not in summary_cells[0].source_text
            footer = app.query_one(FooterWidget)
            assert footer.phase == Phase.ERROR
            assert footer.hint == "Blocked"

    asyncio.run(scenario())
