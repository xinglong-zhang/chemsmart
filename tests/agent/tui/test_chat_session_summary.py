from __future__ import annotations

import asyncio
from pathlib import Path

from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.widgets.cells import (
    AgentMessageCell,
    FinalAnswerCell,
    ToolCallCell,
    ToolChainToggleCell,
)
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


def test_blocked_summary_hides_stale_command_and_remains_visible(
    tmp_path: Path,
):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            transcript = app.query_one(Transcript)
            transcript.clear_cells()
            app.chat_screen._begin_turn()
            transcript.add_cell(
                ToolCallCell(
                    tool="synthesize_command",
                    status="ok",
                    description="Synthesized an intermediate command.",
                    provider_call_id="synth-blocked",
                )
            )
            stale = FinalAnswerCell(
                "```bash\nchemsmart run gaussian opt\n```",
                title="Final Command",
            )
            transcript.add_cell(stale)
            assert transcript.collapse_tool_chain(
                app.chat_screen._active_turn_id
            )

            app.chat_screen._apply_log_entry(
                {
                    "kind": "session_summary",
                    "payload": {
                        "blocked": True,
                        "block_reason": "semantic_reject",
                        "total_steps_executed": 1,
                        "total_steps_planned": 2,
                    },
                }
            )
            await pilot.pause()

            children = list(transcript.query_one("#cells").children)
            visible = [cell for cell in children if cell.display]
            assert len(visible) == 2
            assert isinstance(visible[0], ToolChainToggleCell)
            assert isinstance(visible[1], AgentMessageCell)
            assert "semantic_reject" in visible[1].source_text
            assert not stale.display

    asyncio.run(scenario())
