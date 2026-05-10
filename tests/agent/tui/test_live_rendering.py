from __future__ import annotations

import asyncio
from pathlib import Path

from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.widgets.composer import Composer
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.transcript import Transcript

from .._agent_session_helpers import FakeProvider, planner_plan
from .._loop_helpers import (
    openai_final_response,
    openai_tool_call_response,
    tool_call,
)


def test_live_run_renders_plan_dry_run_and_critic_cells(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    import chemsmart.agent.tools as agent_tools

    plan = planner_plan(single_molecule_xyz_file, "tui_case")
    tool_calls = [
        tool_call(f"call_{index}", step["tool"], step["args"])
        for index, step in enumerate(plan["steps"], start=1)
    ]

    provider = FakeProvider(
        [
            {"__raw_response__": openai_tool_call_response(*tool_calls)},
            {"__raw_response__": openai_final_response("Done.")},
        ]
    )

    def fake_get_provider():
        return provider

    def fake_validate_runtime(job, server=None):
        return {
            "ok": "ok",
            "local_ok": True,
            "local_issues": [],
            "remote_unknown": [],
        }

    monkeypatch.setattr(
        "chemsmart.agent.providers.get_provider", fake_get_provider
    )
    monkeypatch.setattr("chemsmart.agent.core.get_provider", fake_get_provider)
    monkeypatch.setattr(agent_tools, "validate_runtime", fake_validate_runtime)

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            composer = app.query_one(Composer)
            composer.load_text(f"optimize {single_molecule_xyz_file}")
            await pilot.press("enter")
            for _ in range(20):
                await pilot.pause(0.1)
                if (
                    "finished"
                    in str(app.query_one(FooterWidget).renderable).lower()
                ):
                    break
            transcript = app.query_one(Transcript).query_one("#cells")
            child_types = [
                type(child).__name__ for child in transcript.children
            ]
            assert child_types[0] == "UserMessageCell"
            assert "ToolCallCell" in child_types
            assert "AgentMessageCell" in child_types
            footer_text = str(app.query_one(FooterWidget).renderable).lower()
            assert "finished" in footer_text

    asyncio.run(scenario())
