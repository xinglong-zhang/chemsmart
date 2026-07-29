from __future__ import annotations

import asyncio
from pathlib import Path

from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.widgets.composer import Composer
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.transcript import Transcript

from .._agent_session_helpers import FakeProvider
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
    monkeypatch.setenv("CHEMSMART_AGENT_TUI_MODE", "run")
    import chemsmart.agent.tools_command as command_tools

    request = f"optimize {single_molecule_xyz_file}"
    synthesis_payload = {
        "ok": True,
        "status": "ready",
        "command": (
            "chemsmart run gaussian -p test "
            f"-f {single_molecule_xyz_file} -c 0 -m 1 opt"
        ),
        "explanation": "Prepared by the deterministic synthesis tool.",
        "confidence": "high",
        "project": "test",
        "semantic": {"verdict": "ok", "failed_rule_ids": []},
        "decision_trace": {
            "action": "synthesize_command",
            "confidence": "high",
            "evidence": ["The request asks for a Gaussian optimization."],
        },
    }
    tool_calls = [
        tool_call(
            "call_1",
            "synthesize_command",
            {"request": request},
        )
    ]

    provider = FakeProvider(
        [
            {"__raw_response__": openai_tool_call_response(*tool_calls)},
            {"__raw_response__": openai_final_response("Done.")},
        ]
    )

    def fake_get_provider():
        return provider

    monkeypatch.setattr(
        "chemsmart.agent.providers.get_provider", fake_get_provider
    )
    monkeypatch.setattr("chemsmart.agent.core.get_provider", fake_get_provider)
    monkeypatch.setattr(
        command_tools,
        "synthesize_command",
        lambda request: synthesis_payload,
    )

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            composer = app.query_one(Composer)
            composer.load_text(request)
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
