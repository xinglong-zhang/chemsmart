from __future__ import annotations

import asyncio
from pathlib import Path

from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.widgets.cells import (
    AgentMessageCell,
    CriticVerdictCell,
    DryRunInputCell,
    PlanCell,
    RuntimeValidationCell,
    UserMessageCell,
)
from chemsmart.agent.tui.widgets.composer import Composer
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.transcript import Transcript

from .._agent_session_helpers import FakeProvider, critic_ok, planner_plan


def test_chitchat_input_renders_single_reply_cell_only(
    monkeypatch,
    tmp_path: Path,
):
    provider = FakeProvider(
        [
            {
                "steps": [],
                "rationale": "Hello! I can help plan chemsmart workflows.",
                "estimated_cost": "none",
            }
        ]
    )

    def fake_get_provider():
        return provider

    monkeypatch.setattr(
        "chemsmart.agent.providers.get_provider", fake_get_provider
    )
    monkeypatch.setattr("chemsmart.agent.core.get_provider", fake_get_provider)

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            composer = app.query_one(Composer)
            composer.load_text("hello")
            await pilot.press("enter")
            for _ in range(30):
                await pilot.pause(0.1)
                if app.query_one(FooterWidget).phase == Phase.IDLE:
                    break

            transcript = app.query_one(Transcript).query_one("#cells")
            children = list(transcript.children)
            assert [type(child) for child in children] == [
                UserMessageCell,
                AgentMessageCell,
            ]
            reply_cell = children[1]
            assert reply_cell.border_title == "Reply"
            assert reply_cell.source_text == (
                "Hello! I can help plan chemsmart workflows."
            )
            assert not any(
                isinstance(child, (PlanCell, CriticVerdictCell))
                for child in children
            )
            footer = app.query_one(FooterWidget)
            assert footer.phase == Phase.IDLE
            assert footer.hint == "Ready"
            assert len(provider.calls) == 1

    asyncio.run(scenario())


def test_chemistry_request_keeps_full_workflow_rendering(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    import chemsmart.agent.tools as agent_tools

    provider = FakeProvider(
        [planner_plan(single_molecule_xyz_file, "workflow_case"), critic_ok()]
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
            composer.load_text(
                f"single-point on {single_molecule_xyz_file} at B3LYP/6-31G(d) Gaussian"
            )
            await pilot.press("enter")
            for _ in range(30):
                await pilot.pause(0.1)
                transcript = app.query_one(Transcript).query_one("#cells")
                children = list(transcript.children)
                if app.query_one(FooterWidget).phase == Phase.FINISHED and any(
                    isinstance(child, AgentMessageCell)
                    and child.border_title == "Summary"
                    for child in children
                ):
                    break

            transcript = app.query_one(Transcript).query_one("#cells")
            children = list(transcript.children)
            assert any(isinstance(child, PlanCell) for child in children)
            assert any(
                isinstance(child, DryRunInputCell) for child in children
            )
            assert any(
                isinstance(child, RuntimeValidationCell) for child in children
            )
            assert any(
                isinstance(child, CriticVerdictCell) for child in children
            )
            assert any(
                isinstance(child, AgentMessageCell)
                and child.border_title == "Summary"
                for child in children
            )
            footer = app.query_one(FooterWidget)
            assert footer.phase == Phase.FINISHED

    asyncio.run(scenario())
