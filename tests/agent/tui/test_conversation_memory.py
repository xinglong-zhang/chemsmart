from __future__ import annotations

import asyncio
import json
from pathlib import Path

from chemsmart.agent.tui.app import ChemsmartTuiApp
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


def _single_point_plan(filepath: str, label: str) -> dict:
    return {
        "steps": [
            {
                "tool": "build_molecule",
                "args": {"filepath": filepath},
                "rationale": "Reuse the earlier molecule reference.",
            },
            {
                "tool": "build_gaussian_settings",
                "args": {"functional": "b3lyp", "basis": "6-31g*"},
                "rationale": "Keep the follow-up method lightweight.",
            },
            {
                "tool": "build_job",
                "args": {
                    "kind": "gaussian.sp",
                    "molecule": "$step1",
                    "settings": "$step2",
                    "label": label,
                },
                "rationale": "Prepare the single-point follow-up.",
            },
            {
                "tool": "dry_run_input",
                "args": {"job": "$step3"},
                "rationale": "Preview the follow-up input file.",
            },
            {
                "tool": "validate_runtime",
                "args": {"job": "$step3", "server": None},
                "rationale": "Check prerequisites before execution.",
            },
        ],
        "rationale": "Follow the earlier setup with a single-point run.",
        "estimated_cost": "low",
    }


def test_second_tui_request_includes_prior_turn_memory(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    import chemsmart.agent.tools as agent_tools

    first_request = f"build water from {single_molecule_xyz_file}"
    second_request = "now run a single-point on it"
    provider = FakeProvider(
        [
            planner_plan(single_molecule_xyz_file, "tui_memory_opt"),
            critic_ok(),
            _single_point_plan(single_molecule_xyz_file, "tui_memory_sp"),
            critic_ok(),
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
        session_root = tmp_path / "sessions"
        app = ChemsmartTuiApp(session_root=session_root)
        async with app.run_test() as pilot:
            await pilot.pause()
            composer = app.query_one(Composer)

            composer.load_text(first_request)
            await pilot.press("enter")
            for _ in range(30):
                await pilot.pause(0.1)
                footer_text = str(
                    app.query_one(FooterWidget).renderable
                ).lower()
                if "finished" in footer_text:
                    break

            composer.load_text(second_request)
            await pilot.press("enter")
            for _ in range(30):
                await pilot.pause(0.1)
                if len(provider.calls) >= 4:
                    footer_text = str(
                        app.query_one(FooterWidget).renderable
                    ).lower()
                    if "finished" in footer_text:
                        break

            transcript = app.query_one(Transcript).query_one("#cells")
            children = list(transcript.children)
            assert [type(child) for child in children] == [
                UserMessageCell,
                UserMessageCell,
                PlanCell,
                DryRunInputCell,
                RuntimeValidationCell,
                CriticVerdictCell,
                AgentMessageCell,
            ]
            assert children[0].source_text == first_request
            assert children[1].source_text == second_request
            ephemeral_summary_cells = [
                child
                for child in children
                if isinstance(child, AgentMessageCell)
                and getattr(child, "border_title", None) == "Summary"
            ]
            assert len(ephemeral_summary_cells) == 1
            assert app.chat_screen.active_agent_session is not None

        planner_payload = json.loads(
            provider.calls[2]["messages"][1]["content"]
        )
        history = planner_payload["conversation_history"]["recent_turns"]
        assert history[0]["request"] == first_request
        assert any(
            "dry_run_input wrote" in item
            for item in history[0]["reusable_results"]
        )
        assert len(list(session_root.iterdir())) == 1

    asyncio.run(scenario())


def test_clear_resets_transcript_and_session(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    import chemsmart.agent.tools as agent_tools

    provider = FakeProvider(
        [planner_plan(single_molecule_xyz_file, "tui_memory_opt"), critic_ok()]
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
            composer.load_text(f"build water from {single_molecule_xyz_file}")
            await pilot.press("enter")
            for _ in range(30):
                await pilot.pause(0.1)
                footer_text = str(
                    app.query_one(FooterWidget).renderable
                ).lower()
                if "finished" in footer_text:
                    break

            assert app.chat_screen.active_agent_session is not None

            composer.load_text("/clear")
            await pilot.press("enter")
            await pilot.pause()

            transcript = app.query_one(Transcript).query_one("#cells")
            assert list(transcript.children) == []
            assert app.chat_screen.active_agent_session is None
            assert app.chat_screen.active_resume_id is None

    asyncio.run(scenario())
