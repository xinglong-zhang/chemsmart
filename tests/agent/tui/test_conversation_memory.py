from __future__ import annotations

import asyncio
from pathlib import Path

from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.widgets.cells import (
    AgentMessageCell,
    ToolCallCell,
    UserMessageCell,
)
from chemsmart.agent.tui.widgets.composer import Composer
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.transcript import Transcript

from .._agent_session_helpers import FakeProvider, planner_plan
from .._loop_helpers import (
    openai_final_response,
    openai_tool_call_response,
    tool_call,
)


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
    first_plan = planner_plan(single_molecule_xyz_file, "tui_memory_opt")
    second_plan = _single_point_plan(single_molecule_xyz_file, "tui_memory_sp")
    provider = FakeProvider(
        [
            {
                "__raw_response__": openai_tool_call_response(
                    *[
                        tool_call(f"call_a_{i}", step["tool"], step["args"])
                        for i, step in enumerate(first_plan["steps"], start=1)
                    ]
                )
            },
            {"__raw_response__": openai_final_response("First turn done.")},
            {
                "__raw_response__": openai_tool_call_response(
                    *[
                        tool_call(f"call_b_{i}", step["tool"], step["args"])
                        for i, step in enumerate(second_plan["steps"], start=1)
                    ]
                )
            },
            {"__raw_response__": openai_final_response("Second turn done.")},
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
            user_cells = [
                child
                for child in children
                if isinstance(child, UserMessageCell)
            ]
            assert len(user_cells) >= 2
            assert user_cells[0].source_text == first_request
            assert user_cells[1].source_text == second_request
            assert any(isinstance(child, ToolCallCell) for child in children)
            assert any(
                isinstance(child, AgentMessageCell) for child in children
            )
            ephemeral_summary_cells = [
                child
                for child in children
                if isinstance(child, AgentMessageCell)
                and getattr(child, "border_title", None) == "Assistant"
            ]
            assert len(ephemeral_summary_cells) >= 1
            assert app.chat_screen.active_agent_session is not None

        assert len(list(session_root.iterdir())) == 1

    asyncio.run(scenario())


def test_clear_resets_transcript_and_session(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    import chemsmart.agent.tools as agent_tools

    plan = planner_plan(single_molecule_xyz_file, "tui_memory_opt")
    provider = FakeProvider(
        [
            {
                "__raw_response__": openai_tool_call_response(
                    *[
                        tool_call(f"call_c_{i}", step["tool"], step["args"])
                        for i, step in enumerate(plan["steps"], start=1)
                    ]
                )
            },
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
