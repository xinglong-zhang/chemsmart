from __future__ import annotations

import asyncio
from pathlib import Path

from rich.console import Console
from textual.widgets import Static

from chemsmart.agent.core import Plan, _preview_value
from chemsmart.agent.tools import build_molecule
from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.widgets.cells import (
    AgentMessageCell,
    DryRunInputCell,
    ErrorCell,
    PlanCell,
    RuntimeValidationCell,
    UserMessageCell,
)
from chemsmart.agent.tui.widgets.composer import Composer
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.popups import (
    ApprovalOverlay,
    CwdMismatchOverlay,
)
from chemsmart.agent.tui.widgets.transcript import Transcript

from .._agent_session_helpers import planner_plan


def _render_plain(renderable) -> str:
    console = Console(width=240, record=True)
    console.print(renderable)
    return console.export_text()


def test_tool_call_workflow_row_updates_to_done(
    single_molecule_xyz_file,
    tmp_path: Path,
):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            app.chat_screen._apply_log_entry(
                {
                    "kind": "request",
                    "payload": {
                        "request": (
                            f"single-point on {single_molecule_xyz_file}"
                        )
                    },
                }
            )
            app.chat_screen._apply_log_entry(
                {
                    "kind": "plan",
                    "payload": planner_plan(
                        single_molecule_xyz_file,
                        "workflow_case",
                    ),
                }
            )
            await pilot.pause()
            transcript = app.query_one(Transcript).query_one("#cells")
            plan_cell = next(
                child
                for child in transcript.children
                if isinstance(child, PlanCell)
            )

            app.chat_screen._apply_log_entry(
                {
                    "kind": "tool_call",
                    "payload": {
                        "step_index": 1,
                        "tool": "build_molecule",
                        "args": {"filepath": single_molecule_xyz_file},
                    },
                }
            )
            await pilot.pause()
            running_text = _render_plain(plan_cell.renderable)
            assert "build_molecule" in running_text
            assert "Checking" in running_text

            molecule = build_molecule(single_molecule_xyz_file)
            app.chat_screen._apply_log_entry(
                {
                    "kind": "tool_result",
                    "payload": {
                        "step_index": 1,
                        "tool": "build_molecule",
                        "payload": _preview_value(molecule),
                    },
                }
            )
            await pilot.pause()
            done_text = _render_plain(plan_cell.renderable)
            assert "✓ 1. build_molecule" in done_text
            assert "71 atoms" in done_text
            assert Path(single_molecule_xyz_file).stem in done_text
            assert ".xyz" in done_text

    asyncio.run(scenario())


def test_chat_advisory_renders_agent_message_cell(tmp_path: Path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            app.chat_screen._apply_log_entry(
                {
                    "kind": "plan",
                    "payload": {
                        "steps": [],
                        "rationale": "I am chemsmart agent.",
                        "estimated_cost": "low",
                    },
                }
            )
            await pilot.pause()
            transcript = app.query_one(Transcript).query_one("#cells")
            advisory_cells = [
                child
                for child in transcript.children
                if isinstance(child, AgentMessageCell)
                and child.border_title == "Advisory"
            ]
            assert len(advisory_cells) == 1
            assert advisory_cells[0].source_text == "I am chemsmart agent."
            assert not any(
                isinstance(child, PlanCell) for child in transcript.children
            )
            footer = app.query_one(FooterWidget)
            assert footer.phase == Phase.FINISHED
            assert footer.phase != Phase.PLANNING
            assert footer.hint == "Response ready"

    asyncio.run(scenario())


def test_runtime_validation_partial_copy_is_compact():
    cell = RuntimeValidationCell(
        {
            "ok": "partial",
            "local_ok": True,
            "local_issues": [],
            "remote_unknown": [
                "server.queue required",
                "server.account required",
                "server.scratch_dir required",
                "server.modules_or_executable_path required",
            ],
        }
    )

    text = _render_plain(cell.renderable)
    assert "Local OK / 4 remote item(s) needed" in text
    assert "queue / account / scratch_dir / modules_or_executable_path" in text
    assert "- server.queue required" not in text


def test_dry_run_input_cell_shows_path_before_body():
    cell = DryRunInputCell(
        "%chk=water.chk\n# B3LYP/6-31G(d)\n",
        inputfile="/tmp/water.com",
        command="chemsmart run gaussian -f water.xyz -c 0 -m 1 sp",
        cli_grounded=True,
    )

    text = _render_plain(cell.renderable)
    assert "✓ water.com ready" in text
    assert "path: /tmp/water.com" in text
    assert "Generated chemsmart CLI command" in text
    assert "chemsmart run gaussian -f water.xyz -c 0 -m 1 sp" in text
    assert text.index("Generated chemsmart CLI command") < text.index(
        "%chk=water.chk"
    )


def test_plan_cell_renders_rationale_when_no_steps():
    cell = PlanCell(Plan(steps=[], rationale="Answer directly."))

    text = _render_plain(cell.renderable)
    assert "Answer directly." in text
    assert "no data" not in text


def test_approval_popup_includes_mode_and_rollback_labels(tmp_path: Path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            app.push_screen(
                ApprovalOverlay(action="submit_hpc", request="submit water"),
            )
            await pilot.pause()
            summary = app.screen.query_one("#approval-summary", Static)
            text = str(summary.renderable)
            assert "mode:" in text
            assert "rollback:" in text
            assert "remote-effect:" in text

    asyncio.run(scenario())


def test_greeter_copy_mentions_unified_validation_and_execution(tmp_path: Path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            transcript = app.query_one(Transcript).query_one("#cells")
            first = transcript.children[0]
            text = _render_plain(first.renderable)
            assert "/doctor" in text
            assert "same interface" in text
            assert "/run" in text
            assert "/submit" in text

    asyncio.run(scenario())


def test_clear_emits_short_notification(tmp_path: Path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()

            calls: list[tuple[str, float | None]] = []
            original_notify = app.chat_screen.notify

            def spy(message, *args, **kwargs):
                calls.append((message, kwargs.get("timeout")))
                return original_notify(message, *args, **kwargs)

            app.chat_screen.notify = spy

            composer = app.query_one(Composer)
            composer.load_text("/clear")
            await pilot.press("enter")
            await pilot.pause()

            assert any(
                msg == "Transcript cleared." and timeout == 3
                for msg, timeout in calls
            )
            transcript = app.query_one(Transcript).query_one("#cells")
            assert not any(
                isinstance(child, AgentMessageCell)
                and (child.source_text or "").startswith("Transcript cleared")
                for child in transcript.children
            )

    asyncio.run(scenario())


def test_missing_file_error_uses_exact_copy():
    cell = ErrorCell(
        "FileNotFoundError",
        "[Errno 2] No such file or directory: 'missing.xyz'",
        {"tool": "build_molecule"},
    )

    text = _render_plain(cell.renderable)
    assert "Input file not found" in text
    assert (
        "The path in your request does not exist relative to the current working directory."
        in text
    )
    assert "examples/h2o.xyz" in text


def test_cwd_mismatch_popup_uses_designer_copy(tmp_path: Path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            app.push_screen(
                CwdMismatchOverlay(
                    recorded_cwd="/tmp/recorded",
                    current_cwd="/tmp/current",
                )
            )
            await pilot.pause()
            summary = app.screen.query_one(Static)
            text = str(summary.renderable)
            assert "Working directory mismatch" in text
            assert "resume after switching to the recorded directory" in text
            assert "force continue from the current directory" in text

    asyncio.run(scenario())


def test_user_messages_are_right_aligned_without_moving_slash_commands(
    tmp_path: Path,
):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test(size=(100, 24)) as pilot:
            await pilot.pause()
            transcript = app.query_one(Transcript)
            transcript.clear_cells()
            agent = AgentMessageCell("agent " * 20)
            user = UserMessageCell("build water and recommend a method " * 4)
            slash = UserMessageCell("/help")
            transcript.add_cell(agent)
            transcript.add_cell(user)
            transcript.add_cell(slash)
            await pilot.pause()

            transcript_right_edge = (
                transcript.content_region.x + transcript.content_region.width
            )
            user_right_edge = user.bubble.region.x + user.bubble.region.width

            assert user_right_edge >= transcript_right_edge - 1
            assert user.bubble.region.x >= agent.region.x + int(
                transcript.content_region.width * 0.15
            )
            assert slash.bubble.region.x == agent.region.x

    asyncio.run(scenario())
