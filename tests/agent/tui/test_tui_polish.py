from __future__ import annotations

import asyncio
import json
from pathlib import Path
from threading import Event

import pytest
from textual.widgets import TextArea

from chemsmart.agent.model_command_parser import parse_model_command
from chemsmart.agent.provider_adapter import ToolRequest
from chemsmart.agent.registry import ToolRegistry
from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.widgets.cells import (
    AgentMessageCell,
    CommandInterpretationCell,
    FinalAnswerCell,
    ToolCallCell,
    ToolChainToggleCell,
)
from chemsmart.agent.tui.widgets.composer import Composer
from chemsmart.agent.tui.widgets.popups import (
    ProjectWriteOverlay,
    ResponseCopyOverlay,
)
from chemsmart.agent.tui.widgets.transcript import Transcript


class _FinishedWorker:
    is_finished = True


def _write_candidate(
    session_root: Path,
    *,
    project: str = "candidate",
    program: str = "gaussian",
) -> None:
    session_dir = session_root / "candidate-session"
    session_dir.mkdir(parents=True)
    yaml_text = (
        "gas:\n"
        "  functional: b3lyp\n"
        "  basis: def2svp\n"
        "  freq: true\n"
        "solv:\n"
        "  functional: b3lyp\n"
        "  basis: def2svp\n"
        "  freq: false\n"
    )
    entries = [
        {
            "kind": "tool_use_result",
            "payload": {
                "tool": "render_project_yaml",
                "status": "ok",
                "payload": {
                    "ok": True,
                    "project_name": project,
                    "program": program,
                    "yaml_text": yaml_text,
                },
            },
        },
        {
            "kind": "tool_use_result",
            "payload": {
                "tool": "validate_project_yaml",
                "status": "ok",
                "payload": {"ok": True, "verdict": "ok"},
            },
        },
    ]
    (session_dir / "decision_log.jsonl").write_text(
        "\n".join(json.dumps(entry) for entry in entries) + "\n",
        encoding="utf-8",
    )


def test_composer_completion_keeps_caret_at_end(tmp_path: Path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            composer = app.query_one(Composer)
            composer.load_text("/do")
            await pilot.pause()
            await pilot.press("tab")
            await pilot.pause()

            assert composer.text == "/doctor "
            assert composer.cursor_location == (0, len("/doctor "))
            await pilot.press("x")
            assert composer.text == "/doctor x"

    asyncio.run(scenario())


def test_write_project_dialog_overwrites_active_or_writes_new(
    monkeypatch, tmp_path: Path
):
    monkeypatch.chdir(tmp_path)
    session_root = tmp_path / "sessions"
    _write_candidate(session_root)
    active = tmp_path / ".chemsmart" / "gaussian" / "water.yaml"
    active.parent.mkdir(parents=True)
    active.write_text(
        "gas:\n  functional: pbe0\n  basis: def2svp\n",
        encoding="utf-8",
    )

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=session_root)
        async with app.run_test() as pilot:
            await pilot.pause()
            screen = app.chat_screen
            assert screen._activate_workspace_project_path(active)
            captured: list[dict[str, object]] = []
            screen.run_project_yaml_write = (  # type: ignore[method-assign]
                lambda args: captured.append(dict(args)) or _FinishedWorker()
            )

            screen._handle_project_write_command("")
            await pilot.pause()
            assert isinstance(app.screen, ProjectWriteOverlay)
            assert app.screen.target == active
            await pilot.press("y")
            await pilot.pause()
            assert captured[-1]["project_name"] == "water"
            assert captured[-1]["overwrite"] is True
            result = ToolRegistry.default().call(
                "write_project_yaml", captured[-1]
            )
            assert result["ok"] is True
            assert "functional: b3lyp" in active.read_text(encoding="utf-8")

            screen._current_worker = None
            screen._handle_project_write_command("")
            await pilot.pause()
            assert isinstance(app.screen, ProjectWriteOverlay)
            await pilot.press("n")
            await pilot.pause()
            assert captured[-1]["project_name"] == "candidate"
            assert captured[-1]["overwrite"] is False

    asyncio.run(scenario())


def test_agent_project_write_approval_mutates_runtime_arguments(
    monkeypatch, tmp_path: Path
):
    monkeypatch.chdir(tmp_path)
    target = tmp_path / ".chemsmart" / "orca" / "water.yaml"
    target.parent.mkdir(parents=True)
    target.write_text("gas:\n  functional: pbe0\n  basis: def2svp\n")

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            args: dict[str, object] = {
                "project_name": "water",
                "program": "orca",
                "yaml_text": "gas:\n  functional: pbe0\n  basis: def2svp\n",
                "overwrite": False,
            }
            request = ToolRequest(
                request_id="write-1",
                provider="openai",
                provider_call_id="call-write-1",
                name="write_project_yaml",
                arguments_json=json.dumps(args),
                arguments=args,
                raw={},
            )
            screen = app.chat_screen
            waiter = Event()
            screen._pending_approval = True
            screen._pending_tool_request = request
            screen._approval_waiter = waiter
            screen._request_approval("write_project_yaml")
            await pilot.pause()

            assert isinstance(app.screen, ProjectWriteOverlay)
            await pilot.press("y")
            await pilot.pause()
            assert request.arguments["project_name"] == "water"
            assert request.arguments["overwrite"] is True
            assert waiter.is_set()

    asyncio.run(scenario())


def test_completed_tool_chain_collapses_and_toggles(tmp_path: Path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            transcript = app.query_one(Transcript)
            transcript.clear_cells()
            transcript.start_turn("turn-polish")
            tool = ToolCallCell(
                tool="read_project_yaml",
                status="ok",
                description="Read the active project.",
                provider_call_id="read-1",
            )
            answer = FinalAnswerCell("The project is ready.", title="Answer")
            transcript.add_cell(tool)
            transcript.add_cell(answer)

            assert tool.display
            assert transcript.collapse_tool_chain("turn-polish")
            await pilot.pause()
            children = list(transcript.query_one("#cells").children)
            toggle = next(
                child
                for child in children
                if isinstance(child, ToolChainToggleCell)
            )
            assert children[-1] is answer
            assert not tool.display
            assert not toggle.expanded

            toggle.action_toggle()
            await pilot.pause()
            assert tool.display
            assert toggle.expanded

            toggle.action_toggle()
            await pilot.pause()
            assert not tool.display
            assert not toggle.expanded
            transcript.clear_turn_chrome()
            assert tool in transcript.query_one("#cells").children

    asyncio.run(scenario())


def test_completed_tool_chain_keyboard_toggle_can_hide_again(tmp_path: Path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            transcript = app.query_one(Transcript)
            transcript.clear_cells()
            transcript.start_turn("turn-keyboard-toggle")
            tool = ToolCallCell(
                tool="validate_project_yaml",
                status="ok",
                description="Validate the workspace project.",
                provider_call_id="validate-1",
            )
            answer = FinalAnswerCell("The project is valid.", title="Answer")
            transcript.add_cell(tool)
            transcript.add_cell(answer)
            assert transcript.collapse_tool_chain("turn-keyboard-toggle")
            await pilot.pause()

            toggle = transcript.query_one(ToolChainToggleCell)
            toggle.focus()
            await pilot.press("enter")
            await pilot.pause()
            assert toggle.expanded
            assert tool.display

            await pilot.press("enter")
            await pilot.pause()
            assert not toggle.expanded
            assert not tool.display
            assert answer.display

    asyncio.run(scenario())


def test_completed_repair_chain_keeps_only_last_command_visible(
    tmp_path: Path,
):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            transcript = app.query_one(Transcript)
            transcript.clear_cells()
            transcript.start_turn("turn-repair")

            first_tool = ToolCallCell(
                tool="synthesize_command",
                status="error",
                description="Initial command failed the semantic gate.",
                provider_call_id="synth-1",
            )
            first_parser = CommandInterpretationCell(
                parse_model_command(
                    "chemsmart run gaussian -p water -f inputs/h2o.xyz opt"
                )
            )
            first_answer = FinalAnswerCell(
                "```bash\nchemsmart run gaussian -p water "
                "-f inputs/h2o.xyz opt\n```",
                title="Final Command",
            )
            repair_tool = ToolCallCell(
                tool="repair_command",
                status="ok",
                description="Repaired and revalidated the command.",
                provider_call_id="repair-1",
            )
            repaired_parser = CommandInterpretationCell(
                parse_model_command(
                    "chemsmart run --fake gaussian -p water "
                    "-f inputs/h2o.xyz -c 0 -m 1 opt"
                )
            )
            repaired_answer = FinalAnswerCell(
                "```bash\nchemsmart run --fake gaussian -p water "
                "-f inputs/h2o.xyz -c 0 -m 1 opt\n```",
                title="Repaired Command",
            )
            summary = AgentMessageCell(
                "Session finished (2/2 steps).",
                title="Summary",
            )
            for cell in (
                first_tool,
                first_parser,
                first_answer,
                repair_tool,
                repaired_parser,
                repaired_answer,
                summary,
            ):
                transcript.add_cell(cell)

            assert transcript.collapse_tool_chain("turn-repair")
            await pilot.pause()

            children = list(transcript.query_one("#cells").children)
            toggle = next(
                cell
                for cell in children
                if isinstance(cell, ToolChainToggleCell)
            )
            visible = [cell for cell in children if cell.display]
            assert visible == [toggle, repaired_answer]
            assert not first_answer.display
            assert not first_parser.display
            assert not repaired_parser.display
            assert not summary.display

            toggle.action_toggle()
            await pilot.pause()
            assert all(
                cell.display
                for cell in (
                    first_tool,
                    first_parser,
                    first_answer,
                    repair_tool,
                    repaired_parser,
                    summary,
                )
            )

    asyncio.run(scenario())


def test_completed_non_command_turn_keeps_rich_answer_not_summary(
    tmp_path: Path,
):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            transcript = app.query_one(Transcript)
            transcript.clear_cells()
            transcript.start_turn("turn-yaml")
            tool = ToolCallCell(
                tool="validate_project_yaml",
                status="ok",
                description="Validated the candidate YAML.",
                provider_call_id="yaml-1",
            )
            answer = AgentMessageCell(
                "The PBE0 project YAML is valid and ready for approval.",
                title="Assistant",
            )
            summary = AgentMessageCell(
                "Session finished (1/1 steps).",
                title="Summary",
            )
            for cell in (tool, answer, summary):
                transcript.add_cell(cell)

            assert transcript.collapse_tool_chain("turn-yaml")
            await pilot.pause()

            children = list(transcript.query_one("#cells").children)
            visible = [cell for cell in children if cell.display]
            assert len(visible) == 2
            assert isinstance(visible[0], ToolChainToggleCell)
            assert visible[1] is answer
            assert not tool.display
            assert not summary.display

    asyncio.run(scenario())


def test_response_copy_view_supports_selection(monkeypatch, tmp_path: Path):
    copied: list[str] = []

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        monkeypatch.setattr(app, "copy_to_clipboard", copied.append)
        async with app.run_test() as pilot:
            await pilot.pause()
            response = FinalAnswerCell(
                "Energy: -76.269 Eh", title="ORCA result"
            )
            transcript = app.query_one(Transcript)
            transcript.add_cell(response)
            transcript.scroll_end(animate=False)
            await pilot.pause()
            assert await pilot.click(response, offset=(2, 0))
            await pilot.pause()
            assert isinstance(app.screen, ResponseCopyOverlay)
            area = app.screen.query_one("#response-copy-text", TextArea)
            assert area.read_only
            await pilot.press("a")
            await pilot.pause()
            assert copied
            assert copied[-1] == "Energy: -76.269 Eh"

    asyncio.run(scenario())


@pytest.mark.parametrize("viewport", [(80, 24), (120, 36), (160, 45)])
def test_project_write_overlay_fits_viewport(
    viewport: tuple[int, int], tmp_path: Path
):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test(size=viewport) as pilot:
            await pilot.pause()
            target = tmp_path / ".chemsmart" / "gaussian" / "water.yaml"
            app.push_screen(
                ProjectWriteOverlay(
                    target=target,
                    overwrite_project="water",
                    new_project="water-2",
                    target_exists=True,
                )
            )
            await pilot.pause()
            modal = app.screen.query_one("#project-write-modal")
            assert modal.region.width <= viewport[0]
            assert modal.region.height <= viewport[1]

    asyncio.run(scenario())
