from __future__ import annotations

import asyncio
from statistics import quantiles
from time import perf_counter

import pytest
from textual.widgets import OptionList
from textual.worker import Worker, WorkerState

from chemsmart.agent.harness.workflow_state import (
    current_workflow_state,
    reset_workflow_state,
)
from chemsmart.agent.provider_config import AgentProviderConfig
from chemsmart.agent.runtime.calculations import (
    CalculationEvent,
    CalculationRun,
    CalculationStore,
)
from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.config import DEFAULT_KEYBINDINGS, load_tui_config
from chemsmart.agent.tui.events import ToolUseEvent
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.screens.calculations import CalculationMonitor
from chemsmart.agent.tui.state import TuiState, TuiStateReducer
from chemsmart.agent.tui.widgets.calculation_strip import (
    CalculationStatusStrip,
)
from chemsmart.agent.tui.widgets.cells import (
    CalculationReceiptCell,
    FinalAnswerCell,
    ToolCallCell,
)
from chemsmart.agent.tui.widgets.composer import Composer
from chemsmart.agent.tui.widgets.footer import FooterWidget
from chemsmart.agent.tui.widgets.popups import (
    ApprovalOverlay,
    ApprovalResult,
    FilePickerOverlay,
    HistorySearchOverlay,
    ProjectYamlOverlay,
    ShortcutOverlay,
    ToolActivityOverlay,
)
from chemsmart.agent.tui.widgets.slash_palette import SlashCommandPalette
from chemsmart.agent.tui.widgets.transcript import Transcript


class _BusyWorker:
    is_finished = False


class _FinishedWorker:
    is_finished = True


class _TerminalWorker:
    group = "agent-session"
    name = "agent-unified"
    is_finished = True
    error = None
    result: dict[str, object] = {}


def _tool_event(status: str, *, call_id: str = "call-1") -> ToolUseEvent:
    return ToolUseEvent(
        kind=f"tool_use_{status}",
        step_index=1,
        tool="synthesize_command",
        status=status,
        args={"request": "optimize water"} if status == "pending" else {},
        description="Synthesize a grounded ChemSmart CLI command.",
        provider_call_id=call_id,
        payload=(
            {
                "ok": True,
                "status": "ready",
                "command": "chemsmart run gaussian -p test -f h2o.xyz opt",
            }
            if status == "ok"
            else None
        ),
    )


def test_tui_config_uses_defaults_and_rejects_unsafe_bindings(tmp_path):
    path = tmp_path / "agent.yaml"
    path.write_text(
        """
active: fake
providers: {}
tui:
  tool_detail: expanded
  keybindings:
    show_shortcuts: r
    toggle_transcript: ctrl+x
    show_activity: ctrl+x
    unknown_action: f2
""",
        encoding="utf-8",
    )

    config = load_tui_config(path)

    assert (
        config.keybindings["show_shortcuts"]
        == DEFAULT_KEYBINDINGS["show_shortcuts"]
    )
    assert config.keybindings["toggle_transcript"] == "ctrl+x"
    assert (
        config.keybindings["show_activity"]
        == DEFAULT_KEYBINDINGS["show_activity"]
    )
    assert config.tool_detail == "compact"
    assert any("reserved" in issue for issue in config.issues)
    assert any("unknown" in issue for issue in config.issues)
    assert any("duplicate" in issue for issue in config.issues)


def test_tui_state_reducer_updates_shared_state_and_rejects_unknown_fields():
    state = TuiState()
    reducer = TuiStateReducer(state)

    reducer.set_phase(Phase.VALIDATING, "checking generated input")
    reducer.set_tool("validate_runtime", step=3, total=4)
    reducer.set_jobs(queued=2, running=1)
    reducer.update(project="co2", yaml_loaded=True, yaml_label="YAML OK co2")

    assert reducer.state is state
    assert state.phase == Phase.VALIDATING
    assert state.operation == "checking generated input"
    assert state.tool_progress == "tool 3/4 validate_runtime"
    assert state.job_counts == {"queued": 2, "running": 1, "failed": 0}
    assert state.project == "co2"
    with pytest.raises(KeyError, match="Unknown TUI state field"):
        reducer.update(unknown=True)


def test_api_tui_defaults_to_active_agent_runtime(tmp_path):
    app = ChemsmartTuiApp(session_root=tmp_path / "sessions")

    assert app.chat_screen.runtime_v2 == "active"


def test_keyboard_palette_filters_navigates_and_completes(tmp_path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            composer = app.query_one(Composer)
            palette = app.query_one(SlashCommandPalette)
            composer.load_text("/d")
            await pilot.pause()

            assert palette.is_open
            assert palette.selected_item() is not None
            assert palette.selected_item().command == "/dryrun"
            await pilot.press("down")
            await pilot.press("tab")
            await pilot.pause()

            assert composer.text == "/doctor "
            assert not palette.is_open

    asyncio.run(scenario())


def test_disabled_palette_command_does_not_execute(tmp_path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            composer = app.query_one(Composer)
            composer.load_text("/deny")
            await pilot.pause()
            await pilot.press("enter")
            await pilot.pause()

            assert composer.text == "/deny"
            assert not app.chat_screen._pending_approval

    asyncio.run(scenario())


def test_palette_mouse_selection_matches_keyboard_completion(tmp_path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            composer = app.query_one(Composer)
            palette = app.query_one(SlashCommandPalette)
            composer.load_text("/d")
            await pilot.pause()

            doctor_index = next(
                index
                for index, item in enumerate(palette.items)
                if item.command == "/doctor"
            )
            event = OptionList.OptionSelected(palette, doctor_index)
            app.chat_screen.on_option_list_option_selected(event)

            assert composer.text == "/doctor "
            assert composer.has_focus
            assert not palette.is_open

    asyncio.run(scenario())


def test_shortcut_overlay_and_escape_restore_composer(tmp_path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            await pilot.press("f1")
            await pilot.pause()
            assert isinstance(app.screen, ShortcutOverlay)

            await pilot.press("escape")
            await pilot.pause()
            assert app.screen is app.chat_screen
            assert app.query_one(Composer).has_focus

    asyncio.run(scenario())


def test_shift_tab_opens_workspace_yaml_and_escape_restores_composer(
    tmp_path,
    monkeypatch,
):
    workspace = tmp_path / "workspace"
    project_path = workspace / ".chemsmart" / "gaussian" / "water_demo.yaml"
    project_path.parent.mkdir(parents=True)
    project_path.write_text(
        "gas:\n  functional: b3lyp\n  basis: def2svp\n",
        encoding="utf-8",
    )
    monkeypatch.chdir(workspace)

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            composer = app.query_one(Composer)
            assert composer.has_focus

            await pilot.press("shift+tab")
            await pilot.pause()

            assert isinstance(app.screen, ProjectYamlOverlay)
            assert app.screen.path == project_path.resolve()
            assert "functional: b3lyp" in app.screen.yaml_text

            await pilot.press("escape")
            await pilot.pause()
            assert app.screen is app.chat_screen
            assert app.query_one(Composer).has_focus

    asyncio.run(scenario())


def test_shift_tab_cycles_and_selects_multiple_workspace_yamls(
    tmp_path,
    monkeypatch,
):
    workspace = tmp_path / "workspace"
    gaussian_dir = workspace / ".chemsmart" / "gaussian"
    gaussian_dir.mkdir(parents=True)
    first_path = gaussian_dir / "formaldehyde_uv.yaml"
    second_path = gaussian_dir / "water.yaml"
    first_path.write_text(
        "gas:\n  functional: camb3lyp\n  basis: def2svp\n",
        encoding="utf-8",
    )
    second_path.write_text(
        "gas:\n  functional: pbe0\n  basis: def2svp\n",
        encoding="utf-8",
    )
    monkeypatch.chdir(workspace)
    reset_workflow_state()

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            footer = app.query_one(FooterWidget)
            assert footer.state.yaml_label == "YAML SELECT 2"

            await pilot.press("shift+tab")
            await pilot.pause()
            assert isinstance(app.screen, ProjectYamlOverlay)
            assert app.screen.path == first_path.resolve()

            await pilot.press("down")
            await pilot.pause()
            assert app.screen.path == second_path.resolve()
            assert "functional: pbe0" in app.screen.yaml_text

            await pilot.press("enter")
            await pilot.pause()
            assert app.screen is app.chat_screen
            assert current_workflow_state(workspace).project is not None
            assert current_workflow_state(workspace).project.name == "water"
            assert footer.state.yaml_label == "YAML OK gaussian:water"

            await pilot.press("shift+tab")
            await pilot.pause()
            assert isinstance(app.screen, ProjectYamlOverlay)
            assert app.screen.path == second_path.resolve()

    try:
        asyncio.run(scenario())
    finally:
        reset_workflow_state()


def test_written_project_becomes_active_when_other_yaml_exists(
    tmp_path,
    monkeypatch,
):
    workspace = tmp_path / "workspace"
    gaussian_dir = workspace / ".chemsmart" / "gaussian"
    gaussian_dir.mkdir(parents=True)
    (gaussian_dir / "existing.yaml").write_text(
        "gas:\n  functional: b3lyp\n  basis: def2svp\n",
        encoding="utf-8",
    )
    written_path = gaussian_dir / "water.yaml"
    written_path.write_text(
        "gas:\n  functional: pbe0\n  basis: def2svp\n",
        encoding="utf-8",
    )
    monkeypatch.chdir(workspace)
    reset_workflow_state()

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            app.chat_screen._apply_workspace_project_tool_result(
                "write_project_yaml",
                {
                    "ok": True,
                    "project_name": "water",
                    "program": "gaussian",
                    "written_path": str(written_path),
                },
            )
            await pilot.pause()

            footer = app.query_one(FooterWidget)
            assert footer.state.yaml_label == "YAML OK gaussian:water"
            assert current_workflow_state(workspace).project is not None
            assert current_workflow_state(workspace).project.name == "water"

            await pilot.press("shift+tab")
            await pilot.pause()
            assert isinstance(app.screen, ProjectYamlOverlay)
            assert app.screen.path == written_path.resolve()

    try:
        asyncio.run(scenario())
    finally:
        reset_workflow_state()


def test_at_file_picker_refreshes_cold_nested_workspace_index(
    tmp_path,
    monkeypatch,
):
    workspace = tmp_path / "workspace"
    input_path = workspace / "inputs" / "h2o.xyz"
    input_path.parent.mkdir(parents=True)
    input_path.write_text("3\nwater\nO 0 0 0\nH 1 0 0\nH 0 1 0\n")
    monkeypatch.chdir(workspace)

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            await pilot.press("@")
            await pilot.pause()

            assert isinstance(app.screen, FilePickerOverlay)
            for _ in range(20):
                if app.screen._candidates:
                    break
                await pilot.pause(0.05)

            assert [
                path.relative_to(workspace).as_posix()
                for path in app.screen._candidates
            ] == ["inputs/h2o.xyz"]

            await pilot.press("enter")
            await pilot.pause()

            assert app.screen is app.chat_screen
            assert app.query_one(Composer).text == "inputs/h2o.xyz"

    asyncio.run(scenario())


def test_history_search_filters_and_restores_selected_request(tmp_path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            app.chat_screen._request_history = [
                "Prepare a Gaussian TS optimization",
                "Set up an ORCA NEB calculation",
                "Explain the current project YAML",
            ]
            await pilot.press("ctrl+r")
            await pilot.pause()
            assert isinstance(app.screen, HistorySearchOverlay)

            await pilot.press("o", "r", "c", "a")
            await pilot.pause()
            assert app.screen.filtered == ["Set up an ORCA NEB calculation"]
            await pilot.press("enter")
            await pilot.pause()

            assert app.screen is app.chat_screen
            assert app.query_one(Composer).text == (
                "Set up an ORCA NEB calculation"
            )
            assert app.query_one(Composer).has_focus

    asyncio.run(scenario())


def test_tool_activity_shortcut_opens_grouped_call_inventory(tmp_path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            app.chat_screen._apply_tool_use_event(_tool_event("pending"))
            app.chat_screen._apply_tool_use_event(_tool_event("ok"))
            await pilot.pause()
            await pilot.press("ctrl+t")
            await pilot.pause()

            assert isinstance(app.screen, ToolActivityOverlay)
            assert app.screen.entries == [
                app.chat_screen._tool_cells["call-1"].activity_snapshot()
            ]
            assert app.screen.context == {
                "phase": "tool-running",
                "operation": "DRIVING_AUTO",
                "progress": "synthesize_command",
            }

    asyncio.run(scenario())


def test_approval_revise_shortcut_returns_revise_choice(tmp_path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        results: list[ApprovalResult | None] = []
        async with app.run_test() as pilot:
            await pilot.pause()
            app.push_screen(
                ApprovalOverlay(action="run_local", request="optimize water"),
                results.append,
            )
            await pilot.pause()
            await pilot.press("r")
            await pilot.pause()

            assert results == [ApprovalResult("r")]

    asyncio.run(scenario())


def test_one_slot_follow_up_queue_and_waiting_user_gate(tmp_path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            screen = app.chat_screen
            composer = app.query_one(Composer)
            footer = app.query_one(FooterWidget)
            screen._current_worker = _BusyWorker()  # type: ignore[assignment]
            composer.load_text("then add a frequency calculation")

            screen.action_context_tab()
            assert screen._queued_prompt == "then add a frequency calculation"
            assert composer.text == ""
            assert footer.state.queued_prompt

            screen.action_context_tab()
            assert screen._queued_prompt is None
            assert composer.text == "then add a frequency calculation"

            screen.action_context_tab()
            screen._current_worker = _FinishedWorker()  # type: ignore[assignment]
            screen._waiting_for_user = True
            started: list[str] = []
            screen.start_request = started.append  # type: ignore[method-assign]
            screen._maybe_start_queued_prompt()
            assert started == []
            assert screen._queued_prompt is not None

            screen._waiting_for_user = False
            screen._maybe_start_queued_prompt()
            await pilot.pause()
            assert started == ["then add a frequency calculation"]
            assert not footer.state.queued_prompt

    asyncio.run(scenario())


def test_local_clarification_holds_queue_and_uses_waiting_phase(tmp_path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            screen = app.chat_screen
            screen._current_worker = _FinishedWorker()  # type: ignore[assignment]
            screen._queued_prompt = "use charge minus one"
            app.query_one(FooterWidget).set_queued_prompt(True)
            started: list[str] = []
            screen.start_request = started.append  # type: ignore[method-assign]

            screen._publish_synthesis_result(
                {
                    "synthesis": {
                        "status": "needs_clarification",
                        "missing_info": ["charge"],
                    },
                    "provider_type": "local",
                    "provider_model": "mlx-test",
                }
            )
            screen._maybe_start_queued_prompt()
            await pilot.pause()

            assert screen._waiting_for_user
            assert app.query_one(FooterWidget).phase == Phase.WAITING_USER
            assert screen._queued_prompt == "use charge minus one"
            assert started == []

    asyncio.run(scenario())


def test_interrupted_turn_does_not_start_queued_prompt(tmp_path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            screen = app.chat_screen
            worker = _TerminalWorker()
            screen._current_worker = worker  # type: ignore[assignment]
            screen._queued_prompt = "next request"
            started: list[str] = []
            screen.start_request = started.append  # type: ignore[method-assign]

            screen.on_worker_state_changed(
                Worker.StateChanged(worker, WorkerState.CANCELLED)  # type: ignore[arg-type]
            )
            await pilot.pause()

            assert started == []
            assert screen._queued_prompt == "next request"
            assert app.query_one(FooterWidget).phase == Phase.INTERRUPTED

    asyncio.run(scenario())


def test_tool_lifecycle_coalesces_and_preserves_public_evidence(tmp_path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            screen = app.chat_screen
            screen._apply_tool_use_event(_tool_event("pending"))
            screen._apply_tool_use_event(_tool_event("approved"))
            screen._apply_tool_use_event(_tool_event("ok"))
            await pilot.pause()

            cells = app.query_one(Transcript).query_one("#cells")
            tool_cells = [
                cell
                for cell in cells.children
                if isinstance(cell, ToolCallCell)
            ]
            assert len(tool_cells) == 1
            assert tool_cells[0].status == "ok"
            assert tool_cells[0].arguments == {"request": "optimize water"}
            assert tool_cells[0].result is not None
            assert "reasoning" not in tool_cells[0].result

    asyncio.run(scenario())


def test_error_tool_auto_expands_without_stealing_composer_focus(tmp_path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            composer = app.query_one(Composer)
            composer.focus()
            event = _tool_event("error")
            event.reason = "semantic gate rejected the command"
            app.chat_screen._apply_tool_use_event(event)
            await pilot.pause()

            cell = app.chat_screen._tool_cells["call-1"]
            assert cell.expanded
            assert composer.has_focus

    asyncio.run(scenario())


def test_transcript_preserves_turn_order_when_later_cells_arrive(tmp_path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            transcript = app.query_one(Transcript)
            transcript.clear_cells()
            transcript.add_cell(FinalAnswerCell("done", title="Final"))
            cell = ToolCallCell(
                tool="read",
                status="ok",
                description="Read a file.",
                provider_call_id="read-1",
            )
            transcript.add_cell(cell)
            await pilot.pause()

            children = list(transcript.query_one("#cells").children)
            assert isinstance(children[0], FinalAnswerCell)
            assert children[-1] is cell
            assert not cell.expanded
            app.chat_screen.action_toggle_transcript()
            assert cell.expanded
            app.chat_screen.action_toggle_transcript()
            assert not cell.expanded

    asyncio.run(scenario())


def test_late_calculation_completion_updates_original_turn_in_place(tmp_path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            screen = app.chat_screen
            transcript = app.query_one(Transcript)
            transcript.clear_cells()
            transcript.start_turn("turn-1")
            running = CalculationRun(
                run_id="calc-late",
                command="chemsmart run orca -f h2o.xyz -c 0 -m 1 sp",
                cwd=str(tmp_path),
                turn_id="turn-1",
                program="orca",
                kind="sp",
                label="h2o",
                status="running",
                stage="SCF cycle 4",
            )
            screen._on_calculation_run(running.to_dict(), persist=False)
            receipt = screen._calculation_cells["calc-late"]

            transcript.start_turn("turn-2")
            later_answer = FinalAnswerCell(
                "The project YAML is still loaded.", title="Answer"
            )
            transcript.add_cell(later_answer)
            running.status = "completed"
            running.stage = "Calculation completed"
            running.energy = -76.269459371830
            running.normal_termination = True
            screen._on_calculation_run(running.to_dict(), persist=False)
            await pilot.pause()

            children = list(transcript.query_one("#cells").children)
            assert children == [receipt, later_answer]
            assert receipt.run["status"] == "completed"
            assert receipt.run["energy"] == pytest.approx(-76.269459371830)

    asyncio.run(scenario())


def test_one_hundred_late_completions_never_reorder_turn_receipts(tmp_path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            screen = app.chat_screen
            transcript = app.query_one(Transcript)
            transcript.clear_cells()
            expected: list[object] = []
            runs: list[CalculationRun] = []
            for index in range(100):
                turn_id = f"turn-{index}"
                transcript.start_turn(turn_id)
                run = CalculationRun(
                    run_id=f"calc-{index}",
                    command="chemsmart run orca -f h2o.xyz -c 0 -m 1 sp",
                    cwd=str(tmp_path),
                    turn_id=turn_id,
                    program="orca",
                    kind="sp",
                    label=f"h2o-{index}",
                    status="running",
                    stage="SCF running",
                )
                runs.append(run)
                screen._on_calculation_run(run.to_dict(), persist=False)
                receipt = screen._calculation_cells[run.run_id]
                answer = FinalAnswerCell(f"turn {index}", title="Answer")
                transcript.add_cell(answer)
                expected.extend((receipt, answer))
            await pilot.pause()

            for run in reversed(runs):
                run.status = "completed"
                run.stage = "Calculation completed"
                screen._on_calculation_run(run.to_dict(), persist=False)
            await pilot.pause()

            assert list(transcript.query_one("#cells").children) == expected

    asyncio.run(scenario())


def test_calculation_status_strip_and_ctrl_b_share_run_state(tmp_path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            run = CalculationRun(
                run_id="calc-monitor",
                command="chemsmart run orca -f h2o.xyz -c 0 -m 1 sp",
                cwd=str(tmp_path),
                program="orca",
                kind="sp",
                label="h2o",
                project="water",
                method="pbe0",
                basis="def2-svp",
                status="running",
                stage="SCF cycle 6",
                elapsed_s=0.4,
            )
            app.chat_screen._on_calculation_run(run.to_dict(), persist=False)
            await pilot.pause()

            strip = app.query_one(CalculationStatusStrip)
            assert strip.display
            rendered = str(strip.render())
            assert "ORCA" in rendered
            assert "SCF cycle 6" in rendered

            await pilot.press("ctrl+b")
            await pilot.pause()
            assert isinstance(app.screen, CalculationMonitor)
            assert app.screen.runs["calc-monitor"]["status"] == "running"

    asyncio.run(scenario())


def test_local_result_diagnosis_uses_latest_run_without_path_prompt(
    monkeypatch, tmp_path
):
    monkeypatch.chdir(tmp_path)
    session_root = tmp_path / "sessions"
    output = tmp_path / "h2o.out"
    output.write_text(
        "FINAL SINGLE POINT ENERGY -76.269459371830\n"
        "ORCA TERMINATED NORMALLY\n",
        encoding="utf-8",
    )
    run = CalculationRun(
        run_id="calc-diagnose",
        command="chemsmart run orca -f h2o.xyz -c 0 -m 1 sp",
        cwd=str(tmp_path),
        program="orca",
        kind="sp",
        label="h2o",
        status="completed",
        stage="Calculation completed",
        output_path=str(output),
        energy=-76.269459371830,
        normal_termination=True,
    )
    CalculationStore(session_root / "session-one").write_run(run)
    local_config = AgentProviderConfig(
        name="local",
        type="local",
        api_key="",
        model="mlx-local",
        base_url="",
        extra_headers={},
    )
    monkeypatch.setattr(
        "chemsmart.agent.tui.screens.chat._load_tui_provider_config",
        lambda: local_config,
    )

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=session_root)
        async with app.run_test() as pilot:
            await pilot.pause()
            app.chat_screen.start_request("계산 결과를 진단해줘")
            await pilot.pause()
            screenshot = app.export_screenshot()

            assert "-76.269459371830" in screenshot
            assert "normal" in screenshot.lower()
            assert "project name" not in screenshot.lower()

    asyncio.run(scenario())


def test_footer_uses_provider_state_and_responsive_priority(
    monkeypatch, tmp_path
):
    config = AgentProviderConfig(
        name="deepseek",
        type="openai",
        api_key="test",
        model="deepseek-chat",
        base_url="https://example.invalid/v1",
        extra_headers={},
        project="co2",
    )
    monkeypatch.setattr(
        "chemsmart.agent.tui.screens.chat._load_tui_provider_config",
        lambda: config,
    )

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test(size=(160, 45)) as pilot:
            await pilot.pause()
            footer = app.query_one(FooterWidget)
            footer.set_phase(Phase.TOOL_RUNNING)
            footer.set_hint("validate_runtime")
            footer.set_queued_prompt(True)
            footer.set_usage(input_tokens=120, output_tokens=30)
            await pilot.pause()

            text = str(footer.renderable)
            assert "deepseek-chat" in text
            assert "project co2" in text
            assert "next queued" in text
            assert "usage 120/30" in text

    asyncio.run(scenario())


def test_worker_usage_updates_shared_footer_state(tmp_path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            app.chat_screen._sync_footer_usage(
                {"total_input_tokens": 211, "total_output_tokens": 34}
            )

            footer = app.query_one(FooterWidget)
            assert footer.state is app.chat_screen.tui_state
            assert footer.state.usage_input_tokens == 211
            assert footer.state.usage_output_tokens == 34

    asyncio.run(scenario())


@pytest.mark.parametrize("viewport", [(80, 24), (120, 36), (160, 45)])
def test_responsive_viewports_preserve_vertical_widget_boundaries(
    viewport, tmp_path
):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test(size=viewport) as pilot:
            await pilot.pause()
            app.chat_screen._on_calculation_run(
                CalculationRun(
                    run_id="calc-viewport",
                    command="chemsmart run orca -f h2o.xyz -c 0 -m 1 sp",
                    cwd=str(tmp_path),
                    program="orca",
                    kind="sp",
                    label="h2o",
                    status="running",
                    stage="SCF cycle 5",
                ).to_dict(),
                persist=False,
            )
            await pilot.pause()
            transcript = app.query_one(Transcript)
            calculation_strip = app.query_one(CalculationStatusStrip)
            composer = app.query_one(Composer)
            footer = app.query_one(FooterWidget)
            screenshot = app.export_screenshot()

            assert transcript.region.bottom <= calculation_strip.region.y
            assert calculation_strip.region.bottom <= composer.region.y
            assert composer.region.bottom <= footer.region.y
            assert footer.region.bottom <= viewport[1]
            assert "YAML" in screenshot

    asyncio.run(scenario())


def test_one_hundred_tool_events_update_within_budget(tmp_path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            screen = app.chat_screen
            screen._apply_tool_use_event(_tool_event("pending"))
            timings: list[float] = []
            for index in range(100):
                event = _tool_event("approved" if index % 2 == 0 else "ok")
                start = perf_counter()
                screen._apply_tool_use_event(event)
                timings.append((perf_counter() - start) * 1000)
            await pilot.pause()

            p95 = quantiles(timings, n=20)[18]
            assert p95 < 50
            assert len(screen._tool_cells) == 1

    asyncio.run(scenario())


def _ready_orca_payload(command: str) -> dict[str, object]:
    return {
        "ok": True,
        "status": "ready",
        "command": command,
        "confidence": "high",
        "semantic": {
            "verdict": "ok",
            "failed_rule_ids": [],
            "generated_inputs": [
                {
                    "path": "/tmp/h2o_sp_fake.inp",
                    "route": "! pbe0 def2-svp",
                }
            ],
        },
        "intent": {
            "verdict": "ok",
            "failed_rule_ids": [],
        },
    }


def test_run_executes_latest_validated_command_without_mode_switch(
    monkeypatch, tmp_path
):
    command = (
        "chemsmart run -n 1 --no-scratch orca -p orca_sp_demo "
        "-f inputs/h2o.xyz -c 0 -m 1 sp"
    )
    project_dir = tmp_path / ".chemsmart" / "orca"
    project_dir.mkdir(parents=True)
    (project_dir / "orca_sp_demo.yaml").write_text(
        "gas:\n  functional: pbe0\n  basis: def2-svp\n",
        encoding="utf-8",
    )
    monkeypatch.chdir(tmp_path)

    calls: list[tuple[str, bool]] = []

    def fake_execute(
        value,
        *,
        test,
        calculation_context,
        event_sink,
    ):
        calls.append((value, test))
        run = CalculationRun(
            run_id="calc-orca-sp",
            command=value,
            cwd=str(tmp_path),
            session_id=calculation_context.session_id,
            turn_id=calculation_context.turn_id,
            program="orca",
            kind="sp",
            label="h2o",
            project="orca_sp_demo",
            method="pbe0",
            basis="def2-svp",
            status="running",
            stage="SCF cycle 8",
        )
        event_sink(CalculationEvent("started", run))
        run.status = "completed"
        run.stage = "Calculation completed"
        run.returncode = 0
        run.energy = -76.0
        run.scf_cycles = 8
        run.normal_termination = True
        run.output_path = str(tmp_path / "h2o.out")
        event_sink(CalculationEvent("completed", run))
        return {
            "ok": True,
            "status": "ok",
            "command": value,
            "returncode": 0,
            "stdout_tail": "RAW OUTPUT MUST NOT BE IN TRANSCRIPT",
            "stderr_tail": "",
            "semantic": {"verdict": "ok"},
            "calculation": run.to_dict(),
        }

    monkeypatch.setattr(
        "chemsmart.agent.tui.mixins.calculations.execute_chemsmart_command_observed",
        fake_execute,
    )

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            screen = app.chat_screen
            screen._publish_command_tool_result(
                "synthesize_command",
                _ready_orca_payload(command),
            )
            assert screen._ready_command is not None
            assert screen._interaction_mode in {"synthesis", "unified"}

            screen._handle_slash_command("/run")
            for _ in range(50):
                await pilot.pause(0.05)
                if calls and screen._ready_command is None:
                    break

            assert calls == [(command, False)]
            receipt_cells = [
                cell
                for cell in app.query_one(Transcript)
                .query_one("#cells")
                .children
                if isinstance(cell, CalculationReceiptCell)
            ]
            assert len(receipt_cells) == 1
            assert receipt_cells[0].run["status"] == "completed"
            assert receipt_cells[0].run["energy"] == -76.0
            assert (
                "RAW OUTPUT MUST NOT BE IN TRANSCRIPT"
                not in app.export_screenshot()
            )

    asyncio.run(scenario())


def test_second_local_run_requires_explicit_yes(monkeypatch, tmp_path):
    command = (
        "chemsmart run orca -p orca_sp_demo -f inputs/h2o.xyz -c 0 -m 1 sp"
    )
    project_dir = tmp_path / ".chemsmart" / "orca"
    project_dir.mkdir(parents=True)
    (project_dir / "orca_sp_demo.yaml").write_text(
        "gas:\n  functional: pbe0\n  basis: def2-svp\n",
        encoding="utf-8",
    )
    monkeypatch.chdir(tmp_path)

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            screen = app.chat_screen
            screen._publish_command_tool_result(
                "synthesize_command", _ready_orca_payload(command)
            )
            screen._calculation_runs["calc-active"] = {
                "run_id": "calc-active",
                "status": "running",
                "execution_mode": "local",
            }
            calls: list[tuple[object, dict[str, object]]] = []
            screen.run_calculation_request = (  # type: ignore[method-assign]
                lambda ready, **kwargs: calls.append((ready, kwargs))
            )

            screen._handle_slash_command("/run")
            await pilot.pause()
            assert calls == []
            assert screen._ready_command is not None

            screen._handle_slash_command("/run yes")
            await pilot.pause()
            assert len(calls) == 1
            assert screen._ready_command is None

    asyncio.run(scenario())


def test_run_rejects_stale_project_yaml(monkeypatch, tmp_path):
    command = (
        "chemsmart run orca -p orca_sp_demo -f inputs/h2o.xyz -c 0 -m 1 sp"
    )
    project_dir = tmp_path / ".chemsmart" / "orca"
    project_dir.mkdir(parents=True)
    project_path = project_dir / "orca_sp_demo.yaml"
    project_path.write_text(
        "gas:\n  functional: pbe0\n  basis: def2-svp\n",
        encoding="utf-8",
    )
    monkeypatch.chdir(tmp_path)

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            screen = app.chat_screen
            screen._publish_command_tool_result(
                "synthesize_command",
                _ready_orca_payload(command),
            )
            project_path.write_text(
                "gas:\n  functional: b3lyp\n  basis: def2-svp\n",
                encoding="utf-8",
            )
            calls: list[object] = []
            screen.run_slash_tool_request = (  # type: ignore[method-assign]
                lambda *args, **kwargs: calls.append((args, kwargs))
            )

            screen._handle_slash_command("/run")
            await pilot.pause()

            assert calls == []
            assert "changed after validation" in screen._ready_command_problem(
                "run"
            )

    asyncio.run(scenario())
