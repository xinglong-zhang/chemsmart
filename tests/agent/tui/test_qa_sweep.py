from __future__ import annotations

import asyncio
from pathlib import Path

import pytest

from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.phase import Phase
from chemsmart.agent.tui.services.job_poller import (
    _cluster_running_jobs,
    collect_job_snapshot,
    format_jobs_table,
)
from chemsmart.agent.tui.widgets.popups import ApprovalOverlay
from chemsmart.agent.tui.widgets.transcript import Transcript


def _last_cell(app: ChemsmartTuiApp):
    cells = app.query_one(Transcript).query_one("#cells")
    return list(cells.children)[-1]


def _cell_text(cell) -> str:
    source_text = getattr(cell, "source_text", None)
    if isinstance(source_text, str):
        return source_text
    message = getattr(cell, "message", None)
    error_title = getattr(cell, "error_title", None)
    if isinstance(error_title, str):
        return f"{error_title} {message or ''}"
    if isinstance(message, str):
        return message
    return str(getattr(cell, "renderable", ""))


def test_format_jobs_table_handles_empty_rows():
    assert format_jobs_table([]) == "No jobs found."


def test_collect_job_snapshot_raises_for_malformed_state(tmp_path: Path):
    session_dir = tmp_path / "sessions" / "bad-session"
    session_dir.mkdir(parents=True)
    (session_dir / "state.json").write_text("{", encoding="utf-8")

    with pytest.raises(RuntimeError, match="Malformed session state"):
        collect_job_snapshot(tmp_path / "sessions")


def test_cluster_running_jobs_surfaces_scheduler_failures(monkeypatch):
    monkeypatch.setattr(
        "chemsmart.agent.tui.services.job_poller.ClusterHelper.get_gaussian_running_jobs",
        lambda self: (_ for _ in ()).throw(RuntimeError("scheduler offline")),
    )

    with pytest.raises(RuntimeError, match="Scheduler queue lookup failed"):
        _cluster_running_jobs()


@pytest.mark.parametrize(
    ("command", "phase", "expected"),
    [
        ("/sessions", Phase.IDLE, "No sessions found."),
        ("/resume", Phase.IDLE, "No sessions found."),
        ("/cancel 12345.remote", Phase.IDLE, "Confirmation required"),
    ],
)
def test_plain_mode_slash_commands_avoid_overlays(
    command: str,
    phase: Phase,
    expected: str,
    monkeypatch,
    tmp_path: Path,
):
    monkeypatch.setattr(
        "chemsmart.agent.tui.services.job_poller.collect_job_snapshot",
        lambda _session_root: {},
    )

    async def scenario() -> None:
        app = ChemsmartTuiApp(plain=True, session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            app.chat_screen.query_one("#status-footer").set_phase(phase)
            app.chat_screen._handle_slash_command(command)
            await pilot.pause()
            assert app.screen is app.chat_screen
            assert expected in _cell_text(_last_cell(app))

    asyncio.run(scenario())


def test_plain_mode_submit_and_run_commands_stay_inline(tmp_path: Path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(plain=True, session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            app.chat_screen._pending_approval = True
            app.chat_screen._pending_tool_request = type(
                "Req",
                (),
                {"name": "submit_hpc", "arguments": {}},
            )()
            app.chat_screen._pending_approval_index = 1
            app.chat_screen._pending_approval_total = 1
            app.chat_screen._approval_waiter = None
            app.chat_screen._handle_slash_command("/submit")
            await pilot.pause()
            assert app.screen is app.chat_screen

            app.chat_screen._pending_tool_request = type(
                "Req",
                (),
                {"name": "run_local", "arguments": {}},
            )()
            app.chat_screen._handle_slash_command("/run")
            await pilot.pause()
            assert app.screen is app.chat_screen

    asyncio.run(scenario())


@pytest.mark.parametrize(
    ("command", "phase", "expected"),
    [
        ("/queue", Phase.IDLE, "No jobs found."),
        ("/cancel", Phase.IDLE, "Usage: /cancel <job-id> [yes]"),
        ("/extract", Phase.IDLE, "Usage: /extract <job-id|inputfile>"),
        ("/molecule", Phase.IDLE, "Usage: /molecule <path>"),
        ("/dryrun", Phase.DRY_RUN_READY, "There is no active request."),
        (
            "/critic",
            Phase.PLANNING,
            "The critic has not finished yet.",
        ),
        ("/plan", Phase.PLANNING, "The planner has not finished yet."),
        ("/rationale", Phase.IDLE, "No planner rationale is available."),
    ],
)
def test_slash_commands_report_empty_states(
    command: str,
    phase: Phase,
    expected: str,
    monkeypatch,
    tmp_path: Path,
):
    monkeypatch.setattr(
        "chemsmart.agent.tui.services.job_poller.collect_job_snapshot",
        lambda _session_root: {},
    )

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            app.chat_screen.query_one("#status-footer").set_phase(phase)
            if command == "/submit":
                app.chat_screen._pending_approval = True
                app.chat_screen._pending_risky_tool = "submit_hpc"
            if command == "/run":
                app.chat_screen._pending_approval = True
                app.chat_screen._pending_risky_tool = "run_local"
            app.chat_screen._handle_slash_command(command)
            await pilot.pause()
            assert expected in _cell_text(_last_cell(app))

    asyncio.run(scenario())


def test_job_snapshot_failure_surfaces_error_cell(monkeypatch, tmp_path: Path):
    monkeypatch.setattr(
        "chemsmart.agent.tui.services.job_poller.collect_job_snapshot",
        lambda _session_root: (_ for _ in ()).throw(
            RuntimeError("broken job snapshot")
        ),
    )

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause(0.2)
            cell = _last_cell(app)
            assert cell.error_title == "Job poller failed"
            assert cell.message == "Job poller failed"
            assert "broken job snapshot" in str(cell.details)

    asyncio.run(scenario())


def test_approval_overlay_escape_closes_when_revision_input_focused(
    tmp_path: Path,
):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            app.push_screen(
                ApprovalOverlay(action="run_local", request="optimize h2o"),
            )
            await pilot.pause()
            assert isinstance(app.screen, ApprovalOverlay)
            await pilot.press("escape")
            await pilot.pause()
            assert app.screen is app.chat_screen

    asyncio.run(scenario())


def test_ctrl_c_first_press_arms_then_second_press_exits_when_idle(
    tmp_path: Path,
):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        exit_calls: list[bool] = []
        async with app.run_test() as pilot:
            await pilot.pause()
            app.exit = lambda *a, **kw: exit_calls.append(True)  # type: ignore[method-assign]
            chat = app.chat_screen
            assert chat._current_worker is None
            chat.action_soft_cancel()
            await pilot.pause()
            assert exit_calls == []
            assert chat._quit_armed is True
            chat.action_soft_cancel()
            await pilot.pause()
            assert exit_calls == [True]

    asyncio.run(scenario())


def test_ctrl_c_disarms_after_timer_expires(tmp_path: Path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            chat = app.chat_screen
            chat.action_soft_cancel()
            await pilot.pause()
            assert chat._quit_armed is True
            chat._disarm_soft_cancel()
            await pilot.pause()
            assert chat._quit_armed is False
            assert chat._quit_timer is None

    asyncio.run(scenario())


def test_ctrl_c_keypress_dispatches_via_app_to_screen(tmp_path: Path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        exit_calls: list[bool] = []
        async with app.run_test() as pilot:
            await pilot.pause()
            app.exit = lambda *a, **kw: exit_calls.append(True)  # type: ignore[method-assign]
            chat = app.chat_screen
            await pilot.press("ctrl+c")
            await pilot.pause()
            assert chat._quit_armed is True
            assert exit_calls == []
            await pilot.press("ctrl+c")
            await pilot.pause()
            assert exit_calls == [True]

    asyncio.run(scenario())
