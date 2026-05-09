from __future__ import annotations

import asyncio
from pathlib import Path

from rich.console import Console
from textual.widgets import Static

from chemsmart.agent.core import _preview_value
from chemsmart.agent.tools import build_molecule
from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.widgets.cells import (
    DryRunInputCell,
    ErrorCell,
    PlanCell,
    RuntimeValidationCell,
)
from chemsmart.agent.tui.widgets.composer import Composer
from chemsmart.agent.tui.widgets.popups import (
    ApprovalOverlay,
    CwdMismatchOverlay,
)
from chemsmart.agent.tui.widgets.transcript import Transcript

from .._agent_session_helpers import planner_plan


def _render_plain(renderable) -> str:
    console = Console(width=120, record=True)
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
            assert "확인 중" in running_text

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
            assert single_molecule_xyz_file in done_text

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
    assert "로컬 OK / 원격 정보 4개 필요" in text
    assert "queue / account / scratch_dir / modules_or_executable_path" in text
    assert "- server.queue required" not in text


def test_dry_run_input_cell_shows_path_before_body():
    cell = DryRunInputCell(
        "%chk=water.chk\n# B3LYP/6-31G(d)\n",
        inputfile="/tmp/water.com",
    )

    text = _render_plain(cell.renderable)
    assert "✓ water.com ready" in text
    assert "path: /tmp/water.com" in text
    assert text.index("path: /tmp/water.com") < text.index("%chk=water.chk")


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


def test_greeter_copy_mentions_doctor_ask_and_run(tmp_path: Path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            transcript = app.query_one(Transcript).query_one("#cells")
            first = transcript.children[0]
            text = _render_plain(first.renderable)
            assert "/doctor" in text
            assert "ask 모드" in text
            assert "run 모드" in text

    asyncio.run(scenario())


def test_clear_message_uses_designer_copy(tmp_path: Path):
    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test() as pilot:
            await pilot.pause()
            composer = app.query_one(Composer)
            composer.load_text("/clear")
            await pilot.press("enter")
            await pilot.pause()
            transcript = app.query_one(Transcript).query_one("#cells")
            first = transcript.children[0]
            text = _render_plain(first.renderable)
            assert "대화가 비워졌습니다." in text
            assert "/sessions" in text

    asyncio.run(scenario())


def test_missing_file_error_uses_exact_copy():
    cell = ErrorCell(
        "FileNotFoundError",
        "[Errno 2] No such file or directory: 'missing.xyz'",
        {"tool": "build_molecule"},
    )

    text = _render_plain(cell.renderable)
    assert "입력 파일을 찾을 수 없습니다" in text
    assert (
        "요청에 적힌 경로가 현재 작업 폴더 기준으로 존재하지 않습니다." in text
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
            assert "작업 디렉터리가 다릅니다" in text
            assert "기록된 폴더로 이동 후 재개" in text
            assert "현재 폴더에서 강제로 계속" in text

    asyncio.run(scenario())
