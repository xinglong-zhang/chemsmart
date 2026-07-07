from __future__ import annotations

import asyncio
import threading
import time
from pathlib import Path
from types import SimpleNamespace

from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.services import file_index, job_poller
from chemsmart.agent.tui.widgets.transcript import Transcript


def _last_cell_text(app: ChemsmartTuiApp) -> str:
    cells = app.query_one(Transcript).query_one("#cells")
    cell = list(cells.children)[-1]
    source_text = getattr(cell, "source_text", None)
    if isinstance(source_text, str):
        return source_text
    message = getattr(cell, "message", None)
    if isinstance(message, str):
        return message
    return str(getattr(cell, "renderable", ""))


def _all_cell_texts(app: ChemsmartTuiApp) -> list[str]:
    cells = app.query_one(Transcript).query_one("#cells")
    return [
        getattr(cell, "source_text", None)
        or getattr(cell, "message", None)
        or str(getattr(cell, "renderable", ""))
        for cell in cells.children
    ]


def test_tools_slash_command_returns_immediately_and_off_main_thread(
    monkeypatch,
    tmp_path: Path,
):
    call_threads: list[int] = []

    def slow_invoke(self, command, args, catch_exceptions=False):
        call_threads.append(threading.get_ident())
        time.sleep(0.2)
        return SimpleNamespace(
            output="registered tools: 10\n",
            exit_code=0,
        )

    monkeypatch.setattr(
        "chemsmart.agent.tui.screens.chat.CliRunner.invoke",
        slow_invoke,
    )

    async def scenario() -> None:
        app = ChemsmartTuiApp(
            session_root=tmp_path / "sessions",
            job_poll_interval=60.0,
        )
        main_thread_id = threading.get_ident()
        async with app.run_test() as pilot:
            await pilot.pause()
            started = time.perf_counter()
            app.chat_screen._handle_slash_command("/tools")
            elapsed = time.perf_counter() - started
            assert elapsed < 0.05
            await pilot.pause(0.35)
            assert call_threads == [call_threads[0]]
            assert call_threads[0] != main_thread_id
            assert "registered tools: 10" in _last_cell_text(app)

    asyncio.run(scenario())


def test_inline_slash_commands_use_semantic_assertions_and_strip_debug_hint(
    monkeypatch,
    tmp_path: Path,
):
    def fake_invoke(self, command, args, catch_exceptions=False):
        if args == ["tools"]:
            return SimpleNamespace(
                output=(
                    "registered tools: 10\n"
                    "Debug log saved to: /tmp/chemsmart-agent.log\n"
                ),
                exit_code=0,
            )
        return SimpleNamespace(
            output=(
                "AI_PROVIDER=openai OK\n"
                "api.env: OK (key length=32)\n"
                "tools registered: 10\n"
            ),
            exit_code=0,
        )

    monkeypatch.setattr(
        "chemsmart.agent.tui.screens.chat.CliRunner.invoke",
        fake_invoke,
    )

    async def scenario() -> None:
        app = ChemsmartTuiApp(
            session_root=tmp_path / "sessions",
            job_poll_interval=60.0,
        )
        async with app.run_test() as pilot:
            await pilot.pause()
            app.chat_screen._handle_slash_command("/tools")
            await pilot.pause(0.1)
            tools_text = _last_cell_text(app)
            assert "registered tools: 10" in tools_text
            assert "Debug log saved to:" not in tools_text

            app.chat_screen._handle_slash_command("/doctor")
            await pilot.pause(0.1)
            doctor_text = _last_cell_text(app)
            assert "AI_PROVIDER=openai OK" in doctor_text
            assert "tools registered: 10" in doctor_text

    asyncio.run(scenario())


def test_cancel_confirmation_returns_immediately_and_off_main_thread(
    monkeypatch,
    tmp_path: Path,
):
    call_threads: list[int] = []

    def slow_cancel(_job):
        call_threads.append(threading.get_ident())
        time.sleep(0.2)
        return {
            "command": ["scancel", "12345.remote"],
            "stdout": "",
            "stderr": "",
        }

    monkeypatch.setattr(
        "chemsmart.agent.tui.screens.chat.cancel_job",
        slow_cancel,
    )

    async def scenario() -> None:
        app = ChemsmartTuiApp(
            session_root=tmp_path / "sessions",
            job_poll_interval=60.0,
        )
        main_thread_id = threading.get_ident()
        async with app.run_test() as pilot:
            await pilot.pause()
            app.chat_screen._job_snapshot["12345.remote"] = {
                "job_id": "12345.remote",
                "scheduler": "slurm",
            }
            started = time.perf_counter()
            app.chat_screen._handle_cancel_confirmation(
                "12345.remote",
                "yes",
            )
            elapsed = time.perf_counter() - started
            assert elapsed < 0.05
            await pilot.pause(0.35)
            assert call_threads == [call_threads[0]]
            assert call_threads[0] != main_thread_id
            assert "12345.remote" in app.chat_screen._cancelled_job_ids
            assert any(
                "Cancelled `12345.remote`" in text
                for text in _all_cell_texts(app)
            )

    asyncio.run(scenario())


def test_job_snapshot_cache_returns_immediately_when_cold(
    monkeypatch,
    tmp_path: Path,
):
    call_threads: list[int] = []

    def slow_collect(_session_root):
        call_threads.append(threading.get_ident())
        time.sleep(0.2)
        return {
            "job-1": {
                "job_id": "job-1",
                "status": "queued",
            }
        }

    monkeypatch.setattr(job_poller, "collect_job_snapshot", slow_collect)

    started = time.perf_counter()
    cold_snapshot = job_poller.get_cached_job_snapshot(tmp_path / "sessions")
    future = job_poller.request_job_snapshot_refresh(tmp_path / "sessions")
    elapsed = time.perf_counter() - started

    assert cold_snapshot == {}
    assert elapsed < 0.05
    assert future is not None
    future.result(timeout=1.0)
    warm_snapshot = job_poller.get_cached_job_snapshot(tmp_path / "sessions")
    assert warm_snapshot["job-1"]["status"] == "queued"
    assert call_threads == [call_threads[0]]
    assert call_threads[0] != threading.get_ident()


def test_refresh_job_snapshot_returns_immediately_when_cache_is_cold(
    monkeypatch,
    tmp_path: Path,
):
    call_threads: list[int] = []

    def slow_collect(_session_root):
        call_threads.append(threading.get_ident())
        time.sleep(0.2)
        return {
            "job-2": {
                "job_id": "job-2",
                "status": "running",
            }
        }

    monkeypatch.setattr(job_poller, "collect_job_snapshot", slow_collect)

    async def scenario() -> None:
        app = ChemsmartTuiApp(
            session_root=tmp_path / "sessions",
            job_poll_interval=60.0,
        )
        async with app.run_test() as pilot:
            await pilot.pause()
            started = time.perf_counter()
            app.chat_screen._refresh_job_snapshot()
            elapsed = time.perf_counter() - started
            assert elapsed < 0.05
            assert app.chat_screen._job_snapshot == {}
            await pilot.pause(0.35)
            app.chat_screen._refresh_job_snapshot()
            assert app.chat_screen._job_snapshot["job-2"]["status"] == (
                "running"
            )

    asyncio.run(scenario())
    assert call_threads


def test_iter_candidate_files_returns_immediately_when_cold(
    monkeypatch,
    tmp_path: Path,
):
    (tmp_path / "water.xyz").write_text("3\nwater\n")
    original_rglob = Path.rglob
    call_threads: list[int] = []

    def slow_rglob(self, pattern):
        call_threads.append(threading.get_ident())
        time.sleep(0.2)
        return original_rglob(self, pattern)

    monkeypatch.setattr(Path, "rglob", slow_rglob)

    started = time.perf_counter()
    cold_candidates = file_index.iter_candidate_files(tmp_path)
    future = file_index.request_candidate_file_refresh(tmp_path)
    elapsed = time.perf_counter() - started

    assert cold_candidates == []
    assert elapsed < 0.05
    assert future is not None
    future.result(timeout=1.0)
    warm_candidates = file_index.iter_candidate_files(tmp_path)
    assert warm_candidates[0].name == "water.xyz"
    assert call_threads
    assert call_threads[0] != threading.get_ident()


def test_open_file_picker_does_not_rglob_synchronously(
    monkeypatch,
    tmp_path: Path,
):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "water.xyz").write_text("3\nwater\n")
    original_rglob = Path.rglob
    call_threads: list[int] = []

    def slow_rglob(self, pattern):
        call_threads.append(threading.get_ident())
        time.sleep(0.2)
        return original_rglob(self, pattern)

    monkeypatch.setattr(Path, "rglob", slow_rglob)

    async def scenario() -> None:
        app = ChemsmartTuiApp(
            session_root=tmp_path / "sessions",
            job_poll_interval=60.0,
        )
        async with app.run_test() as pilot:
            await pilot.pause()
            started = time.perf_counter()
            app.chat_screen.open_file_picker()
            elapsed = time.perf_counter() - started
            assert elapsed < 0.05
            await pilot.pause(0.1)

    asyncio.run(scenario())
    assert call_threads
    assert call_threads[0] != threading.get_ident()
