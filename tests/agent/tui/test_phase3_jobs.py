from __future__ import annotations

import asyncio
from pathlib import Path

from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.screens.calculations import CalculationMonitor
from chemsmart.agent.tui.services import job_poller
from chemsmart.agent.tui.widgets.composer import Composer

from ._helpers import normalize_svg, write_session_fixture


def _job_snapshot(status: str = "queued") -> dict[str, dict]:
    return {
        "12345.remote": {
            "job_id": "12345.remote",
            "session_id": "session-001",
            "name": "water_opt",
            "scheduler": "slurm",
            "status": status,
            "started": "2026-05-09 00:10",
            "runtime": "5m 12s",
            "host": "remote-hpc",
            "server_name": "remote-hpc",
            "session_dir": "/tmp/session-001",
            "cwd": str(Path.cwd()),
            "output_path": "/tmp/water_opt.log",
            "raw_started": "2026-05-09T00:10:00+00:00",
        }
    }


def test_jobs_panel_keys_and_polling_diff(monkeypatch, tmp_path: Path):
    session_root = tmp_path / "sessions"
    write_session_fixture(session_root)
    snapshots = [_job_snapshot("queued"), _job_snapshot("running")]

    def fake_collect_job_snapshot(_session_root):
        return snapshots[0]

    monkeypatch.setattr(
        "chemsmart.agent.tui.services.job_poller.collect_job_snapshot",
        fake_collect_job_snapshot,
    )
    job_poller.refresh_job_snapshot_cache(session_root)

    async def scenario() -> None:
        app = ChemsmartTuiApp(
            session_root=session_root,
            job_poll_interval=60.0,
        )
        cancel_calls: list[str] = []
        extract_calls: list[str] = []
        app.chat_screen._confirm_cancel = lambda job_id: cancel_calls.append(job_id)  # type: ignore[method-assign]
        app.chat_screen._extract_job_result = lambda job_id: extract_calls.append(job_id)  # type: ignore[method-assign]
        async with app.run_test() as pilot:
            await pilot.pause()
            app.chat_screen.action_open_jobs_panel()
            await pilot.pause()
            assert isinstance(app.screen, CalculationMonitor)
            assert app.screen.jobs["12345.remote"]["status"] == "queued"

            snapshots.pop(0)
            app.chat_screen._refresh_job_snapshot()
            await pilot.pause(0.2)
            app.chat_screen._refresh_job_snapshot()
            await pilot.pause()
            assert isinstance(app.screen, CalculationMonitor)
            assert app.screen.jobs["12345.remote"]["status"] == "running"

            await pilot.press("c")
            await pilot.pause()
            assert cancel_calls == ["12345.remote"]

            app.chat_screen.action_open_jobs_panel()
            await pilot.pause()
            await pilot.press("e")
            await pilot.pause()
            assert extract_calls == ["12345.remote"]

            app.chat_screen.action_open_jobs_panel()
            await pilot.pause()
            await pilot.press("enter")
            await pilot.pause()
            assert isinstance(app.screen, CalculationMonitor)

    asyncio.run(scenario())


def test_jobs_plain_mode_falls_back_to_inline_table(
    monkeypatch, tmp_path: Path
):
    session_root = tmp_path / "sessions"
    write_session_fixture(session_root)
    snapshot = _job_snapshot("queued")

    monkeypatch.setattr(
        "chemsmart.agent.tui.services.job_poller.collect_job_snapshot",
        lambda _session_root: snapshot,
    )
    job_poller.refresh_job_snapshot_cache(session_root)

    async def scenario() -> None:
        app = ChemsmartTuiApp(
            plain=True,
            session_root=session_root,
            job_poll_interval=60.0,
        )
        async with app.run_test() as pilot:
            await pilot.pause()
            composer = app.query_one(Composer)
            composer.load_text("/jobs")
            await pilot.press("enter")
            await pilot.pause()
            assert not isinstance(app.screen, CalculationMonitor)
            assert "12345.remote" in normalize_svg(app.export_screenshot())

    asyncio.run(scenario())
