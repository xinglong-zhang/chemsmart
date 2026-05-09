from __future__ import annotations

import asyncio
from pathlib import Path

import pytest

from chemsmart.agent.tui.app import ChemsmartTuiApp
from chemsmart.agent.tui.widgets.header import ChemsmartHeader
from chemsmart.agent.tui.widgets.transcript import Transcript


@pytest.mark.parametrize("width", [40, 80])
def test_chat_header_renders_wordmark_and_sits_above_transcript(
    monkeypatch, tmp_path: Path, width: int
):
    monkeypatch.setattr(
        "chemsmart.agent.tui.services.job_poller.collect_job_snapshot",
        lambda _session_root: {},
    )

    async def scenario() -> None:
        app = ChemsmartTuiApp(session_root=tmp_path / "sessions")
        async with app.run_test(size=(width, 24)) as pilot:
            await pilot.pause()
            header = app.query_one(ChemsmartHeader)
            transcript = app.query_one(Transcript)
            assert (
                ChemsmartHeader.normalize_wordmark(header.render().plain)
                == ChemsmartHeader.plain_wordmark
            )
            assert header.region.y < transcript.region.y

    asyncio.run(scenario())


def test_chat_header_plain_mode_falls_back_to_ascii(
    monkeypatch, tmp_path: Path
):
    monkeypatch.setattr(
        "chemsmart.agent.tui.services.job_poller.collect_job_snapshot",
        lambda _session_root: {},
    )

    async def scenario() -> None:
        app = ChemsmartTuiApp(plain=True, session_root=tmp_path / "sessions")
        async with app.run_test(size=(80, 24)) as pilot:
            await pilot.pause()
            header = app.query_one(ChemsmartHeader)
            wordmark = header.render()
            assert wordmark.plain == ChemsmartHeader.plain_wordmark
            assert wordmark.spans == []

    asyncio.run(scenario())
