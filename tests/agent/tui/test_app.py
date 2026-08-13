"""A terminal-level smoke observation for the production Textual shell."""

from __future__ import annotations

import asyncio
from pathlib import Path

import pytest

pytest.importorskip("textual")

from chemsmart.agent.tui.app import ChemSmartAgentApp  # noqa: E402
from chemsmart.agent.tui.controller import (  # noqa: E402
    AgentSessionConfigV1,
    AgentTuiController,
)


def test_tui_launches_and_exposes_explicit_approval_chain(tmp_path: Path):
    secret = tmp_path / "secret.env"
    secret.write_text("PROVIDER_KEY=not-used\n", encoding="utf-8")
    envelope = tmp_path / "execution-envelope.json"
    envelope.write_text("{}\n", encoding="utf-8")
    app = ChemSmartAgentApp(
        AgentTuiController(
            AgentSessionConfigV1(
                workspace=tmp_path,
                secret_file=secret,
                execution_envelope_file=envelope,
                review_file=tmp_path / "review-output.json",
            )
        ),
        plain=True,
    )

    async def observe() -> None:
        async with app.run_test(size=(120, 35)) as pilot:
            composer = app.query_one("#composer")
            composer.value = "/help"
            await pilot.press("enter")
            await pilot.pause()

            transcript = app.query_one("#transcript")
            assert len(transcript.lines) > 10
            assert app.controller.phase.value == "ready"

    asyncio.run(observe())
