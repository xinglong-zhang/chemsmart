"""While a session plans, its tool calls narrate live instead of at the end.

The feed is driven by the same append-only events.jsonl the session writes;
a started row mutates in place when its call settles, and the settled text
replaces the pending phrase.
"""

from __future__ import annotations

import asyncio
import json
from pathlib import Path

import pytest

pytest.importorskip("textual")

from chemsmart.agent.tui.app import ChemSmartAgentApp  # noqa: E402
from chemsmart.agent.tui.controller import (  # noqa: E402
    AgentSessionConfigV1,
    AgentTuiController,
)
from chemsmart.agent.tui.transcript import ToolRow  # noqa: E402


def _line(kind: str, payload: dict) -> str:
    return json.dumps({"kind": kind, "payload": payload}) + "\n"


def test_started_rows_settle_in_place_as_their_calls_finish(tmp_path: Path):
    secret = tmp_path / "secret.env"
    secret.write_text("PROVIDER_KEY=not-used\n", encoding="utf-8")
    app = ChemSmartAgentApp(
        AgentTuiController(
            AgentSessionConfigV1(workspace=tmp_path, secret_file=secret)
        ),
        plain=True,
    )
    run_directory = tmp_path / ".chemsmart-agent" / "runs" / "live-x"
    run_directory.mkdir(parents=True)
    stream = run_directory / "events.jsonl"

    async def observe() -> None:
        async with app.run_test(size=(120, 35)) as pilot:
            stream.write_text(
                _line(
                    "tool_started",
                    {"tool": "preview_command", "request_id": "c1"},
                ),
                encoding="utf-8",
            )
            app._start_planning_tail(run_directory)
            await asyncio.sleep(0.5)
            await pilot.pause()
            rows = list(app.query(ToolRow))
            assert len(rows) == 1
            assert rows[0].state == "running"

            with stream.open("a", encoding="utf-8") as handle:
                handle.write(
                    _line(
                        "tool_succeeded",
                        {
                            "tool": "preview_command",
                            "request_id": "c1",
                            "canonical_result": {"result": {"node_id": "sp1"}},
                        },
                    )
                )
            await asyncio.sleep(0.6)
            await pilot.pause()
            rows = list(app.query(ToolRow))
            assert len(rows) == 1, "the settle mutates, never appends"
            assert rows[0].state == "finished"
            app._tail_stop.set()

    asyncio.run(observe())
