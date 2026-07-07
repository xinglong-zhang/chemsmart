from __future__ import annotations

import json
from pathlib import Path

from chemsmart.agent.core import AgentSession
from chemsmart.agent.registry import ToolRegistry

from ._agent_session_helpers import FakeProvider


def _read_decision_log(session_dir: Path) -> list[dict]:
    return [
        json.loads(line)
        for line in (session_dir / "decision_log.jsonl")
        .read_text()
        .splitlines()
        if line.strip()
    ]


def test_chitchat_request_is_classified_and_skips_critic(tmp_path: Path):
    provider = FakeProvider(
        [
            {
                "steps": [],
                "rationale": "Hello! I can help plan chemsmart workflows.",
                "estimated_cost": "none",
            }
        ]
    )
    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path / "sessions",
    )

    result = session.run("hello", dry_submit=True)

    assert result["plan"].intent == "chitchat"
    assert result["advisory_only"] is True
    assert result["is_chitchat"] is True
    assert len(provider.calls) == 1

    entries = _read_decision_log(Path(result["session_dir"]))
    assert [entry["kind"] for entry in entries] == [
        "request",
        "plan",
        "session_summary",
    ]
    assert entries[1]["payload"]["intent"] == "chitchat"
    assert entries[2]["payload"]["request_intent"] == "chitchat"
    assert entries[2]["payload"]["is_chitchat"] is True


def test_zero_step_chemistry_advice_remains_advisory(tmp_path: Path):
    provider = FakeProvider(
        [
            {
                "steps": [],
                "rationale": "Start with B3LYP/6-31G(d), then refine if needed.",
                "estimated_cost": "low",
            }
        ]
    )
    session = AgentSession(
        provider=provider,
        registry=ToolRegistry.default(),
        session_root=tmp_path / "sessions",
    )

    result = session.run(
        "Which method should I use for a neutral organic single-point?",
        dry_submit=True,
    )

    assert result["plan"].intent == "advisory"
    assert result["is_chitchat"] is False
