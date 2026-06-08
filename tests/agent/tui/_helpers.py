from __future__ import annotations

import json
import re
from pathlib import Path

from chemsmart.agent.core import CriticVerdict

SNAPSHOT_ROOT = Path(__file__).parent / "snapshots"


def normalize_svg(svg: str) -> str:
    svg = re.sub(r"terminal-\d+", "terminal-ID", svg)
    svg = re.sub(r'viewBox="[^"]+"', 'viewBox="NORMALIZED"', svg)
    svg = re.sub(r'width="\d+(?:\.\d+)?"', 'width="NORMALIZED"', svg)
    svg = re.sub(r'height="\d+(?:\.\d+)?"', 'height="NORMALIZED"', svg)
    svg = re.sub(r"fill: #[0-9a-fA-F]{6}", "fill: #COLOR", svg)
    svg = re.sub(r'fill="#[0-9a-fA-F]{6}"', 'fill="#COLOR"', svg)
    return svg.strip()


def assert_matches_snapshot(name: str, actual: str) -> None:
    import os
    path = SNAPSHOT_ROOT / f"{name}.svg"
    if os.environ.get("UPDATE_SNAPSHOTS") == "1":
        path.write_text(actual, encoding="utf-8")
        return
    expected = path.read_text(encoding="utf-8")
    assert normalize_svg(actual) == normalize_svg(expected)


def write_session_fixture(
    session_root: Path,
    session_id: str = "session-001",
    *,
    cwd: str | None = None,
) -> Path:
    session_dir = session_root / session_id
    session_dir.mkdir(parents=True, exist_ok=True)
    recorded_cwd = cwd or str(Path.cwd())
    entries = [
        {
            "kind": "request",
            "payload": {"request": "build water + recommend method"},
            "rationale": "build water + recommend method",
            "ts": "2026-05-08T00:00:00Z",
        },
        {
            "kind": "plan",
            "payload": {
                "steps": [
                    {
                        "tool": "build_molecule",
                        "args": {"filepath": "water.xyz"},
                        "rationale": "Load the starting structure.",
                    },
                    {
                        "tool": "dry_run_input",
                        "args": {"job": "$step1"},
                        "rationale": "Write the input for inspection.",
                    },
                ],
                "rationale": "Prepare a small water optimization.",
                "estimated_cost": "low",
            },
            "rationale": "Prepare a small water optimization.",
            "ts": "2026-05-08T00:00:01Z",
        },
        {
            "kind": "tool_result",
            "payload": {
                "step_index": 1,
                "tool": "dry_run_input",
                "artifact": "step_02.json",
                "payload": {
                    "inputfile": "/tmp/water.com",
                    "content": "#p b3lyp/6-31g* opt\n\nwater\n\n0 1\nO 0 0 0\nH 0 0 1\nH 0 1 0\n",
                },
            },
            "rationale": "Write the input for inspection.",
            "ts": "2026-05-08T00:00:02Z",
        },
        {
            "kind": "critic_verdict",
            "payload": CriticVerdict(
                verdict="ok",
                confidence=0.92,
                issues=[],
                rationale="The generated input looks reasonable.",
            ).model_dump(),
            "rationale": "The generated input looks reasonable.",
            "ts": "2026-05-08T00:00:03Z",
        },
        {
            "kind": "session_summary",
            "payload": {
                "total_steps_executed": 2,
                "total_steps_planned": 2,
                "blocked": False,
                "block_reason": None,
            },
            "rationale": "",
            "ts": "2026-05-08T00:00:04Z",
        },
    ]
    (session_dir / "decision_log.jsonl").write_text(
        "\n".join(json.dumps(entry) for entry in entries) + "\n",
        encoding="utf-8",
    )
    (session_dir / "session.json").write_text(
        json.dumps(
            {
                "session_id": session_id,
                "cwd": recorded_cwd,
                "started_at": "2026-05-08T00:00:00Z",
                "request_intent": "opt",
                "total_steps_planned": 2,
                "current_step_index": 2,
                "request": "build water + recommend method",
                "env_snapshot": {},
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    (session_dir / "state.json").write_text(
        (session_dir / "session.json").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    (session_dir / "session_metadata.json").write_text(
        json.dumps(
            {
                "session_id": session_id,
                "request": "build water + recommend method",
                "intent": "opt",
                "request_intent": "opt",
                "plan_steps": 2,
                "executed_steps": 2,
                "critic_verdict": "ok",
                "blocked": False,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return session_dir
