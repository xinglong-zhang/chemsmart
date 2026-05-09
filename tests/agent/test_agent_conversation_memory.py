from __future__ import annotations

import json
from pathlib import Path

from chemsmart.agent.core import AgentSession
from chemsmart.agent.registry import ToolRegistry

from ._agent_session_helpers import FakeProvider, critic_ok, planner_plan


def _single_point_plan(filepath: str, label: str) -> dict:
    return {
        "steps": [
            {
                "tool": "build_molecule",
                "args": {"filepath": filepath},
                "rationale": "Reload the same molecule for the follow-up.",
            },
            {
                "tool": "build_gaussian_settings",
                "args": {"functional": "b3lyp", "basis": "6-31g*"},
                "rationale": "Reuse the same lightweight DFT level.",
            },
            {
                "tool": "build_job",
                "args": {
                    "kind": "gaussian.sp",
                    "molecule": "$step1",
                    "settings": "$step2",
                    "label": label,
                },
                "rationale": "Prepare the single-point follow-up job.",
            },
            {
                "tool": "dry_run_input",
                "args": {"job": "$step3"},
                "rationale": "Preview the single-point input file.",
            },
            {
                "tool": "validate_runtime",
                "args": {"job": "$step3", "server": None},
                "rationale": "Check local prerequisites before execution.",
            },
        ],
        "rationale": "Follow the earlier setup with a single-point job.",
        "estimated_cost": "low",
    }


def test_resumed_session_rebuilds_memory_from_decision_log(
    monkeypatch,
    single_molecule_xyz_file,
    tmp_path: Path,
):
    import chemsmart.agent.tools as agent_tools

    first_request = f"build water from {single_molecule_xyz_file}"
    second_request = "now run a single-point on it"
    session_root = tmp_path / "sessions"

    def fake_validate_runtime(job, server=None):
        return {
            "ok": "ok",
            "local_ok": True,
            "local_issues": [],
            "remote_unknown": [],
        }

    monkeypatch.setattr(agent_tools, "validate_runtime", fake_validate_runtime)

    first_session = AgentSession(
        provider=FakeProvider(
            [planner_plan(single_molecule_xyz_file, "memory_opt"), critic_ok()]
        ),
        registry=ToolRegistry.default(),
        session_root=session_root,
    )
    first_result = first_session.run(first_request, dry_submit=True)

    resumed_provider = FakeProvider(
        [
            _single_point_plan(single_molecule_xyz_file, "memory_sp"),
            critic_ok(),
        ]
    )
    resumed_session = AgentSession.load(
        first_result["session_id"],
        provider=resumed_provider,
        registry=ToolRegistry.default(),
        session_root=session_root,
    )

    assert [
        turn.request for turn in resumed_session.conversation_history.turns
    ] == [first_request]

    second_result = resumed_session.run(second_request, dry_submit=True)
    planner_payload = json.loads(
        resumed_provider.calls[0]["messages"][1]["content"]
    )
    history = planner_payload["conversation_history"]["recent_turns"]

    assert second_result["session_id"] == first_result["session_id"]
    assert history[0]["request"] == first_request
    assert any(
        "dry_run_input wrote" in item
        for item in history[0]["reusable_results"]
    )
