from __future__ import annotations

from click.testing import CliRunner

from chemsmart.agent.cli import agent
from chemsmart.agent.core import Plan
from chemsmart.agent.permissions import ApprovalDecision, PermissionMode


def _loop_result() -> dict[str, object]:
    return {
        "session_id": "session-123",
        "session_dir": "/tmp/session-123",
        "plan": Plan(steps=[], rationale="Answer directly."),
        "plan_text": "Plan:",
        "critic_verdict": None,
        "completed_steps": 0,
        "blocked": False,
        "dry_run_result": None,
        "dry_run_results": [],
        "runtime_result": None,
        "preview_submit": None,
        "results": [],
        "assistant_output": "Done.",
        "tool_requests": [],
        "tool_outcomes": [],
        "loop_state": {},
        "final_message": "Done.",
        "limit_reason": None,
        "advisory_only": True,
        "is_chitchat": False,
        "approval_mode": "driving",
        "driving_mode": True,
        "yolo": True,
        "denials_count": 0,
        "approvals_count": 0,
    }


def test_agent_run_mode_flags_pass_permission_policy(monkeypatch):
    captured: dict[str, object] = {}

    def fake_run_loop(self, request, **kwargs):
        captured["request"] = request
        captured["policy"] = kwargs["policy"]
        captured["approver"] = kwargs["approver"]
        return _loop_result()

    monkeypatch.setattr(
        "chemsmart.agent.cli.AgentSession.run_loop", fake_run_loop
    )

    runner = CliRunner()
    result = runner.invoke(
        agent,
        ["run", "--mode", "driving", "--yolo", "optimize water"],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    assert captured["request"] == "optimize water"
    assert captured["policy"].mode == PermissionMode.DRIVING
    assert captured["policy"].yolo is True
    assert captured["approver"](None) == ApprovalDecision.DENY
