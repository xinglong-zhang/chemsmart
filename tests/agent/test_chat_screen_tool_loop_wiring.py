from __future__ import annotations

from types import SimpleNamespace

from chemsmart.agent.core import AgentSession, Plan
from chemsmart.agent.permissions import (
    ApprovalDecision,
    PermissionMode,
)
from chemsmart.agent.tui.screens.chat import ChatScreen


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
        "approval_mode": "permission",
        "driving_mode": False,
        "yolo": False,
        "denials_count": 0,
        "approvals_count": 0,
    }


def test_run_agent_session_calls_run_loop_with_screen_policy(
    monkeypatch, tmp_path
):
    captured: dict[str, object] = {}

    def fake_run_loop(self, request, **kwargs):
        captured["request"] = request
        captured["policy"] = kwargs["policy"]
        captured["approver"] = kwargs["approver"]
        kwargs["policy"].record(
            "build_molecule", ApprovalDecision.ALLOW_SESSION
        )
        return _loop_result()

    monkeypatch.setattr(AgentSession, "run_loop", fake_run_loop)

    screen = ChatScreen(session_root=tmp_path)
    screen._permission_mode = PermissionMode.PERMISSION
    screen._yolo_enabled = True
    screen._session_allow_tools = {"dry_run_input"}

    result = ChatScreen.run_agent_session.__wrapped__(screen, "load water")

    policy = captured["policy"]
    assert captured["request"] == "load water"
    assert policy.mode == PermissionMode.PERMISSION
    assert policy.yolo is True
    assert policy.session_allow == {
        "dry_run_input",
        "build_molecule",
    }
    assert callable(captured["approver"])
    assert result["session_allow_tools"] == [
        "build_molecule",
        "dry_run_input",
    ]


def test_resume_agent_session_calls_run_loop_with_loaded_request(
    monkeypatch,
    tmp_path,
):
    captured: dict[str, object] = {}

    class DummySession:
        def __init__(self) -> None:
            self.state = SimpleNamespace(request="resume this")

        def run_loop(self, request, **kwargs):
            captured["request"] = request
            captured["policy"] = kwargs["policy"]
            captured["approver"] = kwargs["approver"]
            return _loop_result()

    monkeypatch.setattr(
        AgentSession,
        "load",
        classmethod(lambda cls, session_id, **kwargs: DummySession()),
    )

    screen = ChatScreen(session_root=tmp_path)
    screen._permission_mode = PermissionMode.DRIVING
    screen._yolo_enabled = True

    ChatScreen.resume_agent_session.__wrapped__(
        screen,
        "session-123",
        cwd_override="/tmp/project",
    )

    policy = captured["policy"]
    assert captured["request"] == "resume this"
    assert policy.mode == PermissionMode.DRIVING
    assert policy.yolo is True
    assert callable(captured["approver"])
