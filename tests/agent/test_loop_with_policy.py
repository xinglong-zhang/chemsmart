from __future__ import annotations

from dataclasses import dataclass, field

from chemsmart.agent.core import DecisionLog
from chemsmart.agent.handles import HandleStore
from chemsmart.agent.loop import ToolLoop
from chemsmart.agent.permissions import (
    ApprovalDecision,
    PermissionMode,
    PermissionPolicy,
)

from ._agent_session_helpers import FakeProvider
from ._loop_helpers import (
    ScriptedRegistry,
    openai_final_response,
    openai_tool_call_response,
    tool_call,
)


@dataclass
class FakeApprover:
    decisions: list[ApprovalDecision]
    calls: list[str] = field(default_factory=list)

    def __call__(self, request):
        self.calls.append(request.name)
        return self.decisions.pop(0)


def test_driving_mode_denies_run_local_without_yolo(tmp_path):
    provider = FakeProvider(
        [
            {
                "__raw_response__": openai_tool_call_response(
                    tool_call("call_1", "run_local", {"job": "job_1"})
                )
            },
            {"__raw_response__": openai_final_response("Denied.")},
        ]
    )
    registry = ScriptedRegistry({"run_local": {"ok": True}})
    loop = ToolLoop(
        provider=provider,
        registry=registry,
        handle_store=HandleStore(tmp_path),
        decision_log=DecisionLog(tmp_path / "decision_log.jsonl"),
        policy=PermissionPolicy(mode=PermissionMode.DRIVING),
    )

    result = loop.run_turn(
        messages=[{"role": "user", "content": "Run it."}],
        tool_defs=registry.openai_tool_defs(),
    )

    assert [outcome.status for outcome in result["tool_outcomes"]] == [
        "denied"
    ]
    assert registry.calls == []
    assert result["denials_count"] == 1
    assert result["approvals_count"] == 0


def test_driving_mode_allows_run_local_with_yolo(tmp_path):
    provider = FakeProvider(
        [
            {
                "__raw_response__": openai_tool_call_response(
                    tool_call("call_1", "run_local", {"job": "job_1"})
                )
            },
            {"__raw_response__": openai_final_response("Ran.")},
        ]
    )
    registry = ScriptedRegistry({"run_local": {"ok": True, "status": "done"}})
    loop = ToolLoop(
        provider=provider,
        registry=registry,
        handle_store=HandleStore(tmp_path),
        decision_log=DecisionLog(tmp_path / "decision_log.jsonl"),
        policy=PermissionPolicy(mode=PermissionMode.DRIVING, yolo=True),
    )

    result = loop.run_turn(
        messages=[{"role": "user", "content": "Run it."}],
        tool_defs=registry.openai_tool_defs(),
    )

    assert [outcome.status for outcome in result["tool_outcomes"]] == ["ok"]
    assert registry.calls == [("run_local", {"job": "job_1"})]
    assert result["denials_count"] == 0
    assert result["approvals_count"] == 1


def test_permission_mode_needs_approver_and_allow_session_persists(tmp_path):
    provider = FakeProvider(
        [
            {
                "__raw_response__": openai_tool_call_response(
                    tool_call("call_1", "recommend_method", {"task": "opt"})
                )
            },
            {
                "__raw_response__": openai_tool_call_response(
                    tool_call("call_2", "recommend_method", {"task": "freq"})
                )
            },
            {"__raw_response__": openai_final_response("Done.")},
        ]
    )
    registry = ScriptedRegistry({"recommend_method": {"method": "wb97xd"}})
    approver = FakeApprover([ApprovalDecision.ALLOW_SESSION])
    policy = PermissionPolicy(mode=PermissionMode.PERMISSION)
    loop = ToolLoop(
        provider=provider,
        registry=registry,
        handle_store=HandleStore(tmp_path),
        decision_log=DecisionLog(tmp_path / "decision_log.jsonl"),
        policy=policy,
        approver=approver,
    )

    result = loop.run_turn(
        messages=[{"role": "user", "content": "Recommend twice."}],
        tool_defs=registry.openai_tool_defs(),
    )

    assert [outcome.status for outcome in result["tool_outcomes"]] == [
        "ok",
        "ok",
    ]
    assert approver.calls == ["recommend_method"]
    assert policy.session_allow == {"recommend_method"}
    assert registry.calls == [
        ("recommend_method", {"task": "opt"}),
        ("recommend_method", {"task": "freq"}),
    ]
    assert result["approvals_count"] == 2


def test_permission_mode_without_approver_auto_denies(tmp_path):
    provider = FakeProvider(
        [
            {
                "__raw_response__": openai_tool_call_response(
                    tool_call("call_1", "recommend_method", {"task": "opt"})
                )
            },
            {"__raw_response__": openai_final_response("Denied.")},
        ]
    )
    registry = ScriptedRegistry({"recommend_method": {"method": "wb97xd"}})
    loop = ToolLoop(
        provider=provider,
        registry=registry,
        handle_store=HandleStore(tmp_path),
        decision_log=DecisionLog(tmp_path / "decision_log.jsonl"),
        policy=PermissionPolicy(mode=PermissionMode.PERMISSION),
    )

    result = loop.run_turn(
        messages=[{"role": "user", "content": "Recommend."}],
        tool_defs=registry.openai_tool_defs(),
    )

    assert [outcome.status for outcome in result["tool_outcomes"]] == [
        "denied"
    ]
    assert result["tool_outcomes"][0].error_message == "no_approver"
    assert registry.calls == []
