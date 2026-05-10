from __future__ import annotations

import pytest

from chemsmart.agent.permissions import (
    ApprovalDecision,
    PermissionMode,
    PermissionPolicy,
    ResolvedDecision,
)
from chemsmart.agent.provider_adapter import ToolRequest


def make_request(name: str) -> ToolRequest:
    return ToolRequest(
        request_id=f"openai:{name}",
        provider="openai",
        provider_call_id=f"call_{name}",
        name=name,
        arguments_json="{}",
        arguments={},
        raw={},
    )


@pytest.mark.parametrize(
    ("policy", "tool", "decision", "reason"),
    [
        (
            PermissionPolicy(mode=PermissionMode.DRIVING),
            "recommend_method",
            ResolvedDecision.AUTO_ALLOW,
            "driving_mode",
        ),
        (
            PermissionPolicy(mode=PermissionMode.DRIVING),
            "run_local",
            ResolvedDecision.AUTO_DENY,
            "missing_yolo",
        ),
        (
            PermissionPolicy(mode=PermissionMode.DRIVING, yolo=True),
            "run_local",
            ResolvedDecision.AUTO_ALLOW,
            "yolo",
        ),
        (
            PermissionPolicy(mode=PermissionMode.PERMISSION),
            "recommend_method",
            ResolvedDecision.NEEDS_USER,
            "needs_user",
        ),
        (
            PermissionPolicy(
                mode=PermissionMode.PERMISSION,
                session_allow={"recommend_method"},
            ),
            "recommend_method",
            ResolvedDecision.AUTO_ALLOW,
            "session_rule",
        ),
        (
            PermissionPolicy(
                mode=PermissionMode.PERMISSION,
                session_allow={"recommend_method"},
            ),
            "run_local",
            ResolvedDecision.NEEDS_USER,
            "needs_user",
        ),
        (
            PermissionPolicy(
                mode=PermissionMode.DRIVING,
                driving_denylist={"validate_runtime"},
            ),
            "validate_runtime",
            ResolvedDecision.AUTO_DENY,
            "missing_yolo",
        ),
        (
            PermissionPolicy(
                mode=PermissionMode.DRIVING,
                yolo=True,
                driving_denylist={"validate_runtime"},
            ),
            "validate_runtime",
            ResolvedDecision.AUTO_ALLOW,
            "yolo",
        ),
    ],
)
def test_permission_policy_resolve_matrix(policy, tool, decision, reason):
    resolved = policy.resolve(make_request(tool))

    assert resolved.decision == decision
    assert resolved.reason == reason


def test_permission_policy_record_persists_session_allow():
    policy = PermissionPolicy(mode=PermissionMode.PERMISSION)

    policy.record("recommend_method", ApprovalDecision.ALLOW_ONCE)
    assert "recommend_method" not in policy.session_allow

    policy.record("recommend_method", ApprovalDecision.ALLOW_SESSION)

    assert "recommend_method" in policy.session_allow
