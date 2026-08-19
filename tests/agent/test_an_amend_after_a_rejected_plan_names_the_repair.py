"""A live session amended straight after a rejected plan call.

Its plan_scientific_workflow had been refused, which records nothing, so
the follow-up amend_scientific_workflow could only fail -- and it failed
with "unknown scientific workflow ID", a sentence that is true and names
neither the state (which IDs exist) nor the action (plan again, don't
amend). The refusal now carries both.
"""

from __future__ import annotations

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1


def _host(tmp_path):
    return CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="s1"
        ),
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
    )


def test_with_no_recorded_workflow_the_refusal_says_plan_again(tmp_path):
    host = _host(tmp_path)

    with pytest.raises(ContractError) as refusal:
        host.dispatch(
            turn_id="t1",
            tool_name="amend_scientific_workflow",
            arguments={
                "workflow_id": "dcca-dmso-abs",
                "support_repairs": [
                    {"node_id": "x", "blocked_reason": ""}
                ],
            },
        )

    message = str(refusal.value)
    assert "'dcca-dmso-abs'" in message
    assert "recorded workflow IDs: none" in message
    assert "records nothing" in message
    assert "fresh plan_scientific_workflow" in message


def test_a_wrong_id_is_shown_beside_the_recorded_ones(tmp_path):
    host = _host(tmp_path)

    class _FakePlan:
        workflow_id = "acetone-vertical-s0s1"

    host.scientific_toolchain_plans["f" * 64] = _FakePlan()

    with pytest.raises(ContractError) as refusal:
        host.dispatch(
            turn_id="t1",
            tool_name="amend_scientific_workflow",
            arguments={
                "workflow_id": "acetone-vertical",
                "support_repairs": [
                    {"node_id": "x", "blocked_reason": ""}
                ],
            },
        )

    message = str(refusal.value)
    assert "acetone-vertical-s0s1" in message
    assert "'acetone-vertical'" in message
