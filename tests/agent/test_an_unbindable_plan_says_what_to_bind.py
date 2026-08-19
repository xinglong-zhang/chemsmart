"""A refusal about the host's own bookkeeping teaches the model nothing.

A live session planned a faithful two-arm workflow -- the paper's Gaussian
platform declared blocked beside an executable ORCA substitution -- whose
inputs all named either a producer or nothing, because the geometry the
observable needs did not exist yet. The identity resolver wrote down exactly
why the scientific plan could not bind (`workflow.input.unresolved`,
`workflow.scientific_identity.unbound`) and returned None; the toolchain
binding then raised its own invariant, "scientific toolchain lacks its
task-bound scientific plan", discarding the diagnosis. The session retried
the identical call three times -- the message gave it nothing to change.

The refusal now happens where the findings still exist, and carries them.
"""

from __future__ import annotations

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1


def _plan_args(inputs):
    return dict(
        plan_id="p",
        workflow_id="w",
        task_spec_id="a" * 64,
        calculation_nodes=[
            {
                "node_id": "opt-x",
                "program": "orca",
                "jobtype": "opt",
                "project_role": "r",
                "dependencies": [],
                "inputs": inputs,
                "expected_outputs": [
                    {"output_id": "o", "artifact_class": "orca_output"}
                ],
                "unresolved_fields": [],
                "produces_observables": ["energy"],
                "support_state": "planned",
                "blocked_reason": "",
            }
        ],
        analysis_nodes=[
            {
                "node_id": "ex",
                "analysis_kind": "result_extraction",
                "dependencies": ["opt-x"],
                "inputs": [
                    {
                        "input_id": "r",
                        "source_kind": "program_output",
                        "producer_node_id": "opt-x",
                        "producer_output_id": "o",
                    }
                ],
                "selectors": [
                    {"quantity_id": "energy", "selector": "energy"}
                ],
                "outputs": [
                    {
                        "output_id": "energy",
                        "quantity_kind": "energy",
                        "unit": "hartree",
                    }
                ],
                "expression_nodes": [],
                "expression_output_node_ids": [],
                "support_state": "planned",
                "blocked_reason": "",
                "validation_rules": [],
            }
        ],
        required_output_ids=["energy"],
    )


def test_the_refusal_carries_the_findings_and_the_actions(tmp_path):
    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="s1"
        ),
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
    )
    named_but_unbound = [
        {
            "binding_id": "filename",
            "artifact_class": "geometry_xyz",
            "artifact_id": "geometry-feedbeef",
            "producer_node_id": "",
            "producer_output_id": "",
        }
    ]

    with pytest.raises(ContractError) as refusal:
        host.dispatch(
            turn_id="t1",
            tool_name="plan_scientific_workflow",
            arguments=_plan_args(named_but_unbound),
        )

    message = str(refusal.value)
    assert "opt-x" in message
    assert "workflow.input.unresolved" in message
    assert "workflow.scientific_identity.unbound" in message
    assert "bind_scientific_identity" in message
    assert "blocked_unsupported" in message
    # The internal invariant no longer leaks for this shape.
    assert "task-bound scientific plan" not in message


def test_a_workflow_that_anchors_nothing_is_told_so(tmp_path):
    """Plan 3 of the same protocol: every stage blocked, identities bound,
    and still refused twice -- the general sentence was satisfied vacuously.
    An all-blocked draft consumes no artifact anywhere, and the refusal must
    say that anchoring, not input hygiene, is what is missing."""

    host = CommandCompiledToolHostV1(
        event_store=RuntimeEventStore(
            tmp_path / "events.jsonl", session_id="s2"
        ),
        task_spec_sha256s=("a" * 64,),
        approved_workspace=tmp_path / "workspace",
    )
    args = _plan_args([])
    args["calculation_nodes"][0]["support_state"] = "blocked_unsupported"
    args["calculation_nodes"][0]["blocked_reason"] = (
        "no complex geometry exists in this workspace"
    )

    with pytest.raises(ContractError) as refusal:
        host.dispatch(
            turn_id="t1",
            tool_name="plan_scientific_workflow",
            arguments=args,
        )

    message = str(refusal.value)
    assert "anchors no molecule" in message
    assert "At least one plannable initial node" in message
    assert "bind_scientific_identity" in message
    # The input-hygiene sentence would mislead here; it must not be the text.
    assert "declared blocked_unsupported with its reason" not in message
