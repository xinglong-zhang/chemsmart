from __future__ import annotations

import pytest

from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import CommandCompiledToolHostV1

TASK_SPEC_SHA256 = "a" * 64


def _host(tmp_path, store):
    return CommandCompiledToolHostV1(
        event_store=store,
        task_spec_sha256s=(TASK_SPEC_SHA256,),
        approved_workspace=tmp_path / "workspace",
    )


def _plan(host, *, expected_count):
    return host.dispatch(
        turn_id="turn-plan",
        tool_name="plan_scientific_workflow",
        arguments={
            "plan_id": "typed-validation-plan",
            "workflow_id": "typed-validation-workflow",
            "task_spec_id": TASK_SPEC_SHA256,
            "calculation_nodes": [],
            "analysis_nodes": [
                {
                    "node_id": "derive-count",
                    "analysis_kind": "quantity_expression",
                    "dependencies": [],
                    "inputs": [],
                    "selectors": [],
                    "outputs": [
                        {
                            "output_id": "imaginary_mode_count",
                            "quantity_kind": "count",
                            "unit": "1",
                        }
                    ],
                    "expression_nodes": [
                        {
                            "node_id": "imaginary_mode_count",
                            "operation": "literal",
                            "literal_value": 0,
                            "literal_unit": "1",
                        }
                    ],
                    "expression_output_node_ids": ["imaginary_mode_count"],
                    "support_state": "planned",
                    "blocked_reason": "",
                    "validation_rules": [],
                },
                {
                    "node_id": "validate-count",
                    "analysis_kind": "scientific_validation",
                    "dependencies": ["derive-count"],
                    "inputs": [
                        {
                            "input_id": "observed_count",
                            "source_kind": "analysis_output",
                            "producer_node_id": "derive-count",
                            "producer_output_id": "imaginary_mode_count",
                        }
                    ],
                    "selectors": [],
                    "outputs": [
                        {
                            "output_id": "minimum_verdict",
                            "quantity_kind": "validation_verdict",
                            "unit": "1",
                        }
                    ],
                    "expression_nodes": [],
                    "expression_output_node_ids": [],
                    "support_state": "planned",
                    "blocked_reason": "",
                    "validation_rules": [
                        {
                            "rule_id": "expected_imaginary_count",
                            "predicate": "count_equals",
                            "input_ids": ["observed_count"],
                            "expected_count": expected_count,
                        }
                    ],
                },
                {
                    "node_id": "render-claims",
                    "analysis_kind": "claim_rendering",
                    "dependencies": ["derive-count", "validate-count"],
                    "inputs": [
                        {
                            "input_id": "count_claim",
                            "source_kind": "analysis_output",
                            "producer_node_id": "derive-count",
                            "producer_output_id": "imaginary_mode_count",
                        },
                        {
                            "input_id": "verdict_claim",
                            "source_kind": "analysis_output",
                            "producer_node_id": "validate-count",
                            "producer_output_id": "minimum_verdict",
                        },
                    ],
                    "selectors": [],
                    "outputs": [
                        {
                            "output_id": "imaginary_mode_count",
                            "quantity_kind": "count",
                            "unit": "1",
                        },
                        {
                            "output_id": "minimum_verdict",
                            "quantity_kind": "validation_verdict",
                            "unit": "1",
                        },
                    ],
                    "expression_nodes": [],
                    "expression_output_node_ids": [],
                    "support_state": "planned",
                    "blocked_reason": "",
                    "validation_rules": [],
                },
            ],
            "required_output_ids": [
                "imaginary_mode_count",
                "minimum_verdict",
            ],
        },
    )["result"]


def _decision(host, *, receipts):
    return host.dispatch(
        turn_id="turn-analysis",
        tool_name="record_scientific_decision",
        arguments={
            "decision_id": "typed-validation-decision",
            "task_spec_sha256": TASK_SPEC_SHA256,
            "assumptions": ["The typed count is the planned observable."],
            "method_rationale": "Evaluate the predeclared count predicate.",
            "alternatives": ["Report the failed predicate when it is false."],
            "uncertainties": [
                "This mechanical fixture contains no chemistry."
            ],
            "diagnostics": ["Inspect the host-owned validation verdict."],
            "stage_order": ["derive-count", "validate-count", "render-claims"],
            "evidence_refs": [],
            "postprocessing_receipt_sha256s": list(receipts),
        },
    )["result"]


@pytest.mark.parametrize(
    ("expected_count", "expected_verdict"),
    ((0, True), (1, False)),
)
def test_planned_validation_is_typed_replayable_and_required_for_completion(
    tmp_path, expected_count, expected_verdict
):
    store = RuntimeEventStore(
        tmp_path / "events.jsonl", session_id="validation-session"
    )
    host = _host(tmp_path, store)
    _plan(host, expected_count=expected_count)

    expression = host.dispatch(
        turn_id="turn-analysis",
        tool_name="evaluate_quantity_expression",
        arguments={
            "expression_id": "observed-imaginary-count",
            "inputs": [],
            "nodes": [
                {
                    "node_id": "imaginary_mode_count",
                    "operation": "literal",
                    "literal_value": 0,
                    "literal_unit": "1",
                }
            ],
            "output_node_ids": ["imaginary_mode_count"],
        },
    )["result"]
    early_claim = host.dispatch(
        turn_id="turn-analysis",
        tool_name="record_analysis_claims",
        arguments={
            "task_spec_sha256": TASK_SPEC_SHA256,
            "claims": [
                {
                    "claim_id": "imaginary_mode_count",
                    "receipt_sha256": expression["receipt_sha256"],
                    "quantity_id": "imaginary_mode_count",
                    "display_unit": "1",
                }
            ],
        },
    )["result"]
    _decision(
        host,
        receipts=(
            expression["receipt_sha256"],
            early_claim["receipt_sha256"],
        ),
    )

    assert host._completion_receipts_for_latest_analysis_toolchain() == ()
    before = host.dispatch(
        turn_id="turn-analysis",
        tool_name="inspect_workflow_frontier",
        arguments={"workflow_id": "typed-validation-workflow"},
    )["result"]["workflow_frontier"]
    validator = next(
        node for node in before["nodes"] if node["node_id"] == "validate-count"
    )
    assert validator["state"] == "actionable"
    assert validator["next_tool"] == "evaluate_scientific_validation"

    validation = host.dispatch(
        turn_id="turn-analysis",
        tool_name="evaluate_scientific_validation",
        arguments={
            "workflow_id": "typed-validation-workflow",
            "node_id": "validate-count",
            "inputs": [
                {
                    "input_id": "observed_count",
                    "receipt_sha256": expression["receipt_sha256"],
                    "quantity_id": "imaginary_mode_count",
                }
            ],
        },
    )["result"]
    assert validation["status"] == "evaluated"
    assert validation["all_rules_passed"] is expected_verdict
    assert validation["outputs"][0]["quantity_id"] == "minimum_verdict"
    assert validation["outputs"][0]["value"] == int(expected_verdict)

    claims = host.dispatch(
        turn_id="turn-analysis",
        tool_name="record_analysis_claims",
        arguments={
            "task_spec_sha256": TASK_SPEC_SHA256,
            "claims": [
                {
                    "claim_id": "imaginary_mode_count",
                    "receipt_sha256": expression["receipt_sha256"],
                    "quantity_id": "imaginary_mode_count",
                    "display_unit": "1",
                },
                {
                    "claim_id": "minimum_verdict",
                    "receipt_sha256": validation["receipt_sha256"],
                    "quantity_id": "minimum_verdict",
                    "display_unit": "1",
                },
            ],
        },
    )["result"]
    _decision(
        host,
        receipts=(
            expression["receipt_sha256"],
            validation["receipt_sha256"],
            claims["receipt_sha256"],
        ),
    )

    after = host.dispatch(
        turn_id="turn-analysis",
        tool_name="inspect_workflow_frontier",
        arguments={"workflow_id": "typed-validation-workflow"},
    )["result"]["workflow_frontier"]
    assert "validate-count" in after["completed_node_ids"]
    assert "render-claims" in after["completed_node_ids"]
    completion = host._completion_receipts_for_latest_analysis_toolchain()
    assert len(completion) == 1

    store.terminate(
        turn_id="turn-analysis",
        terminal_state="complete",
        reason="typed validation toolchain completed",
        required_receipt_sha256s=completion,
    )
    assert validation["receipt_sha256"] in (
        store.state().scientific_validation_receipts
    )
    reopened = _host(tmp_path, store)
    replayed = reopened.scientific_validation_receipts[
        validation["receipt_sha256"]
    ]
    assert replayed.all_rules_passed is expected_verdict
    assert replayed.outputs[0].value == int(expected_verdict)
