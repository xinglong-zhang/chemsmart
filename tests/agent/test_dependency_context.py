from __future__ import annotations

import pytest

from chemsmart.agent._contracts import (
    ContractError,
    canonical_json,
    canonical_sha256,
)
from chemsmart.agent.dependency_context import (
    bind_selected_public_records,
    build_dependency_context_public_projection,
    build_predecessor_evidence_ref,
    build_task_dependency_context_policy,
    select_task_dependency_context,
)
from chemsmart.agent.workflows import (
    ScientificWorkflowEdgeV2,
    ScientificWorkflowNodeV2,
    build_scientific_workflow_plan,
)


def _node(node_id: str, stage: str) -> ScientificWorkflowNodeV2:
    return ScientificWorkflowNodeV2(
        node_id=node_id,
        stage=stage,
        requested_program="pyscf",
        program="pyscf",
        engine="cpu",
        project_role="paper-project",
        unresolved_fields=(),
    )


def _fan_out_fan_in_plan():
    return build_scientific_workflow_plan(
        workflow_id="paper-workflow",
        task_spec_sha256="a" * 64,
        scientific_identity_sha256="b" * 64,
        nodes=(
            _node("opt", "opt"),
            _node("sp-low", "sp"),
            _node("sp-high", "sp"),
            _node("irrelevant-sibling", "sp"),
            _node("derived-value", "analysis"),
            _node("report", "report"),
        ),
        edges=(
            ScientificWorkflowEdgeV2(
                edge_id="a-opt-low",
                source_node_id="opt",
                target_node_id="sp-low",
                edge_kind="data",
                artifact_class="geometry_xyz",
                producer_output_id="optimized-geometry",
                consumer_input_id="geometry",
            ),
            ScientificWorkflowEdgeV2(
                edge_id="b-opt-high",
                source_node_id="opt",
                target_node_id="sp-high",
                edge_kind="data",
                artifact_class="geometry_xyz",
                producer_output_id="optimized-geometry",
                consumer_input_id="geometry",
            ),
            ScientificWorkflowEdgeV2(
                edge_id="c-opt-sibling",
                source_node_id="opt",
                target_node_id="irrelevant-sibling",
                edge_kind="control",
            ),
            ScientificWorkflowEdgeV2(
                edge_id="d-low-derived",
                source_node_id="sp-low",
                target_node_id="derived-value",
                edge_kind="data",
                artifact_class="structured_result",
                producer_output_id="electronic-energy",
                consumer_input_id="low-energy",
            ),
            ScientificWorkflowEdgeV2(
                edge_id="e-high-derived",
                source_node_id="sp-high",
                target_node_id="derived-value",
                edge_kind="data",
                artifact_class="structured_result",
                producer_output_id="electronic-energy",
                consumer_input_id="high-energy",
            ),
            ScientificWorkflowEdgeV2(
                edge_id="f-derived-report",
                source_node_id="derived-value",
                target_node_id="report",
                edge_kind="control",
            ),
        ),
    )


def _evidence():
    rows = (
        ("workflow-method", "decision", "", "1" * 64),
        ("opt-geometry", "artifact", "opt", "2" * 64),
        ("low-result", "validation", "sp-low", "3" * 64),
        ("high-result", "validation", "sp-high", "4" * 64),
        (
            "sibling-result",
            "validation",
            "irrelevant-sibling",
            "5" * 64,
        ),
        ("derived-rule", "decision", "derived-value", "6" * 64),
        ("report-result", "artifact", "report", "7" * 64),
        ("private-transcript", "transcript", "sp-low", "8" * 64),
    )
    return tuple(
        build_predecessor_evidence_ref(
            record_id=record_id,
            record_class=record_class,
            node_id=node_id,
            record_sha256=record_sha256,
            size_bytes=100,
        )
        for record_id, record_class, node_id, record_sha256 in rows
    )


def _policy(mode: str, *, max_public_bytes: int = 100_000):
    return build_task_dependency_context_policy(
        policy_id=f"workflow-{mode}",
        mode=mode,
        record_classes=("artifact", "decision", "validation"),
        max_public_bytes=max_public_bytes,
    )


def _select(mode: str, *, evidence=None, max_public_bytes=100_000):
    return select_task_dependency_context(
        selection_id=f"select-{mode}",
        plan=_fan_out_fan_in_plan(),
        target_node_id="derived-value",
        policy=_policy(mode, max_public_bytes=max_public_bytes),
        evidence_refs=_evidence() if evidence is None else evidence,
    )


def test_dependency_projection_selects_all_fan_in_and_ancestors_only():
    context, receipt = _select("dependency_projected")

    assert context is not None
    assert context.direct_predecessor_node_ids == ("sp-low", "sp-high")
    assert context.transitive_ancestor_node_ids == ("opt",)
    assert context.included_node_ids == (
        "opt",
        "sp-low",
        "sp-high",
        "derived-value",
    )
    assert context.excluded_node_ids == ("irrelevant-sibling", "report")
    assert context.non_dependency_node_ids == ()
    assert tuple(edge.edge_id for edge in context.selected_edges) == (
        "a-opt-low",
        "b-opt-high",
        "d-low-derived",
        "e-high-derived",
    )
    assert {item.record_id for item in context.evidence_refs} == {
        "workflow-method",
        "opt-geometry",
        "low-result",
        "high-result",
        "derived-rule",
    }
    assert receipt.status == "selected"
    assert receipt.context_sha256 == context.context_sha256
    reasons = {
        item.evidence_ref_sha256: item.reason for item in receipt.exclusions
    }
    sibling_ref = next(
        item for item in _evidence() if item.record_id == "sibling-result"
    )
    transcript_ref = next(
        item for item in _evidence() if item.record_id == "private-transcript"
    )
    assert reasons[sibling_ref.evidence_ref_sha256] == "node_not_selected"
    assert (
        reasons[transcript_ref.evidence_ref_sha256]
        == "record_class_not_selected"
    )


def test_ordered_predecessors_include_prior_sibling_not_descendant():
    context, _ = _select("ordered_predecessors")

    assert context is not None
    assert context.included_node_ids == (
        "opt",
        "sp-low",
        "sp-high",
        "irrelevant-sibling",
        "derived-value",
    )
    assert context.excluded_node_ids == ("report",)
    assert context.non_dependency_node_ids == ("irrelevant-sibling",)
    assert "c-opt-sibling" in {
        edge.edge_id for edge in context.selected_edges
    }
    assert "sibling-result" in {
        item.record_id for item in context.evidence_refs
    }
    assert "report-result" not in {
        item.record_id for item in context.evidence_refs
    }


def test_target_only_context_exposes_no_predecessor():
    context, _ = _select("target_only")

    assert context is not None
    assert context.included_node_ids == ("derived-value",)
    assert context.direct_predecessor_node_ids == ()
    assert context.transitive_ancestor_node_ids == ()
    assert context.selected_edges == ()
    assert {item.record_id for item in context.evidence_refs} == {
        "workflow-method",
        "derived-rule",
    }


def test_selection_is_canonical_across_evidence_input_order():
    forward, forward_receipt = _select("dependency_projected")
    reverse, reverse_receipt = _select(
        "dependency_projected", evidence=tuple(reversed(_evidence()))
    )

    assert reverse == forward
    assert reverse_receipt == forward_receipt


def test_context_budget_blocks_instead_of_truncating_required_evidence():
    context, receipt = _select("dependency_projected", max_public_bytes=1)

    assert context is None
    assert receipt.status == "blocked_context_budget"
    assert receipt.context_sha256 == ""
    assert receipt.included_evidence_ref_sha256s == ()
    assert "policy permits 1" in receipt.reason
    assert {
        item.reason for item in receipt.exclusions
    }.issuperset({"context_budget_exceeded"})


def test_unknown_node_evidence_fails_closed():
    unknown = build_predecessor_evidence_ref(
        record_id="unknown-record",
        record_class="artifact",
        node_id="not-in-plan",
        record_sha256="9" * 64,
        size_bytes=10,
    )
    with pytest.raises(ContractError, match="unknown node"):
        _select("dependency_projected", evidence=(*_evidence(), unknown))




def _payload_context():
    payloads = {
        "workflow-method": {
            "schema_version": "chemsmart.public-decision.v1",
            "method": "b3lyp",
        },
        "opt-geometry": {
            "schema_version": "chemsmart.public-artifact.v1",
            "artifact_id": "optimized-geometry",
        },
        "low-result": {
            "schema_version": "chemsmart.public-validation.v1",
            "energy_hartree": -40.1,
        },
        "high-result": {
            "schema_version": "chemsmart.public-validation.v1",
            "energy_hartree": -40.2,
        },
        "derived-rule": {
            "schema_version": "chemsmart.public-decision.v1",
            "expression": "high-low",
        },
    }
    node_ids = {
        "workflow-method": "",
        "opt-geometry": "opt",
        "low-result": "sp-low",
        "high-result": "sp-high",
        "derived-rule": "derived-value",
    }
    record_classes = {
        "workflow-method": "decision",
        "opt-geometry": "artifact",
        "low-result": "validation",
        "high-result": "validation",
        "derived-rule": "decision",
    }
    refs = tuple(
        build_predecessor_evidence_ref(
            record_id=record_id,
            record_class=record_classes[record_id],
            node_id=node_ids[record_id],
            record_sha256=canonical_sha256(payload),
            size_bytes=len(canonical_json(payload).encode("utf-8")),
        )
        for record_id, payload in payloads.items()
    )
    context, _ = select_task_dependency_context(
        selection_id="payload-selection",
        plan=_fan_out_fan_in_plan(),
        target_node_id="derived-value",
        policy=_policy("dependency_projected"),
        evidence_refs=refs,
    )
    assert context is not None
    return context, payloads


def test_selected_public_record_binding_is_exact_hash_verified_and_canonical():
    context, payloads = _payload_context()

    bound = bind_selected_public_records(
        context=context,
        records=dict(reversed(tuple(payloads.items()))),
    )

    assert tuple(item["record_id"] for item in bound) == tuple(
        item.record_id for item in context.evidence_refs
    )
    assert tuple(item["record_sha256"] for item in bound) == tuple(
        item.record_sha256 for item in context.evidence_refs
    )
    assert all(
        canonical_sha256(item["public_record"]) == item["record_sha256"]
        for item in bound
    )


def test_public_projection_exposes_only_selected_hash_bound_payloads():
    context, payloads = _payload_context()
    _, receipt = select_task_dependency_context(
        selection_id="payload-selection",
        plan=_fan_out_fan_in_plan(),
        target_node_id="derived-value",
        policy=_policy("dependency_projected"),
        evidence_refs=tuple(
            build_predecessor_evidence_ref(
                record_id=item.record_id,
                record_class=item.record_class,
                node_id=item.node_id,
                record_sha256=item.record_sha256,
                size_bytes=item.size_bytes,
            )
            for item in context.evidence_refs
        ),
    )

    projection = build_dependency_context_public_projection(
        context=context,
        selection_receipt=receipt,
        records=payloads,
    )

    assert set(projection) == {
        "task_dependency_context",
        "context_selection_receipt",
        "selected_predecessor_records",
    }
    assert {
        item["record_id"]
        for item in projection["selected_predecessor_records"]
    } == set(payloads)
    assert "unselected-sibling" not in canonical_json(projection)


def test_public_projection_rejects_a_receipt_from_another_context_mode():
    context, payloads = _payload_context()
    _, wrong_receipt = _select("ordered_predecessors")

    with pytest.raises(ContractError, match="selection receipt .* mismatch"):
        build_dependency_context_public_projection(
            context=context,
            selection_receipt=wrong_receipt,
            records=payloads,
        )


@pytest.mark.parametrize("mutation", ["missing", "unexpected", "changed"])
def test_selected_public_record_binding_rejects_inexact_payload_sets(mutation):
    context, payloads = _payload_context()
    supplied = dict(payloads)
    if mutation == "missing":
        supplied.pop("high-result")
        match = "exactly match selected records"
    elif mutation == "unexpected":
        supplied["unselected-sibling"] = {"value": "must-not-leak"}
        match = "exactly match selected records"
    else:
        supplied["high-result"] = {
            **supplied["high-result"],
            "energy_hartree": -39.0,
        }
        match = "digest mismatch"

    with pytest.raises(ContractError, match=match):
        bind_selected_public_records(context=context, records=supplied)
