"""A finished plan must be reviewable into an approval without re-planning.

``agent run`` consumes an approval binding exact nodes and an exact workflow
digest.  Nothing produced one.  The documented flow is "plan, then run with an
approval file", and the only route to that file was to re-plan and hope a
fresh session emitted a byte-identical workflow -- which is why no numerical
case in the observed sessions had ever been executed.

These tests pin the missing half and, more importantly, pin that it stays
inert: a request is not an approval, and cannot become one by accident.
"""

import json

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.execution import (
    ApprovedNodeBindingV1,
    WorkflowApprovalRequestV1,
    approve_workflow_request,
    build_execution_resource_spec,
    build_producer_edge_rule,
    build_workflow_approval_request,
    workflow_approval_request_json,
    workflow_execution_approval_json,
)

_WORKSPACE = "/tmp/chemsmart-approval-request-probe"


@pytest.fixture
def resources():
    return build_execution_resource_spec(
        execution_target="run",
        cores=2,
        memory_gb=4.0,
        gpu_count=0,
        scratch_policy="none",
        node_timeout_seconds=600,
    )


def _binding(node_id, jobtype, **overrides):
    fields = dict(
        node_id=node_id,
        program="pyscf",
        engine="cpu",
        jobtype=jobtype,
        project_artifact_sha256="c" * 64,
        settings_sha256="d" * 64,
        charge=0,
        multiplicity=1,
        input_mode="initial",
        initial_artifact_id="butane-anti",
        initial_artifact_sha256="e" * 64,
        scientific_identity_sha256="f" * 64,
        producer_edge_sha256="",
    )
    fields.update(overrides)
    return ApprovedNodeBindingV1(**fields)


@pytest.fixture
def request_document(resources):
    return build_workflow_approval_request(
        request_id="butane-populations",
        workflow_id="butane-populations",
        workflow_sha256="a" * 64,
        task_spec_sha256="b" * 64,
        workspace=_WORKSPACE,
        resources=resources,
        node_bindings=(_binding("opt-anti", "opt"),),
    )


def test_a_request_carries_the_body_a_reviewer_must_sign(request_document):
    assert request_document.status == "unapproved"
    assert request_document.workflow_sha256 == "a" * 64
    assert request_document.node_bindings[0].node_id == "opt-anti"
    assert request_document.workspace.endswith(
        "chemsmart-approval-request-probe"
    )


def test_a_request_is_not_the_envelope_agent_run_reads(request_document):
    """Inertness is the safety property, not a naming convention."""

    envelope = json.loads(workflow_approval_request_json(request_document))
    assert "workflow_approval" not in envelope
    assert envelope["workflow_approval_request"]["status"] == "unapproved"


def test_a_request_cannot_be_constructed_already_approved(request_document):
    body = {
        key: value
        for key, value in request_document.__dict__.items()
        if key != "request_sha256"
    }
    body["status"] = "approved"
    with pytest.raises(ContractError, match="never approved"):
        WorkflowApprovalRequestV1(
            **body, request_sha256=request_document.request_sha256
        )


def test_approving_a_reviewed_request_yields_the_run_envelope(
    request_document, resources
):
    approval = approve_workflow_request(
        request_document,
        approval_id="butane-run-1",
        approved_request_sha256=request_document.request_sha256,
        resources=resources,
    )
    assert approval.status == "approved"
    assert approval.workflow_sha256 == request_document.workflow_sha256
    assert approval.node_bindings == request_document.node_bindings
    envelope = json.loads(workflow_execution_approval_json(approval))
    assert envelope["workflow_approval"]["status"] == "approved"


def test_approval_is_refused_when_the_body_changed_after_review(
    request_document, resources
):
    with pytest.raises(ContractError, match="does not match the request"):
        approve_workflow_request(
            request_document,
            approval_id="butane-run-1",
            approved_request_sha256="9" * 64,
            resources=resources,
        )


def test_approval_is_refused_when_the_compute_allocation_moved(
    request_document,
):
    """A plan previewed against two cores was not reviewed for sixteen."""

    larger = build_execution_resource_spec(
        execution_target="run",
        cores=16,
        memory_gb=4.0,
        gpu_count=0,
        scratch_policy="none",
        node_timeout_seconds=600,
    )
    with pytest.raises(ContractError, match="resources have changed"):
        approve_workflow_request(
            request_document,
            approval_id="butane-run-1",
            approved_request_sha256=request_document.request_sha256,
            resources=larger,
        )


def test_a_chained_workflow_keeps_its_producer_edge_through_approval(
    resources,
):
    """The opt -> hess handoff is the shape every thermochemistry case has."""

    edge = build_producer_edge_rule(
        producer_node_id="opt-anti",
        consumer_node_id="hess-anti",
        artifact_kind="geometry_xyz",
        selection_rule="final_geometry",
    )
    document = build_workflow_approval_request(
        request_id="butane-chain",
        workflow_id="butane-chain",
        workflow_sha256="a" * 64,
        task_spec_sha256="b" * 64,
        workspace=_WORKSPACE,
        resources=resources,
        node_bindings=(
            _binding("opt-anti", "opt"),
            _binding(
                "hess-anti",
                "hess",
                input_mode="producer",
                initial_artifact_id="",
                initial_artifact_sha256="",
                scientific_identity_sha256="",
                producer_edge_sha256=edge.edge_sha256,
            ),
        ),
        producer_edges=(edge,),
    )
    approval = approve_workflow_request(
        document,
        approval_id="butane-chain-run",
        approved_request_sha256=document.request_sha256,
        resources=resources,
    )
    assert approval.node("hess-anti").input_mode == "producer"
    assert approval.producer_edges[0].producer_node_id == "opt-anti"


def test_a_request_with_no_nodes_is_refused(resources):
    with pytest.raises(ContractError, match="at least one node"):
        build_workflow_approval_request(
            request_id="empty",
            workflow_id="empty",
            workflow_sha256="a" * 64,
            task_spec_sha256="b" * 64,
            workspace=_WORKSPACE,
            resources=resources,
            node_bindings=(),
        )


def test_a_tampered_request_body_fails_its_own_digest(request_document):
    body = {
        key: value
        for key, value in request_document.__dict__.items()
        if key != "request_sha256"
    }
    body["node_bindings"] = (_binding("opt-anti", "sp"),)
    with pytest.raises(ContractError, match="digest mismatch"):
        WorkflowApprovalRequestV1(
            **body, request_sha256=request_document.request_sha256
        )
