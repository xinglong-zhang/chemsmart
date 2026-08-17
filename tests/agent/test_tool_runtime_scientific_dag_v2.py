from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from chemsmart.agent.execution import derive_ready_node_ids
from chemsmart.agent._contracts import ContractError, canonical_sha256
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.projects import ProjectValidationReceiptV1
from chemsmart.agent.tool_runtime import (
    CommandCompiledToolHostV1,
    _project_v1_execution_run_state,
    _scientific_plan_from_v1_approval,
)
from tests.agent.neutral_workflow_fixture import (
    build_neutral_workflow_fixture,
)


def test_plan_tool_projects_explicit_data_edge_into_scientific_dag(tmp_path):
    fixture = build_neutral_workflow_fixture(tmp_path / "fixture")
    store = RuntimeEventStore(
        tmp_path / "events" / "runtime.jsonl", session_id="session"
    )
    host = CommandCompiledToolHostV1(
        event_store=store,
        task_spec_sha256s=(fixture.public_context.task_spec_sha256,),
        approved_workspace=tmp_path / "workspace",
        **fixture.host_inputs,
    )
    action = next(
        item
        for item in fixture.public_context.next_actions
        if item.tool_name == "plan_command_workflow"
    )

    observed = host.dispatch(
        turn_id="turn-1",
        tool_name=action.tool_name,
        arguments=dict(action.fields),
    )["result"]

    plan = next(iter(host.scientific_workflow_plans.values()))
    assert observed["scientific_workflow_plan"]["plan_sha256"] == (
        plan.plan_sha256
    )
    assert tuple(edge.edge_kind for edge in plan.edges) == ("data",)
    assert plan.edges[0].source_node_id == "node.opt"
    assert plan.edges[0].target_node_id == "node.hess"
    event = store.read_events()[-1]
    assert event.payload["scientific_plan_record"]["plan_sha256"] == (
        plan.plan_sha256
    )


def test_rejected_replan_cannot_displace_the_latest_observed_workflow(
    tmp_path,
):
    """Planned termination remains bound to the accepted immutable plan.

    An execution session may try to replan after loading a frozen approval.
    The approval correctly rejects a different scientific plan.  That failed
    tool call must not make its unobserved command draft the host's latest
    workflow receipt, or Runtime V2 cannot terminate honestly as ``planned``.
    """

    fixture = build_neutral_workflow_fixture(tmp_path / "fixture")
    store = RuntimeEventStore(
        tmp_path / "events" / "runtime.jsonl", session_id="session"
    )
    host = CommandCompiledToolHostV1(
        event_store=store,
        task_spec_sha256s=(fixture.public_context.task_spec_sha256,),
        approved_workspace=tmp_path / "workspace",
        **fixture.host_inputs,
    )
    action = next(
        item
        for item in fixture.public_context.next_actions
        if item.tool_name == "plan_command_workflow"
    )
    accepted_arguments = dict(action.fields)
    accepted = host.dispatch(
        turn_id="turn-1",
        tool_name=action.tool_name,
        arguments=accepted_arguments,
    )["result"]
    accepted_draft_sha256 = accepted["workflow_draft"]["draft_sha256"]
    accepted_plan_sha256 = accepted["scientific_workflow_plan"][
        "plan_sha256"
    ]
    frozen_approval = SimpleNamespace(plan_sha256=accepted_plan_sha256)
    host.frozen_workflow_approval = frozen_approval

    revised_nodes = [dict(node) for node in accepted_arguments["nodes"]]
    revised_nodes[0]["project_role"] = "different-project-role"
    rejected_arguments = {**accepted_arguments, "nodes": revised_nodes}

    with pytest.raises(
        ContractError,
        match="planned workflow differs from frozen execution approval",
    ):
        host.dispatch(
            turn_id="turn-1",
            tool_name=action.tool_name,
            arguments=rejected_arguments,
        )

    assert host.frozen_workflow_approval is frozen_approval
    assert tuple(host.workflow_drafts) == (accepted_draft_sha256,)
    assert host.latest_workflow_draft_receipt() == accepted_draft_sha256
    assert store.state().workflow_receipts == [accepted_draft_sha256]

    store.terminate(
        turn_id="turn-1",
        terminal_state="planned",
        reason="tool budget exhausted after rejected replan",
        required_receipt_sha256s=(accepted_draft_sha256,),
    )
    assert store.state().terminal_state == "planned"


def test_v1_approval_projection_does_not_invent_tuple_order_dependencies():
    bindings = (
        SimpleNamespace(
            node_id="sp-initial",
            program="pyscf",
            engine="cpu",
            jobtype="sp",
            scientific_identity_sha256="b" * 64,
        ),
        SimpleNamespace(
            node_id="opt-initial",
            program="pyscf",
            engine="cpu",
            jobtype="opt",
            scientific_identity_sha256="b" * 64,
        ),
        SimpleNamespace(
            node_id="hess-optimized",
            program="pyscf",
            engine="cpu",
            jobtype="hess",
            scientific_identity_sha256="",
        ),
    )
    edges = (
        SimpleNamespace(
            producer_node_id="opt-initial",
            consumer_node_id="hess-optimized",
            artifact_kind="geometry_xyz",
            selection_rule="validated_final_geometry",
        ),
    )
    approval = SimpleNamespace(
        approval_id="water-approval",
        workflow_id="water-workflow",
        task_spec_sha256="a" * 64,
        approval_sha256="c" * 64,
        node_bindings=bindings,
        producer_edges=edges,
    )

    plan = _scientific_plan_from_v1_approval(approval)
    run = _project_v1_execution_run_state(plan, approval, {})

    assert derive_ready_node_ids(plan, run) == (
        "sp-initial",
        "opt-initial",
    )
    assert tuple(
        (edge.source_node_id, edge.target_node_id) for edge in plan.edges
    ) == (("opt-initial", "hess-optimized"),)


def test_validated_project_repair_can_rebind_without_rewriting_analysis_dag(
    tmp_path,
):
    fixture = build_neutral_workflow_fixture(tmp_path / "fixture")
    store = RuntimeEventStore(
        tmp_path / "events" / "runtime.jsonl", session_id="session"
    )
    host = CommandCompiledToolHostV1(
        event_store=store,
        task_spec_sha256s=(fixture.public_context.task_spec_sha256,),
        approved_workspace=tmp_path / "workspace",
        **fixture.host_inputs,
    )
    action = next(
        item
        for item in fixture.public_context.next_actions
        if item.tool_name == "plan_command_workflow"
    )
    calculation_nodes = []
    for node in dict(action.fields)["nodes"]:
        calculation_nodes.append(
            {
                **node,
                "project_role": "project.pyscf.cpu-hf",
                "produces_observables": ["vibrational_frequencies"],
                "support_state": "planned",
                "blocked_reason": "",
            }
        )
    planned = host.dispatch(
        turn_id="turn-1",
        tool_name="plan_scientific_workflow",
        arguments={
            "plan_id": "water-frequency-plan",
            "workflow_id": "workflow.ds-pm-003",
            "task_spec_id": fixture.public_context.task_spec_sha256,
            "calculation_nodes": calculation_nodes,
            "analysis_nodes": [
                {
                    "node_id": "extract-frequencies",
                    "analysis_kind": "result_extraction",
                    "dependencies": ["node.hess"],
                    "inputs": [
                        {
                            "input_id": "hessian-result",
                            "source_kind": "program_output",
                            "producer_node_id": "node.hess",
                            "producer_output_id": "structured_result",
                        }
                    ],
                    "selectors": [
                        {
                            "quantity_id": "frequency-values",
                            "selector": "frequencies",
                        }
                    ],
                    "outputs": [
                        {
                            "output_id": "frequency-values",
                            "quantity_kind": "vibrational_frequencies",
                            "unit": "cm^-1",
                        }
                    ],
                    "expression_nodes": [],
                    "expression_output_node_ids": [],
                    "support_state": "planned",
                    "blocked_reason": "",
                    "validation_rules": [],
                }
            ],
            "required_output_ids": ["frequency-values"],
        },
    )["result"]
    original_analysis = planned["scientific_toolchain_plan"]["analysis_nodes"]

    original_project = host.artifacts["project.pyscf.cpu-hf"]
    sections = yaml.safe_load(
        Path(original_project.path).read_text(encoding="utf-8")
    )
    rendered = host.dispatch(
        turn_id="turn-1",
        tool_name="render_project_yaml",
        arguments={"program": "pyscf", "sections": sections},
    )["result"]
    promoted = host.dispatch(
        turn_id="turn-1",
        tool_name="promote_project_yaml",
        arguments={
            "render_receipt_sha256": rendered["receipt_sha256"],
            "artifact_id": "primary-v2",
        },
    )["result"]
    replacement_id = promoted["artifact"]["artifact_id"]
    for jobtype in ("opt", "hess"):
        capability = next(
            receipt
            for receipt in host.capabilities.values()
            if receipt.query.program == "pyscf"
            and receipt.query.jobtype == jobtype
            and receipt.query.engine == "cpu"
        )
        body = {
            "schema_version": "chemsmart.project-validation-receipt.v1",
            "project_artifact_id": replacement_id,
            "project_sha256": host.artifacts[replacement_id].sha256,
            "capability_receipt_sha256": capability.receipt_sha256,
            "program": "pyscf",
            "jobtype": jobtype,
            "loader_id": "test.loader",
            "settings_sha256": canonical_sha256({}),
            "settings": (),
            "status": "valid",
            "error_class": "",
            "diagnostic": "",
            "rule_ids": ("project.loader.valid",),
        }
        receipt = ProjectValidationReceiptV1(
            **body, receipt_sha256=canonical_sha256(body)
        )
        host.project_validations[receipt.receipt_sha256] = receipt
    assert host._project_workflow_binding_observation(replacement_id)[
        "status"
    ] == "unbound"

    revised = host.dispatch(
        turn_id="turn-1",
        tool_name="amend_scientific_workflow",
        arguments={
            "workflow_id": "workflow.ds-pm-003",
            "project_replacements": [
                {"node_id": "node.opt", "project_role": replacement_id},
                {"node_id": "node.hess", "project_role": replacement_id},
            ],
        },
    )["result"]

    assert {
        node["project_role"]
        for node in revised["calculation_plan"]["workflow_draft"]["nodes"]
    } == {"primary-v2"}
    assert revised["scientific_toolchain_plan"]["analysis_nodes"] == (
        original_analysis
    )
