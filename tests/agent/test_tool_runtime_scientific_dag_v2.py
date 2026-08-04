from __future__ import annotations

from types import SimpleNamespace

from chemsmart.agent.execution import derive_ready_node_ids
from chemsmart.agent.experiments.program_management_fixtures import (
    ProgramManagementHostFixtureFactoryV1,
)
from chemsmart.agent.runtime.event_store import RuntimeEventStore
from chemsmart.agent.tool_runtime import (
    CommandCompiledToolHostV1,
    _project_v1_execution_run_state,
    _scientific_plan_from_v1_approval,
)


def test_plan_tool_projects_explicit_data_edge_into_scientific_dag(tmp_path):
    factory = ProgramManagementHostFixtureFactoryV1(
        source_tree_root=".",
        materialization_root=tmp_path / "fixture",
    )
    fixture = factory(SimpleNamespace(case_id="DS-PM-003"))
    store = RuntimeEventStore(
        tmp_path / "events" / "runtime.jsonl", session_id="session"
    )
    host = CommandCompiledToolHostV1(
        event_store=store,
        task_spec_sha256s=(fixture.public_context.task_spec_sha256,),
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
