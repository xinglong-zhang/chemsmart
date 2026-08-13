"""CPU image qualification must reject preview-only false greens."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from chemsmart.agent._contracts import (
    ContractError,
    canonical_data,
    canonical_sha256,
)
from chemsmart.agent.runtime.events import (
    OPTIMIZED_GEOMETRY_HANDED_OFF,
    PROGRAM_EXECUTED,
    RUNTIME_TERMINATED,
    SCIENTIFIC_WORKFLOW_MATERIALIZED,
    WORKFLOW_PLANNED,
    RuntimeEvent,
)
from chemsmart.agent.workflows import (
    MaterializedNodeV1,
    ScientificWorkflowEdgeV2,
    ScientificWorkflowNodeV2,
    build_materialized_workflow,
    build_scientific_workflow_plan,
)
from containers.cpu.qualification.validate_live_plan import (
    _EXECUTION_EVENT_KINDS,
    validate,
)


ROOT = Path(__file__).resolve().parents[2]


def _plan(
    stages=("sp", "opt"),
    *,
    edges=(),
    requested_program="pyscf",
    program="pyscf",
    engine="cpu",
):
    nodes = tuple(
        ScientificWorkflowNodeV2(
            node_id=f"water-{stage}",
            stage=stage,
            requested_program=requested_program,
            program=program,
            engine=engine,
            project_role="water-project",
            unresolved_fields=(),
        )
        for stage in stages
    )
    return build_scientific_workflow_plan(
        workflow_id="cpu-live-plan",
        task_spec_sha256="1" * 64,
        scientific_identity_sha256="2" * 64,
        nodes=nodes,
        edges=edges,
    )


def _workflow(
    plan, *, status="previewed", unresolved=(), different_inputs=False
):
    unresolved_ids = set(unresolved)
    nodes = tuple(
        MaterializedNodeV1(
            node_id=node.node_id,
            input_artifact_sha256=(
                ("3" if index == 0 or not different_inputs else "d") * 64
            ),
            project_artifact_sha256="4" * 64,
            project_validation_receipt_sha256="5" * 64,
            environment_receipt_sha256="6" * 64,
            invocation_sha256="7" * 64,
            preflight_receipt_sha256="8" * 64,
            state="previewed",
        )
        for index, node in enumerate(plan.nodes)
        if node.node_id not in unresolved_ids
    )
    return build_materialized_workflow(
        plan=plan,
        live_cli_schema_sha256="9" * 64,
        resource_sha256="a" * 64,
        nodes=nodes,
        unresolved_node_ids=tuple(unresolved_ids),
        status=status,
    )


def _write_stream(
    tmp_path,
    *,
    plan,
    workflow,
    earlier_workflow=None,
    execution_kind="",
    tamper=False,
):
    run = tmp_path / ".chemsmart-agent" / "runs" / "run-1"
    run.mkdir(parents=True)
    plan_receipt = "b" * 64
    materialized_record = canonical_data(workflow)
    payloads = [
        (
            WORKFLOW_PLANNED,
            {
                "receipt_sha256": plan_receipt,
                "status": "planned",
                "scientific_plan_sha256": plan.plan_sha256,
                "scientific_plan_record": canonical_data(plan),
            },
        ),
    ]
    if earlier_workflow is not None:
        earlier_record = canonical_data(earlier_workflow)
        payloads.append(
            (
                SCIENTIFIC_WORKFLOW_MATERIALIZED,
                {
                    "receipt_sha256": canonical_sha256(earlier_record),
                    "workflow_id": earlier_workflow.workflow_id,
                    "plan_sha256": earlier_workflow.plan_sha256,
                    "status": earlier_workflow.status,
                    "record": earlier_record,
                },
            )
        )
    payloads.append(
        (
            SCIENTIFIC_WORKFLOW_MATERIALIZED,
            {
                "receipt_sha256": canonical_sha256(materialized_record),
                "workflow_id": workflow.workflow_id,
                "plan_sha256": workflow.plan_sha256,
                "status": workflow.status,
                "record": materialized_record,
            },
        )
    )
    if execution_kind:
        payloads.append(
            (
                execution_kind,
                (
                    {
                        "receipt_sha256": "c" * 64,
                        "execution_state": "not_started",
                    }
                    if execution_kind == PROGRAM_EXECUTED
                    else {
                        "receipt_sha256": "c" * 64,
                        "status": "validated_handoff",
                    }
                ),
            )
        )
    payloads.append(
        (
            RUNTIME_TERMINATED,
            {
                "terminal_state": "planned",
                "reason": "preview transaction complete",
                "plan_receipt_sha256": plan_receipt,
            },
        )
    )
    events = []
    previous_hash = ""
    for sequence, (kind, payload) in enumerate(payloads, start=1):
        event = RuntimeEvent.create(
            sequence=sequence,
            session_id="cpu-live-plan-session",
            turn_id="turn-1",
            kind=kind,
            payload=payload,
            previous_hash=previous_hash,
        )
        events.append(event.to_dict())
        previous_hash = event.event_hash
    if tamper:
        events[1]["payload"]["status"] = "materialized"
    stream = run / "events.jsonl"
    stream.write_text(
        "".join(json.dumps(event, sort_keys=True) + "\n" for event in events),
        encoding="utf-8",
    )
    return stream


def test_live_plan_qualification_accepts_two_fully_previewed_nodes(tmp_path):
    plan = _plan()
    _write_stream(tmp_path, plan=plan, workflow=_workflow(plan))

    result = validate(tmp_path)

    assert result["status"] == "qualified"
    assert result["materialized_status"] == "previewed"
    assert result["unresolved_node_ids"] == ()
    assert result["planned_node_ids"] == ("water-opt", "water-sp")
    assert result["execution_event_count"] == 0


def test_packaged_live_task_and_ci_bind_the_preview_transaction():
    task = (ROOT / "containers/cpu/qualification/agent-task.txt").read_text(
        encoding="utf-8"
    )
    assert "exactly two independent" in task
    assert "single point" in task
    assert "geometry optimization" in task
    assert "both consume the same initial geometry" in task
    assert "neither consumes output from the other" in task
    assert "Hessian" not in task
    assert "frequency" not in task
    assert "optimized geometry" not in task

    workflow = (
        ROOT / ".github/workflows/cpu-container-candidate.yml"
    ).read_text(encoding="utf-8")
    assert "qualification/validate_live_plan.py" in workflow
    assert "live-plan-qualification.json" in workflow


def test_live_plan_qualification_rejects_partial_materialization(tmp_path):
    plan = _plan()
    _write_stream(
        tmp_path,
        plan=plan,
        workflow=_workflow(plan, status="partial", unresolved=("water-opt",)),
    )

    with pytest.raises(RuntimeError, match="latest.*not fully previewed"):
        validate(tmp_path)


def test_live_plan_qualification_rejects_partial_after_earlier_green(tmp_path):
    plan = _plan()
    _write_stream(
        tmp_path,
        plan=plan,
        earlier_workflow=_workflow(plan),
        workflow=_workflow(plan, status="partial", unresolved=("water-opt",)),
    )

    with pytest.raises(RuntimeError, match="latest.*not fully previewed"):
        validate(tmp_path)


@pytest.mark.parametrize("stages", [("sp",), ("sp", "opt", "hess")])
def test_live_plan_qualification_rejects_wrong_stage_set(tmp_path, stages):
    plan = _plan(stages)
    _write_stream(tmp_path, plan=plan, workflow=_workflow(plan))

    with pytest.raises(
        RuntimeError, match="exactly one SP and one optimization"
    ):
        validate(tmp_path)


@pytest.mark.parametrize(
    ("requested_program", "program", "engine"),
    [
        ("orca", "pyscf", "cpu"),
        ("pyscf", "orca", "cpu"),
        ("pyscf", "pyscf", "gpu"),
    ],
)
def test_live_plan_qualification_rejects_wrong_program_or_engine(
    tmp_path, requested_program, program, engine
):
    plan = _plan(
        requested_program=requested_program,
        program=program,
        engine=engine,
    )
    _write_stream(tmp_path, plan=plan, workflow=_workflow(plan))

    with pytest.raises(RuntimeError, match="requires PySCF CPU nodes"):
        validate(tmp_path)


def test_live_plan_qualification_rejects_dependent_stages(tmp_path):
    edge = ScientificWorkflowEdgeV2(
        edge_id="control.water-sp.water-opt",
        source_node_id="water-sp",
        target_node_id="water-opt",
        edge_kind="control",
    )
    plan = _plan(edges=(edge,))
    _write_stream(tmp_path, plan=plan, workflow=_workflow(plan))

    with pytest.raises(RuntimeError, match="independently resolvable"):
        validate(tmp_path)


def test_live_plan_qualification_rejects_different_inputs(tmp_path):
    plan = _plan()
    _write_stream(
        tmp_path,
        plan=plan,
        workflow=_workflow(plan, different_inputs=True),
    )

    with pytest.raises(RuntimeError, match="consume one initial geometry"):
        validate(tmp_path)


@pytest.mark.parametrize(
    "execution_kind", [PROGRAM_EXECUTED, OPTIMIZED_GEOMETRY_HANDED_OFF]
)
def test_live_plan_qualification_rejects_execution_event(
    tmp_path, execution_kind
):
    plan = _plan()
    _write_stream(
        tmp_path,
        plan=plan,
        workflow=_workflow(plan),
        execution_kind=execution_kind,
    )

    with pytest.raises(RuntimeError, match="observed an execution event"):
        validate(tmp_path)


def test_live_plan_execution_fence_includes_handoff_and_data_binding():
    from chemsmart.agent.runtime.events import WORKFLOW_DATA_EDGE_BOUND

    assert OPTIMIZED_GEOMETRY_HANDED_OFF in _EXECUTION_EVENT_KINDS
    assert WORKFLOW_DATA_EDGE_BOUND in _EXECUTION_EVENT_KINDS


def test_live_plan_qualification_rejects_tampered_hash_chain(tmp_path):
    plan = _plan()
    _write_stream(tmp_path, plan=plan, workflow=_workflow(plan), tamper=True)

    with pytest.raises(ContractError, match="hash verification failed"):
        validate(tmp_path)


def test_live_plan_qualification_rejects_unterminated_jsonl(tmp_path):
    plan = _plan()
    stream = _write_stream(tmp_path, plan=plan, workflow=_workflow(plan))
    stream.write_bytes(stream.read_bytes().removesuffix(b"\n"))

    with pytest.raises(RuntimeError, match="unterminated Runtime line"):
        validate(tmp_path)
