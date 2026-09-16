"""The barrier asks about a node and was reading the run's own word.

`workflow_node_state_changed` carries both:

    "node_id": node_id, "node_state": new_state,
    "workflow_state": updated.state, "record": updated_record

`record` is the **run's** state record, and `record["state"]` is the
run's summary word -- ambiguous, running, validated, failed, blocked,
cancelled. Five of those six are also node-terminal words, so reading
`record["state"]` is a type confusion that type-checks.

It fails in both directions, which is why it has to be driven through
the real emitter rather than a hand-made payload:

- **hang**: a member cancelled while a sibling still runs is written
  when the run reads `running`, so it is never counted finished, and the
  wave waits on it forever -- exactly the case `cohort_completion`'s own
  docstring says it fixed.
- **early wake**: once any member ends `ambiguous` the run's word is
  `ambiguous` and outranks the rest, so the *next* member's ordinary
  `pending -> running` transition marks it finished while its engine is
  starting, and the Agent is woken over a wave still running.

Every witness here builds its events with
`RuntimeEventStore.transition_workflow_run_node`. The previous ones
wrote `{"record": {"node_id": n, "state": s}}` -- a payload the product
has never emitted -- and so measured their own scaffolding.
"""

from __future__ import annotations

from chemsmart.agent.cohort import build_cohort_manifest, cohort_completion
from chemsmart.agent.runtime.event_store import RuntimeEventStore


def _wave(tmp_path, node_ids):
    build_cohort_manifest(
        goal_id="g1",
        cycle=1,
        bundle_sha256="e" * 64,
        node_ids=node_ids,
        max_concurrent_tasks=4,
        created_at="2026-09-16T00:00:00+00:00",
    ).write(tmp_path)


def _store(tmp_path):
    """A real store over a real two-member wave.

    Built here rather than borrowed, because the shared fixture plans one
    node and a barrier needs a sibling: the whole bug is what the *run's*
    word says while another member is still going.
    """

    from chemsmart.agent._contracts import canonical_sha256
    from chemsmart.agent.execution import (
        ProgramExecutionInvocationV1,
        build_execution_resource_spec,
        build_frozen_workflow_approval,
    )
    from chemsmart.agent.workflows import (
        MaterializedNodeV1,
        ScientificWorkflowNodeV2,
        build_materialized_workflow,
        build_scientific_workflow_plan,
    )

    nodes = tuple(
        ScientificWorkflowNodeV2(
            node_id=node_id,
            stage="sp",
            requested_program="pyscf",
            program="pyscf",
            engine="cpu",
            project_role="water-project",
            unresolved_fields=(),
        )
        for node_id in ("sp-a", "sp-b")
    )
    plan = build_scientific_workflow_plan(
        workflow_id="wave-workflow",
        task_spec_sha256="a" * 64,
        scientific_identity_sha256="b" * 64,
        nodes=nodes,
    )
    resources = build_execution_resource_spec(
        execution_target="run",
        cores=4,
        memory_gb=4,
        gpu_count=0,
        scratch_policy="none",
        node_timeout_seconds=600,
    )
    body = {
        "schema_version": "chemsmart.program-execution-invocation.v1",
        "node_id": "sp-a",
        "approval_sha256": "4" * 64,
        "program": "pyscf",
        "engine": "cpu",
        "jobtype": "sp",
        "project_sha256": "e" * 64,
        "input_artifact_id": "water-xyz",
        "input_sha256": "d" * 64,
        "scientific_identity_sha256": plan.scientific_identity_sha256,
        "environment_receipt_sha256": "1" * 64,
        "resource_sha256": resources.resource_sha256,
        "workspace": str(tmp_path.resolve()),
        "argv": ("chemsmart", "run", "pyscf", "sp"),
        "idempotency_key": "5" * 64,
        "status": "ready",
    }
    invocation = ProgramExecutionInvocationV1(
        **body, invocation_sha256=canonical_sha256(body)
    )
    materialized = build_materialized_workflow(
        plan=plan,
        live_cli_schema_sha256="c" * 64,
        resource_sha256=resources.resource_sha256,
        nodes=tuple(
            MaterializedNodeV1(
                node_id=node_id,
                input_artifact_sha256="d" * 64,
                project_artifact_sha256="e" * 64,
                project_validation_receipt_sha256="f" * 64,
                environment_receipt_sha256="1" * 64,
                invocation_sha256=invocation.invocation_sha256,
                preflight_receipt_sha256="3" * 64,
                state="previewed",
            )
            for node_id in ("sp-a", "sp-b")
        ),
        unresolved_node_ids=(),
        status="ready_for_approval",
    )
    approval = build_frozen_workflow_approval(
        approval_id="wave-approval",
        plan=plan,
        materialized_workflow=materialized,
        resources=resources,
        environment_identity_sha256s=("1" * 64,),
    )
    store = RuntimeEventStore(tmp_path / "events.jsonl", session_id="s")
    result = store.reserve_workflow_node_launch(
        turn_id="turn-1",
        plan=plan,
        materialized_workflow=materialized,
        approval=approval,
        invocation=invocation,
        run_id="run.wave-approval",
        timestamp="2026-09-16T00:00:00+00:00",
    )
    return store, result.run_state.run_id, plan


def _node_ids(plan):
    return tuple(str(node.node_id) for node in plan.nodes)


def test_a_member_cancelled_while_a_sibling_runs_is_finished(tmp_path):
    """The hang direction, through the real emitter."""

    store, run_id, plan = _store(tmp_path)
    first, second = _node_ids(plan)
    _wave(tmp_path, (first, second))

    # The sibling is still running, so the *run* reads "running" at the
    # moment the cancellation is written.
    _, state = store.transition_workflow_run_node(
        turn_id="t1",
        run_id=run_id,
        node_id=second,
        new_state="cancelled",
        failure_rule_ids=("execution.cancelled.human",),
        timestamp="2026-09-16T00:00:01+00:00",
    )
    assert state.state == "running", (
        "the fixture no longer reproduces the case: the run's own word "
        f"is {state.state!r}"
    )
    complete, pending = cohort_completion(tmp_path)
    assert second not in pending, (
        f"a cancelled member is still pending ({pending}); the barrier "
        "read the run's word instead of the node's, so this wave can "
        "never satisfy itself"
    )


def test_a_member_mid_engine_is_not_finished_by_a_siblings_timeout(
    tmp_path,
):
    """The early-wake direction: the costly one.

    Both members are already running -- which is what a wave is -- when
    one times out. The run's summary word becomes `ambiguous` and
    outranks everything, so the other member's ordinary
    `running -> engine_complete` transition was written under it and
    counted as terminal, although `engine_complete` is precisely the
    state that still owes a validation. The Agent would be woken over a
    member whose result had not been judged, and `derive_run_outcome`
    would read a stream that element was still appending to.
    """

    store, run_id, plan = _store(tmp_path)
    first, second = _node_ids(plan)
    _wave(tmp_path, (first, second))

    # A wave runs its members concurrently: the sibling is already up.
    store.transition_workflow_run_node(
        turn_id="t1",
        run_id=run_id,
        node_id=second,
        new_state="running",
        plan=plan,
        invocation_sha256="7" * 64,
        timestamp="2026-09-16T00:00:01+00:00",
    )
    store.transition_workflow_run_node(
        turn_id="t1",
        run_id=run_id,
        node_id=first,
        new_state="ambiguous",
        failure_rule_ids=("execution.process.timeout",),
        timestamp="2026-09-16T00:00:02+00:00",
    )
    _, state = store.transition_workflow_run_node(
        turn_id="t1",
        run_id=run_id,
        node_id=second,
        new_state="engine_complete",
        execution_receipt_sha256="8" * 64,
        output_artifact_sha256s=("9" * 64,),
        timestamp="2026-09-16T00:00:03+00:00",
    )
    assert state.state == "ambiguous", (
        "the fixture no longer reproduces the case: the run's own word "
        f"is {state.state!r}"
    )

    complete, pending = cohort_completion(tmp_path)
    assert complete is False and second in pending, (
        "a member that has produced an engine result and not yet been "
        "validated was counted as ended, so the wave would wake the "
        "Agent over an unjudged result"
    )
