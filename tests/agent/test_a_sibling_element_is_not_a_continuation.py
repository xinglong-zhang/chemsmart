"""An element that starts a moment later is a sibling, not a re-invocation.

`run()` reads `durable.run_state is not None` as "a continuation: the
durable stream is the starting truth", and admits through
`continue_workflow_execution_approval_bundle`, which requires the
workspace consumption ledger.

A wave is N processes starting at once. Whichever gets there first
creates the run state; every element after it sees a run state and takes
the continuation path -- and a cohort never writes the workspace
consumption ledger, because membership replaces the claim. So the
siblings die with *"continuation requires the approval's consumption
ledger; this bundle was never claimed in this workspace"*.

Measured live on CUHK (goal `butane-wave-2`, array 2135285): elements 1
and 2 both failed in eight seconds while element 0 ran the Hessian. The
two earlier waves of the same goal survived only because all three
elements happened to reach the check before any of them had written the
run state -- the same code, the same wave, decided by milliseconds.

Membership is admission, at entry as at launch. A continuation is a
*later invocation* of an approval, and a member of the wave that is
running right now is not one.
"""

from __future__ import annotations

from types import SimpleNamespace


from chemsmart.agent._contracts import ContractError
from chemsmart.agent.cohort import build_cohort_manifest


def _executor(tmp_path, *, element, run_state):
    from chemsmart.agent.executor import ApprovedWorkflowExecutor

    run_directory = tmp_path / "run"
    run_directory.mkdir(parents=True, exist_ok=True)
    build_cohort_manifest(
        goal_id="g1",
        cycle=1,
        bundle_sha256="e" * 64,
        node_ids=("a1", "a2"),
        max_concurrent_tasks=4,
        created_at="2026-09-16T00:00:00+00:00",
    ).write(run_directory)

    frontier = SimpleNamespace(
        run_state=run_state,
        data_edge_bindings=(),
        plan=None,
        materialized_workflow=None,
        approval=None,
    )
    host = SimpleNamespace(
        event_store=SimpleNamespace(
            workflow_frontier=lambda **_kw: frontier
        ),
    )
    return ApprovedWorkflowExecutor(
        host=host,
        plan=SimpleNamespace(
            workflow_id="w",
            plan_sha256="b" * 64,
            nodes=(
                SimpleNamespace(node_id="a1", program="orca"),
                SimpleNamespace(node_id="a2", program="orca"),
            ),
        ),
        approval=SimpleNamespace(
            approval_id="a1-approval",
            node_bindings=(
                SimpleNamespace(node_id="a1", program="orca", jobtype="opt"),
                SimpleNamespace(node_id="a2", program="orca", jobtype="opt"),
            ),
        ),
        frozen_approval=SimpleNamespace(approval_sha256="c" * 64),
        initial_artifacts={},
        project_artifacts=(),
        task_spec_sha256="a" * 64,
        run_directory=run_directory,
        execution_bundle=SimpleNamespace(
            non_executable_node_ids=(), bundle_sha256="e" * 64
        ),
        approval_workspace=tmp_path / "workspace",
        claim_workspace_bundle=True,
        cohort_element=element,
    )


def test_a_cohort_element_does_not_claim_a_continuation(tmp_path):
    """The sibling arrives second and must not be refused for it."""

    run_state = SimpleNamespace(
        run_id="run.a1-approval",
        nodes=(
            SimpleNamespace(
                node_id="a1",
                state="running",
                failure_rule_ids=(),
                execution_receipt_sha256="",
            ),
        ),
    )
    executor = _executor(tmp_path, element=1, run_state=run_state)

    claimed: list[str] = []

    import chemsmart.agent.live_session as live_session

    original = live_session.continue_workflow_execution_approval_bundle

    def refuse(*_a, **_kw):
        claimed.append("continued")
        raise ContractError(
            "continuation requires the approval's consumption ledger; "
            "this bundle was never claimed in this workspace"
        )

    live_session.continue_workflow_execution_approval_bundle = refuse
    try:
        # Only the admission decision is under test; the walk itself
        # needs far more of a real host than this fixture supplies.
        assert executor._admits_as_continuation() is False, (
            "a member of the running wave was admitted as a later "
            "invocation of its own approval, and a cohort never writes "
            "the workspace consumption ledger a continuation demands"
        )
    finally:
        live_session.continue_workflow_execution_approval_bundle = original
    assert not claimed


def test_a_real_continuation_still_claims(tmp_path):
    """Not a weakening: a re-invocation is still admitted as one."""

    run_state = SimpleNamespace(
        run_id="run.a1-approval",
        nodes=(),
    )
    executor = _executor(tmp_path, element=None, run_state=run_state)
    assert executor._admits_as_continuation() is True


def test_a_first_element_has_no_continuation_to_admit(tmp_path):
    """Whoever creates the run state was never a continuation."""

    executor = _executor(tmp_path, element=0, run_state=None)
    assert executor._admits_as_continuation() is False
