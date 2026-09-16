"""A wave runs its members and stops. It is not a flowing scheduler.

The executor walks the ready frontier until nothing is ready, which
maximises throughput and is the wrong shape for this contract. The Agent
asked for A1, A2 and A3 as one scientific experiment; their *collective*
evidence is what should trigger the next reasoning turn. If A1 finishes
early and B1's local dependency is satisfied, B1 must not start -- the
Agent has not seen A2 and A3 yet, and B1 is a decision it has not made.

That barrier is scientifically meaningful, not a scheduler detail: it is
what makes the wake a reasoning point rather than an interrupt (owner
ruling, 2026-09-16).
"""

from __future__ import annotations

from chemsmart.agent.cohort import cohort_frontier


def test_only_cohort_members_run():
    ready = ("a1", "a2", "b1")
    assert cohort_frontier(ready, ("a1", "a2", "a3")) == ("a1", "a2")


def test_a_node_whose_dependency_cleared_mid_wave_still_waits():
    """B1 became ready because A1 validated. It is not in the wave, so it
    does not run: the Agent chooses the next wave after seeing all three."""

    assert cohort_frontier(("b1",), ("a1", "a2", "a3")) == ()


def test_an_empty_cohort_frontier_is_how_a_wave_ends():
    assert cohort_frontier((), ("a1",)) == ()


def test_without_a_cohort_every_ready_node_runs():
    """A single-job dispatch and every recorded run keep the old walk."""

    assert cohort_frontier(("a1", "b1"), None) == ("a1", "b1")


def test_the_order_the_frontier_offered_is_kept():
    assert cohort_frontier(("a3", "a1"), ("a1", "a2", "a3")) == ("a3", "a1")


def _executor(tmp_path, *, element):
    """A real executor over a real manifest, as an array element is."""

    from types import SimpleNamespace

    from chemsmart.agent.cohort import build_cohort_manifest
    from chemsmart.agent.executor import ApprovedWorkflowExecutor

    run_directory = tmp_path / "run"
    run_directory.mkdir(parents=True, exist_ok=True)
    build_cohort_manifest(
        goal_id="g1",
        cycle=1,
        bundle_sha256="e" * 64,
        node_ids=("a1", "a2", "a3"),
        max_concurrent_tasks=4,
        created_at="2026-09-16T00:00:00+00:00",
    ).write(run_directory)
    return ApprovedWorkflowExecutor(
        host=SimpleNamespace(),
        plan=SimpleNamespace(
            workflow_id="w",
            plan_sha256="b" * 64,
            nodes=(
                SimpleNamespace(node_id="a1", program="orca"),
                SimpleNamespace(node_id="a2", program="orca"),
                SimpleNamespace(node_id="a3", program="orca"),
            ),
        ),
        approval=SimpleNamespace(node_bindings=()),
        frozen_approval=SimpleNamespace(approval_sha256="c" * 64),
        initial_artifacts={},
        project_artifacts=(),
        task_spec_sha256="a" * 64,
        run_directory=run_directory,
        execution_bundle=SimpleNamespace(non_executable_node_ids=()),
        approval_workspace=tmp_path / "workspace",
        claim_workspace_bundle=False,
        cohort_element=element,
    )


def test_an_array_element_runs_one_calculation_not_the_whole_wave(tmp_path):
    """An array element is one approved scientific calculation node.

    ``authorise_cohort_element`` resolved exactly which one and its
    answer was discarded, so every element bounded its walk to the whole
    cohort: three processes each tried to run all three members. The
    launch fence would have refused the duplicates -- after three
    processes had raced for one reservation and two had done the work of
    finding out they had lost.
    """

    scope = _executor(tmp_path, element=1)._cohort_scope()
    assert scope == ("a2",), (
        f"element 1 may run {scope}, so it walks members that belong to "
        "other elements"
    )
    # Composed with the frontier, which is what the walk actually does.
    assert cohort_frontier(("a1", "a2", "a3"), scope) == ("a2",)


def test_one_process_running_the_whole_wave_keeps_the_whole_wave(tmp_path):
    """A local dispatch of a cohort is one process running every member."""

    executor = _executor(tmp_path, element=None)
    assert executor._cohort_scope() == ("a1", "a2", "a3")


def test_without_a_manifest_the_walk_is_unbounded(tmp_path):
    """Every run recorded before waves existed keeps its own walk."""

    from types import SimpleNamespace

    from chemsmart.agent.executor import ApprovedWorkflowExecutor

    run_directory = tmp_path / "bare"
    run_directory.mkdir()
    executor = ApprovedWorkflowExecutor(
        host=SimpleNamespace(),
        plan=SimpleNamespace(
            workflow_id="w",
            plan_sha256="b" * 64,
            nodes=(SimpleNamespace(node_id="a1", program="orca"),),
        ),
        approval=SimpleNamespace(node_bindings=()),
        frozen_approval=SimpleNamespace(approval_sha256="c" * 64),
        initial_artifacts={},
        project_artifacts=(),
        task_spec_sha256="a" * 64,
        run_directory=run_directory,
        execution_bundle=SimpleNamespace(non_executable_node_ids=()),
        approval_workspace=tmp_path / "workspace",
        claim_workspace_bundle=False,
    )
    assert executor._cohort_scope() is None
