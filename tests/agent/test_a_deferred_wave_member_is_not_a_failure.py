"""A calculation the next wave will run has no ending to answer yet.

The settlement asks which terminal states a revision can answer. A node
that was never launched is not one of them -- correctly, for a launch
that should have happened. But under a wave cohort, *every approved node
outside this wave* is `not_launched` by construction and by design: the
Agent chose to see wave A's evidence before deciding wave B, which is
the entire point of the barrier.

So the first barrier of every multi-wave goal settled
`returned_to_human` with "the run ended in a state no revision can
answer: b1=not_launched, b2=not_launched" -- and the only cohort that
escaped it was one containing the whole approved partition, where the
barrier does nothing. The sequence the round exists for could not run.

The precedent is three lines up in the same block: a stage the plan
declared non-executable is excluded because "it has no ending to
answer." A member of a later wave is the same category, and the cohort
manifest -- digest-bound, written before any element started -- is the
authority that says which nodes this cycle was ever going to run.
"""

from __future__ import annotations

from types import SimpleNamespace


from chemsmart.agent.cohort import build_cohort_manifest


def _driver(tmp_path, *, cohort=None, nodes=()):
    from chemsmart.agent.driver import GoalDriver

    from .test_the_goal_loop_recovers_or_returns import _envelope_file

    workspace = tmp_path / "ws"
    workspace.mkdir(exist_ok=True)
    driver = GoalDriver(
        task="t",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path),
        goal_id="g1",
        granted_by="tester",
    )
    driver.cycles = 1
    driver.run_directory = driver.goal_dir / "runs" / "cycle-1"
    driver.run_directory.mkdir(parents=True, exist_ok=True)
    if cohort is not None:
        build_cohort_manifest(
            goal_id="g1",
            cycle=1,
            bundle_sha256="e" * 64,
            node_ids=tuple(cohort),
            max_concurrent_tasks=4,
            created_at="2026-09-16T00:00:00+00:00",
        ).write(driver.run_directory)
    driver.outcome = SimpleNamespace(
        nodes=tuple(
            SimpleNamespace(node_id=node_id, state=state)
            for node_id, state in nodes
        )
    )
    return driver


def test_a_member_of_a_later_wave_is_not_read_as_a_dead_end(tmp_path):
    """Wave A validated; wave B has not been chosen yet."""

    driver = _driver(
        tmp_path,
        cohort=("a1", "a2"),
        nodes=(
            ("a1", "validated"),
            ("a2", "validated"),
            ("b1", "not_launched"),
            ("b2", "not_launched"),
        ),
    )
    unanswerable = driver._unanswerable_terminal_states()
    assert unanswerable == {}, (
        "the first barrier of a multi-wave goal reads as a dead end: "
        f"{unanswerable}"
    )


def test_a_member_of_this_wave_that_never_launched_is_still_read(tmp_path):
    """Not a weakening: a launch that should have happened still counts."""

    driver = _driver(
        tmp_path,
        cohort=("a1", "a2"),
        nodes=(
            ("a1", "validated"),
            ("a2", "not_launched"),
            ("b1", "not_launched"),
        ),
    )
    assert driver._unanswerable_terminal_states() == {"a2": "not_launched"}


def test_a_failed_member_is_repairable_and_not_a_dead_end(tmp_path):
    driver = _driver(
        tmp_path,
        cohort=("a1", "a2"),
        nodes=(
            ("a1", "failed_wrong_stationary_point"),
            ("a2", "validated"),
            ("b1", "not_launched"),
        ),
    )
    unanswerable = driver._unanswerable_terminal_states()
    assert unanswerable == {"a1": "failed_wrong_stationary_point"}


def test_without_a_cohort_every_unlaunched_node_is_still_read(tmp_path):
    """A single-job dispatch runs the whole partition; a gap is a gap."""

    driver = _driver(
        tmp_path,
        cohort=None,
        nodes=(("a1", "validated"), ("b1", "not_launched")),
    )
    assert driver._unanswerable_terminal_states() == {"b1": "not_launched"}
