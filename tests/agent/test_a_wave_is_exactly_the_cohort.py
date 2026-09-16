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


def test_the_executor_asks_the_cohort_before_it_walks():
    """A barrier nothing reads is a paragraph, not a contract."""

    import inspect

    from chemsmart.agent import executor

    source = inspect.getsource(executor)
    assert "cohort_frontier" in source, (
        "the executor still walks the whole ready frontier, so a wave "
        "runs work the Agent did not ask for in it"
    )
