"""A wave the host cannot dispatch is typed evidence, not an exception.

Readiness is a host/DAG fact and keeps its single authority
(`derive_ready_node_ids`). Which of the ready calculations to run
together is the Agent's scientific strategy. So the host's job is to say,
per proposed member, what it sees -- and to say it as something the Agent
can read and select again from.

It must not raise. A refusal after a goal is under way is an agent-tool
failure: it teaches a session to carry workarounds and control logic for
a question the host owns (owner ruling, 2026-09-16).
"""

from __future__ import annotations

from chemsmart.agent.cohort import validate_wave


def test_a_wave_of_ready_independent_calculations_is_dispatchable():
    verdict = validate_wave(
        proposed=("a1", "a2", "a3"),
        ready=("a1", "a2", "a3", "a4"),
        edges=(),
    )
    assert verdict.dispatchable
    assert verdict.members == ("a1", "a2", "a3")
    assert all(row.status == "ready" for row in verdict.rows)


def test_a_member_that_is_not_ready_says_what_it_waits_on():
    """A producer outside the wave: the member simply is not ready yet."""

    verdict = validate_wave(
        proposed=("a1", "c1"),
        ready=("a1",),
        edges=(("b1", "c1"),),
    )
    assert not verdict.dispatchable
    rows = {row.node_id: row for row in verdict.rows}
    assert rows["a1"].status == "ready"
    assert rows["c1"].status == "not_ready"
    assert "b1" in rows["c1"].detail, (
        "the Agent is told the member is not ready but not what it is "
        "waiting for, so it cannot choose a different wave"
    )


def test_a_producer_inside_the_wave_is_named_as_the_dependency():
    """Both facts are true -- b1 is not ready *and* it consumes a1 -- and
    the one the Agent can act on is the sibling, because the fix is to
    take a1 now and choose b1 after reading its evidence."""

    verdict = validate_wave(
        proposed=("a1", "b1"),
        ready=("a1",),
        edges=(("a1", "b1"),),
    )
    assert not verdict.dispatchable
    rows = {row.node_id: row for row in verdict.rows}
    assert rows["b1"].status == "depends_on"
    assert "a1" in rows["b1"].detail


def test_two_members_of_one_chain_are_named_as_ordered():
    """Independence is the other half: A and B may both be ready in a
    later wave, but if B consumes A they are one experiment, not two."""

    verdict = validate_wave(
        proposed=("a1", "b1"),
        ready=("a1", "b1"),
        edges=(("a1", "b1"),),
    )
    assert not verdict.dispatchable
    rows = {row.node_id: row for row in verdict.rows}
    assert rows["b1"].status == "depends_on"
    assert "a1" in rows["b1"].detail


def test_a_node_outside_the_plan_is_named_rather_than_guessed():
    verdict = validate_wave(proposed=("ghost",), ready=("a1",), edges=())
    rows = {row.node_id: row for row in verdict.rows}
    assert rows["ghost"].status == "not_ready"
    assert not verdict.dispatchable


def test_nothing_raises_whatever_the_agent_proposes():
    """Every shape the model could get wrong returns a verdict."""

    for proposed in ((), ("a1", "a1"), ("",), ("a1", "b1", "ghost")):
        verdict = validate_wave(
            proposed=proposed, ready=("a1",), edges=(("a1", "b1"),)
        )
        assert isinstance(verdict.dispatchable, bool)
        assert verdict.rows is not None


def test_an_empty_wave_is_not_dispatchable_and_says_so():
    verdict = validate_wave(proposed=(), ready=("a1",), edges=())
    assert not verdict.dispatchable
    assert "empty" in verdict.summary.lower()


def test_a_repeated_member_is_one_calculation():
    verdict = validate_wave(proposed=("a1", "a1"), ready=("a1",), edges=())
    assert verdict.members == ("a1",)
    assert verdict.dispatchable


def test_the_verdict_reads_as_a_record():
    verdict = validate_wave(proposed=("a1",), ready=("a1",), edges=())
    record = verdict.public_record()
    assert record["dispatchable"] is True
    assert record["rows"][0]["node_id"] == "a1"
    assert record["rows"][0]["status"] == "ready"
