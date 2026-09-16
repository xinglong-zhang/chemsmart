"""A node that is running right now is not a process that died.

Three sites read ``state == "running"`` with no receipt and concluded the
same thing: a prior invocation was interrupted mid-engine and needs human
reconciliation. That was true while one process walked every node in
turn. Under a wave cohort it is the ordinary appearance of a *healthy
sibling*, and the reservation carried nothing that could tell the two
apart -- no host, no pid, no scheduler id, no expiry. The distinction was
provably underivable from durable state.

The lease closes it without new machinery, because the host already
bounds how long an engine may run: the approved node timeout. A
reservation younger than that bound may be a live engine and the cohort
simply is not terminal yet; one older than it cannot be, and keeps the
interrupted reading it has always had.
"""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path

from chemsmart.agent.runtime.records import reservation_lease_is_live


def _stamp(seconds_ago: float) -> str:
    moment = datetime.now(timezone.utc) - timedelta(seconds=seconds_ago)
    return moment.isoformat()


def test_a_fresh_reservation_is_a_live_sibling():
    assert reservation_lease_is_live(reserved_at=_stamp(30), lease_seconds=900)


def test_a_reservation_past_its_lease_cannot_be_a_live_engine():
    """The node timeout is the host's own bound on an engine's life."""

    assert not reservation_lease_is_live(
        reserved_at=_stamp(1800), lease_seconds=900
    )


def test_a_reservation_that_declares_no_lease_keeps_the_old_reading():
    """Every reservation already on disk predates the lease.

    It carries no bound, so nothing can be concluded about liveness from
    it and the interrupted reading stands -- which is what those streams
    were written under.
    """

    assert not reservation_lease_is_live(
        reserved_at=_stamp(1), lease_seconds=None
    )
    assert not reservation_lease_is_live(
        reserved_at=_stamp(1), lease_seconds=0
    )


def test_an_unreadable_timestamp_is_not_treated_as_live():
    """A lease the host cannot read is not a licence to wait forever."""

    assert not reservation_lease_is_live(
        reserved_at="not-a-timestamp", lease_seconds=900
    )
    assert not reservation_lease_is_live(reserved_at="", lease_seconds=900)


def test_the_lease_is_read_by_every_site_that_calls_running_dead():
    """A lease nothing reads is a field, not a contract."""

    import ast
    from pathlib import Path

    import chemsmart.agent.executor as executor
    import chemsmart.agent.terminal_states as terminal_states

    for module in (terminal_states, executor):
        source = Path(module.__file__).read_text(encoding="utf-8")
        ast.parse(source)
        assert "reservation_lease_is_live" in source, (
            f"{module.__name__} still concludes a running node is dead "
            "without asking whether its lease is live"
        )


def _reserved(node_id, reserved_at, lease_seconds):
    from types import SimpleNamespace

    return SimpleNamespace(
        kind="workflow_node_launch_reserved",
        payload={
            "node_id": node_id,
            "record": {
                "node_id": node_id,
                "reserved_at": reserved_at,
                "lease_seconds": lease_seconds,
            },
        },
    )


def _observed(node_id):
    from types import SimpleNamespace

    return SimpleNamespace(
        kind="program_execution_observed",
        payload={"node_id": node_id, "record": {"node_id": node_id}},
    )


def test_the_barrier_asks_which_members_are_still_running():
    """A wave is over when every member is terminal, and a member inside
    its lease is not terminal -- it is being run right now by another
    array element."""

    from chemsmart.agent.terminal_states import run_live_leases

    events = [
        _reserved("a", _stamp(10), 900),  # running now
        _reserved("b", _stamp(10), 900),
        _observed("b"),  # finished
        _reserved("c", _stamp(9000), 900),  # lease expired: gone
        _reserved("d", _stamp(10), None),  # pre-lease reservation
    ]
    assert run_live_leases(events) == ("a",)


def test_a_live_lease_never_invents_a_terminal_word():
    """`running` is not in NODE_TERMINAL_STATES and NodeTerminalStateV1
    refuses anything outside it, so a liveness answer must not be smuggled
    into the terminal vocabulary."""

    import chemsmart.agent.terminal_states as terminal_states
    from chemsmart.agent.terminal_states import NODE_TERMINAL_STATES

    assert "running" not in NODE_TERMINAL_STATES
    source = Path(terminal_states.__file__).read_text(encoding="utf-8")
    assert 'terminal = "running"' not in source
