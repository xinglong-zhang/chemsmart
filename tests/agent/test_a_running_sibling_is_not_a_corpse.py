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

    terminal_source = Path(terminal_states.__file__).read_text(
        encoding="utf-8"
    )
    ast.parse(terminal_source)
    assert "reservation_lease_is_live" in terminal_source

    executor_source = Path(executor.__file__).read_text(encoding="utf-8")
    ast.parse(executor_source)
    assert "run_live_leases" in executor_source, (
        "the executor still concludes a running node is dead without "
        "asking whether its lease is live"
    )
    # The first wiring read `reserved_at`/`lease_seconds` off
    # WorkflowNodeRunStateV1, which has neither, so the branch could
    # never be taken and a live sibling was still written down as
    # interrupted. A getattr against a missing attribute is how a
    # contract looks wired while being unreachable.
    assert 'getattr(node_state, "lease_seconds"' not in executor_source
    assert 'getattr(node_state, "reserved_at"' not in executor_source


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


def test_the_lease_survives_the_builder_that_mints_the_reservation():
    """Drive the producer, not its spelling.

    The first version of this module asserted only that
    `reservation_lease_is_live` appeared in two files. It did -- and the
    builder that mints reservations accepted `lease_seconds` and
    `reserver` and then never put them in the record body, so every
    reservation carried the defaults and the helper always answered
    False. The whole lease was dead behind a green test. A witness that
    does not drive the producer proves the name exists, not the contract.
    """

    import inspect

    from chemsmart.agent.runtime import records

    source = inspect.getsource(records.build_workflow_node_launch_reservation)
    assert "lease_seconds" in source.split("body = {")[1], (
        "the builder takes a lease and drops it before the record is "
        "constructed, so reservation_lease_is_live can never say yes"
    )


def test_a_reservation_record_round_trips_its_lease():
    from chemsmart.agent.runtime.records import (
        workflow_node_launch_reservation_from_record,
    )

    record = {
        "schema_version": "chemsmart.workflow-node-launch-reservation.v1",
        "reservation_id": "r1",
        "run_id": "run1",
        "workflow_id": "w1",
        "node_id": "n1",
        "plan_sha256": "a" * 64,
        "materialized_workflow_sha256": "b" * 64,
        "approval_id": "ap1",
        "approval_sha256": "c" * 64,
        "invocation_sha256": "d" * 64,
        "consumes_approval": True,
        "state": "running",
        "reserved_at": _stamp(5),
        "admission_sha256": "e" * 64,
        "data_edge_binding_sha256s": (),
        "lease_seconds": 900,
        "reserver": "host 1 SLURM_ARRAY_TASK_ID=2",
    }
    from chemsmart.agent._contracts import canonical_sha256

    record["reservation_sha256"] = canonical_sha256(record)
    reservation = workflow_node_launch_reservation_from_record(record)
    assert reservation.lease_seconds == 900
    assert reservation_lease_is_live(
        reserved_at=reservation.reserved_at,
        lease_seconds=reservation.lease_seconds,
    )
