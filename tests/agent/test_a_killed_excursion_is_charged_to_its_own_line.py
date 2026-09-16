"""An excursion killed before its receipt is still an excursion.

An excursion node investigates an anomaly the host recorded, and the
charter is explicit: it is charged to the grant's own investigation line
"at planning, at launch, in the run outcome and in the goal ledger, never
to the engine-call budget", so the asked observable is never bought with
the grant.

`derive_run_outcome` learns which nodes are excursions from the execution
*receipt*'s own flag. That was harmless while only receipts were counted.
Once a reservation without a receipt became a spent call -- which it is,
and which the launch fence has always agreed -- a controller killed
mid-engine on an excursion node charged it to the engine line, taking the
call out of the budget that buys the science.

The reservation carries the flag now, which is the only place a killed
launch leaves anything at all.
"""

from __future__ import annotations


from chemsmart.agent.runtime.event_store import (
    RuntimeEventStore,
    engine_calls_spent,
)
from chemsmart.agent.terminal_states import derive_run_outcome, read_run_events

from .test_runtime_v2_launch_fence import _frontier


def _reserved(tmp_path, *, excursion):
    store = RuntimeEventStore(tmp_path / "events.jsonl", session_id="s")
    plan, materialized, approval, invocation = _frontier(tmp_path)
    store.reserve_workflow_node_launch(
        turn_id="turn-1",
        plan=plan,
        materialized_workflow=materialized,
        approval=approval,
        invocation=invocation,
        run_id="run.water-approval",
        timestamp="2026-09-16T00:00:00+00:00",
        excursion=excursion,
    )
    return read_run_events(tmp_path / "events.jsonl")


def test_a_killed_excursion_is_not_on_the_engine_line(tmp_path):
    events = _reserved(tmp_path, excursion="d" * 64)
    outcome = derive_run_outcome(events)
    assert outcome.engine_calls_consumed == 0, (
        "an excursion killed before its receipt was charged to the "
        "engine-call budget, so the anomaly investigation was paid for "
        "out of the science"
    )
    assert outcome.excursion_calls_consumed == 1


def test_a_killed_ordinary_node_is_still_on_the_engine_line(tmp_path):
    """Not a weakening: the charge this round added still lands."""

    events = _reserved(tmp_path, excursion="")
    outcome = derive_run_outcome(events)
    assert outcome.engine_calls_consumed == 1
    assert outcome.excursion_calls_consumed == 0


def test_the_fence_and_the_outcome_agree_on_an_excursion(tmp_path):
    """Two authorities, one answer -- which is why this was found."""

    events = _reserved(tmp_path, excursion="d" * 64)
    outcome = derive_run_outcome(events)
    assert engine_calls_spent(
        events, excursion_node_ids=frozenset({"sp-initial"})
    ) == outcome.engine_calls_consumed
    assert engine_calls_spent(
        events,
        excursion_node_ids=frozenset({"sp-initial"}),
        excursions=True,
    ) == outcome.excursion_calls_consumed
