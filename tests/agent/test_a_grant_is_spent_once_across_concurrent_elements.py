"""Four elements cannot each spend the whole grant.

`_require_bounded_launch_budget` counted `self.execution_receipts`, an
instance dict populated only by the nodes *this* process touched, and the
episode clock started at `time.monotonic()` in `__init__`. That was exact
while one process walked every node in turn. With one process per array
element it is exact per element and meaningless per approval: four
elements each see zero prior receipts and a fresh clock, so a grant of
three engine calls authorises four launches and an episode window of six
hours is handed out four times over.

Afterwards the ledger clamps the overrun to zero (`max(calls, 0)`), so
the human's next cycle reads "0 remaining" and never "you were charged
four times".

The count must come from the durable stream, and it must be taken under
the same exclusive lock that already serialises every launch -- otherwise
two elements both read "two spent of three" and both proceed.
"""

from __future__ import annotations

from chemsmart.agent.runtime.event_store import engine_calls_spent


def _reserved(node_id):
    from types import SimpleNamespace

    return SimpleNamespace(
        kind="workflow_node_launch_reserved",
        payload={"node_id": node_id, "record": {"node_id": node_id}},
    )


def _observed(node_id):
    from types import SimpleNamespace

    return SimpleNamespace(
        kind="program_execution_observed",
        payload={"node_id": node_id, "record": {"node_id": node_id}},
    )


def test_a_finished_engine_is_spent():
    assert engine_calls_spent([_observed("a"), _observed("b")]) == 2


def test_a_reserved_but_unfinished_engine_is_spent_too():
    """The call is spent when it is taken, not when it returns.

    Counting only receipts lets four elements reserve simultaneously
    against a grant of three: every one of them reads zero spent,
    because none has finished yet.
    """

    assert engine_calls_spent([_reserved("a"), _reserved("b")]) == 2


def test_one_node_reserved_and_finished_is_one_call():
    assert engine_calls_spent([_reserved("a"), _observed("a")]) == 1


def test_excursion_nodes_are_charged_to_their_own_line():
    events = [_reserved("a"), _observed("a"), _reserved("x"), _observed("x")]
    assert engine_calls_spent(events, excursion_node_ids=frozenset({"x"})) == 1
    assert (
        engine_calls_spent(
            events, excursion_node_ids=frozenset({"x"}), excursions=True
        )
        == 1
    )


def test_the_count_is_taken_where_the_launches_are_serialised():
    """A check outside the lock is two elements both reading 'two of
    three' and both proceeding."""

    import inspect

    from chemsmart.agent.runtime.event_store import RuntimeEventStore

    source = inspect.getsource(RuntimeEventStore.reserve_workflow_node_launch)
    assert "engine_calls_spent" in source, (
        "the budget is still counted outside the exclusive lock that "
        "serialises launches, so concurrent elements can both pass it"
    )
    assert "max_engine_calls" in source
