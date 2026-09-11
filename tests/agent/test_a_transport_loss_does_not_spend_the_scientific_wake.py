"""A cycle the provider ended is not a cycle the science ended.

SUFFICIENCY-1's cycle 1 died on three inter-event timeouts after seven
good turns -- `provider transport or protocol failed`, no workflow, no
claim, no decision. The goal was re-woken to answer that, correctly.
But the guard allows one re-wake per goal and counts them all alike, so
the transport continuation consumed the goal's only opportunity and the
requirement wake this round exists to fire could never have fired,
however short the delivery turned out to be.

Transport continuation and scientific revision are different budgets.
Both remain bounded by the revision and wall grants, which is what
those grants are for.
"""

from __future__ import annotations

import pytest

from chemsmart.agent.terminal_states import (
    PROVIDER_TRANSPORT_TERMINAL_REASON,
    is_provider_transport_terminal,
)

pytestmark = pytest.mark.capability("rule:wake.refusal_is_a_deliverable")


def test_the_transport_terminal_is_one_word_both_organs_read():
    """The loop wrote the sentence and nothing compared against it."""

    assert is_provider_transport_terminal(PROVIDER_TRANSPORT_TERMINAL_REASON)
    assert not is_provider_transport_terminal("host readiness gates passed")
    assert not is_provider_transport_terminal("")

    import inspect

    from chemsmart.agent import loop

    source = inspect.getsource(loop)
    assert "PROVIDER_TRANSPORT_TERMINAL_REASON" in source
    assert f'"{PROVIDER_TRANSPORT_TERMINAL_REASON}"' not in source


def test_a_transport_continuation_leaves_the_scientific_wake_unspent(
    tmp_path,
):
    from chemsmart.agent.goal import GoalLedger

    from .test_an_analysis_only_cycle_does_not_freeze_an_empty_scope import (
        _goal,
    )

    ledger = GoalLedger(tmp_path / "goal")
    ledger.create(_goal())

    def scientific_wakes_spent():
        return [
            entry
            for entry in ledger.entries()
            if entry["kind"] == "rewake_opened"
            and not entry["payload"].get("transport_continuation")
        ]

    ledger.append(
        "rewake_opened",
        {"cycle": 1, "transport_continuation": True, "failure_report": {}},
    )
    assert scientific_wakes_spent() == []

    ledger.append(
        "rewake_opened",
        {"cycle": 2, "transport_continuation": False, "failure_report": {}},
    )
    assert len(scientific_wakes_spent()) == 1
