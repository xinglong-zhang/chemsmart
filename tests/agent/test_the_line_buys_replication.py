"""A recovery wake with a live excursion line says what the line buys:
replication under a named perturbation, never the asked observable.

Two E4 windows measured a line with nothing legal to buy; the design
note names replication as the class payable by construction.
"""

from __future__ import annotations

import pytest

from .test_the_goal_loop_recovers_or_returns import (
    _execute,
    _loop,
    _planning_session,
    _review_payload,
)

pytestmark = pytest.mark.capability("rule:wake.excursion_buys_replication")


def _second_context(tmp_path, excursions):
    contexts = []

    def capture(inner):
        def step(workspace, kwargs):
            contexts.append(kwargs["goal_context"])
            return inner(workspace, kwargs)

        return step

    tmp_path.mkdir(parents=True, exist_ok=True)
    _loop(
        tmp_path,
        excursions=excursions,
        sessions=[
            capture(_planning_session("live-1", review=_review_payload())),
            capture(_planning_session("live-2", review=_review_payload())),
        ],
        executes=[
            _execute(
                tmp_path, failed=True, status="partial", analysis="partial"
            ),
            _execute(
                tmp_path, failed=True, status="partial", analysis="partial"
            ),
        ],
        max_revisions=1,
    )
    return contexts[1]


def test_a_live_line_names_replication_as_its_purchase(tmp_path):
    context = _second_context(tmp_path / "live", excursions=2)
    assert context["budgets"]["excursion_calls_remaining"] == 2
    assert "buys replication before belief" in context["authority"]
    assert "identical input is not a perturbation" in context["authority"]


def test_a_dead_line_says_nothing_about_buying(tmp_path):
    context = _second_context(tmp_path / "dead", excursions=0)
    assert context["budgets"]["excursion_calls_remaining"] == 0
    assert "buys replication" not in context["authority"]
