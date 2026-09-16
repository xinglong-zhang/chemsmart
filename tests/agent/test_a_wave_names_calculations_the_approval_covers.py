"""A wave is a subset of the approved partition, or the host says so.

`_select_execution_wave` judges a proposal against the *planned* draft's
ready frontier -- which is what exists when the model selects, before any
review is resolved. So a node the review later retained as non-executable,
or an id the model mistyped into the plan and then named, can reach the
manifest.

The failure is silent and total. Each element's `cohort_frontier(ready,
(my_node,))` is empty, so every element runs nothing; `derive_run_outcome`
raises "found 0"; the driver's own analysis-only branch records the cycle
with zero engine calls and settles `complete`. An approved partition that
never ran settles as a delivery.

This is a host-internal consistency check between two things the host
owns -- the approved partition and the manifest it is about to write --
not a refusal shown to the Agent. What it drops, it records.
"""

from __future__ import annotations

import json
from types import SimpleNamespace



def _driver(tmp_path, *, non_executable=(), wave=()):
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
    reviews = driver.goal_dir / "reviews"
    reviews.mkdir(parents=True, exist_ok=True)
    (reviews / "cycle-1.json").write_text(
        json.dumps(
            {
                "workflow_execution_review": {
                    "non_executable_node_ids": list(non_executable),
                }
            }
        ),
        encoding="utf-8",
    )
    driver.session = SimpleNamespace(selected_execution_wave=tuple(wave))
    return driver


def test_a_wave_of_approved_calculations_passes_through(tmp_path):
    driver = _driver(tmp_path, wave=("a1", "a2"))
    assert driver._dispatchable_wave() == ("a1", "a2")


def test_a_retained_stage_is_dropped_from_the_wave(tmp_path):
    """Displayed as intent, never approved, and never launchable."""

    driver = _driver(
        tmp_path, non_executable=("irc-blocked",), wave=("a1", "irc-blocked")
    )
    assert driver._dispatchable_wave() == ("a1",)


def test_a_wave_of_nothing_approved_is_not_a_cohort(tmp_path):
    """Empty is the single-job path, not an empty array.

    An array of zero elements dispatches nothing and the goal settles a
    partition that never ran as a complete delivery.
    """

    driver = _driver(
        tmp_path, non_executable=("irc-blocked",), wave=("irc-blocked",)
    )
    assert driver._dispatchable_wave() == ()


def test_a_session_that_selected_nothing_stays_a_single_job(tmp_path):
    driver = _driver(tmp_path, wave=())
    assert driver._dispatchable_wave() == ()
