"""A woken cycle must settle on the stream *it* planned in.

`_project_before_settling` resolves a missing `events_path` through
`_session_events_path(None, workspace)`, whose fallback globs
`live-*/events.jsonl` across the whole workspace and takes the newest.
That is right for the case it was written for -- a planning session that
*raised* seconds ago in this very process, where the guess cannot be
wrong -- and it was silently inherited by the woken cycle, where the
window is minutes to hours by design, because the whole point of
`--dispatch scheduler` is that the scientist does something else while
the array runs.

A session id carries a task-spec digest, not a goal id, so streams from
different goals in one workspace are indistinguishable by name and there
is no goal filter. So: goal A's cycle 2 dispatches, the wake fires forty
minutes later, the scientist has since started goal B in the same
workspace, and goal A's settlement resolves **goal B's stream** --
projecting goal B's declarations and claims into goal A's ledger. One
workspace with two goals is the ordinary way a scientist works on one
system.

The driver already holds the exact path at plan time. Recording the run
id it planned in lets the wake resolve it instead of guessing, and a
guess is only right when it cannot be wrong.
"""

from __future__ import annotations

import json




def _workspace_with_two_streams(tmp_path):
    runs = tmp_path / "ws" / ".chemsmart-agent" / "runs"
    for name in ("live-20260916T100000000000Z-aaaa-1111",
                 "live-20260916T114000000000Z-bbbb-2222"):
        (runs / name).mkdir(parents=True)
        (runs / name / "events.jsonl").write_text(
            json.dumps({"kind": "session_started", "payload": {}}) + "\n",
            encoding="utf-8",
        )
    return tmp_path / "ws"


def test_the_fallback_still_finds_the_newest_when_nothing_was_recorded(
    tmp_path,
):
    """Unchanged for the case it was written for."""

    from chemsmart.agent.driver import _session_events_path

    workspace = _workspace_with_two_streams(tmp_path)
    found = _session_events_path(None, workspace)
    assert found.parent.name.endswith("bbbb-2222")


def test_a_cycle_settles_on_the_stream_it_recorded(tmp_path):
    """The earlier stream is this cycle's; the newer belongs to another
    goal the scientist started while the array ran."""

    from chemsmart.agent.driver import GoalDriver

    from .test_the_goal_loop_recovers_or_returns import _envelope_file

    workspace = _workspace_with_two_streams(tmp_path)
    driver = GoalDriver(
        task="t",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path),
        goal_id="goal-a",
        granted_by="tester",
    )
    driver.cycles = 2
    driver.ledger.directory.mkdir(parents=True, exist_ok=True)
    driver.ledger.append(
        "session_stream_recorded",
        {
            "cycle": 2,
            "run_id": "live-20260916T100000000000Z-aaaa-1111",
        },
    )

    resolved = driver._planned_events_path()
    assert resolved is not None
    assert resolved.parent.name.endswith("aaaa-1111"), (
        "the woken cycle settled on whatever stream was newest in the "
        "workspace, which is another goal's whenever the scientist "
        "started one while the array ran"
    )


def test_a_cycle_that_recorded_nothing_keeps_the_old_behaviour(tmp_path):
    from chemsmart.agent.driver import GoalDriver

    from .test_the_goal_loop_recovers_or_returns import _envelope_file

    workspace = _workspace_with_two_streams(tmp_path)
    driver = GoalDriver(
        task="t",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path),
        goal_id="goal-a",
        granted_by="tester",
    )
    driver.cycles = 2
    assert driver._planned_events_path() is None


def test_a_recorded_stream_that_is_gone_is_not_substituted(tmp_path):
    """A named stream that vanished is absence, not licence to guess."""

    from chemsmart.agent.driver import GoalDriver

    from .test_the_goal_loop_recovers_or_returns import _envelope_file

    workspace = _workspace_with_two_streams(tmp_path)
    driver = GoalDriver(
        task="t",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path),
        goal_id="goal-a",
        granted_by="tester",
    )
    driver.cycles = 2
    driver.ledger.directory.mkdir(parents=True, exist_ok=True)
    driver.ledger.append(
        "session_stream_recorded",
        {"cycle": 2, "run_id": "live-20260916T090000000000Z-cccc-3333"},
    )
    assert driver._planned_events_path() is None


def test_the_plan_records_the_stream_it_planned_in(tmp_path):
    """Composition: what `_plan` writes is what the wake reads back.

    Cycle 1 produces this row before `goal_created` exists, so the row
    waits and the flush after `ledger.create` enters it -- the same
    deferral every other cycle-ending row uses. Both halves are driven:
    deferred with no goal, appended once there is one.
    """

    from types import SimpleNamespace

    from chemsmart.agent.driver import GoalDriver

    from .test_the_goal_loop_recovers_or_returns import _envelope_file

    workspace = _workspace_with_two_streams(tmp_path)
    driver = GoalDriver(
        task="t",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path),
        goal_id="goal-a",
        granted_by="tester",
    )
    driver.cycles = 2
    driver.ledger.directory.mkdir(parents=True, exist_ok=True)

    # No goal record yet: the row waits rather than being dropped.
    driver._record_session_stream(
        SimpleNamespace(run_id="live-20260916T100000000000Z-aaaa-1111")
    )
    assert driver._planned_events_path() is None
    assert driver._pending_ledger_rows

    # With a goal, it is entered and the wake reads it back.
    driver.goal = SimpleNamespace(actor="tester")
    driver._record_session_stream(
        SimpleNamespace(run_id="live-20260916T100000000000Z-aaaa-1111")
    )
    resolved = driver._planned_events_path()
    assert resolved is not None
    assert resolved.parent.name.endswith("aaaa-1111")


def test_a_named_stream_that_is_gone_never_reaches_the_glob(tmp_path):
    """Two states, one `None`, and the caller could not tell them apart.

    `_planned_events_path` returns `None` both when this cycle recorded
    no stream -- correct, and what every goal written before this needs
    -- and when it recorded one whose directory has since gone. The
    second fell through to the same workspace-wide glob the fix exists
    to prevent, while the function's own docstring said it did not:
    "absence rather than licence to substitute another goal's". The
    sentence was right and the return type could not carry it.

    Reachable without anything unusual: a `live-*` directory removed by
    workspace tidying, a scratch-retention policy, or a workspace
    rsynced without the session streams -- which look like logs and are
    the first thing someone drops.
    """

    from chemsmart.agent.driver import GoalDriver

    from .test_the_goal_loop_recovers_or_returns import _envelope_file

    workspace = _workspace_with_two_streams(tmp_path)
    driver = GoalDriver(
        task="t",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path),
        goal_id="goal-a",
        granted_by="tester",
    )
    driver.cycles = 2
    driver.ledger.directory.mkdir(parents=True, exist_ok=True)
    driver.ledger.append(
        "session_stream_recorded",
        {"cycle": 2, "run_id": "live-20260916T090000000000Z-cccc-3333"},
    )

    assert driver._named_its_own_stream() is True
    assert driver._planned_events_path() is None

    driver._project_before_settling()
    assert driver.events_path is None, (
        "a goal that named its stream fell back to whatever was newest "
        "in the workspace, which is the substitution this exists to stop"
    )


def test_a_cycle_that_named_nothing_still_uses_the_fallback(tmp_path):
    """Unchanged for every goal written before the row existed."""

    from chemsmart.agent.driver import GoalDriver

    from .test_the_goal_loop_recovers_or_returns import _envelope_file

    workspace = _workspace_with_two_streams(tmp_path)
    driver = GoalDriver(
        task="t",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path),
        goal_id="goal-a",
        granted_by="tester",
    )
    driver.cycles = 2
    assert driver._named_its_own_stream() is False
    driver._project_before_settling()
    assert driver.events_path is not None
    assert driver.events_path.parent.name.endswith("bbbb-2222")
