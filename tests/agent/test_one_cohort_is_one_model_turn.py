"""One wake per cohort means one model turn, not one ledger row.

Three separate defects sat behind that sentence.

A keyed `run_recorded` deduplicates the *row* and then returns silently,
so two wakes that both passed `resume` before either wrote both carried
on into workspace recording, settlement, and -- for a repairable outcome
-- recovery and a fresh planning turn. The lock protected the row; it
never transferred ownership of the continuation. Two invocations of the
public `agent wake` command are enough; no dependent wake job is needed.

`run_dispatched` is keyed too, but the key is taken *after* `sbatch`. Two
execution attempts on one cycle therefore both reach the scheduler and
get different job ids; only then does the conflicting payload raise. The
second job already exists, and the ledger and the receipt sidecar can
name different jobs.

And a call spent at the launch fence but never receipted -- a controller
killed mid-engine -- was counted by the fence and not by the outcome, so
the goal's grant silently regained a call that had been spent.
"""

from __future__ import annotations

import json
from types import SimpleNamespace

from chemsmart.agent.goal import GOAL_SCHEMA_VERSION, GoalLedger, GoalRecordV1


def _ledger(tmp_path):
    ledger = GoalLedger(tmp_path / "g")
    ledger.directory.mkdir(parents=True, exist_ok=True)
    ledger.create(
        GoalRecordV1(
            schema_version=GOAL_SCHEMA_VERSION,
            goal_id="g1",
            task_spec_sha256="a" * 64,
            scientific_identity_sha256="b" * 64,
            conditions={},
            envelope={"max_engine_calls": 4},
            max_revisions=3,
            granted_by="tester",
            initial_review_sha256="c" * 64,
            created_at="2026-09-16T00:00:00+00:00",
        )
    )
    return ledger


def test_the_ledger_tells_the_caller_whether_it_owns_the_turn(tmp_path):
    """A silent no-op is why both wakes continued."""

    ledger = _ledger(tmp_path)
    payload = {"cycle": 1, "engine_calls_consumed": 2}
    key = "run-recorded:g1:1"
    assert ledger.append("run_recorded", payload, idempotency_key=key) is True
    assert ledger.append("run_recorded", payload, idempotency_key=key) is False


def test_an_unkeyed_append_always_owns_its_row(tmp_path):
    ledger = _ledger(tmp_path)
    assert ledger.append("observables_declared", {"n": 1}) is True
    assert ledger.append("observables_declared", {"n": 1}) is True


def _park(tmp_path, *, dispatch_run, turns=None):
    """A goal parked on a scheduler job, ready for a wake to resume.

    ``turns`` collects one entry per planning session, which is what a
    model turn costs: a provider call over this cycle's evidence.
    """

    from chemsmart.agent.driver import GoalDriver

    from .test_the_goal_loop_recovers_or_returns import (
        _envelope_file,
        _planning_session,
        _review_payload,
    )

    workspace = tmp_path / "ws"
    workspace.mkdir(exist_ok=True)
    driver = GoalDriver(
        task="the goal task",
        workspace=workspace,
        execution_envelope_file=_envelope_file(tmp_path),
        goal_id="goal-w1",
        granted_by="claude-owner-delegated-reviewer",
        plan_session=lambda **kw: (
            turns.append(kw) if turns is not None else None,
            _planning_session("live-1", review=_review_payload())(
                workspace, kw
            ),
        )[1],
        resolve_review=lambda **_kw: ("d" * 64, tmp_path / "bundle.json"),
        dispatch_run=dispatch_run,
        dispatch="scheduler",
        server="canned-slurm",
    )
    return workspace, driver


def _dispatch_recording(calls):
    def dispatch_run(**kwargs):
        calls.append(kwargs["cycle"])
        return SimpleNamespace(
            scheduler="SLURM",
            job_id=f"19{len(calls)}",
            submitted_at="2026-09-01T00:00:00+00:00",
            submit_script=str(kwargs["run_directory"] / "sub.sh"),
        )

    return dispatch_run


def test_two_wakes_on_one_cycle_are_one_model_turn(tmp_path):
    """Drive the public command twice, as a duplicate notification does.

    Both wakes pass ``resume`` before either writes ``run_recorded``,
    which is the whole race: the keyed append deduplicated the *row* and
    returned silently, so the loser was never told it had lost and
    carried on into workspace recording, settlement and -- the cost that
    matters -- a second planning turn against the provider for one
    cycle's evidence. Counting ledger rows cannot see this, because the
    rows are keyed; counting model turns can.
    """

    from click.testing import CliRunner

    from chemsmart.agent.dispatch import EXECUTION_RESULT_FILE
    from chemsmart.cli.agent import agent

    from .test_the_goal_loop_recovers_or_returns import _engine_stream

    calls: list[int] = []
    turns: list[dict] = []
    workspace, driver = _park(
        tmp_path, dispatch_run=_dispatch_recording(calls), turns=turns
    )
    assert driver.run().settlement == "parked"
    planned_once = len(turns)

    # A repairable ending, so a continuing process would plan again --
    # which is exactly the second provider call the barrier exists to
    # prevent.
    _engine_stream(tmp_path, driver.run_directory, failed=True)
    (driver.run_directory / EXECUTION_RESULT_FILE).write_text(
        json.dumps({"status": "partial", "analysis_status": "partial"}),
        encoding="utf-8",
    )

    # The woken driver is rebuilt by `resume`, which supplies the real
    # planning session, so the seam a model turn crosses is the module
    # default rather than this driver's own hook.
    import chemsmart.agent.driver as driver_module

    from .test_the_goal_loop_recovers_or_returns import (
        _planning_session,
        _review_payload,
    )

    def counted_plan_session(**kw):
        turns.append(kw)
        return _planning_session("live-2", review=_review_payload())(
            workspace, kw
        )

    monkey = driver_module._default_plan_session
    driver_module._default_plan_session = counted_plan_session
    import chemsmart.agent.driver as _dm

    _orig_init = _dm.GoalDriver.__init__

    def _init(self, *a, **kw):
        kw.setdefault("plan_session", counted_plan_session)
        return _orig_init(self, *a, **kw)

    _dm.GoalDriver.__init__ = _init

    runner = CliRunner()
    released: list[str] = []
    original = type(driver.ledger).append

    # The second wake is released from inside the first one's outcome
    # phase, after it has derived and before it has written: the only
    # window the race lives in. Neither side is a mock of the thing
    # under test -- both are the real `agent wake`.
    def append_then_race(self, kind, payload, *, idempotency_key=None):
        if kind == "run_recorded" and not released:
            released.append("released")
            type(self).append = original
            try:
                second = runner.invoke(
                    agent,
                    [
                        "wake",
                        "--workspace",
                        str(workspace),
                        "--goal",
                        "goal-w1",
                    ],
                )
                released.append(second.output)
            finally:
                type(self).append = append_then_race
        return original(self, kind, payload, idempotency_key=idempotency_key)

    type(driver.ledger).append = append_then_race
    try:
        first = runner.invoke(
            agent,
            ["wake", "--workspace", str(workspace), "--goal", "goal-w1"],
        )
    finally:
        type(driver.ledger).append = original
        driver_module._default_plan_session = monkey
        _dm.GoalDriver.__init__ = _orig_init

    assert released, "the race never fired; the witness proves nothing"
    rows = [
        json.loads(line)
        for line in driver.ledger.ledger_path.read_text().splitlines()
        if line.strip()
    ]
    assert len([r for r in rows if r["kind"] == "run_recorded"]) == 1

    # One cohort, one wake, one turn. Two wakes that both got past
    # `resume` planned twice against one cycle's evidence.
    assert len(turns) == planned_once + 1, (
        f"one cycle's evidence bought {len(turns) - planned_once} model "
        "turns: the wake that lost the keyed append was not told, and "
        "reasoned from the same evidence beside the winner"
    )
    assert first.exit_code == 0 or "is settled" in first.output, first.output


def test_a_cycle_reaches_the_scheduler_once(tmp_path):
    """The irreversible act is ``sbatch``, so the claim precedes it.

    ``run_dispatched`` is keyed, but the key was taken *after* the
    submission: two attempts on one cycle both reached the scheduler,
    got different job ids, and only then did the conflicting payload
    raise -- with a second job already queued and the ledger and the
    receipt sidecar naming different jobs.
    """

    calls: list[int] = []
    _workspace, driver = _park(
        tmp_path, dispatch_run=_dispatch_recording(calls)
    )
    assert driver.run().settlement == "parked"
    assert calls == [1]

    # A retried execute on the same cycle: a controller that crashed
    # after submitting, or a wake racing the tail.
    driver.phase = "execute"
    driver._execute()

    assert calls == [1], (
        "the cycle reached the scheduler twice; the first job is now "
        "orphaned and two jobs write one run directory"
    )
    rows = [
        json.loads(line)
        for line in driver.ledger.ledger_path.read_text().splitlines()
        if line.strip()
    ]
    dispatched = [row for row in rows if row["kind"] == "run_dispatched"]
    assert len(dispatched) == 1
    assert dispatched[0]["payload"]["job_id"] == "191"


def test_a_claimed_cycle_with_no_job_is_ambiguous_not_resubmitted(tmp_path):
    """Claimed and never dispatched is the one state nobody can resolve.

    A controller killed between the claim and the submission leaves a
    cycle that may or may not have a job. Submitting again is how one
    grant buys two jobs; the honest answer is the word CHEMSMART already
    uses for an interruption it cannot reconcile.
    """

    calls: list[int] = []
    _workspace, driver = _park(
        tmp_path, dispatch_run=_dispatch_recording(calls)
    )
    assert driver.run().settlement == "parked"

    # Erase the dispatch row, keeping the claim: the crash window.
    rows = [
        line
        for line in driver.ledger.ledger_path.read_text().splitlines()
        if line.strip() and json.loads(line)["kind"] != "run_dispatched"
    ]
    driver.ledger.ledger_path.write_text("\n".join(rows) + "\n")

    driver.phase = "execute"
    driver._execute()

    assert calls == [1], "a cycle that may already hold a job was resubmitted"
    assert driver.result is not None
    reasons = " ".join(driver.result.reasons)
    assert "ambiguous" in reasons.lower()


def test_claiming_a_cycle_twice_is_refused(tmp_path):
    ledger = _ledger(tmp_path)
    key = "run-dispatch-claimed:g1:1"
    assert (
        ledger.append(
            "run_dispatch_claimed", {"cycle": 1}, idempotency_key=key
        )
        is True
    )
    assert (
        ledger.append(
            "run_dispatch_claimed", {"cycle": 1}, idempotency_key=key
        )
        is False
    )


def test_a_reserved_call_that_never_receipted_is_still_charged(tmp_path):
    """A controller killed mid-engine spent the call.

    The launch fence counts a call when it is *taken*
    (``engine_calls_spent``); the run outcome counted only receipts. So
    the goal's grant silently regained a call the engine had already
    begun burning, and the next cycle was handed budget that omitted it.
    """

    from chemsmart.agent.runtime.event_store import (
        RuntimeEventStore,
        engine_calls_spent,
    )
    from chemsmart.agent.terminal_states import (
        derive_run_outcome,
        read_run_events,
    )

    from .test_runtime_v2_launch_fence import _reserve

    store = RuntimeEventStore(tmp_path / "events.jsonl", session_id="s")
    _reserve(store, tmp_path)

    events = read_run_events(tmp_path / "events.jsonl")
    fence = engine_calls_spent(events)
    outcome = derive_run_outcome(events)

    assert fence == 1, "the fence charges a taken call"
    assert outcome.engine_calls_consumed == fence, (
        "the two authorities disagree about one node: the fence has "
        f"spent {fence} and the run outcome reports "
        f"{outcome.engine_calls_consumed}, so a crash between the "
        "reservation and the receipt returns a spent call to the grant"
    )
