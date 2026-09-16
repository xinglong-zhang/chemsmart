"""A wave that cannot progress ends typed, and a goal settles once.

Two ways a cohort disappeared without a word.

The poll had no bound. `max_polls` exists and is test-only: the one
production caller passes `poll_seconds` and nothing else, so a wake job
stuck `PENDING` with `Reason=DependencyNeverSatisfied` -- which is a
`Reason`, not a `JobState`, and `PENDING` is not terminal -- polls every
thirty seconds forever. Nothing in the tree reads `Reason` at all, and
with `--mail-type=END,FAIL` on a job that never runs, no human is told
either.

And `goal_settled` -- the one row that ends a goal, and the strongest
charge in the ledger -- was the only charging row with no idempotency
key. Its sole protection was `resume`'s unlocked read-then-act, which is
exactly the pattern the keyed append was introduced to replace.
"""

from __future__ import annotations

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.goal import GOAL_SCHEMA_VERSION, GoalLedger, GoalRecordV1


def _record():
    return GoalRecordV1(
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


def _ledger(tmp_path):
    ledger = GoalLedger(tmp_path / "g")
    ledger.directory.mkdir(parents=True, exist_ok=True)
    ledger.create(_record())
    return ledger


def test_a_goal_settles_once_however_many_wakes_arrive(tmp_path):
    ledger = _ledger(tmp_path)
    for _ in range(3):
        ledger.settle("achieved", reasons=("done",))
    settled = [e for e in ledger.entries() if e["kind"] == "goal_settled"]
    assert len(settled) == 1, (
        f"the goal ended {len(settled)} times; a second settlement can "
        "carry a different word than the first"
    )


def test_two_wakes_cannot_settle_a_goal_two_different_ways(tmp_path):
    """The dangerous half: not a duplicate, a contradiction.

    Two elements' tails racing could both pass `resume`'s check and
    append different words, and a later reader takes whichever it scans
    first.
    """

    ledger = _ledger(tmp_path)
    ledger.settle("achieved", reasons=("done",))
    with pytest.raises(ContractError, match="idempotency key"):
        ledger.settle("returned_to_human", reasons=("something else",))
    settled = [e for e in ledger.entries() if e["kind"] == "goal_settled"]
    assert len(settled) == 1
    assert settled[0]["payload"]["state"] == "achieved"


def test_a_dependency_that_can_never_be_satisfied_is_read(tmp_path):
    """`Reason=DependencyNeverSatisfied` on a PENDING job is the shape
    Slurm gives a wake that will never run."""

    from chemsmart.settings.probe.scheduler_job import parse_scontrol_job

    stdout = (
        "JobId=4242 JobName=goal-g1-cycle-1-wake\n"
        "   JobState=PENDING Reason=DependencyNeverSatisfied "
        "Dependency=afterany:4241(failed)\n"
        "   RunTime=00:00:00 TimeLimit=00:10:00\n"
    )
    state = parse_scontrol_job(0, stdout, "", job_id="4242")
    assert state.known
    assert state.state == "PENDING"
    assert not state.terminal
    assert state.reason == "DependencyNeverSatisfied"
    assert state.dependency_unsatisfiable, (
        "the host cannot see that this wake will never run, so the goal "
        "parks forever with nothing to read"
    )


def test_an_ordinary_pending_reason_is_not_unsatisfiable(tmp_path):
    from chemsmart.settings.probe.scheduler_job import parse_scontrol_job

    stdout = (
        "JobId=4242 JobName=w\n"
        "   JobState=PENDING Reason=Resources Dependency=(null)\n"
    )
    state = parse_scontrol_job(0, stdout, "", job_id="4242")
    assert state.reason == "Resources"
    assert not state.dependency_unsatisfiable


def test_the_wait_gives_up_rather_than_polling_forever(tmp_path):
    """A deadline, not an unbounded loop."""

    import json
    import subprocess

    from chemsmart.agent.dispatch import (
        DISPATCH_RECEIPT_FILE,
        DispatchReceiptV1,
        wait_for_dispatched_run,
    )

    run_directory = tmp_path / "run"
    run_directory.mkdir()
    (run_directory / DISPATCH_RECEIPT_FILE).write_text(
        json.dumps(
            DispatchReceiptV1(
                scheduler="SLURM",
                job_id="4242",
                submitted_at="2026-09-16T00:00:00+00:00",
                submit_command="sbatch x.sh",
                submit_script=str(run_directory / "x.sh"),
                run_directory=str(run_directory),
                approval_file=str(tmp_path / "bundle.json"),
                goal_id="g1",
                cycle=1,
                wake_command="python -m chemsmart agent wake",
            ).public_record()
        ),
        encoding="utf-8",
    )

    def stuck(argv, **_kw):
        return subprocess.CompletedProcess(
            argv,
            0,
            stdout=(
                "JobId=4242\n   JobState=PENDING "
                "Reason=DependencyNeverSatisfied\n"
            ),
        )

    # The sleep gives up, so an unbounded loop fails the test instead of
    # hanging it. Without this the production behaviour is an infinite
    # spin, which is the defect.
    slept: list[float] = []

    def bounded_sleep(seconds):
        slept.append(seconds)
        if len(slept) > 3:
            raise AssertionError(
                "wait_for_dispatched_run polled past its own deadline on a "
                "dependency that can never be satisfied"
            )

    why = wait_for_dispatched_run(
        run_directory, poll_seconds=5, runner=stuck, sleep=bounded_sleep
    )
    assert "never" in why.lower() or "satisfied" in why.lower(), why
    assert len(slept) <= 1, (
        f"polled {len(slept)} times on a dependency that can never be "
        "satisfied; the host can read that on the first look"
    )
