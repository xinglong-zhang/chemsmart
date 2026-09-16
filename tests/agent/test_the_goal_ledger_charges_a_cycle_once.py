"""A cycle is charged once, whatever races to record it.

`GoalLedger.budgets` subtracts `engine_calls_consumed` from every
`run_recorded` row, so a second row for one cycle spends the grant twice.
The only thing standing between a goal and that double charge was
`GoalDriver.resume` re-reading the ledger and filtering out cycles it had
already seen -- a read-then-act with no mutual exclusion, on an
append-only file opened with a bare `open("a")`, while the runtime event
store beside it had held `flock` and an idempotency key for its own
appends all along.

Nothing had to race for this to matter: `chemsmart agent wake` is a
public command a human can run twice, and the plan's dependent wake job
makes a second wake an ordinary scheduler event rather than an accident.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.goal import GoalLedger


def _ledger(tmp_path: Path) -> GoalLedger:
    ledger = GoalLedger(tmp_path / "goal-dir")
    ledger.directory.mkdir(parents=True, exist_ok=True)
    return ledger


def _kinds(ledger: GoalLedger) -> list[str]:
    return [entry["kind"] for entry in ledger.entries()]


def test_one_cycle_recorded_twice_is_charged_once(tmp_path):
    ledger = _ledger(tmp_path)
    payload = {"cycle": 1, "engine_calls_consumed": 3}
    key = "run-recorded:g1:1:abc"

    ledger.append("run_recorded", payload, idempotency_key=key)
    ledger.append("run_recorded", payload, idempotency_key=key)

    assert _kinds(ledger) == ["run_recorded"], (
        "the same cycle was recorded twice, so its engine calls are "
        f"subtracted twice from the grant: {ledger.entries()}"
    )


def test_a_key_that_means_something_else_is_refused_not_merged(tmp_path):
    """Two different facts under one key is a contradiction, not a repeat."""

    ledger = _ledger(tmp_path)
    key = "run-recorded:g1:1:abc"
    ledger.append(
        "run_recorded",
        {"cycle": 1, "engine_calls_consumed": 3},
        idempotency_key=key,
    )
    with pytest.raises(ContractError, match="idempotency key"):
        ledger.append(
            "run_recorded",
            {"cycle": 1, "engine_calls_consumed": 9},
            idempotency_key=key,
        )
    assert len(ledger.entries()) == 1


def test_unkeyed_appends_are_untouched(tmp_path):
    """Most rows are a running account and repeat legitimately."""

    ledger = _ledger(tmp_path)
    ledger.append("observables_declared", {"n": 1})
    ledger.append("observables_declared", {"n": 1})
    assert _kinds(ledger) == [
        "observables_declared",
        "observables_declared",
    ]


def test_concurrent_processes_cannot_interleave_or_duplicate_a_row(
    tmp_path,
):
    """Real processes, real file descriptors, one row.

    Threads share a process and `flock` is per-process, so a thread test
    would pass against no lock at all. This forks interpreters.
    """

    import subprocess

    directory = tmp_path / "goal-dir"
    directory.mkdir(parents=True, exist_ok=True)
    program = (
        "import sys;"
        "sys.path.insert(0, %r);"
        "from chemsmart.agent.goal import GoalLedger;"
        "ledger = GoalLedger(%r);"
        "ledger.append('run_recorded', {'cycle': 1, "
        "'engine_calls_consumed': 3}, idempotency_key='k')"
    ) % (str(Path.cwd()), str(directory))

    workers = [
        subprocess.Popen(
            [sys.executable, "-c", program],
            env={**os.environ, "PYTHONPATH": str(Path.cwd())},
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        for _ in range(8)
    ]
    failures = []
    for worker in workers:
        _out, err = worker.communicate(timeout=120)
        if worker.returncode != 0:
            failures.append(err.decode()[-400:])

    rows = [
        line
        for line in (directory / "ledger.jsonl").read_text().splitlines()
        if line.strip()
    ]
    # Every line is whole: a torn write is the other thing a missing lock
    # produces, and it is invisible until something parses the file.
    for line in rows:
        json.loads(line)
    assert len(rows) == 1, (
        f"{len(rows)} rows for one cycle from 8 processes; "
        f"worker errors: {failures}"
    )


def test_the_grant_itself_is_spent_once(tmp_path):
    """The consumer, not just the row count.

    `budgets` is what the plan gate and the launch gate read. A second
    `run_recorded` row for one cycle subtracts its engine calls, its
    excursion calls and its engine seconds from the human's grant again,
    which is how a goal could run more work than was approved while every
    individual receipt looked correct.
    """

    from chemsmart.agent.goal import GOAL_SCHEMA_VERSION, GoalRecordV1

    ledger = _ledger(tmp_path)
    record = GoalRecordV1(
        schema_version=GOAL_SCHEMA_VERSION,
        goal_id="g1",
        task_spec_sha256="a" * 64,
        scientific_identity_sha256="b" * 64,
        conditions={},
        envelope={
            "max_engine_calls": 10,
            "max_excursion_calls": 2,
            "episode_wall_time_seconds": 1000.0,
        },
        max_revisions=3,
        granted_by="tester",
        initial_review_sha256="c" * 64,
        created_at="2026-09-16T00:00:00+00:00",
    )
    payload = {
        "cycle": 1,
        "engine_calls_consumed": 4,
        "excursion_calls_consumed": 1,
        "engine_wall_seconds": 250.0,
    }
    key = "run-recorded:g1:1"
    for _ in range(3):
        ledger.append("run_recorded", payload, idempotency_key=key)

    budgets = ledger.budgets(record)
    assert (
        budgets.engine_calls_remaining == 6
    ), "one cycle's engine calls were subtracted more than once"
    assert budgets.excursion_calls_remaining == 1
    assert budgets.wall_seconds_remaining == pytest.approx(750.0)
