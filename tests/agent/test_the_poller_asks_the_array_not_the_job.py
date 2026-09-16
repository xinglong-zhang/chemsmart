"""Polling a cohort asks about its elements, not about a scalar job.

`scontrol show job <N>` on an array prints one record per element, whose
`JobId=` is `N_0`, `N_1`, ... and never the bare `N` the reader compares
against; the field harvest is first-occurrence-wins, so a three-element
array collapses into whichever element scontrol happened to print first.
A wave of three where one element is still RUNNING and the first printed
element has COMPLETED reads as "job N COMPLETED" -- the poll returns, the
goal re-enters, and the Agent reasons over a wave that is still running.

`squeue_array_command`, `parse_squeue_array` and `SchedulerArrayStateV1`
were written for exactly this and read by nothing.
"""

from __future__ import annotations

import json
import subprocess

from chemsmart.agent.cohort import build_cohort_manifest
from chemsmart.agent.dispatch import (
    DISPATCH_RECEIPT_FILE,
    DispatchReceiptV1,
    wait_for_dispatched_run,
)

_SQUEUE_RUNNING = (
    "4242_0|4242|0|COMPLETED|00:05:00|2026-09-16T12:00:00|"
    "2026-09-16T12:00:05|2026-09-16T12:05:05\n"
    "4242_1|4242|1|RUNNING|00:02:00|2026-09-16T12:00:00|"
    "2026-09-16T12:00:05|N/A\n"
    "4242_2|4242|2|RUNNING|00:02:00|2026-09-16T12:00:00|"
    "2026-09-16T12:00:05|N/A\n"
)
_SQUEUE_DONE = (
    "4242_0|4242|0|COMPLETED|00:05:00|2026-09-16T12:00:00|"
    "2026-09-16T12:00:05|2026-09-16T12:05:05\n"
    "4242_1|4242|1|CANCELLED|00:03:00|2026-09-16T12:00:00|"
    "2026-09-16T12:00:05|2026-09-16T12:03:05\n"
    "4242_2|4242|2|FAILED|00:04:00|2026-09-16T12:00:00|"
    "2026-09-16T12:00:05|2026-09-16T12:04:05\n"
)


def _cohort_run(tmp_path, *, wake_job="4243"):
    run_directory = tmp_path / "run"
    run_directory.mkdir(parents=True, exist_ok=True)
    build_cohort_manifest(
        goal_id="g1",
        cycle=1,
        bundle_sha256="e" * 64,
        node_ids=("a1", "a2", "a3"),
        max_concurrent_tasks=4,
        created_at="2026-09-16T00:00:00+00:00",
    ).write(run_directory)
    (run_directory / DISPATCH_RECEIPT_FILE).write_text(
        json.dumps(
            DispatchReceiptV1(
                scheduler="SLURM",
                job_id="4242",
                submitted_at="2026-09-16T12:00:00+00:00",
                submit_command="sbatch x.sh",
                submit_script="x.sh",
                run_directory=str(run_directory),
                approval_file=str(tmp_path / "bundle.json"),
                goal_id="g1",
                cycle=1,
                wake_command="wake",
                wake_job_id=wake_job,
                cohort_node_ids=("a1", "a2", "a3"),
            ).public_record()
        ),
        encoding="utf-8",
    )
    return run_directory


def _runner(outputs):
    seen = []

    def runner(argv, **_kw):
        seen.append(list(argv))
        return subprocess.CompletedProcess(argv, 0, stdout=outputs.pop(0))

    runner.seen = seen
    return runner


def test_a_running_wave_is_not_reported_finished(tmp_path):
    """Two members still running; the first printed element completed."""

    run_directory = _cohort_run(tmp_path)
    runner = _runner([_SQUEUE_RUNNING, _SQUEUE_DONE])
    why = wait_for_dispatched_run(
        run_directory, poll_seconds=0, runner=runner, sleep=lambda _s: None
    )
    assert "a1" not in why or "still" in why
    assert why.startswith("array 4242"), why
    assert "3" in why or "ended" in why, why
    # It asked squeue about the array, expanded, not scontrol about a
    # scalar job id the array's records never carry.
    first = runner.seen[0]
    assert first[0] == "squeue", first
    assert "-r" in first, ("compressed rows hide membership", first)
    assert any("4242" in str(part) for part in first)


def test_a_wave_whose_members_all_ended_ends_the_wait(tmp_path):
    """Terminality, not success: cancelled and failed members ended."""

    run_directory = _cohort_run(tmp_path)
    runner = _runner([_SQUEUE_DONE])
    why = wait_for_dispatched_run(
        run_directory, poll_seconds=0, runner=runner, sleep=lambda _s: None
    )
    assert why == "array 4242: every element ended", why


def test_an_array_the_scheduler_forgot_is_not_a_finished_array(tmp_path):
    """'I cannot find it' is not 'every member finished'."""

    run_directory = _cohort_run(tmp_path)
    runner = _runner([""])
    why = wait_for_dispatched_run(
        run_directory, poll_seconds=0, runner=runner, sleep=lambda _s: None
    )
    assert "no longer knows" in why, why


def test_the_single_job_path_still_asks_scontrol(tmp_path):
    """A run of one is not a degenerate array."""

    run_directory = tmp_path / "solo"
    run_directory.mkdir()
    (run_directory / DISPATCH_RECEIPT_FILE).write_text(
        json.dumps(
            DispatchReceiptV1(
                scheduler="SLURM",
                job_id="191",
                submitted_at="2026-09-16T12:00:00+00:00",
                submit_command="sbatch x.sh",
                submit_script="x.sh",
                run_directory=str(run_directory),
                approval_file=str(tmp_path / "bundle.json"),
                goal_id="g1",
                cycle=1,
                wake_command="wake",
            ).public_record()
        ),
        encoding="utf-8",
    )
    runner = _runner(["JobId=191 JobState=COMPLETED RunTime=00:00:02\n"])
    why = wait_for_dispatched_run(
        run_directory, poll_seconds=0, runner=runner, sleep=lambda _s: None
    )
    assert why == "job 191 COMPLETED", why
    assert runner.seen[0][0] == "scontrol", runner.seen[0]
