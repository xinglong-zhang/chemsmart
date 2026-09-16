"""An array's state is a multiset, and one scalar cannot carry it.

`SchedulerJobStateV1` has one `state`, one `exit_code`, one `end_time`.
An array job has one of each *per element*, and the barrier asks a
question no scalar can answer: has **every** member reached a terminal
state?

Measured on CUHK (job 2135107, a two-element array at %1), verbatim:

    squeue -j 2135107 -o '%i|%T|...'
        2135107_1|PENDING|...
        2135107_0|RUNNING|...

so the existing filter, which compares `%i` against the bare array id,
matches nothing. And `scontrol show job 2135107` opens with

    JobId=2135107 ArrayJobId=2135107 ArrayTaskId=1 ... JobState=PENDING

while the element that already started appears further down under its
own `JobId=2135108`. The scontrol field harvest is first-occurrence-wins
across the whole output, so it does not fail loudly -- it reports the
*pending aggregate's* state as though it were the array's, hiding the
running element entirely. squeue with `%F` is the honest authority: one
line per element, in a pinned format, joined on ArrayJobId.
"""

from __future__ import annotations

from chemsmart.settings.probe.scheduler_job import (
    parse_squeue_array,
    squeue_array_command,
)

# Captured verbatim from CUHK job 2135107.
_RUNNING = (
    "2135107_1|2135107|1|PENDING|0:00|2026-09-16T17:46:50|N/A|N/A\n"
    "2135107_0|2135107|0|RUNNING|0:11|2026-09-16T17:46:50|"
    "2026-09-16T17:46:51|2026-09-16T17:48:51\n"
)
_DONE = (
    "2135107_1|2135107|1|COMPLETED|0:20|2026-09-16T17:46:50|"
    "2026-09-16T17:47:10|2026-09-16T17:47:30\n"
    "2135107_0|2135107|0|FAILED|0:11|2026-09-16T17:46:50|"
    "2026-09-16T17:46:51|2026-09-16T17:47:02\n"
)


def test_the_command_asks_for_the_array_join_key():
    command = squeue_array_command("2135107")
    assert "%F" in " ".join(command), (
        "without ArrayJobId there is nothing to join an element to its "
        "cohort by"
    )
    assert "2135107" in command


def test_a_cohort_mid_flight_is_not_terminal():
    cohort = parse_squeue_array(0, _RUNNING, "", array_job_id="2135107")
    assert cohort.known
    assert len(cohort.elements) == 2
    assert {e.task_id: e.state for e in cohort.elements} == {
        "0": "RUNNING",
        "1": "PENDING",
    }
    assert not cohort.terminal, (
        "a cohort with a running element was called finished, which is "
        "the barrier firing over a half-done wave"
    )


def test_the_barrier_is_terminality_not_success():
    """One element failed and one completed: the wave is over, and the
    failure is evidence the Agent must see, not a reason to wait."""

    cohort = parse_squeue_array(0, _DONE, "", array_job_id="2135107")
    assert cohort.terminal
    assert sorted(e.state for e in cohort.elements) == ["COMPLETED", "FAILED"]
    assert cohort.unfinished == ()


def test_another_array_s_elements_are_not_this_cohort():
    mixed = _DONE + "2199999_0|2199999|0|RUNNING|0:01|x|y|z\n"
    cohort = parse_squeue_array(0, mixed, "", array_job_id="2135107")
    assert len(cohort.elements) == 2
    assert cohort.terminal


def test_an_array_the_scheduler_has_forgotten_is_not_terminal():
    """MinJobAge purges a finished array, and 'I cannot see it' is not
    'every member finished'. Saying otherwise would wake the model over
    a cohort whose evidence the host never checked."""

    cohort = parse_squeue_array(0, "", "", array_job_id="2135107")
    assert not cohort.known
    assert not cohort.terminal
    assert cohort.elements == ()


def test_a_failed_query_is_not_terminal_either():
    cohort = parse_squeue_array(
        1, "", "slurm_load_jobs error", array_job_id="1"
    )
    assert not cohort.known
    assert not cohort.terminal


def test_the_query_expands_compressed_array_rows():
    """squeue compresses pending elements by default.

    `squeue -j 4242` reads `4242_[1-6%4]` for six pending calculations,
    so a barrier asking "is every member terminal" would count one row
    for six and could report a range string as an unfinished task id.
    """

    assert "-r" in squeue_array_command("4242")


def test_a_compressed_row_makes_the_cohort_unknown_not_undercounted():
    """A site that returns a range anyway must not be read as one
    element: unknown is honest, six-counted-as-one is not."""

    compressed = "4242_[1-6%4]|4242|1-6%4|PENDING|0:00|x|N/A|N/A\n"
    cohort = parse_squeue_array(0, compressed, "", array_job_id="4242")
    assert not cohort.known
    assert not cohort.terminal
    assert cohort.unfinished == ()
