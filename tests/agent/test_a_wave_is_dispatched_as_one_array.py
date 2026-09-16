"""A wave reaches the scheduler as one throttled array and one wake.

The dispatcher submitted one job per cycle that walked every approved
node serially. A wave is N independent calculations the Agent wants to
see together, so it is one array of N elements, `%N`-throttled by the
host, plus exactly one job whose only work is to wake the goal once the
array is over.

The wake is a separate job with `--dependency=afterany`, not a tail on
each element: N tails would wake the model N times, and `afterany` fires
whatever the elements' exit statuses are, because the barrier is
terminality and not success. `--kill-on-invalid-dep=yes` turns a
dependency that can never be satisfied into a CANCELLED job -- a state
this tree already classifies as terminal -- instead of a PENDING job
nobody can see.
"""

from __future__ import annotations

from types import SimpleNamespace

from chemsmart.agent.dispatch import (
    build_cohort_dispatch_script,
    build_wake_dispatch_script,
)
from chemsmart.settings.server import Server


def _server(**kw):
    base = dict(
        SCHEDULER="SLURM",
        SUBMIT_COMMAND="sbatch",
        NUM_CORES=64,
        MEM_GB=160,
        NUM_HOURS=24,
        NUM_GPUS=0,
        QUEUE_NAME="chpc",
    )
    base.update(kw)
    return Server("canned", **base)


def _submitter(server, label="goal-g1-cycle-1"):
    job = SimpleNamespace(label=label, folder="/tmp", PROGRAM=None)
    return server.get_submitter(job)


def _script(tmp_path, size=3, server=None):
    server = server or _server()
    return build_cohort_dispatch_script(
        submitter=_submitter(server),
        python="/env/bin/python",
        approval_file=tmp_path / "bundle.json",
        workspace=tmp_path / "ws",
        run_directory=tmp_path / "ws" / "run",
        cohort_size=size,
    )


def test_a_wave_of_three_is_one_array_of_three(tmp_path):
    lines = _script(tmp_path, size=3).splitlines()
    assert any(
        line.startswith("#SBATCH --array=0-2%") for line in lines
    ), lines


def test_the_host_throttles_it(tmp_path):
    directive = next(
        line
        for line in _script(tmp_path, size=7).splitlines()
        if line.startswith("#SBATCH --array=")
    )
    assert directive == "#SBATCH --array=0-6%4", directive
    tighter = next(
        line
        for line in _script(
            tmp_path, size=7, server=_server(MAX_CONCURRENT_TASKS=2)
        ).splitlines()
        if line.startswith("#SBATCH --array=")
    )
    assert tighter == "#SBATCH --array=0-6%2"


def test_each_element_runs_its_own_calculation(tmp_path):
    body = _script(tmp_path, size=3)
    assert "--cohort-element ${SLURM_ARRAY_TASK_ID}" in body, body
    # Its own result file, so N writers do not define the cycle by who
    # finished last.
    assert "execution-result.${SLURM_ARRAY_TASK_ID}.json" in body


def test_no_element_wakes_the_goal(tmp_path):
    """N tails would wake the model N times."""

    assert " agent wake " not in _script(tmp_path, size=3)


def test_the_wake_is_one_job_that_fires_however_the_array_ended(tmp_path):
    script = build_wake_dispatch_script(
        submitter=_submitter(_server(), label="goal-g1-cycle-1-wake"),
        python="/env/bin/python",
        workspace=tmp_path / "ws",
        goal_id="g1",
        array_job_id="4242",
    )
    lines = script.splitlines()
    assert "#SBATCH --dependency=afterany:4242" in lines, lines
    assert "#SBATCH --kill-on-invalid-dep=yes" in lines, (
        "a dependency that can never be satisfied stays PENDING forever "
        "and nothing in the tree can see it"
    )
    wake = next(line for line in lines if " agent wake " in line)
    assert wake.endswith("--goal g1")
    # One element's array is not this job's work.
    assert "--array" not in script


def test_the_wake_is_small_enough_to_schedule_immediately(tmp_path):
    """It runs one host process, so it must not inherit the cohort's own
    allocation: a wake asking 64 cores queues behind real science."""

    script = build_wake_dispatch_script(
        submitter=_submitter(_server(), label="w"),
        python="/env/bin/python",
        workspace=tmp_path / "ws",
        goal_id="g1",
        array_job_id="4242",
    )
    cores = next(
        line for line in script.splitlines() if "--ntasks-per-node=" in line
    )
    assert "--ntasks-per-node=1" in cores, cores
