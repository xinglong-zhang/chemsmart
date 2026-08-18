"""Stopping an engine and destroying its evidence are different acts.

ChemSmart's runner stages work in scratch and copies results back in
`_postrun`. An external SIGTERM -- the scheduler's timeout, an operator's
scancel, the agent host's node deadline -- killed the wrapper before that
copy-back, so everything the engine had computed was stranded in scratch and
the job folder received nothing.

Measured on a real approved run: a 73-point relaxed scan of dimercaprol timed
out with 55 points converged, and the typed evidence recorded only
`orca.result.output_count` -- no surface, no geometries -- while the data sat
intact in scratch. The parsers deliberately read a truncated surface as
partial evidence; the runner's death was the only reason they never got it.

The runner now traps SIGTERM, terminates the child, and lets the ordinary
postrun path run. Completion state stays honest -- the job is still incomplete,
scratch deletion is already skipped for incomplete jobs, and the wrapper still
exits nonzero. This drives the real `run()` with a real child process and a
real signal, because a fake would test the handler's spelling rather than its
timing.
"""

from __future__ import annotations

import os
import signal
import subprocess
import sys
import threading

import pytest

from chemsmart.jobs.runner import JobRunner


class _RecordingRunner(JobRunner):
    """The smallest concrete runner that records its lifecycle."""

    JOBTYPES = ["probe"]
    PROGRAM = "probe"
    FAKE = False
    SCRATCH = False

    def __init__(self):
        # Deliberately not calling super().__init__: the base constructor
        # resolves servers and executables this probe does not need. Only the
        # attributes run() touches are set.
        self.scratch = False
        self.delete_scratch = False
        self.server = None
        self.calls: list[str] = []

    def _prerun(self, job):
        self.calls.append("prerun")

    def _write_input(self, job):
        self.calls.append("write_input")

    def _get_command(self, job):
        return [sys.executable, "-c", "import time; time.sleep(600)"]

    def _create_process(self, job, command, env):
        return subprocess.Popen(command, env=env)

    def _postrun(self, job, **kwargs):
        self.calls.append("postrun")

    def _postrun_cleanup(self, job):
        self.calls.append("postrun_cleanup")

    def _update_os_environ(self, job):
        return os.environ.copy()


class _Job:
    def is_complete(self):
        return False

    def __repr__(self):
        return "probe-job"


def test_sigterm_stops_the_child_and_still_runs_postrun():
    runner = _RecordingRunner()
    # Deliver the signal the way a scheduler does: from outside the run call,
    # while the child is alive.
    timer = threading.Timer(
        0.5, lambda: os.kill(os.getpid(), signal.SIGTERM)
    )
    timer.start()
    try:
        with pytest.raises(subprocess.CalledProcessError) as failure:
            runner.run(_Job())
    finally:
        timer.cancel()

    # The child died from the forwarded terminate, not from a full sleep.
    assert failure.value.returncode == -signal.SIGTERM
    # The copy-back path ran even though the job was externally stopped.
    assert "postrun" in runner.calls
    assert runner.calls.index("postrun") > runner.calls.index("prerun")


def test_the_previous_handler_is_restored_afterwards():
    original = signal.getsignal(signal.SIGTERM)
    runner = _RecordingRunner()
    timer = threading.Timer(
        0.3, lambda: os.kill(os.getpid(), signal.SIGTERM)
    )
    timer.start()
    try:
        with pytest.raises(subprocess.CalledProcessError):
            runner.run(_Job())
    finally:
        timer.cancel()

    assert signal.getsignal(signal.SIGTERM) is original


def test_an_untouched_run_behaves_exactly_as_before():
    """No signal, quick child: the ordinary lifecycle, in order."""

    runner = _RecordingRunner()
    runner._get_command = lambda job: [sys.executable, "-c", "pass"]

    returncode = runner.run(_Job())

    assert returncode == 0
    assert runner.calls == [
        "prerun",
        "write_input",
        "postrun",
        "postrun_cleanup",
    ]
