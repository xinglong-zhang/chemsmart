"""Physical concurrency is the host's, and it is on by default.

`%N` in `--array=0-9%4` bounds *simultaneously running array tasks of one
array job*. It is not nodes, not jobs, not allocations, and not the
scientific width of a wave: the Agent says how many independent
calculations the science needs, and the host decides how many of them run
at once (owner ruling, 2026-09-16).

Two defects this holds shut. The throttle was opt-in -- `num_nodes=None`
wrote a bare `--array=0-N` with no cap -- and every test passed it
explicitly, so the suite proved the throttle was *writable* and never
that it was *on*. And the parameter was called `num_nodes` while meaning
"max concurrent tasks", so a caller passing an honest node count of 1
would silently serialise a whole cohort.
"""

from __future__ import annotations

import re

from chemsmart.settings.server import Server
from chemsmart.settings.submitters import DEFAULT_MAX_CONCURRENT_TASKS

_DIRECTIVE = re.compile(r"^#SBATCH --array=(\d+)-(\d+)(?:%(\d+))?\s*$", re.M)


def _server(**overrides):
    base = dict(
        SCHEDULER="SLURM",
        SUBMIT_COMMAND="sbatch",
        NUM_CORES=64,
        MEM_GB=160,
        NUM_HOURS=24,
        NUM_GPUS=0,
        QUEUE_NAME="chpc",
    )
    base.update(overrides)
    return Server("canned", **base)


def _directive(server, count):
    import io
    from types import SimpleNamespace

    job = SimpleNamespace(label="wave", folder="/tmp", PROGRAM="XTB")
    submitter = server.get_submitter(job)
    submitter.jobs = [job] * count
    buffer = io.StringIO()
    submitter._write_array_scheduler_options(buffer, None)
    match = _DIRECTIVE.search(buffer.getvalue())
    assert match, buffer.getvalue()
    first, last, throttle = match.groups()
    return int(first), int(last), (None if throttle is None else int(throttle))


def test_the_default_is_four_not_unthrottled():
    assert DEFAULT_MAX_CONCURRENT_TASKS == 4
    first, last, throttle = _directive(_server(), 7)
    assert (first, last) == (0, 6), "a wave of seven is seven elements"
    assert throttle == 4, (
        "the cohort was submitted unthrottled; the default is a safety and "
        "fairness default, not an opt-in"
    )


def test_the_operator_owns_it_through_the_server_profile():
    _, _, throttle = _directive(_server(MAX_CONCURRENT_TASKS=2), 7)
    assert throttle == 2


def test_scientific_width_is_not_reduced_by_the_throttle():
    """Seven calculations are one cohort and one wake, never 'four then
    three'. The host queues the remainder inside the same array."""

    first, last, throttle = _directive(_server(MAX_CONCURRENT_TASKS=4), 7)
    assert last - first + 1 == 7
    assert throttle == 4


def test_a_cohort_smaller_than_the_throttle_is_not_padded():
    first, last, throttle = _directive(_server(), 2)
    assert (first, last) == (0, 1)
    assert throttle == 4
