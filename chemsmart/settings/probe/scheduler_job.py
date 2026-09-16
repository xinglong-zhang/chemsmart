"""Scheduler job identity and state: what a submission said, what the
scheduler still knows.

``sbatch`` answers a submission with one line naming the job; ``scontrol``
and ``squeue`` answer for a job the scheduler still remembers, which on a
site without accounting storage is only until ``MinJobAge`` has passed.
Every parser here is a pure function over ``(returncode, stdout, stderr)``
in the probe convention, so scheduler behaviour is provable from canned
real outputs without a cluster. A job the scheduler no longer knows is
reported as unknown, never guessed at: after that window the run's own
durable record is the completion evidence.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from .units import ProbeUnitError

#: Scheduler states after which a job will not run again.
TERMINAL_JOB_STATES = frozenset(
    {
        "BOOT_FAIL",
        "CANCELLED",
        "COMPLETED",
        "DEADLINE",
        "FAILED",
        "NODE_FAIL",
        "OUT_OF_MEMORY",
        "PREEMPTED",
        "TIMEOUT",
    }
)

#: job id|state|time used|submit|start|end -- pinned, never human default.

#: element|array job id|task id|state|time used|submit|start|end.
#: `%F` is the join: `%i` prints ``<arrayid>_<task>`` and never the bare
#: array id, so an element cannot be matched to its cohort without it.
_SQUEUE_ARRAY_FORMAT = "%i|%F|%K|%T|%M|%V|%S|%e"

_SBATCH_LINE = re.compile(r"Submitted batch job (\d+)")
_BARE_JOB_ID = re.compile(r"^\d+(?:\.\S+)?$")


@dataclass(frozen=True)
class SchedulerJobStateV1:
    """What the scheduler currently says about one job.

    ``known`` is False when the scheduler no longer has the job at all;
    every other field is then empty, and the caller must read the run's
    own record instead of inferring an outcome from absence.
    """

    job_id: str
    known: bool
    state: str = ""
    exit_code: str = ""
    submit_time: str = ""
    start_time: str = ""
    end_time: str = ""
    run_seconds: int | None = None
    #: Slurm's own word for *why* a job is where it is. A dependency that
    #: can never be satisfied is a `Reason`, never a `JobState`, so a
    #: reader that looks only at the state waits on it forever.
    reason: str = ""

    @property
    def dependency_unsatisfiable(self) -> bool:
        """Whether this job is waiting on something that cannot happen.

        `PENDING` is correctly not terminal, so nothing else in the tree
        would ever stop waiting. Slurm says `DependencyNeverSatisfied`
        once it knows, and that is the whole signal.
        """

        return self.reason.strip().lower() == "dependencyneversatisfied"

    @property
    def terminal(self) -> bool:
        return self.known and self.state.split()[0] in TERMINAL_JOB_STATES


def parse_submission(returncode: int, stdout: str, stderr: str) -> str:
    """The job id a submission command printed.

    Accepts ``Submitted batch job N`` (sbatch), a bare ``N`` (sbatch
    ``--parsable``) and ``N.server`` (qsub). A submission that exited
    nonzero or printed no id is refused with the scheduler's own words.
    """

    stdout, stderr = stdout or "", stderr or ""
    if returncode != 0:
        raise ProbeUnitError(
            f"submission exited {returncode}: {stderr.strip() or 'no output'}"
        )
    match = _SBATCH_LINE.search(stdout)
    if match:
        return match.group(1)
    for line in stdout.splitlines():
        token = line.strip().split(";")[0]
        if _BARE_JOB_ID.match(token):
            return token
    raise ProbeUnitError(
        "submission printed no job id: "
        f"{stdout.strip() or stderr.strip() or 'no output'}"
    )


def scontrol_job_command(job_id: str) -> tuple[str, ...]:
    return ("scontrol", "show", "job", str(job_id))


def parse_elapsed_seconds(text: str) -> int | None:
    """Seconds from a Slurm elapsed field: ``[D-]HH:MM:SS`` or the shortened
    ``M:SS`` / ``H:MM:SS`` that ``squeue %M`` prints."""

    cleaned = text.strip()
    if not cleaned or cleaned in {"N/A", "INVALID", "UNLIMITED"}:
        return None
    days = 0
    if "-" in cleaned:
        day_text, cleaned = cleaned.split("-", 1)
        days = int(day_text)
    parts = cleaned.split(":")
    try:
        numbers = [int(part) for part in parts]
    except ValueError:
        return None
    if len(numbers) == 3:
        hours, minutes, seconds = numbers
    elif len(numbers) == 2:
        hours, (minutes, seconds) = 0, numbers
    elif len(numbers) == 1:
        hours, minutes, seconds = 0, 0, numbers[0]
    else:
        return None
    return ((days * 24 + hours) * 60 + minutes) * 60 + seconds


def _unknown(job_id: str) -> SchedulerJobStateV1:
    return SchedulerJobStateV1(job_id=str(job_id), known=False)


def parse_scontrol_job(
    returncode: int, stdout: str, stderr: str, *, job_id: str
) -> SchedulerJobStateV1:
    """State from ``scontrol show job``; unknown when the scheduler has
    forgotten the job (``Invalid job id specified``)."""

    stdout, stderr = stdout or "", stderr or ""
    if returncode != 0 or "Invalid job id" in stderr:
        return _unknown(job_id)
    fields: dict[str, str] = {}
    for token in stdout.split():
        key, sep, value = token.partition("=")
        if sep and key not in fields:
            fields[key] = value
    if fields.get("JobId", "") != str(job_id):
        return _unknown(job_id)
    return SchedulerJobStateV1(
        job_id=str(job_id),
        known=True,
        state=fields.get("JobState", ""),
        exit_code=fields.get("ExitCode", ""),
        submit_time=fields.get("SubmitTime", ""),
        start_time=fields.get("StartTime", ""),
        end_time=fields.get("EndTime", ""),
        run_seconds=parse_elapsed_seconds(fields.get("RunTime", "")),
        # Harvested all along and then discarded. A dependency that can
        # never be satisfied says so here and nowhere else.
        reason=fields.get("Reason", ""),
    )


@dataclass(frozen=True)
class SchedulerArrayElementV1:
    """One array element's own state."""

    element_id: str
    task_id: str
    state: str
    submit_time: str = ""
    start_time: str = ""
    end_time: str = ""
    run_seconds: int | None = None

    @property
    def terminal(self) -> bool:
        return self.state in TERMINAL_JOB_STATES


@dataclass(frozen=True)
class SchedulerArrayStateV1:
    """What one array job's elements are doing, as a multiset.

    The barrier's question is "has every member reached a terminal
    state", which a scalar job state cannot answer: an array has one
    state per element, and squashing them loses exactly the fact the
    barrier needs.
    """

    array_job_id: str
    known: bool
    elements: tuple[SchedulerArrayElementV1, ...] = ()

    @property
    def unfinished(self) -> tuple[str, ...]:
        """Task ids not yet in a terminal state."""

        return tuple(
            element.task_id
            for element in self.elements
            if not element.terminal
        )

    @property
    def terminal(self) -> bool:
        """Whether every member has ended, however it ended.

        An array the scheduler cannot see is **not** terminal: "I cannot
        find it" is not "every member finished", and treating it as such
        would wake the model over a cohort whose evidence nobody read.
        Terminality, not success -- a failed or cancelled element has
        ended, and its outcome is evidence the Agent must see.
        """

        if not self.known or not self.elements:
            return False
        return all(element.terminal for element in self.elements)


def squeue_array_command(array_job_id: str) -> tuple[str, ...]:
    """One row per element, expanded.

    By default squeue *compresses* an array's pending elements into a
    single range row -- `%i` reads `4242_[1-6%4]` -- so a barrier asking
    "is every member terminal" would count one row for six calculations
    and could report a range string as an unfinished task id. `-r`
    (`--array`) expands them.
    """

    return (
        "squeue",
        "-h",
        "-r",
        "-t",
        "all",
        "-j",
        str(array_job_id),
        "-o",
        _SQUEUE_ARRAY_FORMAT,
    )


def parse_squeue_array(
    returncode: int, stdout: str, stderr: str, *, array_job_id: str
) -> SchedulerArrayStateV1:
    """Every element of one array job, joined on ArrayJobId.

    Measured on CUHK (job 2135107): ``squeue -j <arrayid>`` prints
    ``%i`` as ``<arrayid>_<task>``, so a filter comparing it against the
    bare array id matches nothing. ``scontrol`` is worse than useless
    here -- it opens with the pending aggregate under the bare array id
    and lists each started element under its *own* job id, while the
    field harvest is first-occurrence-wins across the whole output, so
    it reports one arbitrary record's state as the array's.
    """

    if returncode != 0:
        return SchedulerArrayStateV1(
            array_job_id=str(array_job_id), known=False
        )
    elements = []
    for line in (stdout or "").splitlines():
        parts = [part.strip() for part in line.split("|")]
        if len(parts) != 8 or parts[1] != str(array_job_id):
            continue
        # A compressed range row would make one line stand for several
        # calculations, so it is not read as an element. Asking with -r
        # should prevent it; a site that returns one anyway makes the
        # cohort unknown rather than silently undercounted.
        if not parts[2].isdigit():
            return SchedulerArrayStateV1(
                array_job_id=str(array_job_id), known=False
            )
        elements.append(
            SchedulerArrayElementV1(
                element_id=parts[0],
                task_id=parts[2],
                state=parts[3],
                run_seconds=parse_elapsed_seconds(parts[4]),
                submit_time=parts[5],
                start_time=parts[6],
                end_time=parts[7],
            )
        )
    return SchedulerArrayStateV1(
        array_job_id=str(array_job_id),
        known=bool(elements),
        elements=tuple(elements),
    )


__all__ = [
    "TERMINAL_JOB_STATES",
    "SchedulerJobStateV1",
    "parse_elapsed_seconds",
    "parse_scontrol_job",
    "parse_squeue_array",
    "squeue_array_command",
    "SchedulerArrayElementV1",
    "SchedulerArrayStateV1",
    "parse_submission",
    "scontrol_job_command",
]
