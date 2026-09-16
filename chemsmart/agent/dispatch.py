"""Scheduler dispatch of one approved run: submit, park, and let the job
wake the goal.

The approved partition is not re-planned or re-approved here. The job
script runs the same provider-free executor a local run uses, inside the
allocation, and its own tail re-enters the goal with ``chemsmart agent
wake`` -- so no poller, no scheduler accounting, and no second decision
are needed. The scheduler directives come from the operator's server
profile through the existing submitters; what this module adds is the
script body and the receipt that names the job.
"""

from __future__ import annotations

import io
import json
import shlex
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable, Optional

from chemsmart.agent._contracts import ContractError

#: What a wake job asks for. It runs one host process -- reading the
#: run's own receipts and re-entering the goal -- and nothing about the
#: cohort's allocation applies to it.
_WAKE_MEMORY_GB = 4

#: Written in the run directory when a run is handed to a scheduler.
DISPATCH_RECEIPT_FILE = "dispatch.receipt.json"

#: Where the job writes the executor's own result record; the driver that
#: resumes the goal reads the same status words a local run returns.
EXECUTION_RESULT_FILE = "execution-result.json"


@dataclass(frozen=True)
class DispatchReceiptV1:
    """What was submitted, where, and how the goal will be woken."""

    scheduler: str
    job_id: str
    submitted_at: str
    submit_command: str
    submit_script: str
    run_directory: str
    approval_file: str
    goal_id: str
    cycle: int
    wake_command: str
    schema_version: str = "chemsmart.goal-dispatch-receipt.v1"
    #: What the scheduler was asked for, what the envelope requested, and
    #: the profile ceiling that bounded it -- kept because the approved
    #: review renders the envelope's numbers and a reader of this receipt
    #: is entitled to see whether the two agree. Absent on receipts minted
    #: before the envelope had any route to a scheduler directive.
    scheduler_request: Optional[dict[str, Any]] = None

    def public_record(self) -> dict[str, Any]:
        return asdict(self)


class _AgentRunJob:
    """The duck a submitter needs: a label, a folder, no program."""

    PROGRAM = None

    def __init__(self, *, label: str, folder: str) -> None:
        self.label = label
        self.folder = folder

    def is_complete(self) -> bool:
        return False


def _wake_command(*, python: str, workspace: Path, goal_id: str) -> str:
    return (
        f"{shlex.quote(python)} -m chemsmart agent wake "
        f"--workspace {shlex.quote(str(workspace))} "
        f"--goal {shlex.quote(goal_id)}"
    )


def build_dispatch_script(
    *,
    submitter: Any,
    python: str,
    approval_file: Path,
    workspace: Path,
    run_directory: Path,
    goal_id: str,
    wake: bool = True,
) -> str:
    """The job script: the scheduler's directives, then the executor's
    continuation, then the wake."""

    buffer = io.StringIO()
    submitter._write_bash_header(buffer)
    submitter._write_scheduler_options(buffer)
    submitter._write_extra_commands(buffer)
    submitter._write_change_to_job_directory(buffer)
    result_file = run_directory / EXECUTION_RESULT_FILE
    buffer.write(
        "# The approved partition, executed provider-free inside the "
        "allocation.\n"
    )
    buffer.write(
        f"{shlex.quote(python)} -m chemsmart agent run "
        f"--approval-file {shlex.quote(str(approval_file))} "
        f"--workspace {shlex.quote(str(workspace))} "
        f"--run-directory {shlex.quote(str(run_directory))} "
        f"--json > {shlex.quote(str(result_file))}\n"
    )
    if wake:
        buffer.write(
            "# The job's own tail re-enters the goal: no poller and no "
            "accounting needed.\n"
        )
        buffer.write(
            _wake_command(python=python, workspace=workspace, goal_id=goal_id)
            + "\n"
        )
    return buffer.getvalue()


def build_cohort_dispatch_script(
    *,
    submitter: Any,
    python: str,
    approval_file: Path,
    workspace: Path,
    run_directory: Path,
    cohort_size: int,
) -> str:
    """One array: each element runs one approved calculation of the wave.

    No element wakes the goal. N tails would wake the model N times, and
    a wave is one reasoning point; the wake is a separate job that fires
    once the whole array is over.

    Each element writes its own result file, so N writers do not define
    the cycle by whoever finished last.
    """

    buffer = io.StringIO()
    submitter._write_bash_header(buffer)
    submitter._write_array_scheduler_options(buffer, None, count=cohort_size)
    submitter._write_extra_commands(buffer)
    submitter._write_change_to_job_directory(buffer)
    element = "${SLURM_ARRAY_TASK_ID}"
    result = run_directory / f"execution-result.{element}.json"
    buffer.write(
        "# One approved calculation of this wave, executed provider-free.\n"
    )
    buffer.write(
        f"{shlex.quote(python)} -m chemsmart agent run "
        f"--approval-file {shlex.quote(str(approval_file))} "
        f"--workspace {shlex.quote(str(workspace))} "
        f"--run-directory {shlex.quote(str(run_directory))} "
        f"--cohort-element {element} "
        f"--json > {result}\n"
    )
    return buffer.getvalue()


def build_wake_dispatch_script(
    *,
    submitter: Any,
    python: str,
    workspace: Path,
    goal_id: str,
    array_job_id: str,
) -> str:
    """The cohort's single re-entry, gated on the whole array.

    ``afterany`` because the barrier is terminality and not success: a
    failed or cancelled element has ended, and its outcome is evidence
    the Agent must see. ``--kill-on-invalid-dep=yes`` because a
    dependency that can never be satisfied otherwise stays PENDING
    forever, which is a `Reason` and not a `JobState` and which nothing
    would ever stop waiting on; cancelled is a state this tree already
    classifies as terminal.

    It asks for one task: a wake runs one host process, and a wake
    inheriting the cohort's own allocation would queue behind real
    science.
    """

    from chemsmart.settings.scheduler_request import SchedulerRequestV1

    # A wake runs one host process. Left on the cohort's own allocation it
    # would ask for the whole node and queue behind real science, for a
    # job whose entire work is re-entering the goal.
    submitter.kwargs["scheduler_request"] = SchedulerRequestV1(
        cores=1,
        memory_gb=_WAKE_MEMORY_GB,
        gpu_count=0,
        hours=int(getattr(submitter.server, "num_hours", 0) or 1),
        sealed=False,
    )
    buffer = io.StringIO()
    submitter._write_bash_header(buffer)
    submitter._write_scheduler_options(buffer)
    buffer.write(f"#SBATCH --dependency=afterany:{array_job_id}\n")
    buffer.write("#SBATCH --kill-on-invalid-dep=yes\n")
    submitter._write_extra_commands(buffer)
    submitter._write_change_to_job_directory(buffer)
    buffer.write(
        "# The wave is over: the goal re-enters once, with all of it.\n"
    )
    buffer.write(
        _wake_command(python=python, workspace=workspace, goal_id=goal_id)
        + "\n"
    )
    return buffer.getvalue()


def dispatch_run_to_scheduler(
    *,
    approval_file: Path,
    workspace: Path,
    run_directory: Path,
    goal_id: str,
    cycle: int,
    server: str | None = None,
    python: str | None = None,
    wake: bool = True,
    resources: Any = None,
    sealed: bool = True,
    envelope: Any = None,
) -> DispatchReceiptV1:
    """Submit one approved run and return the receipt naming its job.

    Args:
        resources: The approved ``ExecutionResourceSpecV1``, which becomes
            the scheduler request. ``None`` falls back to the server
            profile's own numbers, which is what this function did before
            the envelope had any route here -- kept so an older caller
            behaves as it did rather than silently changing allocation.
        sealed: Whether the sealed-job memory ceiling applies. Sealed is
            the default for an Agent dispatch; ``chemsmart agent goal
            --unsealed`` drops the headroom and holds the request to the
            profile alone, which is the escape for a profile too small to
            leave anything above it.
        envelope: The approved execution envelope, whose episode window
            and postprocessing reserve are the wall clock asked for. The
            operator's NUM_HOURS caps it and no longer sets it.
    """

    from chemsmart.settings.server import Server

    resolved = Server.from_servername(server) if server else Server.current()
    scheduler = str(getattr(resolved, "scheduler", "") or "")
    if scheduler.upper() in {"", "LOCAL", "NONE", "UNKNOWN SCHEDULER"}:
        raise ContractError(
            "scheduler dispatch needs a server profile with a scheduler; "
            f"{resolved} has none. Run locally, or name a server with "
            "--server."
        )
    run_directory = Path(run_directory)
    run_directory.mkdir(parents=True, exist_ok=True)
    interpreter = python or sys.executable
    job = _AgentRunJob(
        label=f"goal-{goal_id}-cycle-{cycle}", folder=str(run_directory)
    )
    from chemsmart.settings.scheduler_request import (
        resolve_scheduler_request,
    )

    request = resolve_scheduler_request(
        resources=resources,
        server=resolved,
        sealed=sealed,
        envelope=envelope,
    )
    submitter = resolved.get_submitter(job, scheduler_request=request)
    script = build_dispatch_script(
        submitter=submitter,
        python=interpreter,
        approval_file=Path(approval_file).resolve(),
        workspace=Path(workspace).resolve(),
        run_directory=run_directory.resolve(),
        goal_id=goal_id,
        wake=wake,
    )
    (run_directory / submitter.submit_script).write_text(
        script, encoding="utf-8"
    )
    submission = resolved.submit_prepared(job)
    receipt = DispatchReceiptV1(
        scheduler=submission.scheduler,
        job_id=submission.job_id,
        submitted_at=submission.submitted_at,
        submit_command=submission.submit_command,
        submit_script=submission.submit_script,
        run_directory=str(run_directory.resolve()),
        approval_file=str(Path(approval_file).resolve()),
        goal_id=goal_id,
        cycle=int(cycle),
        wake_command=_wake_command(
            python=interpreter,
            workspace=Path(workspace).resolve(),
            goal_id=goal_id,
        ),
        scheduler_request=request.public_record(),
    )
    (run_directory / DISPATCH_RECEIPT_FILE).write_text(
        json.dumps(receipt.public_record(), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return receipt


def read_dispatch_receipt(run_directory: Path) -> DispatchReceiptV1 | None:
    try:
        record = json.loads(
            (Path(run_directory) / DISPATCH_RECEIPT_FILE).read_text(
                encoding="utf-8"
            )
        )
    except (OSError, json.JSONDecodeError):
        return None
    fields = {
        name: record.get(name)
        for name in DispatchReceiptV1.__dataclass_fields__
        if name in record
    }
    return DispatchReceiptV1(**fields)


def wait_for_dispatched_run(
    run_directory: Path,
    *,
    poll_seconds: float = 30.0,
    runner: Callable[..., "subprocess.CompletedProcess[str]"] = subprocess.run,
    sleep: Callable[[float], None] = time.sleep,
    max_polls: int | None = None,
) -> str:
    """Block until the job has written its result or the scheduler says it
    is over. Returns why the wait ended."""

    from chemsmart.settings.probe.scheduler_job import (
        parse_scontrol_job,
        scontrol_job_command,
    )

    run_directory = Path(run_directory)
    receipt = read_dispatch_receipt(run_directory)
    polls = 0
    while True:
        if (run_directory / EXECUTION_RESULT_FILE).is_file():
            return "result recorded"
        if receipt is None:
            return "no dispatch receipt to wait on"
        completed = runner(
            list(scontrol_job_command(receipt.job_id)),
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        state = parse_scontrol_job(
            completed.returncode,
            completed.stdout,
            completed.stderr,
            job_id=receipt.job_id,
        )
        if not state.known:
            return f"scheduler no longer knows job {receipt.job_id}"
        if state.terminal:
            return f"job {receipt.job_id} {state.state}"
        if state.dependency_unsatisfiable:
            # PENDING is correctly not terminal, so without this the loop
            # waits on a job the scheduler has already decided can never
            # run -- every thirty seconds, forever, with no mail and
            # nothing on disk to read.
            return (
                f"job {receipt.job_id} can never run: "
                f"{state.reason or 'dependency never satisfied'}"
            )
        polls += 1
        if max_polls is not None and polls >= max_polls:
            return f"job {receipt.job_id} still {state.state}"
        sleep(poll_seconds)


__all__ = [
    "DISPATCH_RECEIPT_FILE",
    "EXECUTION_RESULT_FILE",
    "DispatchReceiptV1",
    "build_cohort_dispatch_script",
    "build_dispatch_script",
    "build_wake_dispatch_script",
    "dispatch_run_to_scheduler",
    "read_dispatch_receipt",
    "wait_for_dispatched_run",
]
