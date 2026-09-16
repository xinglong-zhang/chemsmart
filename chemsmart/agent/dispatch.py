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
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Optional

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.cohort import (
    build_cohort_manifest,
    cohort_completion,
    execution_result_file,
)

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
    #: The job that will re-enter the goal, when it is not this job's own
    #: tail. A cohort's elements deliberately wake nothing -- N tails
    #: would wake the model N times -- so the goal is parked on the array
    #: and re-entered by a separate dependent job, and a reader of this
    #: receipt is entitled to know which job that is. Empty for the
    #: single-job path, whose own tail wakes.
    wake_job_id: str = ""
    #: The wave this run is, in the manifest's own order. Empty for a run
    #: that is not a cohort.
    cohort_node_ids: tuple[str, ...] = ()

    def public_record(self) -> dict[str, Any]:
        record = asdict(self)
        record["cohort_node_ids"] = list(self.cohort_node_ids)
        return record


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


def _write_dispatch_receipt(
    directory: Path,
    *,
    scheduler: str,
    job_id: str,
    submitted_at: str,
    submit_command: str,
    submit_script: str,
    run_directory: Path,
    approval_file: Path | str,
    goal_id: str,
    cycle: int,
    wake_command: str,
    request: Any,
    wake_job_id: str,
    cohort_node_ids: tuple[str, ...],
) -> DispatchReceiptV1:
    """Write the receipt, twice: once for the allocation, once for the
    jobs. One writer, so the two spellings cannot drift."""

    receipt = DispatchReceiptV1(
        scheduler=str(scheduler),
        job_id=str(job_id),
        submitted_at=str(submitted_at),
        submit_command=str(submit_command),
        submit_script=str(submit_script),
        run_directory=str(Path(run_directory).resolve()),
        approval_file=str(Path(approval_file).resolve()),
        goal_id=str(goal_id),
        cycle=int(cycle),
        wake_command=str(wake_command),
        scheduler_request=request.public_record(),
        wake_job_id=str(wake_job_id),
        cohort_node_ids=tuple(str(item) for item in cohort_node_ids),
    )
    (Path(directory) / DISPATCH_RECEIPT_FILE).write_text(
        json.dumps(receipt.public_record(), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return receipt


def _approved_bundle_digest(approval_file: Path | str) -> str:
    """The bundle's own declared digest, which is what admission reads.

    Hashing the file's bytes looked stronger and was a second authority
    for one question: a bundle's digest is `canonical_sha256` over its
    *content*, so the two never agree. Measured live on CUHK (goal
    `butane-wave-1`, array 2135192): the manifest carried the file hash,
    and the one element that got as far as running was refused by
    `authorise_cohort_element` with "this cohort was dispatched for
    another approval" -- against its own approval.
    """

    try:
        record = json.loads(Path(approval_file).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ContractError(
            f"the approval bundle at {approval_file} cannot be read, so "
            "a cohort cannot be bound to it"
        ) from exc
    body = record.get("workflow_execution_approval_bundle") or record
    digest = str((body or {}).get("bundle_sha256") or "")
    if not digest:
        raise ContractError(
            f"the approval bundle at {approval_file} declares no "
            "bundle_sha256, so a cohort manifest cannot name the "
            "approval it belongs to"
        )
    return digest


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

    # The guard the submitter layer already owns: without an array
    # directive the base writer falls back to an ordinary single-job
    # header, so a wave of N becomes one job running one calculation with
    # an empty element index. "A scheduler failure must not wear a job's
    # failure word" is that guard's own sentence, and this builder
    # called the writer directly and walked past it.
    submitter._require_array_support()
    buffer = io.StringIO()
    submitter._write_bash_header(buffer)
    submitter._write_array_scheduler_options(buffer, None, count=cohort_size)
    submitter._write_extra_commands(buffer)
    submitter._write_change_to_job_directory(buffer)
    # The element index is a shell variable here, so the host's own
    # naming is applied to a placeholder and then substituted: two
    # authors for one path is how a writer and a reader stop agreeing.
    # The submitter declares its own index variable; spelling Slurm's
    # here was a second authority that agrees only on Slurm.
    element = "${" + str(submitter.ARRAY_TASK_ID_VARIABLE) + "}"
    # The directory is quoted and the element index expands inside the
    # quotes. Every other path in this script is quoted and this one was
    # not: with a workspace under "My Drive", bash redirects to `/My` and
    # hands the rest to `chemsmart agent run` as stray arguments, so
    # every element dies before doing anything. Double quotes, because
    # single quotes would redirect to a literal `${SLURM_ARRAY_TASK_ID}`.
    directory = str(execution_result_file(run_directory, element=0).parent)
    result = '"{}/execution-result.{}.json"'.format(
        directory.replace("\\", "\\\\").replace('"', '\\"'), element
    )
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
    cohort_node_ids: tuple[str, ...] = (),
    bundle_sha256: str = "",
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
        cohort_node_ids: The wave the Agent selected, in the order it
            selected it. Given, the run is submitted as one throttled
            array -- one element per approved calculation -- plus one
            dependent job that wakes the goal when every element has
            reached a terminal state. Empty keeps the single-job path,
            whose own tail wakes: a run of one is not a degenerate array,
            it is the path every goal used before waves existed.
        bundle_sha256: The approved bundle the cohort belongs to, bound
            into the manifest digest so a manifest cannot be read against
            a different approval. Left empty it is read from the
            bundle's own declared digest, which is what admission
            compares against.
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
    # The allocation is decided here and the job ids are not, so the
    # allocation is written here. An element resolves what it may use
    # from this receipt, and the receipt used to be written after both
    # submissions -- while the manifest is written first precisely
    # because elements start promptly, and four of them started in the
    # same second on CUHK. An element inside that window found no
    # receipt, fell back to the *approval's* resources, and published
    # them: one element running outside its allocation and the rest
    # dying on bytes that disagree. The job ids are appended when the
    # scheduler names them.
    _write_dispatch_receipt(
        run_directory,
        scheduler=scheduler,
        job_id="",
        submitted_at="",
        submit_command="",
        submit_script=str(submitter.submit_script),
        run_directory=run_directory,
        approval_file=approval_file,
        goal_id=goal_id,
        cycle=cycle,
        wake_command=_wake_command(
            python=interpreter,
            workspace=Path(workspace).resolve(),
            goal_id=goal_id,
        ),
        request=request,
        wake_job_id="",
        cohort_node_ids=cohort_node_ids,
    )
    wake_job_id = ""
    if cohort_node_ids:
        # The manifest is written first and never after: it is the only
        # index-to-node mapping, and an element that starts promptly
        # resolves itself through it before its author would otherwise
        # have got there.
        manifest = build_cohort_manifest(
            goal_id=goal_id,
            cycle=int(cycle),
            bundle_sha256=(
                str(bundle_sha256) or _approved_bundle_digest(approval_file)
            ),
            node_ids=tuple(str(item) for item in cohort_node_ids),
            max_concurrent_tasks=submitter.max_concurrent_tasks(),
            created_at=datetime.now(timezone.utc)
            .isoformat()
            .replace("+00:00", "+00:00"),
        )
        manifest.write(run_directory)
        script = build_cohort_dispatch_script(
            submitter=submitter,
            python=interpreter,
            approval_file=Path(approval_file).resolve(),
            workspace=Path(workspace).resolve(),
            run_directory=run_directory.resolve(),
            cohort_size=len(manifest.node_ids),
        )
    else:
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
    if cohort_node_ids and wake:
        # A separate job, because the barrier is the wave's and not any
        # element's. It is submitted after the array so it can name it;
        # `--kill-on-invalid-dep=yes` is what stops an unsatisfiable
        # dependency from parking the goal on a PENDING job forever.
        wake_job = _AgentRunJob(
            label=f"goal-{goal_id}-cycle-{cycle}-wake",
            folder=str(run_directory),
        )
        wake_submitter = resolved.get_submitter(
            wake_job, scheduler_request=request
        )
        wake_script = build_wake_dispatch_script(
            submitter=wake_submitter,
            python=interpreter,
            workspace=Path(workspace).resolve(),
            goal_id=goal_id,
            array_job_id=str(submission.job_id),
        )
        (run_directory / wake_submitter.submit_script).write_text(
            wake_script, encoding="utf-8"
        )
        wake_job_id = str(resolved.submit_prepared(wake_job).job_id)
    return _write_dispatch_receipt(
        run_directory,
        scheduler=submission.scheduler,
        job_id=submission.job_id,
        submitted_at=submission.submitted_at,
        submit_command=submission.submit_command,
        submit_script=submission.submit_script,
        run_directory=run_directory,
        approval_file=approval_file,
        goal_id=goal_id,
        cycle=cycle,
        wake_command=_wake_command(
            python=interpreter,
            workspace=Path(workspace).resolve(),
            goal_id=goal_id,
        ),
        request=request,
        wake_job_id=wake_job_id,
        cohort_node_ids=cohort_node_ids,
    )


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
    # JSON has one sequence and this receipt has two: a cohort read back
    # as a list would not compare equal to the one that was written, and
    # the cohort is membership, which is the barrier's own question.
    if "cohort_node_ids" in fields:
        fields["cohort_node_ids"] = tuple(
            str(item) for item in fields["cohort_node_ids"] or ()
        )
    return DispatchReceiptV1(**fields)


def cohort_elements_may_still_run(
    run_directory: Path,
    *,
    runner: Callable[..., "subprocess.CompletedProcess[str]"] | None = None,
) -> bool | None:
    """Whether any element of this wave's array could still be working.

    The one question the scheduler is allowed to answer about a wave, and
    it is not a scientific one: the durable stream owns what every
    calculation *means*, and this owns only whether a process can still
    add to it.

    It exists because the stream cannot tell "still running" from
    "process gone". A launch lease answers that for a node that reserved;
    an element that died before reserving leaves nothing at all, so the
    barrier waits on it forever -- measured live on CUHK, where three
    elements of a wave of four died in their first second and the goal
    could never be woken.

    Returns:
        bool | None: ``None`` when there is nothing to ask -- no cohort,
        no dispatch receipt, or a scheduler that no longer knows the job
        -- in which case the stream's own answer stands unchanged.
    """

    from chemsmart.settings.probe.scheduler_job import (
        parse_squeue_array,
        squeue_array_command,
    )

    # Resolved here rather than as a default argument: a default binds
    # the function object at import, so a caller that replaces
    # `subprocess.run` -- a test, or a host that routes its own
    # processes -- would be silently ignored.
    runner = runner or subprocess.run
    receipt = read_dispatch_receipt(Path(run_directory))
    if receipt is None or not receipt.cohort_node_ids:
        return None
    completed = runner(
        list(squeue_array_command(receipt.job_id)),
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    array = parse_squeue_array(
        completed.returncode,
        completed.stdout,
        completed.stderr,
        array_job_id=receipt.job_id,
    )
    if not array.known:
        # An array the scheduler has forgotten is one whose elements have
        # all ended -- it is purged when nothing of it is left. Saying
        # "unknown" here would reinstate the permanent park.
        return False
    return not array.terminal


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
        parse_squeue_array,
        scontrol_job_command,
        squeue_array_command,
    )

    run_directory = Path(run_directory)
    receipt = read_dispatch_receipt(run_directory)
    polls = 0
    while True:
        complete, _pending = cohort_completion(run_directory)
        if complete is True:
            return "wave complete"
        if (
            complete is None
            and (run_directory / EXECUTION_RESULT_FILE).is_file()
        ):
            return "result recorded"
        if receipt is None:
            return "no dispatch receipt to wait on"
        if receipt.cohort_node_ids:
            # An array has one state per element and `scontrol show job
            # N` prints one record per element, each with `JobId=N_0`,
            # `N_1`, ... and never the bare `N` the scalar reader
            # compares against -- while its field harvest is
            # first-occurrence-wins. A wave of three whose first printed
            # element completed while two still ran therefore read as
            # "job N COMPLETED", returned, and let the Agent reason over
            # a running wave. The barrier's question is a multiset, so
            # this asks squeue for the elements, expanded (`-r`, because
            # a queued array prints one compressed row that names no
            # membership at all).
            completed = runner(
                list(squeue_array_command(receipt.job_id)),
                capture_output=True,
                text=True,
                timeout=30,
                check=False,
            )
            array = parse_squeue_array(
                completed.returncode,
                completed.stdout,
                completed.stderr,
                array_job_id=receipt.job_id,
            )
            if not array.known:
                return f"scheduler no longer knows job {receipt.job_id}"
            if array.terminal:
                return f"array {receipt.job_id}: every element ended"
            polls += 1
            if max_polls is not None and polls >= max_polls:
                return (
                    f"array {receipt.job_id}: "
                    f"{len(array.unfinished)} element(s) still running"
                )
            sleep(poll_seconds)
            continue
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
    "cohort_elements_may_still_run",
    "build_dispatch_script",
    "build_wake_dispatch_script",
    "dispatch_run_to_scheduler",
    "read_dispatch_receipt",
    "wait_for_dispatched_run",
]
