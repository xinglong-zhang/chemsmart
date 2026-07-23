"""
Shared batch job infrastructure.

Provides the abstract ``BatchJob`` base class for orchestrating collections
of engine-specific jobs, plus submit-time helpers for scheduler arrays:

- serial local child-job execution with full resources per child
- array-task execution of a single child (scheduler array env)
- fault-tolerant execution with aggregated failures
- ``batch_entry`` attach and per-task CLI rewrite for array runscripts
- nestable single-child invoke via ``--child-index`` / lazy child factories
- ``NestableJobMixin`` shared ``get_array_child_job`` / ``get_array_child_jobs`` API

Cluster concurrency is via ``chemsmart sub`` scheduler arrays, not
in-process multi-node fan-out.
"""

import logging
import os
from abc import ABCMeta
from contextlib import contextmanager
from enum import Enum
from typing import (
    Any,
    Callable,
    Iterator,
    Mapping,
    Optional,
    Sequence,
    Type,
    TypeVar,
)

from chemsmart.jobs.job import Job
from chemsmart.utils.mixins import RegistryMeta

logger = logging.getLogger(__name__)

BatchJobT = TypeVar("BatchJobT", bound="BatchJob")
RewriteCliFn = Callable[
    [Sequence[str], Optional[Mapping[str, Any]]],
    list[str],
]


class BatchExecutionMode(str, Enum):
    """How ``BatchJob.run()`` executes children.

    ``LOCAL_BATCH``
        Run all children in-process (serial policy by default).
    ``ARRAY_TASK``
        Run one child selected by the scheduler array task id.
    """

    LOCAL_BATCH = "local_batch"
    ARRAY_TASK = "array_task"


_ARRAY_TASK_ID_ENV_VARS = (
    "SLURM_ARRAY_TASK_ID",
    "PBS_ARRAYID",
    "LSB_JOBINDEX",
)


def resolve_array_task_id() -> Optional[int]:
    """Return the 1-based scheduler array task id, or ``None`` if unset.

    Checks ``SLURM_ARRAY_TASK_ID``, ``PBS_ARRAYID``, then ``LSB_JOBINDEX``.
    """
    for key in _ARRAY_TASK_ID_ENV_VARS:
        value = os.environ.get(key)
        if value is None or value == "":
            continue
        try:
            return int(value)
        except ValueError as exc:
            raise ValueError(
                f"Invalid {key}={value!r}; expected an integer task id."
            ) from exc
    return None


@contextmanager
def cleared_array_task_env() -> Iterator[None]:
    """Temporarily remove scheduler array task id environment variables.

    Top-level ``BatchJob._run_array_task`` uses this around ``child.run()`` so
    nestable parents (crest/QRC/dias/traj) do not treat the outer array task
    id as a nested child index. Nestable array submit still re-invokes the
    parent CLI with these variables set (no outer ``BatchJob`` array wrapper).
    """
    saved: dict[str, str] = {}
    for key in _ARRAY_TASK_ID_ENV_VARS:
        if key in os.environ:
            saved[key] = os.environ.pop(key)
    try:
        yield
    finally:
        os.environ.update(saved)


def resolve_batch_execution_mode() -> BatchExecutionMode:
    """Return ``ARRAY_TASK`` when a scheduler array task id is set.

    Otherwise return ``LOCAL_BATCH``.
    """
    if resolve_array_task_id() is not None:
        return BatchExecutionMode.ARRAY_TASK
    return BatchExecutionMode.LOCAL_BATCH


class BatchExecutionError(RuntimeError):
    """Raised when one or more child jobs fail in a batch run."""


_LEGACY_JOB_LIST_DEPRECATION = (
    "Returning a bare list of Job instances from CLI is deprecated; "
    "return a BatchJob instead."
)


def warn_legacy_job_list(*, stacklevel: int = 2) -> None:
    """Emit the deprecation warning for bare ``list[Job]`` CLI results."""
    import warnings

    warnings.warn(
        _LEGACY_JOB_LIST_DEPRECATION,
        DeprecationWarning,
        stacklevel=stacklevel,
    )


class BatchJobMeta(RegistryMeta, ABCMeta):
    """Metaclass combining registry support with abstract-base semantics."""


class NestableJobMixin:
    """Shared array-child API for crest/QRC/dias/traj nestable parents.

    Subclasses set ``child_index`` (1-based, optional) in ``__init__`` and
    implement:

    - ``num_array_children`` — child count without necessarily building them
    - ``get_array_child_job(index)`` — build one 0-based child lazily

    ``get_array_child_jobs`` defaults to mapping those primitives. Override it
    only when a bulk build is cheaper (e.g. DI-AS phase concatenation).
    """

    child_index: Optional[int] = None

    @property
    def num_array_children(self) -> int:
        """Return the number of nestable array children."""
        raise NotImplementedError(
            f"{type(self).__name__} must implement num_array_children"
        )

    def get_array_child_job(self, index: int) -> Job:
        """Build the nestable child at 0-based *index*."""
        raise NotImplementedError(
            f"{type(self).__name__} must implement get_array_child_job"
        )

    def validate_array_child_index(self, index: int) -> int:
        """Validate 0-based *index*; return ``num_array_children``."""
        total = self.num_array_children
        if index < 0 or index >= total:
            raise ValueError(
                f"Array child index {index} out of range for {total} "
                f"child job(s) of {self.label!r}; expected 0..{total - 1}."
            )
        return total

    def get_array_child_jobs(self) -> list[Job]:
        """Return all nestable children for scheduler array submission."""
        return [
            self.get_array_child_job(i) for i in range(self.num_array_children)
        ]


class BatchJob(Job, metaclass=BatchJobMeta):
    """
    Abstract controller for running a batch of child jobs.

    Orchestration is engine-agnostic. Subclasses set ``PROGRAM`` and hold
    the child job list.

    ``run()`` selects execution mode from the environment:

    - ``array_task`` — one child at the 1-based scheduler array task id,
      with full resources
    - ``local_batch`` — all children serially with full resources

    ``fail_fast`` stops serial local submission after the first
    unsuccessful outcome.
    ``nested_serial`` marks crest/QRC/dias/traj nested batches for
    ``policy=serial_nested`` logging.
    ``rewrite_cli`` is the optional per-task CLI rewriter used by
    ``chemsmart sub`` when submitting this batch as a scheduler array.
    """

    PROGRAM: Optional[str] = None
    REGISTERABLE: bool = False

    def __init__(
        self,
        jobs: Optional[Sequence[Job]],
        fail_fast: bool = False,
        write_outcome_logs: bool = False,
        label: str = "batch_job",
        jobrunner: Any = None,
        nested_serial: bool = False,
        rewrite_cli: Optional[RewriteCliFn] = None,
        **kwargs,
    ) -> None:
        super().__init__(
            molecule=None,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
        self.jobs: list[Job] = list(jobs) if jobs is not None else []
        self.fail_fast = bool(fail_fast)
        self.write_outcome_logs = bool(write_outcome_logs)
        self.nested_serial = bool(nested_serial)
        self.rewrite_cli = rewrite_cli
        self._last_batch_outcomes: list[dict[str, Any]] = []

    def run(self, **kwargs: Any) -> None:
        """Run this batch in ``local_batch`` or ``array_task`` mode."""
        mode = resolve_batch_execution_mode()
        if mode is BatchExecutionMode.ARRAY_TASK:
            self._run_array_task(**kwargs)
            return
        self._run_local_batch(**kwargs)

    def _run_array_task(self, **kwargs: Any) -> None:
        """Run the single child selected by the scheduler array task id.

        ``SLURM_ARRAY_TASK_ID`` / ``PBS_ARRAYID`` / ``LSB_JOBINDEX`` are
        treated as 1-based indexes into ``self.jobs``.

        Array-task env vars are cleared around ``child.run()`` so nestable
        children do not reuse the outer task id for nested selection.
        """
        task_id = resolve_array_task_id()
        if task_id is None:
            raise RuntimeError(
                "array_task mode requires SLURM_ARRAY_TASK_ID, "
                "PBS_ARRAYID, or LSB_JOBINDEX."
            )
        total_jobs = len(self.jobs)
        if total_jobs == 0:
            raise ValueError(f"BatchJob {self} has no child jobs to run.")
        child_index = task_id - 1
        if child_index < 0 or child_index >= total_jobs:
            raise ValueError(
                f"Array task id {task_id} out of range for {total_jobs} "
                f"child job(s); expected 1..{total_jobs}."
            )

        child = self.jobs[child_index]
        runner = self.jobrunner
        if runner is not None:
            cores = runner.num_cores
            mem_gb = runner.mem_gb
        else:
            cores = None
            mem_gb = None

        logger.info(
            "BatchJob %r: execution=%s, task=%s/%s, cores=%s, "
            "mem_gb=%s, child=%s",
            self.label,
            BatchExecutionMode.ARRAY_TASK.value,
            task_id,
            total_jobs,
            cores,
            mem_gb,
            child.label,
        )
        # Isolate nestable selection from this outer array task id.
        with cleared_array_task_env():
            outcome = self._submit_job(child, **kwargs)
        outcomes = [outcome]
        self._last_batch_outcomes = outcomes
        if self.write_outcome_logs:
            self._write_outcome_logs(outcomes)
        self._raise_if_failures(outcomes, total_jobs=1)

    def _run_local_batch(self, **kwargs: Any) -> None:
        """Run all children serially with full resources per child."""
        total_jobs = len(self.jobs)
        runner = self.jobrunner
        if runner is not None:
            cores = runner.num_cores
            mem_gb = runner.mem_gb
        else:
            cores = None
            mem_gb = None
        policy = "serial_nested" if self.nested_serial else "serial"
        logger.info(
            "BatchJob %r: execution=%s, children=%s, policy=%s, "
            "cores=%s, mem_gb=%s",
            self.label,
            BatchExecutionMode.LOCAL_BATCH.value,
            total_jobs,
            policy,
            cores,
            mem_gb,
        )
        outcomes = self._run_jobs_serially(self.jobs, **kwargs)

        self._last_batch_outcomes = outcomes
        if self.write_outcome_logs:
            self._write_outcome_logs(outcomes)
        self._raise_if_failures(outcomes, total_jobs=total_jobs)

    def _raise_if_failures(
        self,
        outcomes: Sequence[dict[str, Any]],
        *,
        total_jobs: int,
    ) -> None:
        """Raise ``BatchExecutionError`` summarizing attempted failures."""
        failures = [item for item in outcomes if not item["success"]]
        if not failures:
            return

        attempted = len(outcomes)
        not_started = max(0, total_jobs - attempted)
        lines = [f"- {item['label']}: {item['error']}" for item in failures]

        if not_started:
            summary = (
                f"{attempted} attempted, {len(failures)} failed, "
                f"{not_started} not started"
            )
        else:
            summary = f"{len(failures)} of {attempted} batch job(s) failed"

        raise BatchExecutionError(summary + ":\n" + "\n".join(lines))

    def _run_jobs_serially(
        self,
        jobs: Sequence[Job],
        **kwargs: Any,
    ) -> list[dict[str, Any]]:
        """Submit child jobs one-by-one.

        When ``fail_fast`` is enabled, stop after the first unsuccessful
        outcome. Jobs not started are omitted from the returned outcomes;
        ``_raise_if_failures`` derives the not-started count from
        ``total_jobs - len(outcomes)``.
        """
        outcomes: list[dict[str, Any]] = []
        for index, job in enumerate(jobs):
            outcome = self._submit_job(job, **kwargs)
            outcomes.append(outcome)
            if self.fail_fast and not outcome["success"]:
                remaining = len(jobs) - (index + 1)
                if remaining:
                    logger.warning(
                        "fail_fast enabled: stopping serial batch after "
                        f"{job.label} failed; {remaining} job(s) not started."
                    )
                break
        return outcomes

    def _submit_job(
        self,
        job: Job,
        **kwargs: Any,
    ) -> dict[str, Any]:
        """Configure a runner copy and execute a single child job."""
        try:
            self._build_jobrunner(job)
            job.run(**kwargs)
            if not job.is_complete():
                msg = "job incomplete after execution"
                logger.error(f"Job {job.label}: {msg}")
                return {
                    "label": job.label,
                    "success": False,
                    "error": msg,
                }
            return {
                "label": job.label,
                "success": True,
                "error": "",
            }
        except Exception as e:
            logger.error(
                f"Job {job.label} failed during batch execution: {e}",
                exc_info=True,
            )
            return {
                "label": job.label,
                "success": False,
                "error": str(e),
            }

    def _build_jobrunner(self, job: Job) -> Any:
        """Copy the batch runner onto *job* with full resource allocation.

        Serial local batches do not split ``num_cores`` or ``mem_gb`` across
        children.
        """
        return Job._propagate_runner(self.jobrunner, job)

    def _write_outcome_logs(self, outcomes: Sequence[dict[str, Any]]) -> None:
        """Write batch outcomes to success.log and failed.log."""
        success_path = os.path.join(self.folder, "success.log")
        failed_path = os.path.join(self.folder, "failed.log")

        successes = [item for item in outcomes if item["success"]]
        failures = [item for item in outcomes if not item["success"]]

        with open(success_path, "w") as fh:
            for item in successes:
                fh.write(f"{item['label']}\n")

        with open(failed_path, "w") as fh:
            for item in failures:
                error = item["error"] or "unknown error"
                fh.write(f"{item['label']}\t{error}\n")

    def _backup_files(self) -> None:
        pass

    def is_complete(self) -> bool:
        """Return True when all child jobs are complete."""
        if not self.jobs:
            return True
        return all(job.is_complete() for job in self.jobs)


def run_nestable_job(parent: Job, run_local: Callable[[], None]) -> None:
    """Run a nestable parent in single-child or local serial mode.

    Nestable parents (crest/QRC/dias/traj) call this from ``_run``.

    Selection order:

    1. Explicit ``parent.child_index`` (1-based), typically from
       ``--child-index`` rewritten into array runscripts
    2. Scheduler array task id (``SLURM_ARRAY_TASK_ID`` / ``PBS_ARRAYID`` /
       ``LSB_JOBINDEX``) as a legacy fallback
    3. Otherwise *run_local* runs the full serial nested workflow

    Selected children are built via ``parent.get_array_child_job`` so siblings
    are not constructed.
    """
    child_index = parent.child_index
    if child_index is not None:
        _run_one_nestable_child(parent, int(child_index))
        return
    task_id = resolve_array_task_id()
    if task_id is not None:
        _run_one_nestable_child(parent, task_id)
        return
    run_local()


def run_child_jobs_as_batch(
    *,
    batch_cls: Type[BatchJobT],
    jobs: Sequence[Job],
    parent: Job,
    label_suffix: str = "_batch",
    fail_fast: bool = False,
) -> BatchJobT:
    """Run sibling jobs through an engine ``BatchJob``.

    Nested children always run serially, each with the parent jobrunner's
    full ``num_cores`` and ``mem_gb``. Concurrent nested children are not
    used. Independent parent jobs may still run concurrently when submitted
    as a top-level batch via ``chemsmart sub``.

    Scheduler array / ``--child-index`` selection is handled at the nestable
    parent ``_run`` boundary via ``run_nestable_job``, not in this helper.

    ``fail_fast`` controls whether execution stops after the first
    unsuccessful child (default: run all children, then raise on failures).

    Returns the completed ``BatchJob`` instance.
    """
    runner = parent.jobrunner
    if runner is not None:
        cores = runner.num_cores
        mem_gb = runner.mem_gb
    else:
        cores = None
        mem_gb = None
    logger.info(
        "Running nested batch of %s child job(s) serially "
        "(full parent resources per child: cores=%s, mem_gb=%s).",
        len(jobs),
        cores,
        mem_gb,
    )
    batch_job = batch_cls(
        jobs=jobs,
        fail_fast=fail_fast,
        label=f"{parent.label}{label_suffix}",
        jobrunner=parent.jobrunner,
        nested_serial=True,
    )
    with cleared_array_task_env():
        batch_job.run()
    return batch_job


def _run_one_nestable_child(parent: Job, child_index: int) -> None:
    """Build and run one nestable child by 1-based *child_index*."""
    total = parent.num_array_children
    if total == 0:
        raise ValueError(
            f"Nestable job {parent.label!r} has no child jobs for array task."
        )
    zero_based = child_index - 1
    if zero_based < 0 or zero_based >= total:
        raise ValueError(
            f"Child index {child_index} out of range for {total} "
            f"nested child job(s) of {parent.label!r}; expected 1..{total}."
        )
    child = parent.get_array_child_job(zero_based)
    child_runner = Job._propagate_runner(parent.jobrunner, child)
    if child_runner is not None:
        cores = child_runner.num_cores
        mem_gb = child_runner.mem_gb
    else:
        cores = None
        mem_gb = None
    logger.info(
        "Nestable job %r: execution=%s, task=%s/%s, cores=%s, "
        "mem_gb=%s, child=%s",
        parent.label,
        BatchExecutionMode.ARRAY_TASK.value,
        child_index,
        total,
        cores,
        mem_gb,
        child.label,
    )
    child.run()
    if not child.is_complete():
        raise BatchExecutionError(
            f"Nestable array child {child.label!r} of {parent.label!r} "
            f"(task {child_index}/{total}) is incomplete after execution."
        )


def get_nestable_array_children(job: Any) -> Optional[list[Job]]:
    """Return nestable children for array submit, or ``None`` if not nestable."""
    try:
        get_children = job.get_array_child_jobs
    except AttributeError:
        return None
    children = list(get_children())
    if not children:
        return None
    return children


def prepare_nestable_batch_jobs(jobs: Sequence[Any]) -> RewriteCliFn:
    """Attach 1-based ``child_index`` entries for nestable array submit.

    Returns ``rewrite_nestable_cli_args`` for ``BatchJob.rewrite_cli``.
    """
    if not jobs:
        raise ValueError("Cannot prepare nestable batch jobs: empty job list.")
    entries = [
        {"child_index": task_id, "label": job.label}
        for task_id, job in enumerate(jobs, start=1)
    ]
    attach_batch_entries(jobs, entries)
    return rewrite_nestable_cli_args


def rewrite_nestable_cli_args(
    cli_args: Sequence[str],
    batch_entry: Optional[Mapping[str, Any]],
) -> list[str]:
    """Inject ``--child-index`` for a nestable array child."""
    if not batch_entry:
        return list(cli_args)

    args = list(cli_args)
    child_index = batch_entry.get("child_index")
    if child_index is not None:
        set_cli_option(
            args,
            long_opt="--child-index",
            value=str(child_index),
            insert_before=find_job_subcommand_token(args),
        )
    return args


# ---------------------------------------------------------------------------
# Submit-time batch_entry / per-task CLI helpers
# ---------------------------------------------------------------------------

_FILENAME_OPTIONS = frozenset({"-f", "--filename"})
_INDEX_OPTIONS = frozenset({"-i", "--index", "--si", "--structure-index"})
_PROGRAM_TOKENS = frozenset({"gaussian", "orca", "run", "sub"})


def get_job_batch_entry(job: Any) -> Optional[dict[str, Any]]:
    """Return ``job.batch_entry`` when it is a mapping, else ``None``."""
    try:
        entry = job.batch_entry
    except AttributeError:
        return None
    if not isinstance(entry, Mapping):
        return None
    return dict(entry)


def set_job_batch_entry(job: Any, entry: Mapping[str, Any]) -> None:
    """Attach an explicit batch-entry mapping on *job*."""
    job.batch_entry = dict(entry)


def attach_batch_entries(
    jobs: Sequence[Any],
    entries: Sequence[Mapping[str, Any]],
) -> None:
    """Attach one ``batch_entry`` mapping per child job."""
    if len(jobs) != len(entries):
        raise ValueError(
            f"Cannot attach batch entries: {len(jobs)} job(s) vs "
            f"{len(entries)} entr(y/ies)."
        )
    for job, entry in zip(jobs, entries):
        set_job_batch_entry(job, entry)


def prepare_batch_jobs(
    jobs: Sequence[Any],
    molecule_indices: Optional[Sequence[int]],
    *,
    filepath: Optional[str] = None,
) -> Optional[RewriteCliFn]:
    """Attach per-task entries that narrow shared multi-molecule CLI args.

    Used by homogeneous fan-out (opt/sp/ts/…). Each entry keeps the shared
    ``filepath`` and a single ``molecule_index`` for ``-i`` narrowing.

    Returns ``rewrite_batch_cli_args`` when entries were attached, else
    ``None``. Callers should pass the return value as ``BatchJob.rewrite_cli``.
    """
    if molecule_indices is None or len(jobs) <= 1:
        return None
    if len(jobs) != len(molecule_indices):
        raise ValueError(
            f"Cannot prepare batch jobs: {len(jobs)} job(s) vs "
            f"{len(molecule_indices)} molecule index(es)."
        )
    paired_jobs: list[Any] = []
    entries: list[dict[str, Any]] = []
    for job, index in zip(jobs, molecule_indices):
        entry: dict[str, Any] = {"molecule_index": int(index)}
        if filepath is not None:
            entry["filepath"] = str(filepath)
        paired_jobs.append(job)
        entries.append(entry)
    attach_batch_entries(paired_jobs, entries)
    return rewrite_batch_cli_args


def drop_cli_option(
    tokens: list[str], option_names: set[str] | frozenset[str]
) -> None:
    """Remove matching options and their following value tokens."""
    idx = 0
    while idx < len(tokens):
        if tokens[idx] in option_names:
            del tokens[idx]
            if idx < len(tokens):
                del tokens[idx]
            continue
        idx += 1


def set_cli_option(
    tokens: list[str],
    *,
    long_opt: str,
    value: str,
    short_opt: Optional[str] = None,
    insert_before: Optional[str] = None,
    prefer_short: bool = False,
) -> None:
    """Set or insert a long/short option pair to *value*."""
    if long_opt in tokens:
        pos = tokens.index(long_opt)
        if pos + 1 < len(tokens):
            tokens[pos + 1] = value
        return
    if short_opt is not None and short_opt in tokens:
        pos = tokens.index(short_opt)
        if pos + 1 < len(tokens):
            tokens[pos + 1] = value
        return

    insert_idx = len(tokens)
    if insert_before is not None and insert_before in tokens:
        insert_idx = tokens.index(insert_before)
    opt = short_opt if prefer_short and short_opt is not None else long_opt
    tokens[insert_idx:insert_idx] = [opt, value]


def set_cli_option_after(
    tokens: list[str],
    *,
    long_opt: str,
    value: str,
    short_opt: Optional[str] = None,
    insert_after: Optional[str] = None,
) -> None:
    """Replace an option pair, inserting after *insert_after* when absent."""
    names = {long_opt}
    if short_opt is not None:
        names.add(short_opt)
    drop_cli_option(tokens, names)
    insert_idx = len(tokens)
    if insert_after is not None and insert_after in tokens:
        insert_idx = tokens.index(insert_after) + 1
    tokens[insert_idx:insert_idx] = [long_opt, value]


def replace_filename_option(tokens: list[str], filepath: str) -> None:
    """Replace the value of ``-f``/``--filename`` with *filepath*."""
    for idx in range(len(tokens) - 1):
        if tokens[idx] in _FILENAME_OPTIONS:
            tokens[idx + 1] = str(filepath)
            return


def find_job_subcommand_token(tokens: Sequence[str]) -> Optional[str]:
    """Return the job subcommand token (e.g. ``opt``/``pka``), not a path."""
    for token in reversed(tokens):
        if token.startswith("-"):
            continue
        if token in _PROGRAM_TOKENS:
            continue
        return token
    return None


def _resolve_index_value(batch_entry: Mapping[str, Any]) -> Optional[str]:
    for key in ("index", "molecule_index", "fragment_index"):
        if key in batch_entry and batch_entry[key] is not None:
            return str(batch_entry[key])
    return None


def apply_shared_batch_cli_options(
    args: list[str],
    batch_entry: Mapping[str, Any],
    *,
    insert_before: Optional[str] = None,
) -> None:
    """Apply CLI options shared by homogeneous and heterogeneous batch entries."""
    filepath = batch_entry.get("filepath")
    if filepath is not None:
        replace_filename_option(args, str(filepath))

    index_value = _resolve_index_value(batch_entry)
    if index_value is not None:
        prefer_short = any(token in {"-i", "--si"} for token in args)
        drop_cli_option(args, _INDEX_OPTIONS)
        set_cli_option(
            args,
            long_opt="--index",
            short_opt="-i",
            value=index_value,
            insert_before=insert_before,
            prefer_short=prefer_short,
        )

    label = batch_entry.get("label")
    if label is not None:
        set_cli_option(
            args,
            long_opt="--label",
            short_opt="-l",
            value=str(label),
            insert_before=insert_before,
        )

    if "charge" in batch_entry and batch_entry["charge"] is not None:
        set_cli_option(
            args,
            long_opt="--charge",
            short_opt="-c",
            value=str(batch_entry["charge"]),
            insert_before=insert_before,
        )

    if (
        "multiplicity" in batch_entry
        and batch_entry["multiplicity"] is not None
    ):
        set_cli_option(
            args,
            long_opt="--multiplicity",
            short_opt="-m",
            value=str(batch_entry["multiplicity"]),
            insert_before=insert_before,
        )


def rewrite_batch_cli_args(
    cli_args: Sequence[str],
    batch_entry: Optional[Mapping[str, Any]],
) -> list[str]:
    """Rewrite shared submit CLI args from a child ``batch_entry``.

    Applies shared fields (filepath, index, label, charge, multiplicity).
    Domain-specific overlays (e.g. pKa ``batch`` → ``submit``) belong in
    their own rewriter that calls this first.
    """
    if not batch_entry:
        return list(cli_args)

    args = list(cli_args)
    apply_shared_batch_cli_options(
        args,
        batch_entry,
        insert_before=find_job_subcommand_token(args),
    )
    return args


def resolve_array_cli_args(
    jobs: Sequence[Any],
    shared_cli_args: Sequence[str],
    rewrite_cli: Optional[RewriteCliFn] = None,
) -> list[str] | list[list[str]]:
    """Return per-task CLI args when children carry ``batch_entry``.

    Batches without ``batch_entry`` keep a single shared CLI list.
    Batches with ``batch_entry`` require *rewrite_cli* and return one CLI
    list per child for array runscripts.
    """
    entries = [get_job_batch_entry(job) for job in jobs]
    if not any(entry is not None for entry in entries):
        return list(shared_cli_args)
    if rewrite_cli is None:
        raise ValueError(
            "BatchJob children have batch_entry but no "
            "rewrite_cli callback was provided for per-task CLI args."
        )
    return [rewrite_cli(shared_cli_args, entry) for entry in entries]
