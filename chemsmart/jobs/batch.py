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
from contextlib import contextmanager
from enum import Enum
from typing import (
    Any,
    Callable,
    Collection,
    Iterator,
    Mapping,
    Optional,
    Sequence,
    Type,
    TypeVar,
)

from chemsmart.jobs.job import Job

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

    Used so nested serial ``BatchJob.run()`` stays in ``local_batch`` mode when
    an outer scheduler array task id is present.
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


class BatchJob(Job):
    """
    Abstract controller for running a batch of child jobs.

    Orchestration is engine-agnostic. Subclasses set ``PROGRAM`` and hold
    the child job list.

    ``run()`` selects execution mode from the environment:

    - ``array_task`` — one child at the 1-based scheduler array task id,
      with full resources
    - ``local_batch`` — all children serially with full resources

    ``rewrite_cli`` is the optional per-task CLI rewriter used by
    ``chemsmart sub`` when submitting this batch as a scheduler array.
    """

    PROGRAM: Optional[str] = None
    REGISTERABLE: bool = False

    def __init__(
        self,
        jobs: Optional[Sequence[Job]],
        write_outcome_logs: bool = False,
        label: str = "batch_job",
        jobrunner: Any = None,
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
        self.write_outcome_logs = bool(write_outcome_logs)
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

        Array-task env vars are cleared around ``child.run()`` so nested
        serial batches and legacy nestable env fallback do not reuse the
        outer task id.
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
        child_index = task_id - 1  # use 1-based task id
        if child_index < 0 or child_index >= total_jobs:
            raise ValueError(
                f"Array task id {task_id} out of range for {total_jobs} "
                f"child job(s); expected 1..{total_jobs}."
            )

        child = self.jobs[child_index]

        logger.info(
            "BatchJob %r: execution=%s, task=%s/%s, cores=%s, "
            "mem_gb=%s, child=%s",
            self.label,
            BatchExecutionMode.ARRAY_TASK.value,
            task_id,
            total_jobs,
            *self._runner_resources(),
            child.label,
        )
        # Keep nested BatchJob.run() / legacy nestable fallback off this task id.
        with cleared_array_task_env():
            outcome = self._submit_job(child, **kwargs)
        self._finish_batch_run([outcome])

    def _run_local_batch(self, **kwargs: Any) -> None:
        """Run all children serially with full resources per child."""
        total_jobs = len(self.jobs)
        logger.info(
            "BatchJob %r: execution=%s, children=%s, policy=serial, "
            "cores=%s, mem_gb=%s",
            self.label,
            BatchExecutionMode.LOCAL_BATCH.value,
            total_jobs,
            *self._runner_resources(),
        )
        outcomes = [self._submit_job(job, **kwargs) for job in self.jobs]
        self._finish_batch_run(outcomes)

    def _runner_resources(self) -> tuple[Optional[Any], Optional[Any]]:
        """Return ``(num_cores, mem_gb)`` from the batch runner, if set."""
        runner = self.jobrunner
        if runner is None:
            return None, None
        return runner.num_cores, runner.mem_gb

    def _finish_batch_run(
        self,
        outcomes: Sequence[dict[str, Any]],
    ) -> None:
        """Store outcomes, optionally log them, and raise on failures."""
        self._last_batch_outcomes = list(outcomes)
        if self.write_outcome_logs:
            self._write_outcome_logs(outcomes)
        self._raise_if_failures(outcomes)

    def _raise_if_failures(self, outcomes: Sequence[dict[str, Any]]) -> None:
        """Raise ``BatchExecutionError`` summarizing failed outcomes."""
        failures = [item for item in outcomes if not item["success"]]
        if not failures:
            return

        lines = [f"- {item['label']}: {item['error']}" for item in failures]
        summary = f"{len(failures)} of {len(outcomes)} batch job(s) failed"
        raise BatchExecutionError(summary + ":\n" + "\n".join(lines))

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

    When ``parent.child_index`` is set (1-based, typically from ``--child-index``
    rewritten into array runscripts), only that child runs. Otherwise *run_local*
    runs the full serial nested workflow.

    Selected children are built via ``parent.get_array_child_job`` so siblings
    are not constructed.
    """
    child_index = parent.child_index
    if child_index is not None:
        _run_one_nestable_child(parent, int(child_index))
        return
    run_local()


def run_child_jobs_as_batch(
    *,
    batch_cls: Type[BatchJobT],
    jobs: Sequence[Job],
    parent: Job,
    label_suffix: str = "_batch",
) -> BatchJobT:
    """Run sibling jobs through an engine ``BatchJob``.

    Nested children always run serially, each with the parent jobrunner's
    full ``num_cores`` and ``mem_gb``. Concurrent nested children are not
    used. Independent parent jobs may still run concurrently when submitted
    as a top-level batch via ``chemsmart sub``.

    Scheduler array / ``--child-index`` selection is handled at the nestable
    parent ``_run`` boundary via ``run_nestable_job``, not in this helper.

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
        label=f"{parent.label}{label_suffix}",
        jobrunner=parent.jobrunner,
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

    Returns ``rewrite_batch_cli_args`` for ``BatchJob.rewrite_cli``.
    """
    if not jobs:
        raise ValueError("Cannot prepare nestable batch jobs: empty job list.")
    entries = [
        {"child_index": task_id, "label": job.label}
        for task_id, job in enumerate(jobs, start=1)
    ]
    attach_batch_entries(jobs, entries)
    return rewrite_batch_cli_args


# ---------------------------------------------------------------------------
# Submit-time batch_entry / per-task CLI helpers
# ---------------------------------------------------------------------------

_FILENAME_OPTIONS = frozenset({"-f", "--filename"})
_INDEX_OPTIONS = frozenset({"-i", "--index", "--si", "--structure-index"})
# Engine job-group tokens (gaussian/orca subcommands). Nested actions such as
# ``batch``/``submit`` are intentionally excluded so shared options insert
# under the engine group, not under the nested command.
_ENGINE_JOB_TOKENS = frozenset(
    {
        "com",
        "crest",
        "dias",
        "inp",
        "irc",
        "link",
        "modred",
        "nci",
        "neb",
        "opt",
        "pka",
        "qrc",
        "resp",
        "scan",
        "sp",
        "td",
        "traj",
        "ts",
        "userjob",
        "wbi",
    }
)

# Known ``batch_entry`` keys mapped to CLI flags.
# Engine-group options insert under ``gaussian``/``orca`` (before the job
# token). Nestable-command options insert after that token (``dias``/``qrc``/…).
_BATCH_ENTRY_ENGINE_OPTION_FIELDS = (("label", "--label", None),)
_BATCH_ENTRY_NESTABLE_OPTION_FIELDS = (("child_index", "--child-index", None),)


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
    ``filepath``, a single ``molecule_index`` for ``-i`` narrowing, and the
    child ``label`` for per-task array run scripts.

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
        job_label = getattr(job, "label", None)
        if job_label is not None:
            entry["label"] = job_label
        paired_jobs.append(job)
        entries.append(entry)
    attach_batch_entries(paired_jobs, entries)
    return rewrite_batch_cli_args


def _patch_cli_option(
    tokens: list[str],
    *,
    long_opt: Optional[str] = None,
    short_opt: Optional[str] = None,
    value: Optional[str] = None,
    drop: Optional[Collection[str]] = None,
    insert_before: Optional[str] = None,
    insert_after: Optional[str] = None,
    prefer_short: bool = False,
) -> None:
    """Remove, update, or insert one CLI option pair in *tokens*."""
    if drop:
        idx = 0
        drop_names = set(drop)
        while idx < len(tokens):
            if tokens[idx] in drop_names:
                del tokens[idx]
                if idx < len(tokens):
                    del tokens[idx]
                continue
            idx += 1

    if value is None:
        return
    if long_opt is None:
        raise ValueError("long_opt is required when value is set")

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
    if insert_after is not None and insert_after in tokens:
        insert_idx = tokens.index(insert_after) + 1
    elif insert_before is not None and insert_before in tokens:
        insert_idx = tokens.index(insert_before)

    opt = short_opt if prefer_short and short_opt is not None else long_opt
    tokens[insert_idx:insert_idx] = [opt, value]


def _find_job_subcommand_token(tokens: Sequence[str]) -> Optional[str]:
    """Return the engine job-group token (e.g. ``opt``/``pka``).

    Matches known gaussian/orca job groups so nested actions
    (``batch``/``submit``) and option values (e.g. ``direct``) are not
    treated as insertion anchors. Shared engine options then insert under
    ``gaussian``/``orca``.
    """
    for token in tokens:
        if token in _ENGINE_JOB_TOKENS:
            return token
    return None


def _resolve_index_value(batch_entry: Mapping[str, Any]) -> Optional[str]:
    for key in ("index", "molecule_index", "fragment_index"):
        if key in batch_entry and batch_entry[key] is not None:
            return str(batch_entry[key])
    return None


def _apply_batch_entry_to_cli(
    args: list[str],
    batch_entry: Mapping[str, Any],
    *,
    insert_before: Optional[str] = None,
) -> None:
    """Map known ``batch_entry`` fields onto shared submit CLI tokens."""
    filepath = batch_entry.get("filepath")
    if filepath is not None:
        _patch_cli_option(
            args,
            long_opt="--filename",
            short_opt="-f",
            value=str(filepath),
            insert_before=insert_before,
            prefer_short=any(token in _FILENAME_OPTIONS for token in args),
        )

    index_value = _resolve_index_value(batch_entry)
    if index_value is not None:
        prefer_short = any(token in {"-i", "--si"} for token in args)
        _patch_cli_option(
            args,
            long_opt="--index",
            short_opt="-i",
            value=index_value,
            drop=_INDEX_OPTIONS,
            insert_before=insert_before,
            prefer_short=prefer_short,
        )

    for entry_key, long_opt, short_opt in _BATCH_ENTRY_ENGINE_OPTION_FIELDS:
        if entry_key not in batch_entry or batch_entry[entry_key] is None:
            continue
        _patch_cli_option(
            args,
            long_opt=long_opt,
            short_opt=short_opt,
            value=str(batch_entry[entry_key]),
            insert_before=insert_before,
        )

    # Nestable flags belong on the job command (after ``dias``/``qrc``/…),
    # not on the engine group where Click would reject them.
    for entry_key, long_opt, short_opt in _BATCH_ENTRY_NESTABLE_OPTION_FIELDS:
        if entry_key not in batch_entry or batch_entry[entry_key] is None:
            continue
        _patch_cli_option(
            args,
            long_opt=long_opt,
            short_opt=short_opt,
            value=str(batch_entry[entry_key]),
            insert_after=insert_before,
        )


def rewrite_batch_cli_args(
    cli_args: Sequence[str],
    batch_entry: Optional[Mapping[str, Any]],
) -> list[str]:
    """Rewrite shared submit CLI args from a child ``batch_entry``.

    Applies shared ``batch_entry`` fields (filepath, index, label,
    child_index). Domain-specific overlays (e.g. pKa charge/multiplicity
    and ``batch`` → ``submit``) belong in their own rewriter that calls
    this first.
    """
    if not batch_entry:
        return list(cli_args)

    args = list(cli_args)
    _apply_batch_entry_to_cli(
        args,
        batch_entry,
        insert_before=_find_job_subcommand_token(args),
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
