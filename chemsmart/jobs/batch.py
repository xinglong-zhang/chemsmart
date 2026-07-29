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
"""

import logging
import os
from contextlib import contextmanager
from enum import Enum
from functools import partial
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
        treated as 1-based indexes into ``self.jobs``. Array-task environment
        variables are cleared around ``child.run()`` so nested serial
        ``BatchJob.run()`` stays in ``local_batch`` mode.
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
        """Propagate the batch runner onto *job* with full resources."""
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
    """Run sibling jobs serially through an engine ``BatchJob``.

    Each child uses the parent jobrunner's full ``num_cores`` and ``mem_gb``.
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


def prepare_nestable_batch_jobs(
    jobs: Sequence[Any],
    *,
    job_token: str,
) -> RewriteCliFn:
    """Attach 1-based ``child_index`` entries for nestable array submit.

    Each array task rebuilds the nestable parent and selects one child via
    ``--child-index``. *job_token* is the engine job-group token of the
    nestable command (e.g. ``qrc``/``dias``), used to place the rewritten
    options. Returns the rewriter for ``BatchJob.rewrite_cli``.
    """
    if not jobs:
        raise ValueError("Cannot prepare nestable batch jobs: empty job list.")
    entries = [
        {"child_index": task_id} for task_id, _ in enumerate(jobs, start=1)
    ]
    attach_batch_entries(jobs, entries)
    return make_batch_cli_rewriter(job_token)


# ---------------------------------------------------------------------------
# Submit-time batch_entry / per-task CLI helpers
# ---------------------------------------------------------------------------

_FILENAME_OPTIONS = frozenset({"-f", "--filename"})
_INDEX_OPTIONS = frozenset({"-i", "--index", "--si", "--structure-index"})

# Known ``batch_entry`` keys mapped to CLI flags.
# Engine-group options insert under ``gaussian``/``orca`` (before the job
# token). Nestable-command options insert after that token (``dias``/``qrc``/…).
_BATCH_ENTRY_ENGINE_OPTION_FIELDS = (("label", "--label", "-l"),)
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
    job_token: str,
    filepath: Optional[str] = None,
    labels: Optional[Sequence[str]] = None,
) -> Optional[RewriteCliFn]:
    """Attach per-task entries that narrow shared multi-molecule CLI args.

    Used by homogeneous fan-out (opt/sp/ts/…). Each entry keeps the shared
    ``filepath``, a single ``molecule_index`` for ``-i`` narrowing, and the
    ``--label`` value to replay in the per-task array run script.

    Args:
        jobs: Child jobs of the batch, one per molecule index.
        molecule_indices: 1-based molecule indices of the children.
        job_token: Engine job-group token of the invoked command (e.g.
            ``opt``/``sp``), used to place the rewritten options.
        filepath: Shared input file replayed as ``-f``/``--filename``.
        labels: ``--label`` value to replay per task. Required when the
            command derives the child label from settings (e.g. single-point
            solvent suffixes); it must be the label *before* that derivation
            so replaying the task CLI reproduces the child label exactly.
            Defaults to ``job.label``.

    Returns:
        The bound per-task CLI rewriter when entries were attached, else
        ``None``. Callers pass the return value as ``BatchJob.rewrite_cli``.
    """
    if molecule_indices is None or len(jobs) <= 1:
        return None
    if len(jobs) != len(molecule_indices):
        raise ValueError(
            f"Cannot prepare batch jobs: {len(jobs)} job(s) vs "
            f"{len(molecule_indices)} molecule index(es)."
        )
    if labels is not None and len(jobs) != len(labels):
        raise ValueError(
            f"Cannot prepare batch jobs: {len(jobs)} job(s) vs "
            f"{len(labels)} label(s)."
        )
    replay_labels = [job.label for job in jobs] if labels is None else labels
    paired_jobs: list[Any] = []
    entries: list[dict[str, Any]] = []
    for job, index, replay_label in zip(jobs, molecule_indices, replay_labels):
        entry: dict[str, Any] = {"molecule_index": int(index)}
        if filepath is not None:
            entry["filepath"] = str(filepath)
        entry["label"] = str(replay_label)
        paired_jobs.append(job)
        entries.append(entry)
    attach_batch_entries(paired_jobs, entries)
    return make_batch_cli_rewriter(job_token)


# Boolean long flags in reconstructed CLI appear alone (no following value).
# Needed so ``--forces opt`` still treats ``opt`` as the job-group token, while
# ``--label opt`` / ``-l opt`` / ``-a opt`` / ``-r opt`` treat ``opt`` as a value.
_BOOLEAN_LONG_OPTIONS = frozenset(
    {
        "--run-in-parallel",
        "--fake",
        "--debug",
        "--stream",
        "--test",
        "--forces",
        "--remove-solvent",
        "--skip-completed",
        "--scratch",
        "--delete-scratch",
    }
)


def _find_job_token_index(tokens: Sequence[str], job_token: str) -> int:
    """Return the job-group command index, skipping option-value collisions.

    The first bare ``job_token`` wins. A match is skipped when it is the value
    of the preceding option (``-l/-a/-r/--label opt … opt``).
    """
    for idx, token in enumerate(tokens):
        if token != job_token:
            continue
        prev = tokens[idx - 1] if idx else ""
        if prev.startswith("-") and not prev.startswith("--no-"):
            if prev.startswith("--"):
                if prev not in _BOOLEAN_LONG_OPTIONS:
                    continue
            elif prev not in {"-R", "-S"}:
                continue
        return idx
    raise ValueError(
        f"Cannot locate job token {job_token!r} in CLI args {list(tokens)!r}."
    )


def patch_cli_option(
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
    if insert_after is not None:
        insert_idx = _find_job_token_index(tokens, insert_after) + 1
    elif insert_before is not None:
        insert_idx = _find_job_token_index(tokens, insert_before)

    opt = short_opt if prefer_short and short_opt is not None else long_opt
    tokens[insert_idx:insert_idx] = [opt, value]


def _resolve_index_value(batch_entry: Mapping[str, Any]) -> Optional[str]:
    for key in ("index", "molecule_index", "fragment_index"):
        if key in batch_entry and batch_entry[key] is not None:
            return str(batch_entry[key])
    return None


def _apply_batch_entry_to_cli(
    args: list[str],
    batch_entry: Mapping[str, Any],
    *,
    job_token: str,
) -> None:
    """Map known ``batch_entry`` fields onto shared submit CLI tokens."""
    filepath = batch_entry.get("filepath")
    if filepath is not None:
        patch_cli_option(
            args,
            long_opt="--filename",
            short_opt="-f",
            value=str(filepath),
            insert_before=job_token,
            prefer_short=any(token in _FILENAME_OPTIONS for token in args),
        )

    index_value = _resolve_index_value(batch_entry)
    if index_value is not None:
        prefer_short = any(token in {"-i", "--si"} for token in args)
        patch_cli_option(
            args,
            long_opt="--index",
            short_opt="-i",
            value=index_value,
            drop=_INDEX_OPTIONS,
            insert_before=job_token,
            prefer_short=prefer_short,
        )

    for entry_key, long_opt, short_opt in _BATCH_ENTRY_ENGINE_OPTION_FIELDS:
        if entry_key not in batch_entry or batch_entry[entry_key] is None:
            continue
        patch_cli_option(
            args,
            long_opt=long_opt,
            short_opt=short_opt,
            value=str(batch_entry[entry_key]),
            insert_before=job_token,
        )

    # Nestable options insert after the job-group token (``dias``/``qrc``/…).
    for entry_key, long_opt, short_opt in _BATCH_ENTRY_NESTABLE_OPTION_FIELDS:
        if entry_key not in batch_entry or batch_entry[entry_key] is None:
            continue
        patch_cli_option(
            args,
            long_opt=long_opt,
            short_opt=short_opt,
            value=str(batch_entry[entry_key]),
            insert_after=job_token,
        )


def rewrite_batch_cli_args(
    cli_args: Sequence[str],
    batch_entry: Optional[Mapping[str, Any]],
    *,
    job_token: str,
) -> list[str]:
    """Rewrite shared submit CLI args from a child ``batch_entry``.

    Applies shared ``batch_entry`` fields (filepath, index, label,
    child_index). Engine-group options are placed before *job_token*, the
    engine job-group token of the invoked command (e.g. ``opt``/``qrc``), and
    nestable options right after it. Domain-specific overlays (e.g. pKa
    charge/multiplicity and ``batch`` → ``submit``) belong in their own
    rewriter that calls this first.

    Raises:
        ValueError: If *job_token* is absent from *cli_args*, which would
            place the rewritten options under the wrong command.
    """
    if not batch_entry:
        return list(cli_args)

    args = list(cli_args)
    if job_token not in args:
        raise ValueError(
            f"Cannot rewrite batch CLI args: job token {job_token!r} is "
            f"not present in {args!r}."
        )
    _apply_batch_entry_to_cli(args, batch_entry, job_token=job_token)
    return args


def make_batch_cli_rewriter(job_token: str) -> RewriteCliFn:
    """Return ``rewrite_batch_cli_args`` bound to *job_token*."""
    return partial(rewrite_batch_cli_args, job_token=job_token)


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
