"""
Shared batch job infrastructure.

Provides the abstract ``BatchJob`` base class for orchestrating collections
of engine-specific jobs. Engine-agnostic behavior includes:

- serial and parallel child-job scheduling
- fault-tolerant execution with aggregated failures
- scheduler/node allocation detection (SLURM/PBS)
- distribution of child jobs across allocated nodes

Subclasses implement engine-specific runner adaptation (for example,
pinning a copied runner to a scheduler node).
"""

import logging
import os
import subprocess
import threading
import time
import types
from abc import ABCMeta, abstractmethod
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import suppress
from typing import Any, Optional, Sequence, Type, TypeVar

from chemsmart.jobs.job import Job
from chemsmart.jobs.runner import get_serial_mode, get_submitter_worker_count
from chemsmart.utils.mixins import RegistryMeta

logger = logging.getLogger(__name__)

BatchJobT = TypeVar("BatchJobT", bound="BatchJob")


class BatchExecutionError(RuntimeError):
    """Raised when one or more child jobs fail in a batch run."""


class BatchJobMeta(RegistryMeta, ABCMeta):
    """Metaclass combining registry support with abstract-base semantics."""


class BatchJob(Job, metaclass=BatchJobMeta):
    """
    Abstract controller for running a batch of child jobs.

    Orchestration is engine-agnostic. Subclasses implement
    ``_configure_runner_for_node`` for engine-specific node adaptation.

    ``no_run_in_parallel`` selects serial vs parallel child scheduling.
    ``fail_fast`` stops serial submission after the first unsuccessful
    outcome (ignored in parallel mode).
    """

    PROGRAM: Optional[str] = None
    REGISTERABLE: bool = False

    def __init__(
        self,
        jobs: Optional[Sequence[Job]],
        no_run_in_parallel: Optional[bool] = None,
        fail_fast: bool = False,
        write_outcome_logs: bool = False,
        label: str = "batch_job",
        jobrunner: Any = None,
        **kwargs,
    ) -> None:
        super().__init__(
            molecule=None,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
        self.jobs: list[Job] = list(jobs) if jobs is not None else []
        runner_serial_mode = get_serial_mode(jobrunner)
        # - explicit True/False at call site wins
        # - None defers to runner policy
        if no_run_in_parallel is None:
            self.no_run_in_parallel = runner_serial_mode.no_run_in_parallel
        else:
            self.no_run_in_parallel = bool(no_run_in_parallel)
        self.fail_fast = bool(fail_fast)
        self.write_outcome_logs = bool(write_outcome_logs)

        # Cache completion checks to avoid repeatedly reparsing output files
        # from the head-node monitoring loop.
        self._status_cache: dict[
            int,
            tuple[Optional[tuple[tuple[str, float], ...]], float, bool],
        ] = {}
        self._status_cache_ttl_seconds: float = 2.0
        self._status_cache_lock = threading.Lock()
        self._last_batch_outcomes: list[dict[str, Any]] = []
        self._jobs_not_started: int = 0

    def run(self, **kwargs: Any) -> None:
        """Run all child jobs in this batch."""
        self._invalidate_status_cache()
        self._jobs_not_started = 0
        self._run(**kwargs)

    def enable_serial_local_execution(self) -> None:
        """Configure this batch for serial local execution with full resources.

        Children run one at a time. Each child receives the batch jobrunner's
        full ``num_cores`` and ``mem_gb`` (no core/memory splitting).
        """
        if not self.no_run_in_parallel:
            logger.info(
                "BatchJob local execution is serial with full resources; "
                "concurrent children are disabled. Use chemsmart sub for "
                "cluster concurrency."
            )
        self.no_run_in_parallel = True

    def _run(self, **kwargs: Any) -> None:
        """Dispatch execution to multi-node, serial, or parallel mode."""
        nodes = self._get_allocated_nodes()
        total_jobs = len(self.jobs)

        if nodes and len(nodes) > 1:
            outcomes = self._run_multi_node(nodes, **kwargs)
        elif self.no_run_in_parallel:
            runner = self.jobrunner
            if runner is not None:
                cores = runner.num_cores
                mem_gb = runner.mem_gb
            else:
                cores = None
                mem_gb = None
            logger.info(
                "Running batch of %s jobs serially "
                "(full resources per child: cores=%s, mem_gb=%s).",
                total_jobs,
                cores,
                mem_gb,
            )
            outcomes = self._run_jobs_serially(self.jobs, **kwargs)
        else:
            if self.fail_fast:
                logger.info(
                    "fail_fast is set but ignored in parallel mode; "
                    "all child jobs will be submitted."
                )
            logger.info(f"Running batch of {total_jobs} jobs in parallel.")
            outcomes = self._run_jobs_in_parallel(self.jobs, **kwargs)

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
        node: Optional[str] = None,
        **kwargs: Any,
    ) -> list[dict[str, Any]]:
        """Submit child jobs one-by-one.

        When ``fail_fast`` is enabled, stop after the first unsuccessful
        outcome. Jobs not started are omitted from the returned outcomes and
        counted via ``_jobs_not_started``.
        """
        outcomes: list[dict[str, Any]] = []
        for index, job in enumerate(jobs):
            outcome = self._submit_job(job, node=node, **kwargs)
            outcomes.append(outcome)
            if self.fail_fast and not outcome["success"]:
                remaining = len(jobs) - (index + 1)
                self._jobs_not_started += remaining
                if remaining:
                    logger.warning(
                        "fail_fast enabled: stopping serial batch after "
                        f"{job.label} failed; {remaining} job(s) not started."
                    )
                break
        return outcomes

    def _run_jobs_in_parallel(
        self,
        jobs: Sequence[Job],
        **kwargs: Any,
    ) -> list[dict[str, Any]]:
        """Submit child jobs concurrently using a thread pool."""
        outcomes: list[dict[str, Any]] = []
        max_workers = get_submitter_worker_count(self.jobrunner, len(jobs))
        logger.info(f"Using up to {max_workers} parallel submitter workers")
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_job = {
                executor.submit(self._submit_job, job, None, **kwargs): job
                for job in jobs
            }

            for future in as_completed(future_to_job):
                job = future_to_job[future]
                try:
                    outcomes.append(future.result())
                except Exception as e:
                    logger.error(
                        f"Unexpected failure while coordinating job {job.label}: {e}",
                        exc_info=True,
                    )
                    outcomes.append(
                        {
                            "label": job.label,
                            "success": False,
                            "error": str(e),
                            "node": None,
                        }
                    )
        return outcomes

    def _submit_job(
        self,
        job: Job,
        node: Optional[str] = None,
        **kwargs: Any,
    ) -> dict[str, Any]:
        """Configure a runner copy and execute a single child job."""
        try:
            self._build_jobrunner(job, node=node)
            job.run(**kwargs)
            is_complete = self._job_is_complete_cached(
                job,
                force_refresh=True,
            )
            if not is_complete:
                location = f" on node {node}" if node else ""
                msg = "job incomplete after execution"
                logger.error(f"Job {job.label}{location}: {msg}")
                return {
                    "label": job.label,
                    "success": False,
                    "error": msg,
                    "node": node,
                }
            return {
                "label": job.label,
                "success": True,
                "error": "",
                "node": node,
            }
        except Exception as e:
            location = f" on node {node}" if node else ""
            logger.error(
                f"Job {job.label} failed during batch execution{location}: {e}",
                exc_info=True,
            )
            self._invalidate_status_cache(job)
            return {
                "label": job.label,
                "success": False,
                "error": str(e),
                "node": node,
            }

    def _build_jobrunner(
        self,
        job: Job,
        node: Optional[str] = None,
    ) -> Any:
        """Copy the batch runner onto *job* with full resource allocation.

        Serial local batches do not split ``num_cores`` or ``mem_gb`` across
        children. Optional *node* adaptation is applied after the copy.
        """
        child_runner = Job._propagate_runner(self.jobrunner, job)
        if child_runner is not None and node is not None:
            child_runner = self._configure_runner_for_node(
                runner=child_runner,
                node=node,
                job=job,
            )
            job.jobrunner = child_runner
        return child_runner

    def _get_allocated_nodes(self) -> Optional[list[str]]:
        """
        Detect allocated compute nodes from environment variables.

        Supports SLURM and PBS and returns a list of unique node names, or
        ``None`` when no multi-node allocation is detected.
        """
        nodelist = os.environ.get("SLURM_JOB_NODELIST")
        if nodelist:
            try:
                output = subprocess.check_output(
                    ["scontrol", "show", "hostnames", nodelist],
                    universal_newlines=True,
                )
                nodes = [
                    node.strip()
                    for node in output.strip().split("\n")
                    if node.strip()
                ]
                return nodes or None
            except (subprocess.SubprocessError, FileNotFoundError):
                return [
                    node.strip()
                    for node in nodelist.split(",")
                    if node.strip()
                ] or None

        nodefile = os.environ.get("PBS_NODEFILE")
        if nodefile and os.path.exists(nodefile):
            with open(nodefile, "r") as f:
                nodes = [line.strip() for line in f if line.strip()]
            return sorted(set(nodes)) or None

        return None

    def _run_multi_node(
        self, nodes: Sequence[str], **kwargs: Any
    ) -> list[dict[str, Any]]:
        """Distribute and execute child jobs across scheduler-allocated nodes."""
        num_nodes = len(nodes)
        logger.info(
            f"Distributing {len(self.jobs)} jobs across {num_nodes} nodes: {nodes}"
        )

        job_chunks = self._split_jobs_across_nodes(self.jobs, num_nodes)

        outcomes: list[dict[str, Any]] = []
        max_workers = get_submitter_worker_count(self.jobrunner, num_nodes)
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_node_chunk: dict[Any, tuple[str, Sequence[Job]]] = {}
            for node, jobs_chunk in zip(nodes, job_chunks):
                if jobs_chunk:
                    future = executor.submit(
                        self._run_chunk_on_node,
                        jobs_chunk,
                        node,
                        **kwargs,
                    )
                    future_to_node_chunk[future] = (node, jobs_chunk)

            for future in as_completed(future_to_node_chunk):
                node, jobs_chunk = future_to_node_chunk[future]
                try:
                    outcomes.extend(future.result())
                except Exception as e:
                    logger.error(f"Node execution failed: {e}", exc_info=True)
                    first_job_label = (
                        jobs_chunk[0].label if jobs_chunk else "unknown_chunk"
                    )
                    outcomes.append(
                        {
                            "label": f"node:{node}:{first_job_label}",
                            "success": False,
                            "error": str(e),
                            "node": node,
                        }
                    )
        return outcomes

    def _run_chunk_on_node(
        self,
        jobs: Sequence[Job],
        node: str,
        **kwargs: Any,
    ) -> list[dict[str, Any]]:
        """Run one chunk of jobs serially on a single node."""
        logger.info(f"Node {node} starting to process {len(jobs)} jobs.")
        outcomes = self._run_jobs_serially(jobs, node=node, **kwargs)
        logger.info(f"Node {node} finished processing.")
        return outcomes

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

    @staticmethod
    def _split_jobs_across_nodes(
        jobs: Sequence[Job],
        num_nodes: int,
    ) -> list[list[Job]]:
        """Split jobs into near-even chunks while preserving input order."""
        if num_nodes <= 0:
            return []

        base_size, remainder = divmod(len(jobs), num_nodes)
        chunks = []
        start = 0
        for idx in range(num_nodes):
            stop = start + base_size + (1 if idx < remainder else 0)
            chunks.append(jobs[start:stop])
            start = stop
        return chunks

    @abstractmethod
    def _configure_runner_for_node(
        self,
        runner: Any,
        node: str,
        job: Job,
    ) -> Any:
        """
        Adapt a copied jobrunner for execution on a specific node.

        Subclasses can patch engine-specific command generation here.
        """
        raise NotImplementedError

    def _backup_files(self) -> None:
        pass

    @staticmethod
    def _append_str_path(candidates: list[str], path: Any) -> None:
        """Append *path* to *candidates* when it is a non-empty string."""
        if isinstance(path, str) and path:
            candidates.append(path)

    def _job_status_signature(
        self,
        job: Job,
    ) -> Optional[tuple[tuple[str, float], ...]]:
        """Return file-mtime signature used to validate cached status."""
        candidates: list[str] = []
        try:
            self._append_str_path(candidates, job.outputfile)
        except AttributeError:
            pass
        try:
            self._append_str_path(candidates, job.joblog)
        except AttributeError:
            pass
        try:
            self._append_str_path(candidates, job.errfile)
        except AttributeError:
            pass

        signature: list[tuple[str, float]] = []
        for path in candidates:
            if os.path.exists(path):
                with suppress(OSError):
                    signature.append((path, os.path.getmtime(path)))

        if not signature:
            return None
        signature.sort(key=lambda item: item[0])
        return tuple(signature)

    def _invalidate_status_cache(self, job: Optional[Job] = None) -> None:
        """Invalidate status cache for all jobs or a single job."""
        with self._status_cache_lock:
            if job is None:
                self._status_cache.clear()
                return
            self._status_cache.pop(id(job), None)

    def _job_is_complete_cached(
        self,
        job: Job,
        *,
        force_refresh: bool = False,
    ) -> bool:
        """Return child completion status using a short-lived cache."""
        now = time.monotonic()
        signature = self._job_status_signature(job)
        key = id(job)

        if not force_refresh:
            with self._status_cache_lock:
                cached = self._status_cache.get(key)
            if cached is not None:
                cached_signature, cached_at, cached_value = cached
                signature_match = (
                    signature is not None and cached_signature == signature
                )
                ttl_match = (
                    signature is None
                    and (now - cached_at) <= self._status_cache_ttl_seconds
                )
                if signature_match or ttl_match:
                    return cached_value

        value = job.is_complete()
        with self._status_cache_lock:
            self._status_cache[key] = (signature, now, value)
        return value

    def _wrap_runner_command_for_node(self, runner: Any, node: str) -> Any:
        """Apply SLURM node pinning to runner command generation."""
        if not os.environ.get("SLURM_JOB_NODELIST"):
            return runner

        try:
            original_get_command = runner._get_command
        except AttributeError:
            return runner

        def patched_get_command_slurm(self_runner: Any, job_obj: Job) -> str:
            command = original_get_command(job_obj)
            prefix = f"srun --nodelist={node} --exclusive -N1 -n1 "
            return prefix + command

        runner._get_command = types.MethodType(
            patched_get_command_slurm,
            runner,
        )
        return runner

    def is_complete(self) -> bool:
        """Return True when all child jobs are complete."""
        if not self.jobs:
            return True
        return all(self._job_is_complete_cached(job) for job in self.jobs)


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
        no_run_in_parallel=True,
        fail_fast=fail_fast,
        label=f"{parent.label}{label_suffix}",
        jobrunner=parent.jobrunner,
    )
    batch_job.enable_serial_local_execution()
    batch_job.run()
    return batch_job
