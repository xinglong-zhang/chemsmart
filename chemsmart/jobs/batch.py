"""
Shared batch job infrastructure.

This module provides the abstract ``BatchJob`` base class for orchestrating
collections of engine-specific jobs. It centralizes engine-agnostic behavior
such as:
- serial vs parallel submission
- fault-tolerant execution of child jobs
- scheduler/node allocation detection (SLURM/PBS)
- distribution of child jobs across allocated nodes

Concrete subclasses only need to implement engine-specific runner adaptation,
for example pinning a copied runner to a scheduler node.
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
from typing import Any, Optional, Sequence

from chemsmart.jobs.job import Job
from chemsmart.jobs.runner import get_serial_mode
from chemsmart.utils.mixins import RegistryMeta

logger = logging.getLogger(__name__)


class BatchExecutionError(RuntimeError):
    """Raised when one or more child jobs fail in a batch run."""


class BatchJobMeta(RegistryMeta, ABCMeta):
    """Metaclass combining registry support with abstract-base semantics."""


class BatchJob(Job, metaclass=BatchJobMeta):
    """
    Abstract controller for running a batch of child jobs.

    ``BatchJob`` intentionally keeps orchestration logic engine-agnostic and
    delegates engine-specific execution details to subclasses via
    ``_configure_runner_for_node``.
    """

    PROGRAM: Optional[str] = None
    REGISTERABLE: bool = False

    def __init__(
        self,
        jobs: Optional[Sequence[Job]],
        run_in_serial: Optional[bool] = None,
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
        if run_in_serial is None:
            self.run_in_serial = runner_serial_mode.run_in_serial
        else:
            self.run_in_serial = bool(run_in_serial)
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

    def run(self, **kwargs: Any) -> None:
        """
        Run the batch of jobs.

        This intentionally preserves the historical batch-job behavior of
        executing child jobs directly instead of relying on ``Job.run()``'s
        completion short-circuit at the batch container level.
        """
        self._invalidate_status_cache()
        self._run(**kwargs)

    def _run(self, **kwargs: Any) -> None:
        """Dispatch execution to multi-node, serial, or parallel mode."""
        nodes = self._get_allocated_nodes()

        if nodes and len(nodes) > 1:
            outcomes = self._run_multi_node(nodes, **kwargs)
        elif self.run_in_serial:
            logger.info(f"Running batch of {len(self.jobs)} jobs serially.")
            outcomes = self._run_jobs_serially(self.jobs, **kwargs)
        else:
            logger.info(f"Running batch of {len(self.jobs)} jobs in parallel.")
            outcomes = self._run_jobs_in_parallel(self.jobs, **kwargs)

        self._last_batch_outcomes = outcomes
        if self.write_outcome_logs:
            self._write_outcome_logs(outcomes)
        failures = [item for item in outcomes if not item["success"]]
        if failures:
            lines = [
                f"- {item['label']}: {item['error']}" for item in failures
            ]
            raise BatchExecutionError(
                f"{len(failures)} of {len(outcomes)} batch job(s) failed:\n"
                + "\n".join(lines)
            )

    def _run_jobs_serially(
        self,
        jobs: Sequence[Job],
        node: Optional[str] = None,
        **kwargs: Any,
    ) -> list[dict[str, Any]]:
        """Submit child jobs one-by-one."""
        outcomes: list[dict[str, Any]] = []
        for job in jobs:
            outcomes.append(self._submit_job(job, node=node, **kwargs))
        return outcomes

    def _run_jobs_in_parallel(
        self,
        jobs: Sequence[Job],
        **kwargs: Any,
    ) -> list[dict[str, Any]]:
        """Submit child jobs concurrently using a thread pool."""
        outcomes: list[dict[str, Any]] = []
        with ThreadPoolExecutor() as executor:
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
            self._job_is_complete_cached(job, force_refresh=True)
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
        """Copy the runner onto *job* and apply optional node adaptation.

        Uses ``Job._propagate_runner`` for the safe copy-and-assign step,
        then applies node-specific configuration if a *node* is given.
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
        with ThreadPoolExecutor(max_workers=num_nodes) as executor:
            futures = []
            for node, jobs_chunk in zip(nodes, job_chunks):
                if jobs_chunk:
                    futures.append(
                        executor.submit(
                            self._run_chunk_on_node,
                            jobs_chunk,
                            node,
                            **kwargs,
                        )
                    )

            for future in as_completed(futures):
                try:
                    outcomes.extend(future.result())
                except Exception as e:
                    logger.error(f"Node execution failed: {e}", exc_info=True)
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

    def _job_status_signature(
        self,
        job: Job,
    ) -> Optional[tuple[tuple[str, float], ...]]:
        """Return file-mtime signature used to validate cached status."""
        candidates: list[str] = []
        for attr in ("outputfile", "joblog", "errfile"):
            path = getattr(job, attr, None)
            if isinstance(path, str) and path:
                candidates.append(path)

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

        if not hasattr(runner, "_get_command"):
            return runner

        original_get_command = runner._get_command

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
