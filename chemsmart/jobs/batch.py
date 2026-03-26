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
from abc import ABCMeta, abstractmethod
from concurrent.futures import ThreadPoolExecutor, as_completed

from chemsmart.jobs.job import Job
from chemsmart.utils.mixins import RegistryMeta

logger = logging.getLogger(__name__)


class BatchJobMeta(RegistryMeta, ABCMeta):
    """Metaclass combining registry support with abstract-base semantics."""


class BatchJob(Job, metaclass=BatchJobMeta):
    """
    Abstract controller for running a batch of child jobs.

    ``BatchJob`` intentionally keeps orchestration logic engine-agnostic and
    delegates engine-specific execution details to subclasses via
    ``_configure_runner_for_node``.
    """

    PROGRAM = None
    REGISTERABLE = False

    def __init__(
        self,
        jobs,
        run_in_serial=True,
        label="batch_job",
        jobrunner=None,
        **kwargs,
    ):
        super().__init__(
            molecule=None,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
        self.jobs = list(jobs) if jobs is not None else []
        self.run_in_serial = run_in_serial

    def run(self, **kwargs):
        """
        Run the batch of jobs.

        This intentionally preserves the historical batch-job behavior of
        executing child jobs directly instead of relying on ``Job.run()``'s
        completion short-circuit at the batch container level.
        """
        self._run(**kwargs)

    def _run(self, **kwargs):
        nodes = self._get_allocated_nodes()

        if nodes and len(nodes) > 1:
            self._run_multi_node(nodes, **kwargs)
        elif self.run_in_serial:
            logger.info(f"Running batch of {len(self.jobs)} jobs serially.")
            self._run_jobs_serially(self.jobs, **kwargs)
        else:
            logger.info(f"Running batch of {len(self.jobs)} jobs in parallel.")
            self._run_jobs_in_parallel(self.jobs, **kwargs)

    def _run_jobs_serially(self, jobs, node=None, **kwargs):
        for job in jobs:
            self._submit_job(job, node=node, **kwargs)

    def _run_jobs_in_parallel(self, jobs, **kwargs):
        with ThreadPoolExecutor() as executor:
            future_to_job = {
                executor.submit(self._submit_job, job, None, **kwargs): job
                for job in jobs
            }

            for future in as_completed(future_to_job):
                job = future_to_job[future]
                try:
                    future.result()
                except Exception as e:
                    logger.error(
                        f"Unexpected failure while coordinating job {job.label}: {e}",
                        exc_info=True,
                    )

    def _submit_job(self, job, node=None, **kwargs):
        try:
            runner = self._build_jobrunner(job, node=node)
            if runner is not None:
                job.jobrunner = runner
            job.run(**kwargs)
        except Exception as e:
            location = f" on node {node}" if node else ""
            logger.error(
                f"Job {job.label} failed during batch execution{location}: {e}",
                exc_info=True,
            )

    def _build_jobrunner(self, job, node=None):
        if not self.jobrunner:
            return None

        runner = self.jobrunner.copy()
        if node is not None:
            runner = self._configure_runner_for_node(
                runner=runner,
                node=node,
                job=job,
            )
        return runner

    def _get_allocated_nodes(self):
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

    def _run_multi_node(self, nodes, **kwargs):
        num_nodes = len(nodes)
        logger.info(
            f"Distributing {len(self.jobs)} jobs across {num_nodes} nodes: {nodes}"
        )

        job_chunks = self._split_jobs_across_nodes(self.jobs, num_nodes)

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
                    future.result()
                except Exception as e:
                    logger.error(f"Node execution failed: {e}", exc_info=True)

    def _run_chunk_on_node(self, jobs, node, **kwargs):
        logger.info(f"Node {node} starting to process {len(jobs)} jobs.")
        self._run_jobs_serially(jobs, node=node, **kwargs)
        logger.info(f"Node {node} finished processing.")

    @staticmethod
    def _split_jobs_across_nodes(jobs, num_nodes):
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
    def _configure_runner_for_node(self, runner, node, job):
        """
        Adapt a copied jobrunner for execution on a specific node.

        Subclasses can patch engine-specific command generation here.
        """
        raise NotImplementedError

    def _backup_files(self):
        pass

    def is_complete(self):
        if not self.jobs:
            return True
        return all(job.is_complete() for job in self.jobs)
