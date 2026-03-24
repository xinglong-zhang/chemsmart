"""
Gaussian batch job implementation.

This module provides the GaussianBatchJob class for executing a collection
of Gaussian jobs either serially or in parallel, with support for
multi-node execution on SLURM clusters.
"""

import logging
import os
import subprocess
import types
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

from chemsmart.jobs.job import Job

logger = logging.getLogger(__name__)


class GaussianBatchJob(Job):
    """
    Generic Gaussian job class for running a batch of calculations.

    Executes a list of job instances either serially or in parallel.
    Inherits from Job directly as it acts as a controller/container.
    """

    PROGRAM = "gaussian"

    def __init__(
        self,
        jobs,
        run_in_serial=True,
        label="batch_job",
        jobrunner=None,
        **kwargs,
    ):
        """
        Initialize a Gaussian batch job.

        Args:
            jobs (list[Job]): List of jobs to execute.
            run_in_serial (bool): If True, execute jobs one by one.
            label (str): Label for the batch job.
            jobrunner (JobRunner): Execution backend.
            **kwargs: Additional arguments for the base Job class.
        """
        super().__init__(
            molecule=None, label=label, jobrunner=jobrunner, **kwargs
        )
        self.jobs = jobs
        self.run_in_serial = run_in_serial

    def run(self):
        """
        Run the batch of jobs.

        Implements fault tolerance: if one job fails, the batch continues.
        Supports parallel execution using ThreadPoolExecutor.
        Supports node-aware parallel execution if SLURM/PBS allocation is detected.
        """
        # Check for node allocation
        nodes = self._get_allocated_nodes()

        if nodes and len(nodes) > 1:
            self._run_multi_node(nodes)
        elif self.run_in_serial:
            logger.info(f"Running batch of {len(self.jobs)} jobs serially.")
            for job in self.jobs:
                try:
                    # Provide specific job runner to child jobs
                    if self.jobrunner:
                        job.jobrunner = self.jobrunner.copy()
                    job.run()
                except Exception as e:
                    logger.error(
                        f"Job {job.label} failed during serial batch execution: {e}",
                        exc_info=True,
                    )
        else:
            logger.info(f"Running batch of {len(self.jobs)} jobs in parallel.")
            # Use ThreadPoolExecutor because job.run() typically waits on subprocesses
            # or IO, releasing the GIL.
            with ThreadPoolExecutor() as executor:
                future_to_job = {}
                for job in self.jobs:
                    # Provide specific job runner to child jobs
                    if self.jobrunner:
                        job.jobrunner = self.jobrunner.copy()
                    future_to_job[executor.submit(job.run)] = job

                for future in as_completed(future_to_job):
                    job = future_to_job[future]
                    try:
                        future.result()
                    except Exception as e:
                        logger.error(
                            f"Job {job.label} failed during parallel batch execution: {e}",
                            exc_info=True,
                        )

    def _get_allocated_nodes(self):
        """
        Detect allocated nodes from environment variables.
        Supports SLURM and PBS.
        Returns a list of unique node names, or None if not found/single node.
        """
        # SLURM
        nodelist = os.environ.get("SLURM_JOB_NODELIST")
        if nodelist:
            try:
                # Use scontrol to expand the nodelist
                output = subprocess.check_output(
                    ["scontrol", "show", "hostnames", nodelist],
                    universal_newlines=True,
                )
                nodes = [
                    n.strip() for n in output.strip().split("\n") if n.strip()
                ]
                return nodes
            except (subprocess.SubprocessError, FileNotFoundError):
                # Fallback if scontrol fails or not found (e.g. testing)
                # This might happen if SLURM_JOB_NODELIST is set but scontrol is not in path
                return nodelist.split(",")

        # PBS
        nodefile = os.environ.get("PBS_NODEFILE")
        if nodefile and os.path.exists(nodefile):
            with open(nodefile, "r") as f:
                nodes = [line.strip() for line in f if line.strip()]
            return sorted(list(set(nodes)))  # Unique nodes

        return None

    def _run_multi_node(self, nodes):
        """
        Distribute jobs across multiple nodes.
        Each node runs its assigned chunk of jobs sequentially (if run_in_serial is True)
        or manages them (if parallel).
        However, per requirements: 'ensure each node processes its assigned subset of molecules sequentially'
        implies we should run sequentially on each node.
        """
        num_nodes = len(nodes)
        logger.info(
            f"Distributing {len(self.jobs)} jobs across {num_nodes} nodes: {nodes}"
        )

        # Split jobs into chunks for each node
        job_chunks = np.array_split(self.jobs, num_nodes)

        # We use a ThreadPool to manage the remote executions
        # Each thread manages one node
        with ThreadPoolExecutor(max_workers=num_nodes) as executor:
            futures = []
            for i, node in enumerate(nodes):
                jobs_chunk = job_chunks[i]
                if len(jobs_chunk) > 0:
                    futures.append(
                        executor.submit(
                            self._run_chunk_on_node, jobs_chunk, node
                        )
                    )

            for future in as_completed(futures):
                try:
                    future.result()
                except Exception as e:
                    logger.error(f"Node execution failed: {e}", exc_info=True)

    def _run_chunk_on_node(self, jobs, node):
        """
        Run a list of jobs on a specific node.
        Uses srun (for SLURM) to pin the execution.
        """
        logger.info(f"Node {node} starting to process {len(jobs)} jobs.")
        for job in jobs:
            try:
                # Prepare a node-specific runner
                if self.jobrunner:
                    runner = self.jobrunner.copy()

                    # Patch _get_command to inject srun/execution pinning
                    # We capture the original method
                    original_get_command = runner._get_command

                    def patched_get_command_slurm(self_runner, job_obj):
                        cmd = original_get_command(job_obj)
                        # cmd is typically a string: "g16 input.com"
                        # We want: "srun --nodelist=node --exclusive -N1 -n1 g16 input.com"
                        # Note: srun flags:
                        # --nodelist: target node
                        # --exclusive: don't share CPU with other steps?
                        # -N1 -n1: 1 node, 1 task (Gaussian itself handles threading via OMP_NUM_THREADS)

                        prefix = f"srun --nodelist={node} --exclusive -N1 -n1 "
                        return prefix + cmd

                    # Detect if we should use SLURM srun
                    if os.environ.get("SLURM_JOB_NODELIST"):
                        runner._get_command = types.MethodType(
                            patched_get_command_slurm, runner
                        )
                    # For PBS, pinning is harder without pbsdsh, but let's assume SLURM for now as requested

                    job.jobrunner = runner

                job.run()
            except Exception as e:
                logger.error(
                    f"Job {job.label} failed on node {node}: {e}",
                    exc_info=True,
                )
        logger.info(f"Node {node} finished processing.")

    def _run(self, **kwargs):
        # Required implementation of abstract method, though run() is overridden.
        self.run()

    def _backup_files(self):
        # Implementation required by Job abstracts, but no files to backup for batch container
        pass

    def is_complete(self):
        """
        Check if all jobs in the batch are complete.
        """
        if not self.jobs:
            return True
        return all(job.is_complete() for job in self.jobs)
