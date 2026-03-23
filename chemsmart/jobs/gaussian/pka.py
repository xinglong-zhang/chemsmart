"""
Gaussian pKa calculation job implementation.

This module provides the GaussianpKaJob class for performing pKa
calculations using Gaussian with a proper thermodynamic cycle:
1. Gas phase optimization + frequency for both HA and A-
2. Solution phase single point for both HA and A- at the same level of theory

Using the same level of theory ensures proper error cancellation for
solvation free energy calculations.
"""

import logging
import os
import subprocess
import types
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

from chemsmart.jobs.job import Job

logger = logging.getLogger(__name__)


class GaussianpKaBatchJob(Job):
    """
    Gaussian job class for running a batch of pKa calculations.

    Executes a list of GaussianpKaJob instances either serially or in parallel.
    Inherits from Job directly as it acts as a controller/container.
    """

    TYPE = "g16pka"
    PROGRAM = "gaussian"

    def __init__(
        self,
        jobs,
        run_in_serial=True,
        label="batch_pka",
        jobrunner=None,
        **kwargs,
    ):
        """
        Initialize a Gaussian pKa batch job.

        Args:
            jobs (list[GaussianpKaJob]): List of pKa jobs to execute.
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
        Run the batch of pKa jobs.

        Implements fault tolerance: if one job fails, the batch continues.
        Supports parallel execution using ThreadPoolExecutor.
        Supports node-aware parallel execution if SLURM/PBS allocation is detected.
        """
        # Check for node allocation
        nodes = self._get_allocated_nodes()

        if nodes and len(nodes) > 1:
            self._run_multi_node(nodes)
        elif self.run_in_serial:
            logger.info(
                f"Running batch of {len(self.jobs)} pKa jobs serially."
            )
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
            logger.info(
                f"Running batch of {len(self.jobs)} pKa jobs in parallel."
            )
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
        """
        Execute the pKa calculation.

        Runs gas phase optimization jobs sequentially (default), then solution phase SP jobs.
        When `self.parallel` is True, run per-species opt->SP pipelines concurrently.
        """
        if self.parallel:
            logger.info("Running pKa calculation in parallel mode")
            # Run target and reference chains in parallel while preserving per-species dependency
            self._run_parallel()
            return

        # Default sequential behaviour preserved
        run_in_serial = self.jobrunner and getattr(
            self.jobrunner, "run_in_serial", False
        )

        # Run gas phase optimization jobs for target acid (HA, A-)
        self._run_opt_jobs()

        if run_in_serial and not self._opt_jobs_are_complete():
            logger.info("Opt jobs incomplete, halting serial execution.")
            return

        # Run gas phase optimization jobs for reference acid (HB, B-) if provided
        if self.has_reference_jobs:
            self._run_ref_opt_jobs()
            if run_in_serial and not self._ref_opt_jobs_are_complete():
                logger.info(
                    "Ref Opt jobs incomplete, halting serial execution."
                )
                return

        # Run solution phase SP jobs for target acid
        self._run_sp_jobs()

        if run_in_serial and not self._sp_jobs_are_complete():
            logger.info("SP jobs incomplete, halting serial execution.")
            return

        # Run solution phase SP jobs for reference acid if provided
        if self.has_reference_jobs:
            self._run_ref_sp_jobs()

    def is_complete(self):
        """
        Check if all pKa jobs are complete.

        Returns:
            bool: True if all optimization jobs and SP jobs
                have completed successfully (including reference jobs if provided).
        """
        # Check target acid optimization jobs
        if not self._opt_jobs_are_complete():
            return False

        # Check target acid SP jobs
        if not self._sp_jobs_are_complete():
            return False

        # Check reference acid jobs if provided
        if self.has_reference_jobs:
            if not self._ref_opt_jobs_are_complete():
                return False
            if not self._ref_sp_jobs_are_complete():
                return False

        return True

    def _opt_jobs_are_complete(self):
        """
        Verify completion status of both gas phase optimization jobs.

        Returns:
            bool: True if all optimization jobs are complete.
        """
        return all(job.is_complete() for job in self.opt_jobs)

    def _ref_opt_jobs_are_complete(self):
        """
        Verify completion status of both reference gas phase optimization jobs.

        Returns:
            bool: True if all reference optimization jobs are complete,
                or True if no reference jobs are configured.
        """
        if not self.has_reference_jobs:
            return True
        return all(job.is_complete() for job in self.ref_opt_jobs)

    def _sp_jobs_are_complete(self):
        """
        Verify completion status of both solution phase SP jobs.

        Returns:
            bool: True if all SP jobs are complete.
        """
        if self.sp_jobs is None:
            return False
        return all(job.is_complete() for job in self.sp_jobs)

    def _ref_sp_jobs_are_complete(self):
        """
        Verify completion status of both reference solution phase SP jobs.

        Returns:
            bool: True if all reference SP jobs are complete,
                or True if no reference jobs are configured.
        """
        if not self.has_reference_jobs:
            return True
        if self.ref_sp_jobs is None:
            return False
        return all(job.is_complete() for job in self.ref_sp_jobs)

    @property
    def protonated_output(self):
        """
        Get the output of the protonated gas phase optimization job.

        Returns:
            Gaussian16Output or None: Parsed output for the protonated job.
        """
        return self.protonated_job._output()

    @property
    def conjugate_base_output(self):
        """
        Get the output of the conjugate base gas phase optimization job.

        Returns:
            Gaussian16Output or None: Parsed output for the conjugate base job.
        """
        return self.conjugate_base_job._output()

    @property
    def protonated_sp_output(self):
        """
        Get the output of the protonated solution phase SP job.

        Returns:
            Gaussian16Output or None: Parsed output for the protonated SP job.
        """
        return self.protonated_sp_job._output()

    @property
    def conjugate_base_sp_output(self):
        """
        Get the output of the conjugate base solution phase SP job.

        Returns:
            Gaussian16Output or None: Parsed output for the conjugate base SP job.
        """
        return self.conjugate_base_sp_job._output()

    # =========================================================================
    # Reference acid output properties
    # =========================================================================

    @property
    def ref_acid_output(self):
        """
        Get the output of the reference acid gas phase optimization job.

        Returns:
            Gaussian16Output or None: Parsed output for the reference acid job.
        """
        if not self.has_reference_jobs:
            return None
        return self.ref_acid_job._output()

    @property
    def ref_conjugate_base_output(self):
        """
        Get the output of the reference conjugate base gas phase optimization job.

        Returns:
            Gaussian16Output or None: Parsed output for the reference conjugate base job.
        """
        if not self.has_reference_jobs:
            return None
        return self.ref_conjugate_base_job._output()

    @property
    def ref_acid_sp_output(self):
        """
        Get the output of the reference acid solution phase SP job.

        Returns:
            Gaussian16Output or None: Parsed output for the reference acid SP job.
        """
        if not self.has_reference_jobs:
            return None
        return self.ref_acid_sp_job._output()

    @property
    def ref_conjugate_base_sp_output(self):
        """
        Get the output of the reference conjugate base solution phase SP job.

        Returns:
            Gaussian16Output or None: Parsed output for the reference conjugate base SP job.
        """
        if not self.has_reference_jobs:
            return None
        return self.ref_conjugate_base_sp_job._output()

    # =========================================================================
    # Thermochemistry extraction
    # =========================================================================

    def get_pka_outputs(self):
        """
        Get Gaussian16pKaOutput objects for completed pKa jobs.

        Creates Gaussian16pKaOutput objects using the output files from
        completed optimization jobs. This provides access to electronic
        energies (E) and quasi-harmonic Gibbs free energies (qh-G(T)) for
        all species.

        Returns:
            dict: Dictionary with Gaussian16pKaOutput objects:
                - 'HA': Output for protonated acid
                - 'A': Output for conjugate base
                - 'HB': Output for reference acid (if available)
                - 'B': Output for reference conjugate base (if available)

        Raises:
            ValueError: If optimization jobs are not complete.

        Example:
            job = GaussianpKaJob(...)
            job.run()  # Run all jobs

            outputs = job.get_pka_outputs()
            print(f"E(HA) = {outputs['HA'].electronic_energy_in_units}")
            print(f"qh-G(HA) = {outputs['HA'].qh_gibbs_free_energy}")
            print(f"E(A-) = {outputs['A'].electronic_energy_in_units}")
            print(f"qh-G(A-) = {outputs['A'].qh_gibbs_free_energy}")
        """
        from chemsmart.io.gaussian.output import Gaussian16pKaOutput

        if not self._opt_jobs_are_complete():
            raise ValueError(
                "Cannot get thermochemistry: optimization jobs are not complete. "
                "Run the pKa jobs first using job.run()."
            )

        # Get output file paths from completed opt jobs
        ha_file = (
            self.protonated_job.outputfile if self.protonated_job else None
        )
        a_file = (
            self.conjugate_base_job.outputfile
            if self.conjugate_base_job
            else None
        )

        # Get reference files if available
        hb_file = None
        b_file = None
        if self.has_reference_jobs and self._ref_opt_jobs_are_complete():
            hb_file = (
                self.ref_acid_job.outputfile if self.ref_acid_job else None
            )
            b_file = (
                self.ref_conjugate_base_job.outputfile
                if self.ref_conjugate_base_job
                else None
            )

        return Gaussian16pKaOutput.from_pka_settings(
            settings=self.settings,
            ha_file=ha_file,
            a_file=a_file,
            hb_file=hb_file,
            b_file=b_file,
        )

    def compute_thermochemistry(self):
        """
        Compute and return thermochemistry results for all species.

        Convenience method that computes thermochemistry for all pKa species
        and returns the results dictionary.

        Returns:
            dict: Dictionary with thermochemistry data for each species.
                See Gaussian16pKaOutput.compute_pka_thermochemistry() for details.
        """
        from chemsmart.io.gaussian.output import Gaussian16pKaOutput

        if not self._opt_jobs_are_complete():
            raise ValueError(
                "Cannot compute thermochemistry: optimization jobs are not complete. "
                "Run the pKa jobs first using job.run()."
            )

        # Get output file paths from completed opt jobs
        ha_file = (
            self.protonated_job.outputfile if self.protonated_job else None
        )
        a_file = (
            self.conjugate_base_job.outputfile
            if self.conjugate_base_job
            else None
        )

        # Get reference files if available
        hb_file = None
        b_file = None
        if self.has_reference_jobs and self._ref_opt_jobs_are_complete():
            hb_file = (
                self.ref_acid_job.outputfile if self.ref_acid_job else None
            )
            b_file = (
                self.ref_conjugate_base_job.outputfile
                if self.ref_conjugate_base_job
                else None
            )

        return Gaussian16pKaOutput.compute_pka_thermochemistry(
            ha_file=ha_file,
            a_file=a_file,
            hb_file=hb_file,
            b_file=b_file,
            temperature=self.settings.temperature,
            concentration=self.settings.concentration,
            pressure=self.settings.pressure,
            cutoff_entropy_grimme=self.settings.cutoff_entropy_grimme,
            cutoff_enthalpy=self.settings.cutoff_enthalpy,
            energy_units=self.settings.energy_units,
        )

    def print_thermochemistry(self):
        """Print formatted thermochemistry summary to stdout."""
        from chemsmart.io.gaussian.output import Gaussian16pKaOutput

        if not self._opt_jobs_are_complete():
            raise ValueError(
                "Cannot print thermochemistry: optimization jobs are not complete. "
                "Run the pKa jobs first using job.run()."
            )

        # Get output file paths from completed opt jobs
        ha_file = (
            self.protonated_job.outputfile if self.protonated_job else None
        )
        a_file = (
            self.conjugate_base_job.outputfile
            if self.conjugate_base_job
            else None
        )

        # Get reference files if available
        hb_file = None
        b_file = None
        if self.has_reference_jobs and self._ref_opt_jobs_are_complete():
            hb_file = (
                self.ref_acid_job.outputfile if self.ref_acid_job else None
            )
            b_file = (
                self.ref_conjugate_base_job.outputfile
                if self.ref_conjugate_base_job
                else None
            )

        Gaussian16pKaOutput.print_pka_summary(
            ha_file=ha_file,
            a_file=a_file,
            hb_file=hb_file,
            b_file=b_file,
            temperature=self.settings.temperature,
            concentration=self.settings.concentration,
            pressure=self.settings.pressure,
            cutoff_entropy_grimme=self.settings.cutoff_entropy_grimme,
            cutoff_enthalpy=self.settings.cutoff_enthalpy,
            energy_units=self.settings.energy_units,
        )
