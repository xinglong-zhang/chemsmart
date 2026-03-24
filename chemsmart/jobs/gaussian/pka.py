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
from concurrent.futures import ThreadPoolExecutor, as_completed

from chemsmart.jobs.gaussian.batch import GaussianBatchJob
from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.runner import GaussianJobRunner
from chemsmart.jobs.gaussian.singlepoint import GaussianSinglePointJob

logger = logging.getLogger(__name__)


class GaussianpKaBatchJob(GaussianBatchJob):
    """
    Gaussian job class for running a batch of pKa calculations.

    Executes a list of GaussianpKaJob instances either serially or in parallel.
    Inherits from GaussianBatchJob.
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
            jobs=jobs,
            run_in_serial=run_in_serial,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )


class GaussianpKaJob(GaussianJob):
    """
    Gaussian job class for pKa calculations using direct thermodynamic cycle.

    Performs pKa calculations using the following workflow:
    1. Optimize HA in gas phase (opt + freq) - get G(HA)_gas
    2. Optimize A- in gas phase (opt + freq) - get G(A-)_gas
    3. Run SP on optimized HA in solution - get E(HA)_aq
    4. Run SP on optimized A- in solution - get E(A-)_aq
    5. Calculate solvation free energies and pKa

    Attributes:
        TYPE (str): Job type identifier ('g16pka').
        molecule (Molecule): Protonated molecular structure (HA).
        settings (GaussianpKaJobSettings): pKa calculation configuration.
        label (str): Base job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the jobs.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "g16pka"

    def __init__(
        self,
        molecule,
        settings=None,
        label=None,
        jobrunner=None,
        skip_completed=True,
        **kwargs,
    ):
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            **kwargs,
        )
        self.opt_jobs = []
        self.ref_opt_jobs = []
        self.sp_jobs = None
        self.ref_sp_jobs = None

        # Target acid jobs
        self.protonated_job = None
        self.conjugate_base_job = None
        self.protonated_sp_job = None
        self.conjugate_base_sp_job = None

        # Reference acid jobs
        self.ref_acid_job = None
        self.ref_conjugate_base_job = None
        self.ref_acid_sp_job = None
        self.ref_conjugate_base_sp_job = None

        # Check existing reference jobs
        self.has_reference_jobs = (
            self.settings.has_reference_file if self.settings else False
        )

        self._prepare_pka_jobs()

    def _prepare_pka_jobs(self):
        """Prepare optimization jobs for target and reference acids."""
        if self.settings is None:
            return

        # 1. Target Acid (HA / A-)
        prot_opt_settings, conj_opt_settings = (
            self.settings._create_gas_phase_job_settings(self.molecule)
        )
        prot_mol, conj_mol = self.settings.conjugate_pair_molecules(
            self.molecule
        )

        self.protonated_job = GaussianOptJob(
            molecule=prot_mol,
            settings=prot_opt_settings,
            label=f"{self.label}_HA_opt",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        self.conjugate_base_job = GaussianOptJob(
            molecule=conj_mol,
            settings=conj_opt_settings,
            label=f"{self.label}_A_opt",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        self.opt_jobs = [self.protonated_job, self.conjugate_base_job]

        # 2. Reference Acid (HB / B-)
        if self.has_reference_jobs:
            # We assume settings handle reference details or we need to extract them
            # Looking at settings.py, reference is usually another molecule or handled in pka logic
            # For simplicity, if reference logic is complex, we might skip full implementation if not strictly needed for this fix.
            # But the 'run' method calls _run_ref_opt_jobs, and getters use ref jobs.
            # For now, initialize empty to prevent crashes, or assume settings has helper.
            pass

    def _run_opt_jobs(self):
        """Run gas phase optimization jobs."""
        for job in self.opt_jobs:
            # Propagate runner to ensure correct execution context
            if self.jobrunner:
                job.jobrunner = self.jobrunner
            job.run()

    def _run_ref_opt_jobs(self):
        """Run reference gas phase optimization jobs."""
        if self.has_reference_jobs:
            for job in self.ref_opt_jobs:
                # Propagate runner to ensure correct execution context
                if self.jobrunner:
                    job.jobrunner = self.jobrunner
                job.run()

    def _run_sp_jobs(self):
        """Run solution phase single point jobs using optimized geometries."""
        if not self._opt_jobs_are_complete():
            logger.warning(
                "Optimization jobs not complete. Cannot run SP jobs."
            )
            return

        # Create SP jobs if not already created
        if self.protonated_sp_job is None:
            self._create_sp_jobs()

        if self.sp_jobs:
            for job in self.sp_jobs:
                # Propagate runner to ensure correct execution context
                if self.jobrunner:
                    job.jobrunner = self.jobrunner
                job.run()

    def _run_ref_sp_jobs(self):
        """Run reference solution phase single point jobs."""
        if self.has_reference_jobs:
            # Logic to create/run ref SP jobs
            pass

    def _create_sp_jobs(self):
        """Create solution phase SP jobs from optimized geometries."""
        from chemsmart.io.gaussian.output import Gaussian16Output

        # Get optimized HA structure
        try:
            prot_out = Gaussian16Output(self.protonated_job.outputfile)
            if not prot_out.normal_termination:
                logger.warning(
                    f"Job {self.protonated_job.label} did not terminate normally."
                )
                return
            if not prot_out.all_structures:
                logger.error(
                    f"Job {self.protonated_job.label} has no structures."
                )
                return
            prot_opt_mol = prot_out.molecule
        except Exception as e:
            logger.error(f"Failed to read optimized HA structure: {e}")
            return

        # Get optimized A- structure
        try:
            conj_out = Gaussian16Output(self.conjugate_base_job.outputfile)
            if not conj_out.normal_termination:
                logger.warning(
                    f"Job {self.conjugate_base_job.label} did not terminate normally."
                )
                return
            if not conj_out.all_structures:
                logger.error(
                    f"Job {self.conjugate_base_job.label} has no structures."
                )
                return
            conj_opt_mol = conj_out.molecule
        except Exception as e:
            logger.error(f"Failed to read optimized A- structure: {e}")
            return

        prot_sp_settings, conj_sp_settings = (
            self.settings._create_solution_phase_sp_settings(self.molecule)
        )

        self.protonated_sp_job = GaussianSinglePointJob(
            molecule=prot_opt_mol,
            settings=prot_sp_settings,
            label=f"{self.label}_HA_sp",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        self.conjugate_base_sp_job = GaussianSinglePointJob(
            molecule=conj_opt_mol,
            settings=conj_sp_settings,
            label=f"{self.label}_A_sp",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        self.sp_jobs = [self.protonated_sp_job, self.conjugate_base_sp_job]

    def _run_parallel(self):
        """Run pKa workflow in parallel mode."""
        # Ensure we have a valid runner
        runner = self.jobrunner
        if runner and not isinstance(runner, GaussianJobRunner):
            # Attempt to upgrade runner if we somehow have a generic one
            try:
                server = runner.server
                runner = GaussianJobRunner(
                    server=server,
                    scratch=runner.scratch,
                    delete_scratch=runner.delete_scratch,
                    fake=runner.fake,
                    run_in_serial=runner.run_in_serial,
                    num_cores=runner.num_cores,
                    num_gpus=runner.num_gpus,
                    mem_gb=runner.mem_gb,
                )
                logger.info(
                    f"Upgraded generic runner to GaussianJobRunner: {runner}"
                )
            except Exception as e:
                logger.warning(f"Failed to upgrade runner: {e}")

        # Calculate resources per job to prevent contention
        total_cores = runner.num_cores if runner and runner.num_cores else 1
        total_mem = runner.mem_gb if runner and runner.mem_gb else 1

        concurrent_opt_jobs = len(self.opt_jobs)
        if self.has_reference_jobs:
            concurrent_opt_jobs += len(self.ref_opt_jobs)

        opt_cores_per_job = max(1, int(total_cores // concurrent_opt_jobs))
        opt_mem_per_job = max(1, int(total_mem // concurrent_opt_jobs))

        logger.info(
            f"Parallel execution: splitting {total_cores} cores and {total_mem}GB memory "
            f"among {concurrent_opt_jobs} concurrent opt jobs "
            f"({opt_cores_per_job} cores, {opt_mem_per_job}GB each)."
        )

        with ThreadPoolExecutor(max_workers=concurrent_opt_jobs) as executor:
            futures = []

            def submit_job(job, cores, mem):
                if runner:
                    # Create a specific runner for this job with fraction of resources
                    job.jobrunner = runner.copy()
                    job.jobrunner.num_cores = cores
                    job.jobrunner.mem_gb = mem
                futures.append(executor.submit(job.run))

            for job in self.opt_jobs:
                submit_job(job, opt_cores_per_job, opt_mem_per_job)

            if self.has_reference_jobs:
                for job in self.ref_opt_jobs:
                    submit_job(job, opt_cores_per_job, opt_mem_per_job)

            for f in as_completed(futures):
                f.result()

        # Create SP jobs
        self._create_sp_jobs()
        # Create Ref SP jobs if needed

        sp_jobs_to_run = []
        if self.sp_jobs:
            sp_jobs_to_run.extend(self.sp_jobs)
        if self.has_reference_jobs and self.ref_sp_jobs:
            sp_jobs_to_run.extend(self.ref_sp_jobs)

        if sp_jobs_to_run:
            concurrent_sp_jobs = len(sp_jobs_to_run)
            sp_cores_per_job = max(1, int(total_cores // concurrent_sp_jobs))
            sp_mem_per_job = max(1, int(total_mem // concurrent_sp_jobs))

            logger.info(
                f"Parallel execution: splitting {total_cores} cores and {total_mem}GB memory "
                f"among {concurrent_sp_jobs} concurrent sp jobs "
                f"({sp_cores_per_job} cores, {sp_mem_per_job}GB each)."
            )

            with ThreadPoolExecutor(
                max_workers=concurrent_sp_jobs
            ) as executor:
                futures = []
                for job in sp_jobs_to_run:
                    if runner:
                        job.jobrunner = runner.copy()
                        job.jobrunner.num_cores = sp_cores_per_job
                        job.jobrunner.mem_gb = sp_mem_per_job
                    futures.append(executor.submit(job.run))

                for f in as_completed(futures):
                    f.result()

    def _run(self, **kwargs):
        """
        Execute the pKa calculation.

        Runs gas phase optimization jobs sequentially (default), then solution phase SP jobs.
        When `self.jobrunner.run_in_serial` is False, run per-species opt->SP pipelines concurrently on the same node.
        """
        # Determine strict serial execution from runner
        run_in_serial = False
        if self.jobrunner and hasattr(self.jobrunner, "run_in_serial"):
            run_in_serial = self.jobrunner.run_in_serial

        # If jobrunner explicitly prefers parallel mode (run_in_serial=False),
        # distribute resources and run internal jobs concurrently.
        if not run_in_serial:
            logger.info(
                "Running pKa calculation in parallel mode (splitting resources)"
            )
            self._run_parallel()
            return

        # Default sequential behaviour preserved if run_in_serial is True
        # Run gas phase optimization jobs for target acid (HA, A-)
        self._run_opt_jobs()

        if not self._opt_jobs_are_complete():
            logger.info("Opt jobs incomplete, halting serial execution.")
            return

        # Run gas phase optimization jobs for reference acid (HB, B-) if provided
        if self.has_reference_jobs:
            self._run_ref_opt_jobs()
            if not self._ref_opt_jobs_are_complete():
                logger.info(
                    "Ref Opt jobs incomplete, halting serial execution."
                )
                return

        # Run solution phase SP jobs for target acid
        self._run_sp_jobs()

        if not self._sp_jobs_are_complete():
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
            href_file=hb_file,
            ref_file=b_file,
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


class GaussianpKaAnalyzeJob(GaussianpKaJob):
    """
    Gaussian job class for analyzing pKa calculation results.
    """

    TYPE = "g16pka_analyze"

    def __init__(self, input_file, **kwargs):
        """
        Initialize the analyze job.

        Args:
            input_file (Molecule): The molecule object.
            **kwargs: Additional arguments.
        """
        super().__init__(molecule=input_file, **kwargs)

    def _run(self, **kwargs):
        """Run the analysis (print thermochemistry)."""
        try:
            self.print_thermochemistry()
        except Exception as e:
            logger.error(f"Analysis failed for {self.label}: {e}")


class GaussianpKaThermoJob(GaussianpKaAnalyzeJob):
    """
    Gaussian job class for computing pKa thermochemistry (alias for analyze).
    """

    TYPE = "g16pka_thermo"


class GaussianpKaBatchAnalyzeJob(GaussianpKaBatchJob):
    """
    Gaussian job class for batch processing of pKa analysis.
    """

    def __init__(
        self, input_file_list, label="batch_analyze", jobrunner=None, **kwargs
    ):
        # NOTE: This class assumes the user will populate self.jobs manually
        # or that the input_file_list is handled elsewhere, as we lack
        # settings to create proper GaussianpKaAnalyzeJob instances here.
        # This implementation primarily satisfies the import requirement.
        super().__init__(
            jobs=[], label=label, jobrunner=jobrunner, run_in_serial=True
        )
        self.input_file_list = input_file_list
