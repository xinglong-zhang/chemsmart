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
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed

from chemsmart.jobs.gaussian.batch import GaussianBatchJob
from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.runner import GaussianJobRunner
from chemsmart.jobs.gaussian.settings import GaussianpKaJobSettings
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
        parallel=False,
        **kwargs,
    ):
        if not isinstance(settings, GaussianpKaJobSettings):
            raise ValueError(
                f"Settings must be instance of GaussianpKaJobSettings for {self.__class__.__name__}, but is {settings} instead!"
            )

        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            **kwargs,
        )

        # parallel support – opt-in only, matching ORCApKaJob behaviour.
        # When False (default), jobs run serially regardless of the jobrunner's
        # run_in_serial flag, so software that does not support concurrent
        # processes (e.g. single-licence Gaussian) is never affected.
        self.parallel = bool(parallel)
        if self.jobrunner and getattr(self.jobrunner, "run_in_serial", False):
            if self.parallel:
                logger.info(
                    "Parallel execution disabled due to run_in_serial=True in JobRunner"
                )
            self.parallel = False

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

        # 2. Reference Acid (HRef / Ref-)
        if self.has_reference_jobs:
            self.ref_opt_jobs = self._prepare_ref_opt_jobs()
            self.ref_acid_job, self.ref_conjugate_base_job = self.ref_opt_jobs

    def _prepare_ref_opt_jobs(self):
        """Prepare gas phase optimization jobs for HRef and Ref-."""
        ref_acid_mol, ref_conjugate_base_mol = (
            self.settings.reference_pair_molecules()
        )
        ref_acid_settings, ref_conjugate_base_settings = (
            self.settings.reference_pair_job_settings()
        )

        ref_acid_job = GaussianOptJob(
            molecule=ref_acid_mol,
            settings=ref_acid_settings,
            label=f"{self.label}_HRef_opt",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        ref_conjugate_base_job = GaussianOptJob(
            molecule=ref_conjugate_base_mol,
            settings=ref_conjugate_base_settings,
            label=f"{self.label}_Ref_opt",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        return [ref_acid_job, ref_conjugate_base_job]

    def _optimized_molecule_from_job(self, job, fallback_molecule):
        """Return optimized geometry for a finished job, or a fallback molecule."""
        out = job._output()
        if out is not None and out.normal_termination:
            return out.molecule
        return fallback_molecule

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
        if self.sp_jobs is None:
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
            if not self._ref_opt_jobs_are_complete():
                logger.warning(
                    "Reference optimization jobs not complete. Cannot run reference SP jobs."
                )
                return

            if self.ref_sp_jobs is None:
                self._create_ref_sp_jobs()

            if self.ref_sp_jobs:
                for job in self.ref_sp_jobs:
                    if self.jobrunner:
                        job.jobrunner = self.jobrunner
                    job.run()

    def _create_sp_jobs(self):
        """Create solution phase SP jobs from optimized geometries."""
        _, conj_fallback_mol = self.settings.conjugate_pair_molecules(
            self.molecule
        )
        prot_opt_mol = self._optimized_molecule_from_job(
            self.protonated_job, self.molecule
        )
        conj_opt_mol = self._optimized_molecule_from_job(
            self.conjugate_base_job, conj_fallback_mol
        )

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

    def _create_ref_sp_jobs(self):
        """Create reference solution phase SP jobs from optimized geometries."""
        ref_acid_fallback_mol, ref_conjugate_base_fallback_mol = (
            self.settings.reference_pair_molecules()
        )
        ref_acid_opt_mol = self._optimized_molecule_from_job(
            self.ref_acid_job, ref_acid_fallback_mol
        )
        ref_conjugate_base_opt_mol = self._optimized_molecule_from_job(
            self.ref_conjugate_base_job, ref_conjugate_base_fallback_mol
        )
        ref_acid_sp_settings, ref_conjugate_base_sp_settings = (
            self.settings.reference_pair_sp_job_settings()
        )

        self.ref_acid_sp_job = GaussianSinglePointJob(
            molecule=ref_acid_opt_mol,
            settings=ref_acid_sp_settings,
            label=f"{self.label}_HRef_sp",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        self.ref_conjugate_base_sp_job = GaussianSinglePointJob(
            molecule=ref_conjugate_base_opt_mol,
            settings=ref_conjugate_base_sp_settings,
            label=f"{self.label}_Ref_sp",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        self.ref_sp_jobs = [
            self.ref_acid_sp_job,
            self.ref_conjugate_base_sp_job,
        ]

    # ------------------------------------------------------------------
    # Imaginary frequency validation
    # ------------------------------------------------------------------

    @staticmethod
    def _validate_imaginary_frequencies(job, role):
        """Check for imaginary frequencies after an optimisation job.

        Returns:
            None if validation passes; an error message string otherwise.
        """
        out = job._output()
        if out is None:
            return (
                f"[{role}] Optimisation produced no output – "
                f"cannot validate frequencies for {job.label}."
            )
        if not out.normal_termination:
            return None  # abnormal termination handled separately

        freqs = getattr(out, "vibrational_frequencies", None)
        if freqs is None:
            return None  # no frequency data – skip

        imaginary = [f for f in freqs if f < 0.0]
        if imaginary:
            return (
                f"[{role}] Imaginary frequency check FAILED for "
                f"{job.label}: found {len(imaginary)} imaginary "
                f"mode(s) {imaginary}. The optimised geometry is not a "
                f"true minimum – please re-optimise."
            )
        return None

    # ------------------------------------------------------------------
    # Parallel execution helpers
    # ------------------------------------------------------------------

    def _run_opt_worker(self, job, role, runner, cores, mem):
        """Run a single opt job and validate its frequencies.

        Returns a dict ``{"role", "label", "success", "error"}``.
        """
        result = {
            "role": role,
            "label": job.label,
            "success": False,
            "error": None,
        }

        if runner:
            job.jobrunner = runner.copy()
            job.jobrunner.num_cores = cores
            job.jobrunner.mem_gb = mem

        job.run()

        freq_err = self._validate_imaginary_frequencies(job, role)
        if freq_err:
            result["error"] = freq_err
            logger.error(freq_err)
            return result

        result["success"] = True
        return result

    @staticmethod
    def _run_sp_worker(job, role, runner, cores, mem):
        """Run a single SP job.

        Returns a dict ``{"role", "label", "success", "error"}``.
        """
        result = {
            "role": role,
            "label": job.label,
            "success": False,
            "error": None,
        }

        if runner:
            job.jobrunner = runner.copy()
            job.jobrunner.num_cores = cores
            job.jobrunner.mem_gb = mem

        job.run()
        result["success"] = True
        return result

    @staticmethod
    def _collect_futures(future_to_role):
        """Wait for all futures, collecting successes and failures.

        Returns ``(successes, failures)`` where each element is a list
        of label strings or error message strings respectively.
        """
        successes = []
        failures = []

        for fut in as_completed(future_to_role):
            role = future_to_role[fut]
            try:
                res = fut.result()
            except Exception:
                tb = traceback.format_exc()
                msg = (
                    f"[{role}] Worker crashed with unhandled exception:\n{tb}"
                )
                logger.error(msg)
                failures.append(msg)
            else:
                if res["success"]:
                    successes.append(res["label"])
                else:
                    failures.append(res["error"])

        return successes, failures

    def _run_parallel(self):
        """Run pKa workflow in parallel mode.

        Uses ``ThreadPoolExecutor`` so that a crash or validation failure
        in one worker does **not** prevent other workers from finishing.
        After all workers complete, a summary is logged and – if any
        worker failed – a ``RuntimeError`` is raised with every failure
        message.
        """
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

        # ── Phase 1: optimisation jobs ──────────────────────────────────
        opt_job_specs = list(zip(self.opt_jobs, ("HA_opt", "A_opt")))
        if self.has_reference_jobs:
            opt_job_specs.extend(
                zip(self.ref_opt_jobs, ("HRef_opt", "Ref_opt"))
            )

        concurrent_opt_jobs = len(opt_job_specs)
        opt_cores_per_job = max(1, int(total_cores // concurrent_opt_jobs))
        opt_mem_per_job = max(1, int(total_mem // concurrent_opt_jobs))

        logger.info(
            f"Parallel execution: splitting {total_cores} cores and {total_mem}GB memory "
            f"among {concurrent_opt_jobs} concurrent opt jobs "
            f"({opt_cores_per_job} cores, {opt_mem_per_job}GB each)."
        )

        all_successes = []
        all_failures = []

        with ThreadPoolExecutor(max_workers=concurrent_opt_jobs) as executor:
            future_to_role = {}
            for job, role in opt_job_specs:
                fut = executor.submit(
                    self._run_opt_worker,
                    job,
                    role,
                    runner,
                    opt_cores_per_job,
                    opt_mem_per_job,
                )
                future_to_role[fut] = role

            successes, failures = self._collect_futures(future_to_role)
            all_successes.extend(successes)
            all_failures.extend(failures)

        # ── Phase 2: SP jobs (only if opt phase had no fatal errors) ────
        # Create SP jobs from optimised geometries
        self._create_sp_jobs()
        if self.has_reference_jobs:
            self._create_ref_sp_jobs()

        sp_job_specs = []
        if self.sp_jobs:
            sp_job_specs.extend(zip(self.sp_jobs, ("HA_sp", "A_sp")))
        if self.has_reference_jobs and self.ref_sp_jobs:
            sp_job_specs.extend(zip(self.ref_sp_jobs, ("HRef_sp", "Ref_sp")))
        sp_jobs_to_run = [job for job, _ in sp_job_specs]
        sp_roles = [role for _, role in sp_job_specs]

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
                future_to_role = {}
                for job, role in zip(sp_jobs_to_run, sp_roles):
                    fut = executor.submit(
                        self._run_sp_worker,
                        job,
                        role,
                        runner,
                        sp_cores_per_job,
                        sp_mem_per_job,
                    )
                    future_to_role[fut] = role

                successes, failures = self._collect_futures(future_to_role)
                all_successes.extend(successes)
                all_failures.extend(failures)

        # ── Final reporting ─────────────────────────────────────────────
        total = len(all_successes) + len(all_failures)
        logger.info(
            f"Parallel pKa run complete: "
            f"{len(all_successes)}/{total} succeeded, "
            f"{len(all_failures)}/{total} failed."
        )
        for label in all_successes:
            logger.info(f"  ✓ {label}")
        for err in all_failures:
            logger.error(f"  ✗ {err}")

        if all_failures:
            summary = (
                f"{len(all_failures)} of {total} parallel pKa worker(s) failed:\n"
                + "\n".join(f"  - {e}" for e in all_failures)
            )
            raise RuntimeError(summary)

    def _run(self, **kwargs):
        """
        Execute the pKa calculation.

        Runs gas phase optimization jobs sequentially by default, then
        solution phase SP jobs.  Pass ``parallel=True`` to the constructor
        to enable concurrent execution of per-species jobs on the same node.
        """
        if self.parallel:
            logger.info(
                "Running Gaussian pKa calculation in parallel mode (splitting resources)"
            )
            self._run_parallel()
            return

        # Default sequential behaviour
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
        if not self.opt_jobs:
            return False
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
        if not self.ref_opt_jobs:
            return False
        return all(job.is_complete() for job in self.ref_opt_jobs)

    def _sp_jobs_are_complete(self):
        """
        Verify completion status of both solution phase SP jobs.

        Returns:
            bool: True if all SP jobs are complete.
        """
        if not self.sp_jobs:
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
        if not self.ref_sp_jobs:
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

    @property
    def original_mol(self):
        """
        Original molecule used to initialize the job (usually HA).
        """
        return self.molecule

    @property
    def conjugate_base_mol(self):
        """
        Conjugate base molecule (A-).
        """
        if self.conjugate_base_job:
            return self.conjugate_base_job.molecule
        # Fallback if jobs not prepared yet (unlikely given __init__)
        _, conj_mol = self.settings.conjugate_pair_molecules(self.molecule)
        return conj_mol


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
