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

from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.settings import GaussianpKaJobSettings
from chemsmart.jobs.gaussian.singlepoint import GaussianSinglePointJob

logger = logging.getLogger(__name__)


class GaussianpKaJob(GaussianJob):
    """
    Gaussian job class for pKa calculations using direct thermodynamic cycle.

    Performs pKa calculations using the following workflow:
    1. Optimize HA in gas phase (opt + freq) - get G(HA)_gas
    2. Optimize A- in gas phase (opt + freq) - get G(A-)_gas
    3. Run SP on optimized HA in solution - get E(HA)_aq
    4. Run SP on optimized A- in solution - get E(A-)_aq
    5. Calculate solvation free energies and pKa

    The thermodynamic cycle:
        HA(g) → A-(g) + H+(g)     [gas phase]
         ↓        ↓       ↓       [solvation]
        HA(aq) → A-(aq) + H+(aq)  [aqueous phase]

    Important: Uses the SAME level of theory (functional/basis) for both
    gas and solution phases to ensure proper error cancellation in the
    solvation free energy calculations.

    Attributes:
        TYPE (str): Job type identifier ('g16pka').
        molecule (Molecule): Protonated molecular structure (HA).
        settings (GaussianpKaJobSettings): pKa calculation configuration.
        label (str): Base job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the jobs.
        skip_completed (bool): If True, completed jobs are not rerun.

    Example:
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.gaussian.pka import GaussianpKaJob
        from chemsmart.jobs.gaussian.settings import GaussianpKaJobSettings

        mol = Molecule.from_filepath("acetic_acid.xyz")
        mol.charge = 0
        mol.multiplicity = 1

        settings = GaussianpKaJobSettings(
            proton_index=10,
            functional="B3LYP",
            basis="6-311+G(d,p)",
            solvent_model="SMD",
            solvent_id="water"
        )

        job = GaussianpKaJob(
            molecule=mol,
            settings=settings,
            label="acetic_acid_pka",
            jobrunner=jobrunner
        )

        # Run all jobs sequentially (gas opt + solution sp)
        job.run()

        # Access individual jobs
        protonated_opt_job, conjugate_base_opt_job = job.opt_jobs  # Gas phase
        protonated_sp_job, conjugate_base_sp_job = job.sp_jobs     # Solution phase
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
        """
        Initialize a Gaussian pKa calculation job.

        Args:
            molecule (Molecule): Protonated molecular structure (HA).
            settings (GaussianpKaJobSettings): pKa calculation configuration.
            label (str, optional): Base label for the jobs.
            jobrunner (JobRunner, optional): Job execution handler.
            skip_completed (bool): Skip completed jobs. Default is True.
            **kwargs: Additional keyword arguments for parent class.

        Raises:
            ValueError: If settings is not a GaussianpKaJobSettings instance.
            ValueError: If proton_index is not specified in settings.
        """
        if not isinstance(settings, GaussianpKaJobSettings):
            raise ValueError(
                f"Settings must be instance of GaussianpKaJobSettings, "
                f"but got {type(settings).__name__} instead!"
            )

        if settings.proton_index is None:
            raise ValueError(
                "proton_index must be specified in GaussianpKaJobSettings "
                "to identify which proton to remove for the conjugate base."
            )

        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            **kwargs,
        )

        # Cache for the created jobs
        self._opt_jobs = None
        self._sp_jobs = None
        self._ref_opt_jobs = None
        self._ref_sp_jobs = None

        # Parallel execution flag
        self.parallel = bool(parallel)
        # Lock to protect shared state when running in threads
        import threading

        self._pka_lock = threading.Lock()

    # =========================================================================
    # Basename helpers for label derivation
    # =========================================================================

    @property
    def _acid_basename(self):
        """Basename for the target acid (HA).

        Returns ``self.label`` unchanged so that the protonated job
        carries the same base name as the parent pKa job.
        """
        return self.label

    @property
    def _conjugate_base_label(self):
        """Label for the conjugate base (A⁻) derived from the acid basename."""
        return f"{self._acid_basename}_cb"

    @property
    def _ref_basename(self):
        """Basename for the reference acid (HB).

        Derived from the reference geometry filename stem to keep it
        unique when multiple HA molecules share the same HB.
        """
        import os

        if not self.settings.has_reference_file:
            return None
        return os.path.splitext(
            os.path.basename(self.settings.reference_file)
        )[0]

    @property
    def _ref_conjugate_base_label(self):
        """Label for the reference conjugate base (B⁻)."""
        ref = self._ref_basename
        if ref is None:
            return None
        return f"{ref}_cb"

    @classmethod
    def settings_class(cls):
        """
        Get the settings class used by this job type.

        Returns:
            Type[GaussianpKaJobSettings]: Settings class for pKa jobs.
        """
        return GaussianpKaJobSettings

    @property
    def protonated_molecule(self):
        """
        Get the protonated molecule (HA).

        Returns:
            Molecule: The protonated form with appropriate charge/multiplicity.
        """
        protonated_mol, _ = self.settings.conjugate_pair_molecules(
            self.molecule
        )
        return protonated_mol

    @property
    def conjugate_base_molecule(self):
        """
        Get the conjugate base molecule (A-).

        Returns:
            Molecule: The deprotonated form with proton removed.
        """
        _, conjugate_base_mol = self.settings.conjugate_pair_molecules(
            self.molecule
        )
        return conjugate_base_mol

    @property
    def opt_jobs(self):
        """
        Get both gas phase optimization jobs (protonated and conjugate base).

        Returns:
            tuple: A tuple of (protonated_opt_job, conjugate_base_opt_job),
                where each is a GaussianOptJob configured for
                gas phase optimization with frequency calculations.
        """
        if self._opt_jobs is None:
            self._opt_jobs = self._prepare_opt_jobs()
        return self._opt_jobs

    @property
    def protonated_job(self):
        """
        Get the gas phase optimization job for the protonated form (HA).

        Returns:
            GaussianOptJob: Gas phase optimization job for the protonated molecule.
        """
        return self.opt_jobs[0]

    @property
    def conjugate_base_job(self):
        """
        Get the gas phase optimization job for the conjugate base (A-).

        Returns:
            GaussianOptJob: Gas phase optimization job for the conjugate base molecule.
        """
        return self.opt_jobs[1]

    @property
    def sp_jobs(self):
        """
        Get both solution phase SP jobs (protonated and conjugate base).

        Returns:
            tuple: A tuple of (protonated_sp_job, conjugate_base_sp_job),
                where each is a GaussianSinglePointJob configured for
                solution phase single point calculations.
        """
        if self._sp_jobs is None:
            self._sp_jobs = self._prepare_sp_jobs()
        return self._sp_jobs

    @property
    def protonated_sp_job(self):
        """
        Get the solution phase SP job for the protonated form (HA).

        Returns:
            GaussianSinglePointJob: Solution phase SP job for the protonated molecule.
        """
        return self.sp_jobs[0]

    @property
    def conjugate_base_sp_job(self):
        """
        Get the solution phase SP job for the conjugate base (A-).

        Returns:
            GaussianSinglePointJob: Solution phase SP job for the conjugate base molecule.
        """
        return self.sp_jobs[1]

    # =========================================================================
    # Reference acid (HB) jobs for proton exchange cycle
    # =========================================================================

    @property
    def has_reference_jobs(self):
        """Check if reference acid jobs are configured."""
        return self.settings.has_reference_file

    @property
    def reference_molecule(self):
        """
        Get the reference acid molecule (HB).

        Returns:
            Molecule: The reference acid with appropriate charge/multiplicity.

        Raises:
            ValueError: If reference file is not provided.
        """
        return self.settings.get_reference_molecule()

    @property
    def reference_conjugate_base_molecule(self):
        """
        Get the reference conjugate base molecule (B-).

        Returns:
            Molecule: The reference conjugate base with proton removed.

        Raises:
            ValueError: If reference settings are invalid.
        """
        return self.settings.get_reference_conjugate_base_molecule()

    @property
    def ref_opt_jobs(self):
        """
        Get both gas phase optimization jobs for reference acid (HB and B-).

        Returns:
            tuple or None: A tuple of (ref_acid_opt_job, ref_conjugate_base_opt_job),
                or None if no reference file is provided.
        """
        if not self.has_reference_jobs:
            return None
        if self._ref_opt_jobs is None:
            self._ref_opt_jobs = self._prepare_ref_opt_jobs()
        return self._ref_opt_jobs

    @property
    def ref_acid_job(self):
        """
        Get the gas phase optimization job for the reference acid (HB).

        Returns:
            GaussianOptJob or None: Gas phase optimization job for reference acid.
        """
        if not self.has_reference_jobs:
            return None
        return self.ref_opt_jobs[0]

    @property
    def ref_conjugate_base_job(self):
        """
        Get the gas phase optimization job for the reference conjugate base (B-).

        Returns:
            GaussianOptJob or None: Gas phase optimization job for reference conjugate base.
        """
        if not self.has_reference_jobs:
            return None
        return self.ref_opt_jobs[1]

    @property
    def ref_sp_jobs(self):
        """
        Get both solution phase SP jobs for reference acid (HB and B-).

        Returns:
            tuple or None: A tuple of (ref_acid_sp_job, ref_conjugate_base_sp_job),
                or None if no reference file is provided.
        """
        if not self.has_reference_jobs:
            return None
        if self._ref_sp_jobs is None:
            self._ref_sp_jobs = self._prepare_ref_sp_jobs()
        return self._ref_sp_jobs

    @property
    def ref_acid_sp_job(self):
        """
        Get the solution phase SP job for the reference acid (HB).

        Returns:
            GaussianSinglePointJob or None: Solution phase SP job for reference acid.
        """
        if not self.has_reference_jobs:
            return None
        return self.ref_sp_jobs[0]

    @property
    def ref_conjugate_base_sp_job(self):
        """
        Get the solution phase SP job for the reference conjugate base (B-).

        Returns:
            GaussianSinglePointJob or None: Solution phase SP job for reference conjugate base.
        """
        if not self.has_reference_jobs:
            return None
        return self.ref_sp_jobs[1]

    def _prepare_opt_jobs(self):
        """
        Create GAS PHASE optimization jobs for both protonated and conjugate base forms.

        Generates two GaussianOptJob objects, one for the protonated form
        and one for the conjugate base, both configured with appropriate
        settings for pKa calculations.

        Returns:
            tuple: A tuple of (protonated_job, conjugate_base_job).
        """
        # Get molecules
        protonated_mol, conjugate_base_mol = (
            self.settings.conjugate_pair_molecules(self.molecule)
        )

        # Get job settings
        protonated_settings, conjugate_base_settings = (
            self.settings.conjugate_pair_job_settings(self.molecule)
        )

        # Create job labels (acid basename + _cb for conjugate base)
        protonated_label = self._acid_basename
        conjugate_base_label = self._conjugate_base_label

        # Create protonated job (HA)
        protonated_job = GaussianOptJob(
            molecule=protonated_mol,
            settings=protonated_settings,
            label=protonated_label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

        # Create conjugate base job (A-)
        conjugate_base_job = GaussianOptJob(
            molecule=conjugate_base_mol,
            settings=conjugate_base_settings,
            label=conjugate_base_label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

        logger.debug(
            f"pKa jobs created: {protonated_job}, {conjugate_base_job}"
        )

        return protonated_job, conjugate_base_job

    def _prepare_ref_opt_jobs(self):
        """
        Create GAS PHASE optimization jobs for reference acid (HB) and its conjugate base (B-).

        Returns:
            tuple: A tuple of (ref_acid_job, ref_conjugate_base_job).

        Raises:
            ValueError: If reference file is not provided.
        """
        if not self.has_reference_jobs:
            raise ValueError(
                "Cannot prepare reference opt jobs: no reference file provided."
            )

        # Get reference molecules
        ref_acid_mol, ref_conjugate_base_mol = (
            self.settings.reference_pair_molecules()
        )

        # Get job settings
        ref_acid_settings, ref_conjugate_base_settings = (
            self.settings.reference_pair_job_settings()
        )

        # Create job labels (ref basename + _cb for conjugate base)
        ref_acid_label = self._ref_basename
        ref_conjugate_base_label = self._ref_conjugate_base_label

        # Create reference acid job (HB)
        ref_acid_job = GaussianOptJob(
            molecule=ref_acid_mol,
            settings=ref_acid_settings,
            label=ref_acid_label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

        # Create reference conjugate base job (B-)
        ref_conjugate_base_job = GaussianOptJob(
            molecule=ref_conjugate_base_mol,
            settings=ref_conjugate_base_settings,
            label=ref_conjugate_base_label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

        logger.debug(
            f"Reference pKa jobs created: {ref_acid_job}, {ref_conjugate_base_job}"
        )

        return ref_acid_job, ref_conjugate_base_job

    def _prepare_ref_sp_jobs(self):
        """
        Create SP jobs for reference acid (HB) and its conjugate base (B-).

        SP jobs use the optimized geometries from the reference optimization jobs.

        Returns:
            tuple: A tuple of (ref_acid_sp_job, ref_conjugate_base_sp_job).

        Raises:
            ValueError: If reference file is not provided.
        """
        if not self.has_reference_jobs:
            raise ValueError(
                "Cannot prepare reference SP jobs: no reference file provided."
            )

        # Get SP job settings
        ref_acid_sp_settings, ref_conjugate_base_sp_settings = (
            self.settings.reference_pair_sp_job_settings()
        )

        # Create job labels
        ref_acid_sp_label = f"{self._ref_basename}_sp"
        ref_conjugate_base_sp_label = f"{self._ref_conjugate_base_label}_sp"

        # Get optimized molecules from completed ref opt jobs
        ref_acid_opt_output = self.ref_acid_job._output()
        if (
            ref_acid_opt_output is not None
            and ref_acid_opt_output.normal_termination
        ):
            ref_acid_mol = ref_acid_opt_output.molecule
        else:
            # Fall back to initial molecule if opt not complete
            ref_acid_mol = self.reference_molecule

        ref_conjugate_base_opt_output = self.ref_conjugate_base_job._output()
        if (
            ref_conjugate_base_opt_output is not None
            and ref_conjugate_base_opt_output.normal_termination
        ):
            ref_conjugate_base_mol = ref_conjugate_base_opt_output.molecule
        else:
            # Fall back to initial molecule if opt not complete
            ref_conjugate_base_mol = self.reference_conjugate_base_molecule

        # Create reference acid SP job (HB)
        ref_acid_sp_job = GaussianSinglePointJob(
            molecule=ref_acid_mol,
            settings=ref_acid_sp_settings,
            label=ref_acid_sp_label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

        # Create reference conjugate base SP job (B-)
        ref_conjugate_base_sp_job = GaussianSinglePointJob(
            molecule=ref_conjugate_base_mol,
            settings=ref_conjugate_base_sp_settings,
            label=ref_conjugate_base_sp_label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

        logger.debug(
            f"Reference SP jobs created: {ref_acid_sp_job}, {ref_conjugate_base_sp_job}"
        )

        return ref_acid_sp_job, ref_conjugate_base_sp_job

    def _prepare_sp_jobs(self):
        """
        Create SP jobs for both protonated and conjugate base forms.

        SP jobs use the optimized geometries from the optimization jobs.
        This method should only be called after optimization jobs are complete.

        Returns:
            tuple: A tuple of (protonated_sp_job, conjugate_base_sp_job).
        """
        # Get SP job settings
        protonated_sp_settings, conjugate_base_sp_settings = (
            self.settings._create_solution_phase_sp_settings(self.molecule)
        )

        # Create job labels
        protonated_sp_label = f"{self._acid_basename}_sp"
        conjugate_base_sp_label = f"{self._conjugate_base_label}_sp"

        # Get optimized molecules from completed opt jobs
        # If opt jobs are complete, use their optimized structures
        protonated_opt_output = self.protonated_job._output()
        if (
            protonated_opt_output is not None
            and protonated_opt_output.normal_termination
        ):
            protonated_mol = protonated_opt_output.molecule
        else:
            # Fall back to initial molecule if opt not complete
            protonated_mol = self.protonated_molecule

        conjugate_base_opt_output = self.conjugate_base_job._output()
        if (
            conjugate_base_opt_output is not None
            and conjugate_base_opt_output.normal_termination
        ):
            conjugate_base_mol = conjugate_base_opt_output.molecule
        else:
            # Fall back to initial molecule if opt not complete
            conjugate_base_mol = self.conjugate_base_molecule

        # Create protonated SP job (HA)
        protonated_sp_job = GaussianSinglePointJob(
            molecule=protonated_mol,
            settings=protonated_sp_settings,
            label=protonated_sp_label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

        # Create conjugate base SP job (A-)
        conjugate_base_sp_job = GaussianSinglePointJob(
            molecule=conjugate_base_mol,
            settings=conjugate_base_sp_settings,
            label=conjugate_base_sp_label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

        logger.debug(
            f"SP jobs created: {protonated_sp_job}, {conjugate_base_sp_job}"
        )

        return protonated_sp_job, conjugate_base_sp_job

    def _run_opt_jobs(self):
        """
        Execute both GAS PHASE optimization jobs (protonated and conjugate base) sequentially.

        Runs the optimization + frequency jobs for both the protonated
        form (HA) and the conjugate base (A-) in sequence.
        """
        for job in self.opt_jobs:
            logger.info(f"Running gas phase optimization job: {job}")
            job.run()

    def _run_ref_opt_jobs(self):
        """
        Execute both GAS PHASE optimization jobs for reference acid (HB and B-) sequentially.

        Runs the optimization + frequency jobs for both the reference acid
        form (HB) and its conjugate base (B-) in sequence.
        """
        if not self.has_reference_jobs:
            logger.debug(
                "No reference jobs to run (no reference file provided)"
            )
            return

        for job in self.ref_opt_jobs:
            logger.info(f"Running reference gas phase optimization job: {job}")
            job.run()

    def _run_sp_jobs(self):
        """
        Execute both SOLUTION PHASE SP jobs (protonated and conjugate base) sequentially.

        Runs the SP jobs for both the protonated form (HA) and the
        conjugate base (A-) in sequence. Should only be called after
        optimization jobs are complete.
        """
        # Clear cached SP jobs to get fresh ones with optimized geometries
        self._sp_jobs = None

        for job in self.sp_jobs:
            logger.info(f"Running solution phase SP job: {job}")
            job.run()

    def _run_ref_sp_jobs(self):
        """
        Execute both SOLUTION PHASE SP jobs for reference acid (HB and B-) sequentially.

        Runs the SP jobs for both the reference acid form (HB) and its
        conjugate base (B-) in sequence. Should only be called after
        reference optimization jobs are complete.
        """
        if not self.has_reference_jobs:
            logger.debug(
                "No reference SP jobs to run (no reference file provided)"
            )
            return

        # Clear cached ref SP jobs to get fresh ones with optimized geometries
        self._ref_sp_jobs = None

        for job in self.ref_sp_jobs:
            logger.info(f"Running reference solution phase SP job: {job}")
            job.run()

    def _make_sp_job_for_role(self, role):
        """Create a single SP job for the given role: 'HA' or 'A'.

        This is used by the parallel worker after the corresponding opt job
        finishes so we can pick up the optimized geometry if available.
        """
        # Get SP settings for both species
        protonated_sp_settings, conjugate_base_sp_settings = (
            self.settings._create_solution_phase_sp_settings(self.molecule)
        )

        # Choose appropriate settings and molecule based on role
        if role == "HA":
            sp_settings = protonated_sp_settings
            opt_job = self.protonated_job
            sp_label = f"{self._acid_basename}_sp"
        else:
            sp_settings = conjugate_base_sp_settings
            opt_job = self.conjugate_base_job
            sp_label = f"{self._conjugate_base_label}_sp"

        # Prefer optimized geometry if opt produced one
        out = opt_job._output()
        if out is not None and getattr(out, "normal_termination", False):
            mol = out.molecule
        else:
            # fall back to prepared molecules
            mol = (
                self.protonated_molecule
                if role == "HA"
                else self.conjugate_base_molecule
            )

        sp_job = GaussianSinglePointJob(
            molecule=mol,
            settings=sp_settings,
            label=sp_label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        return sp_job

    def _worker_opt_then_sp(self, opt_job, role, sp_index):
        """Thread worker: run opt_job, then create and run its SP job.

        sp_index indicates where to store the created sp job in self._sp_jobs.
        """
        logger.info(f"[parallel] Running opt job for role={role}: {opt_job}")
        opt_job.run()

        # create SP job using optimized geometry if available
        sp_job = self._make_sp_job_for_role(role)
        logger.info(f"[parallel] Running SP job for role={role}: {sp_job}")
        sp_job.run()

        # store sp_job into cached list
        with self._pka_lock:
            if self._sp_jobs is None:
                self._sp_jobs = [None, None]
            self._sp_jobs[sp_index] = sp_job

    def _run_parallel(self):
        """Run opt->sp for both species in parallel (per-species dependency)."""
        import threading

        # Ensure opt jobs are created
        opt_jobs = self.opt_jobs

        # Prepare sp cache
        with self._pka_lock:
            self._sp_jobs = [None, None]

        threads = []
        # HA is index 0, A is index 1
        threads.append(
            threading.Thread(
                target=self._worker_opt_then_sp, args=(opt_jobs[0], "HA", 0)
            )
        )
        threads.append(
            threading.Thread(
                target=self._worker_opt_then_sp, args=(opt_jobs[1], "A", 1)
            )
        )

        # If reference jobs exist, run them in parallel as well
        if self.has_reference_jobs:
            ref_opt_jobs = self.ref_opt_jobs
            # create ref sp cache
            with self._pka_lock:
                self._ref_sp_jobs = [None, None]
            threads.append(
                threading.Thread(
                    target=self._worker_opt_then_sp,
                    args=(ref_opt_jobs[0], "HB", 0),
                )
            )
            threads.append(
                threading.Thread(
                    target=self._worker_opt_then_sp,
                    args=(ref_opt_jobs[1], "B", 1),
                )
            )

        # start all threads
        for t in threads:
            t.start()
        # join
        for t in threads:
            t.join()

    def _run(self):
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
        # Run gas phase optimization jobs for target acid (HA, A-)
        self._run_opt_jobs()

        # Run gas phase optimization jobs for reference acid (HB, B-) if provided
        if self.has_reference_jobs:
            self._run_ref_opt_jobs()

        # Run solution phase SP jobs for target acid
        self._run_sp_jobs()

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
