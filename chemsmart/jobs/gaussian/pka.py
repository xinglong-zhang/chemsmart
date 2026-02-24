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
    def pka_jobs(self):
        """Alias for opt_jobs for backward compatibility."""
        return self.opt_jobs

    @property
    def protonated_job(self):
        """
        Get the gas phase optimization job for the protonated form (HA).

        Returns:
            GaussianOptJob: Gas phase optimization job for the protonated molecule.
        """
        return self.opt_jobs[0]

    @property
    def protonated_opt_job(self):
        """Alias for protonated_job."""
        return self.protonated_job

    @property
    def conjugate_base_job(self):
        """
        Get the gas phase optimization job for the conjugate base (A-).

        Returns:
            GaussianOptJob: Gas phase optimization job for the conjugate base molecule.
        """
        return self.opt_jobs[1]

    @property
    def conjugate_base_opt_job(self):
        """Alias for conjugate_base_job."""
        return self.conjugate_base_job

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

        # Create job labels
        protonated_label = f"{self.label}_HA"
        conjugate_base_label = f"{self.label}_A"

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
            self.settings._create_sp_job_settings(self.molecule)
        )

        # Create job labels
        protonated_sp_label = f"{self.label}_HA_sp"
        conjugate_base_sp_label = f"{self.label}_A_sp"

        # Get optimized molecules from completed opt jobs
        # If opt jobs are complete, use their optimized structures
        protonated_opt_output = self.protonated_job._output()
        if (
            protonated_opt_output is not None
            and protonated_opt_output.is_complete
        ):
            protonated_mol = protonated_opt_output.molecule
        else:
            # Fall back to initial molecule if opt not complete
            protonated_mol = self.protonated_molecule

        conjugate_base_opt_output = self.conjugate_base_job._output()
        if (
            conjugate_base_opt_output is not None
            and conjugate_base_opt_output.is_complete
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

    def _run_pka_jobs(self):
        """Alias for _run_opt_jobs for backward compatibility."""
        self._run_opt_jobs()

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

    def _run(self):
        """
        Execute the pKa calculation.

        Runs gas phase optimization jobs sequentially, then solution phase SP jobs.
        """
        # Run gas phase optimization jobs first
        self._run_opt_jobs()

        # Run solution phase SP jobs
        self._run_sp_jobs()

    def is_complete(self):
        """
        Check if all pKa jobs are complete.

        Returns:
            bool: True if all optimization jobs and SP jobs
                have completed successfully.
        """
        # Check optimization jobs
        if not self._opt_jobs_are_complete():
            return False

        # Check SP jobs
        return self._sp_jobs_are_complete()

    def _opt_jobs_are_complete(self):
        """
        Verify completion status of both gas phase optimization jobs.

        Returns:
            bool: True if all optimization jobs are complete.
        """
        return all(job.is_complete() for job in self.opt_jobs)

    def _pka_jobs_are_complete(self):
        """Alias for _opt_jobs_are_complete for backward compatibility."""
        return self._opt_jobs_are_complete()

    def _sp_jobs_are_complete(self):
        """
        Verify completion status of both solution phase SP jobs.

        Returns:
            bool: True if all SP jobs are complete.
        """
        if self.sp_jobs is None:
            return False
        return all(job.is_complete() for job in self.sp_jobs)

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
