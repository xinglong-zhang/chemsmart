"""
Gaussian pKa calculation job implementation.

This module provides the GaussianpKaJob class for performing pKa
calculations using Gaussian. It creates and runs two sequential
optimization jobs: one for the protonated form (HA) and one for
the conjugate base (A-).
"""

import logging

from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.settings import GaussianpKaJobSettings

logger = logging.getLogger(__name__)


class GaussianpKaJob(GaussianJob):
    """
    Gaussian job class for pKa calculations.

    Performs pKa calculations by creating and running two sequential
    optimization jobs: one for the protonated form (HA) and one for
    the conjugate base (A-). Both jobs include frequency calculations
    for thermochemistry.

    The pKa calculation workflow:
    1. Create protonated molecule (HA) with charge q
    2. Create conjugate base molecule (A-) with charge q-1
    3. Run optimization + frequency on HA
    4. Run optimization + frequency on A-
    5. Calculate pKa from free energy difference (post-processing)

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
            solvation_model="SMD",
            solvent_id="water"
        )

        job = GaussianpKaJob(
            molecule=mol,
            settings=settings,
            label="acetic_acid_pka",
            jobrunner=jobrunner
        )

        # Run both jobs sequentially
        job.run()

        # Or access individual jobs
        protonated_job, conjugate_base_job = job.pka_jobs
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
        self._pka_jobs = None

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
    def pka_jobs(self):
        """
        Get both pKa jobs (protonated and conjugate base).

        Returns:
            tuple: A tuple of (protonated_job, conjugate_base_job),
                where each is a GaussianOptJob configured for
                optimization with frequency calculations.
        """
        if self._pka_jobs is None:
            self._pka_jobs = self._prepare_pka_jobs()
        return self._pka_jobs

    @property
    def protonated_job(self):
        """
        Get the optimization job for the protonated form (HA).

        Returns:
            GaussianOptJob: Optimization job for the protonated molecule.
        """
        return self.pka_jobs[0]

    @property
    def conjugate_base_job(self):
        """
        Get the optimization job for the conjugate base (A-).

        Returns:
            GaussianOptJob: Optimization job for the conjugate base molecule.
        """
        return self.pka_jobs[1]

    def _prepare_pka_jobs(self):
        """
        Create optimization jobs for both protonated and conjugate base forms.

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

    def _run_pka_jobs(self):
        """
        Execute both pKa jobs (protonated and conjugate base) sequentially.

        Runs the optimization + frequency jobs for both the protonated
        form (HA) and the conjugate base (A-) in sequence.
        """
        for job in self.pka_jobs:
            logger.info(f"Running pKa job: {job}")
            job.run()

    def _run(self):
        """
        Execute the pKa calculation.

        Runs both optimization jobs sequentially.
        """
        self._run_pka_jobs()

    def is_complete(self):
        """
        Check if both pKa jobs are complete.

        Returns:
            bool: True if both protonated and conjugate base jobs
                have completed successfully.
        """
        return self._pka_jobs_are_complete()

    def _pka_jobs_are_complete(self):
        """
        Verify completion status of both pKa jobs.

        Returns:
            bool: True if all pKa jobs are complete.
        """
        return all(job.is_complete() for job in self.pka_jobs)

    @property
    def protonated_output(self):
        """
        Get the output of the protonated job.

        Returns:
            Gaussian16Output or None: Parsed output for the protonated job.
        """
        return self.protonated_job._output()

    @property
    def conjugate_base_output(self):
        """
        Get the output of the conjugate base job.

        Returns:
            Gaussian16Output or None: Parsed output for the conjugate base job.
        """
        return self.conjugate_base_job._output()
