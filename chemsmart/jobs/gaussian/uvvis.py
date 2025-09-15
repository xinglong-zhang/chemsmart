"""
Gaussian UV-Vis spectroscopy calculation job implementation.

This module provides the GaussianUVVISJob class for performing
UV-Vis spectroscopy calculations using Gaussian. These calculations
combine ground state optimization with excited state calculations
to predict absorption spectra including solvent effects.

Note: This implementation is currently incomplete and under development.
"""

from copy import deepcopy

from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.singlepoint import GaussianSinglePointJob


class GaussianUVVISJob(GaussianJob):
    """
    Gaussian job class for UV-Vis spectroscopy calculations.

    Performs multi-step UV-Vis spectroscopy calculations including
    ground state geometry optimization followed by vertical excitation
    calculations. Handles solvent effects through equilibrium and
    non-equilibrium solvation models.

    The workflow includes:
    1. Ground state geometry optimization with equilibrium solvation
    2. Vertical excitation calculations with linear response solvation

    Note: This job type is currently incomplete and under development.

    Attributes:
        TYPE (str): Job type identifier ('g16uvvis').
        folder (str): Working directory for the UV-Vis workflow.
        atoms: Molecular structure used for the calculations (API subject
            to change; intended to be a Molecule instance in the finalized
            implementation).
        settings (GaussianJobSettings): Calculation configuration applied
            to both ground-state optimization and vertical excitation steps.
    """

    # TODO: incomplete job
    TYPE = "g16uvvis"

    def __init__(self, folder, atoms, settings):
        """
        Initialize a Gaussian UV-Vis spectroscopy calculation job.

        Note: This initialization signature differs from the standard
        GaussianJob pattern and needs to be updated.

        Args:
            folder (str): Working directory for the calculation.
            atoms: Molecular structure (needs update to use Molecule).
            settings (GaussianJobSettings): Calculation configuration.
        """
        super().__init__(folder=folder, atoms=atoms, settings=settings)
        self.settings = settings

    def _gs_geom_opt(self):
        """
        Create ground state geometry optimization job.

        Generates a GaussianOptJob for ground state geometry optimization
        with equilibrium solvation effects. This is the first step in
        the UV-Vis calculation workflow.

        Returns:
            GaussianOptJob: Ground state optimization job with appropriate
                solvent model settings for equilibrium solvation.
        """
        # Ground state geometry optimization and frequencies
        # add equilibrium solvation if necessary
        gs_geom_opt_settings = deepcopy(self.settings)
        return GaussianOptJob(
            folder=self.folder,
            atoms=self.molecule,
            settings=gs_geom_opt_settings,
        )

    def _run_gs_geom_opt(self, jobrunner, queue_manager=None):
        """
        Execute the ground state geometry optimization job.

        Runs the ground state optimization job using the provided
        job runner and optional queue manager. This must complete
        before vertical excitation calculations can proceed.

        Args:
            jobrunner: Job execution handler.
            queue_manager: Optional queue management system.
        """
        self._gs_geom_opt().run(
            jobrunner=jobrunner, queue_manager=queue_manager
        )

    def _gs_geom_opt_complete(self):
        """
        Check if ground state geometry optimization is complete.

        Verifies that the ground state optimization job has
        finished successfully before proceeding to vertical
        excitation calculations.

        Returns:
            bool: True if ground state optimization is complete,
                False otherwise.
        """
        return self._gs_geom_opt().is_complete()

    def _vertical_excitation(self):
        """
        Create vertical excitation calculation job.

        Generates a single-point TD-DFT calculation for vertical
        excitation with linear response (non-equilibrium) solvation
        at the ground state equilibrium geometry. This is the second
        step in the UV-Vis calculation workflow.

        Returns:
            GaussianSinglePointJob or None: Vertical excitation job
                if ground state optimization is complete, None otherwise.
        """
        # Vertical excitation with linear response solvation
        # single-point TD-DFT of vertical excitation (with default solvation),
        # at ground state equilibrium geometry
        # i.e., linear response, non-equilibrium.
        if self._gs_geom_opt_complete():
            sp_settings = deepcopy(self.settings)
            sp_settings.job_type = "sp"
            sp_settings.freq = False
            return GaussianSinglePointJob(
                folder=self.folder, atoms=self.molecule, settings=sp_settings
            )
        return None

    def _run_vertical_excitation(self, jobrunner, queue_manager=None):
        """
        Execute the vertical excitation calculation job.

        Runs the vertical excitation TD-DFT calculation using the
        provided job runner and optional queue manager. This step
        calculates the electronic absorption spectrum.

        Args:
            jobrunner: Job execution handler.
            queue_manager: Optional queue management system.

        Returns:
            bool: True if vertical excitation calculation completes
                successfully, False otherwise.
        """
        self._vertical_excitation().run(
            jobrunner=jobrunner, queue_manager=queue_manager
        )
        return self._vertical_excitation().is_complete()

    def _run(self, jobrunner, queue_manager=None, **kwargs):
        """
        Execute the complete UV-Vis spectroscopy calculation workflow.

        Performs the full two-step UV-Vis calculation including:
        1. Ground state geometry optimization with equilibrium solvation
        2. Vertical excitation calculation with linear response solvation

        Args:
            jobrunner: Job execution handler for running calculations.
            queue_manager: Optional queue management system.
            **kwargs: Additional keyword arguments for job execution.
        """
        self._run_gs_geom_opt(jobrunner=jobrunner, queue_manager=queue_manager)
        self._run_vertical_excitation(
            jobrunner=jobrunner, queue_manager=queue_manager
        )
