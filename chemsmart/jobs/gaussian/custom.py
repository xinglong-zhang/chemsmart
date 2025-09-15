"""
Gaussian custom job implementation for flexible calculations.

This module provides the GaussianCustomJob class for performing
custom Gaussian calculations with user-defined parameters.
Serves as a general-purpose job type for calculations that don't
fit into standard specialized job categories.
"""

from chemsmart.jobs.gaussian.job import GaussianJob


class GaussianCustomJob(GaussianJob):
    """
    Gaussian job class for custom user-defined calculations.

    Provides a flexible job type for running custom Gaussian
    calculations that may not fit into standard specialized
    categories. Allows users to define their own calculation
    parameters and job configurations.

    Attributes:
        TYPE (str): Job type identifier ('g16job').
        molecule (Molecule): Molecular structure used for the calculation.
        settings (GaussianJobSettings): Calculation configuration options.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "g16job"

    def __init__(
        self, molecule, settings=None, label=None, jobrunner=None, **kwargs
    ):
        """
        Initialize a Gaussian custom calculation job.

        Sets up a custom calculation with the specified molecular
        structure and user-defined calculation settings.

        Args:
            molecule (Molecule): Molecular structure for calculation.
            settings (GaussianJobSettings, optional): Calculation configuration.
            label (str, optional): Job identifier for file naming.
            jobrunner (JobRunner, optional): Job execution handler.
            **kwargs: Additional keyword arguments for parent class.
        """
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
