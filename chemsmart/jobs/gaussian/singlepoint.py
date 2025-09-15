"""
Gaussian single point energy calculation job implementation.

This module provides the GaussianSinglePointJob class for performing
single point energy calculations using Gaussian. These calculations
compute the energy and properties of a molecule at a fixed geometry
without optimization.
"""

from chemsmart.jobs.gaussian.job import GaussianJob


class GaussianSinglePointJob(GaussianJob):
    """
    Gaussian job class for single point energy calculations.

    Performs single point energy calculations at fixed molecular
    geometries. Useful for computing energies, electronic properties,
    and other molecular descriptors at specific conformations.

    Attributes:
        TYPE (str): Job type identifier ('g16sp').
        molecule (Molecule): Molecular structure for calculation.
        settings (GaussianJobSettings): Calculation configuration options.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "g16sp"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        """
        Initialize a Gaussian single point calculation job.

        Sets up a single point energy calculation with the
        specified molecular structure and calculation settings.

        Args:
            molecule (Molecule): Molecular structure for calculation.
            settings (GaussianJobSettings): Calculation configuration.
            label (str): Job identifier for file naming.
            jobrunner (JobRunner, optional): Job execution handler.
            **kwargs: Additional keyword arguments for parent class.

        Raises:
            ValueError: If `settings` is not a GaussianJobSettings instance
                or `molecule` is not a Molecule (validated by base class).
        """
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
