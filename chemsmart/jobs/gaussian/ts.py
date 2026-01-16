"""
Gaussian transition state optimization job implementation.

This module provides the GaussianTSJob class for performing
transition state optimization calculations using Gaussian.
These calculations find saddle points on potential energy surfaces
corresponding to chemical reaction transition states.
"""

from chemsmart.jobs.gaussian.job import GaussianJob


class GaussianTSJob(GaussianJob):
    """
    Gaussian job class for transition state optimization calculations.

    Performs optimization to find transition state structures, which
    are first-order saddle points on the potential energy surface.
    Used for studying chemical reaction mechanisms and computing
    activation barriers.

    Attributes:
        TYPE (str): Job type identifier ('g16ts').
        molecule (Molecule): Transition state guess structure.
        settings (GaussianJobSettings): Calculation configuration options.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "g16ts"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        """
        Initialize a Gaussian transition state optimization job.

        Sets up a transition state optimization calculation with
        the specified molecular structure and calculation settings.

        Args:
            molecule (Molecule): Molecular structure (transition state guess).
            settings (GaussianJobSettings): Calculation configuration.
            label (str): Job identifier for file naming.
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
