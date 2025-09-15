"""
Gaussian modified redundant coordinate optimization job implementation.

This module provides the GaussianModredJob class for performing
geometry optimizations with modified redundant internal coordinates
using Gaussian. These calculations allow for constrained optimizations
where specific internal coordinates are frozen or scanned.
"""

from chemsmart.jobs.gaussian.job import GaussianJob


class GaussianModredJob(GaussianJob):
    """
    Gaussian job class for modified redundant coordinate optimizations.

    Performs geometry optimizations using modified redundant internal
    coordinates. Allows for freezing specific bonds, angles, or dihedrals
    during optimization, enabling constrained geometry optimizations
    and potential energy surface scans.

    Attributes:
        TYPE (str): Job type identifier ('g16modred').
        molecule (Molecule): Molecular structure used for optimization.
        settings (GaussianJobSettings): Configuration including modredundant
            coordinate specifications.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "g16modred"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        """
        Initialize a Gaussian modredundant coordinate optimization job.

        Sets up a constrained optimization calculation with the
        specified molecular structure and coordinate constraints.

        Args:
            molecule (Molecule): Molecular structure for optimization.
            settings (GaussianJobSettings): Calculation configuration
                including modredundant coordinate specifications.
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
