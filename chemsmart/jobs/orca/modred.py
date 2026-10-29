"""
ORCA modredundant job implementation.

This module contains the ORCAModredJob class for running geometry
optimizations with modified redundant coordinates using ORCA.
"""

import logging

from chemsmart.jobs.orca.job import ORCAJob

logger = logging.getLogger(__name__)


class ORCAModredJob(ORCAJob):
    """
    ORCA modredundant optimization job.

    This class handles constrained optimizations using modified redundant
    coordinates to fix specific geometric parameters during optimization.

    Attributes:
        TYPE (str): Job type identifier ('orcamodred').
        molecule: Molecule object used for the optimization.
        settings: ORCAJobSettings including modredundant parameters.
        label (str): Job identifier used for file naming.
        jobrunner: Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "orcamodred"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        """
        Initialize ORCAModredJob.

        Args:
            molecule: Molecule object for the optimization
            settings: ORCAJobSettings instance with modred parameters
            label: Job label for identification
            jobrunner: Job runner instance
            **kwargs: Additional keyword arguments
        """
        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )
