"""
ORCA optimization job implementation.

This module contains the ORCAOptJob class for running geometry
optimizations using ORCA.
"""

import logging

from chemsmart.jobs.orca.job import ORCAJob

logger = logging.getLogger(__name__)


class ORCAOptJob(ORCAJob):
    """
    ORCA geometry optimization job.

    This class handles molecular geometry optimizations to find energy
    minima on the potential energy surface.

    Attributes:
        TYPE (str): Job type identifier ('orcaopt').
        molecule: Molecule object used for the optimization.
        settings: ORCAJobSettings configuration for the job.
        label (str): Job identifier used for file naming.
        jobrunner: Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "orcaopt"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        """
        Initialize ORCAOptJob.

        Args:
            molecule: Molecule object for the optimization
            settings: ORCAJobSettings instance
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
