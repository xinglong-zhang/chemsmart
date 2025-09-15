"""
ORCA single point job implementation.

This module contains the ORCASinglePointJob class for running single point
energy calculations using ORCA.
"""

import logging

from chemsmart.jobs.orca.job import ORCAJob

logger = logging.getLogger(__name__)


class ORCASinglePointJob(ORCAJob):
    """
    ORCA single point energy calculation job.

    This class handles single point energy calculations at fixed geometries
    without geometry optimization.

    Attributes:
        TYPE (str): Job type identifier ('orcasp').
        molecule: Molecule object used for the calculation.
        settings: ORCAJobSettings configuration for the job.
        label (str): Job identifier used for file naming.
        jobrunner: Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "orcasp"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        """
        Initialize ORCASinglePointJob.

        Args:
            molecule: Molecule object for the calculation
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
