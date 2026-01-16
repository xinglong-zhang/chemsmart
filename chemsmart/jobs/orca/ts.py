"""
ORCA transition state job implementation.

This module contains the ORCATSJob class for running transition state
optimizations using ORCA.
"""

import logging

from chemsmart.jobs.orca.job import ORCAJob

logger = logging.getLogger(__name__)


class ORCATSJob(ORCAJob):
    """
    ORCA transition state optimization job.

    This class handles transition state optimizations to locate saddle
    points on potential energy surfaces for reaction pathway analysis.

    Attributes:
        TYPE (str): Job type identifier ('orcats').
        molecule: Molecule object used for the TS optimization.
        settings: ORCATSJobSettings configuration for the job.
        label (str): Job identifier used for file naming.
        jobrunner: Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "orcats"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        """
        Initialize ORCATSJob.

        Args:
            molecule: Molecule object for the TS optimization
            settings: ORCATSJobSettings instance
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
