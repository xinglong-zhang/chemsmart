"""
ORCA scan job implementation.

This module contains the ORCAScanJob class for running potential energy
surface scans using ORCA.
"""

import logging

from chemsmart.jobs.orca.job import ORCAJob

logger = logging.getLogger(__name__)


class ORCAScanJob(ORCAJob):
    """
    ORCA potential energy surface scan job.

    This class handles systematic scans of geometric parameters to map
    potential energy surfaces and explore reaction pathways.

    Attributes:
        TYPE (str): Job type identifier ('orcascan').
        molecule: Molecule object used for the scan.
        settings: ORCAJobSettings configuration including scan parameters.
        label (str): Job identifier used for file naming.
        jobrunner: Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "orcascan"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        """
        Initialize ORCAScanJob.

        Args:
            molecule: Molecule object for the scan
            settings: ORCAJobSettings instance with scan parameters
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
