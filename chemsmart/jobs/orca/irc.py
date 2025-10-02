"""
ORCA IRC job implementation.

This module contains the ORCAIRCJob class for running intrinsic reaction
coordinate calculations using ORCA.
"""

import logging

from chemsmart.jobs.orca.job import ORCAJob

logger = logging.getLogger(__name__)


class ORCAIRCJob(ORCAJob):
    """
    ORCA IRC (Intrinsic Reaction Coordinate) job.

    This class handles IRC calculations to trace reaction pathways from
    transition states to reactants and products.

    Attributes:
        TYPE (str): Job type identifier ('orcairc').
        molecule: Molecule object used for the IRC calculation.
        settings: ORCAIRCJobSettings configuration for the job.
        label (str): Job identifier used for file naming.
        jobrunner: Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "orcairc"

    def __init__(self, molecule, settings, label, jobrunner=None, **kwargs):
        """
        Initialize ORCAIRCJob.

        Args:
            molecule: Molecule object for the IRC calculation
            settings: ORCAIRCJobSettings instance
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
