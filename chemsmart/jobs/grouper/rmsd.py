"""
Basic RMSD grouper job.

This module provides the job class for basic RMSD grouping.
"""

from chemsmart.jobs.grouper.job import GrouperJob


class RMSDGrouperJob(GrouperJob):
    """Job for basic RMSD grouping."""

    TYPE = "rmsd_grouper"

    def __init__(self, molecules, threshold=0.5, **kwargs):
        super().__init__(
            molecules=molecules,
            grouping_strategy="rmsd",
            threshold=threshold,
            **kwargs,
        )
