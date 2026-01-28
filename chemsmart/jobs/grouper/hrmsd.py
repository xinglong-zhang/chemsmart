"""
Hungarian RMSD grouper job.

This module provides the job class for Hungarian RMSD grouping.
"""

from chemsmart.jobs.grouper.job import GrouperJob


class HRMSDGrouperJob(GrouperJob):
    """Job for Hungarian RMSD grouping."""

    TYPE = "hrmsd_grouper"

    def __init__(self, molecules, threshold=0.5, **kwargs):
        super().__init__(
            molecules=molecules,
            grouping_strategy="hrmsd",
            threshold=threshold,
            **kwargs,
        )
