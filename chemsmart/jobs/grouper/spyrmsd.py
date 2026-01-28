"""
SpyRMSD grouper job.

This module provides the job class for SpyRMSD grouping.
"""

from chemsmart.jobs.grouper.job import GrouperJob


class SpyRMSDGrouperJob(GrouperJob):
    """Job for SpyRMSD grouping with symmetry correction."""

    TYPE = "spyrmsd_grouper"

    def __init__(self, molecules, threshold=0.5, **kwargs):
        super().__init__(
            molecules=molecules,
            grouping_strategy="spyrmsd",
            threshold=threshold,
            **kwargs,
        )
