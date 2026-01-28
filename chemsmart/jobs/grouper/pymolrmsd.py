"""
PyMOL RMSD grouper job.

This module provides the job class for PyMOL-based RMSD grouping.
"""

from chemsmart.jobs.grouper.job import GrouperJob


class PymolRMSDGrouperJob(GrouperJob):
    """Job for PyMOL RMSD grouping."""

    TYPE = "pymolrmsd_grouper"

    def __init__(self, molecules, threshold=0.5, **kwargs):
        super().__init__(
            molecules=molecules,
            grouping_strategy="pymolrmsd",
            threshold=threshold,
            **kwargs,
        )
