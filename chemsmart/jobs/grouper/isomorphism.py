"""
RDKit isomorphism grouper job.

This module provides the job class for isomorphism grouping.
"""

from chemsmart.jobs.grouper.job import GrouperJob


class IsomorphismGrouperJob(GrouperJob):
    """Job for RDKit isomorphism grouping."""

    TYPE = "isomorphism_grouper"

    def __init__(self, molecules, **kwargs):
        super().__init__(
            molecules=molecules,
            grouping_strategy="isomorphism",
            **kwargs,
        )
