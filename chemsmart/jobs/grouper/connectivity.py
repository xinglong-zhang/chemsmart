"""
Connectivity-based grouper job.

This module provides the job class for connectivity grouping.
"""

from chemsmart.jobs.grouper.job import GrouperJob


class ConnectivityGrouperJob(GrouperJob):
    """Job for connectivity-based grouping."""

    TYPE = "connectivity_grouper"

    def __init__(self, molecules, **kwargs):
        super().__init__(
            molecules=molecules,
            grouping_strategy="connectivity",
            **kwargs,
        )
