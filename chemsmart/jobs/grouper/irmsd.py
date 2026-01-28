"""
Invariant RMSD (iRMSD) grouper job.

This module provides the job class for iRMSD grouping.
"""

from chemsmart.jobs.grouper.job import GrouperJob


class IRMSDGrouperJob(GrouperJob):
    """Job for invariant RMSD grouping."""

    TYPE = "irmsd_grouper"

    def __init__(
        self,
        molecules,
        threshold=0.125,
        check_stereo: str = "auto",
        **kwargs,
    ):
        super().__init__(
            molecules=molecules,
            grouping_strategy="irmsd",
            threshold=threshold,
            check_stereo=check_stereo,
            **kwargs,
        )
        self.check_stereo = check_stereo
