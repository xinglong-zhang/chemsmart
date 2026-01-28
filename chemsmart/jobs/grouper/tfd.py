"""
Torsion Fingerprint Deviation (TFD) grouper job.

This module provides the job class for TFD grouping.
"""

from chemsmart.jobs.grouper.job import GrouperJob


class TFDGrouperJob(GrouperJob):
    """Job for Torsion Fingerprint Deviation grouping."""

    TYPE = "tfd_grouper"

    def __init__(
        self,
        molecules,
        threshold=0.1,
        use_weights: bool = True,
        **kwargs,
    ):
        super().__init__(
            molecules=molecules,
            grouping_strategy="torsion",
            threshold=threshold,
            use_weights=use_weights,
            **kwargs,
        )
        self.use_weights = use_weights
