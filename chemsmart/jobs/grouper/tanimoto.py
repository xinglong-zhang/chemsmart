"""
Tanimoto fingerprint similarity grouper job.

This module provides the job class for Tanimoto similarity grouping.
"""

from chemsmart.jobs.grouper.job import GrouperJob


class TanimotoGrouperJob(GrouperJob):
    """Job for Tanimoto fingerprint similarity grouping."""

    TYPE = "tanimoto_grouper"

    def __init__(
        self,
        molecules,
        threshold=0.9,
        fingerprint_type: str = "rdkit",
        **kwargs,
    ):
        super().__init__(
            molecules=molecules,
            grouping_strategy="tanimoto",
            threshold=threshold,
            fingerprint_type=fingerprint_type,
            **kwargs,
        )
        self.fingerprint_type = fingerprint_type
