"""
Formula-based grouper job.

This module provides the job class for formula grouping.
"""

from chemsmart.jobs.grouper.job import GrouperJob


class FormulaGrouperJob(GrouperJob):
    """Job for formula-based grouping."""

    TYPE = "formula_grouper"

    def __init__(self, molecules, **kwargs):
        super().__init__(
            molecules=molecules,
            grouping_strategy="formula",
            **kwargs,
        )
