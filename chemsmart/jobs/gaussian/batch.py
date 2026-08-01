"""Gaussian ``BatchJob`` for running a collection of Gaussian jobs."""

from chemsmart.jobs.batch import BatchJob


class GaussianBatchJob(BatchJob):
    """Batch controller for a collection of Gaussian child jobs."""

    PROGRAM = "gaussian"
