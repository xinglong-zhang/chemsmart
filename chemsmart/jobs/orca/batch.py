"""ORCA ``BatchJob`` for running a collection of ORCA jobs."""

from chemsmart.jobs.batch import BatchJob


class ORCABatchJob(BatchJob):
    """Batch controller for a collection of ORCA child jobs."""

    PROGRAM = "orca"
