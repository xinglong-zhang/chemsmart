import logging
from abc import abstractmethod

logger = logging.getLogger(__name__)


class InputWriter:
    """Class that writes input files for a job.

    Args:
        job (Job): The job for which the input file is being written.
        jobrunner (JobRunner): The job runner for the job.
    """

    def __init__(self, job, jobrunner):
        self.job = job
        self.jobrunner = jobrunner

    @abstractmethod
    def write(self):
        """Write the input file for the job."""
        raise NotImplementedError
