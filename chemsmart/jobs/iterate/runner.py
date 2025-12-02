import logging
import os
from contextlib import suppress
from shutil import copy

from chemsmart.jobs.runner import JobRunner

logger = logging.getLogger(__name__)


class IterateJobRunner(JobRunner):
    """
    Job runner for Iterate jobs.
    
    Iterate jobs are special in that they don't call external programs.
    They run purely in Python to generate molecular structures.
    """

    JOBTYPES = ["iterate"]

    PROGRAM = "Iterate"
    FAKE = False
    SCRATCH = False

    def __init__(
        self, server, scratch=None, fake=False, scratch_dir=None, **kwargs
    ):
        if scratch is None:
            scratch = self.SCRATCH
        super().__init__(
            server=server,
            scratch=scratch,
            scratch_dir=scratch_dir,
            fake=fake,
            **kwargs,
        )
        logger.debug(f"Jobrunner server: {self.server}")
        logger.debug(f"Jobrunner scratch: {self.scratch}")
        logger.debug(f"Jobrunner fake mode: {self.fake}")
    
    def _prerun(self, job):
        self._assign_variables(job)

    def _assign_variables(self, job):
        """
        Assign necessary variables from the job to the runner.
        """
        self.job_outputfile = job.outputfile
        self.job_errorfile = job.errorfile

        if self.scratch and self.scratch_dir:
            self._set_up_variables_in_scratch(job)
        else:
            self._set_up_variables_in_job_directory(job)

        def _set_up_variables_in_scratch(self, job):
            pass

        def _set_up_variables_in_job_directory(self, job):
            pass

        def _write_input(self, job):
            if self.scratch and self.job_inputfile != job.inputfile:
                logger.info(
                    f"Copying input file {job.inputfile} to {self.job_inputfile}"
                )
                copy(job.inputfile, self.job_inputfile)
        
        def _get_command(self, job):
            return None
        
        