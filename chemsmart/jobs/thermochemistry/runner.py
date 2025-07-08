import logging
import os
from contextlib import suppress
from shutil import copy, rmtree

from chemsmart.jobs.runner import JobRunner

logger = logging.getLogger(__name__)


class ThermochemistryJobRunner(JobRunner):
    """Job runner for thermochemistry analysis jobs."""

    JOBTYPES = ["thermochemistry"]

    PROGRAM = "Thermochemistry"
    FAKE = False
    SCRATCH = False  # Thermochemistry jobs are lightweight, so scratch is not typically needed

    def __init__(
        self, server, scratch=None, fake=False, scratch_dir=None, **kwargs
    ):
        """Initialize the ThermochemistryJobRunner."""
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
        """Prepare the job environment before running."""
        self._assign_variables(job)

    def _assign_variables(self, job):
        """Set up file paths for input, output, and error files."""
        self.job_outputfile = job.outputfile
        self.job_errfile = job.errfile

        if self.scratch and self.scratch_dir:
            self._set_up_variables_in_scratch(job)
        else:
            self._set_up_variables_in_job_directory(job)

    def _set_up_variables_in_scratch(self, job):
        """Set up file paths in a scratch directory."""
        scratch_job_dir = os.path.join(self.scratch_dir, job.label)
        if not os.path.exists(scratch_job_dir):
            with suppress(FileExistsError):
                os.makedirs(scratch_job_dir)
        self.running_directory = scratch_job_dir
        logger.debug(f"Running directory: {self.running_directory}")

        job_inputfile = os.path.basename(job.inputfile)
        scratch_job_inputfile = os.path.join(scratch_job_dir, job_inputfile)
        self.job_inputfile = os.path.abspath(scratch_job_inputfile)

        job_outfile = job.label + ".dat"
        scratch_job_outfile = os.path.join(scratch_job_dir, job_outfile)
        self.job_outputfile = os.path.abspath(scratch_job_outfile)

        job_errfile = job.label + ".err"
        scratch_job_errfile = os.path.join(scratch_job_dir, job_errfile)
        self.job_errfile = os.path.abspath(scratch_job_errfile)

    def _set_up_variables_in_job_directory(self, job):
        """Set up file paths in the job's directory."""
        self.running_directory = job.folder
        logger.debug(f"Running directory: {self.running_directory}")
        self.job_inputfile = os.path.abspath(job.inputfile)
        self.job_outputfile = os.path.abspath(job.outputfile)
        self.job_errfile = os.path.abspath(job.errfile)

    def _write_input(self, job):
        """Copy the input file to the running directory if using scratch."""
        if self.scratch and self.job_inputfile != job.inputfile:
            logger.info(
                f"Copying input file {job.inputfile} to {self.job_inputfile}"
            )
            copy(job.inputfile, self.job_inputfile)

    def _get_command(self, job):
        """No external command is needed for thermochemistry jobs."""
        return None

    def _create_process(self, job, command, env):
        """Run the thermochemistry calculation directly."""
        try:
            job.compute_thermochemistry()
            return 0  # Return 0 to indicate success
        except Exception as e:
            with open(self.job_errfile, "w") as err:
                err.write(
                    f"Error during thermochemistry calculation: {str(e)}\n"
                )
            logger.error(f"Error processing job {job.label}: {str(e)}")
            return 1  # Return 1 to indicate failure

    def _run(self, process, **kwargs):
        """Run the thermochemistry job."""
        pass

    def _get_executable(self):
        """No external executable is needed for thermochemistry jobs."""
        return None

    def _postrun(self, job):
        """Handle post-run tasks, such as copying files from scratch and cleanup."""
        if self.scratch:
            # Copy output and error files to job folder
            for file in [self.job_outputfile, self.job_errfile]:
                if not os.path.exists(file):
                    logger.info(f"Copying file {file} to {job.folder}")
                    copy(file, job.folder)

            # Clean up scratch directory
            if os.path.exists(self.running_directory):
                logger.info(
                    f"Removing scratch directory: {self.running_directory}"
                )
                rmtree(self.running_directory)

        if job.is_complete():
            self._remove_err_files(job)
