"""
Job runner for thermochemistry analysis jobs.

This module provides the ThermochemistryJobRunner class for managing
the execution of thermochemistry calculations, handling file operations,
and managing computational resources for thermochemical analysis.
"""

import logging
import os
from contextlib import suppress
from shutil import copy

from chemsmart.jobs.runner import JobRunner

logger = logging.getLogger(__name__)


class ThermochemistryJobRunner(JobRunner):
    """
    Job runner for thermochemistry analysis calculations.
    
    This class manages the execution environment for thermochemistry jobs,
    handling file operations, scratch directory management, and process
    execution for both single and Boltzmann-averaged calculations.

    Attributes:
        JOBTYPES (list): Supported job types (['thermochemistry', 'boltzmann']).
        PROGRAM (str): Program identifier ('Thermochemistry').
        FAKE (bool): Whether this runner operates in fake/test mode.
        SCRATCH (bool): Whether to use scratch directories by default.
        server: Server configuration used for execution.
        scratch (bool): Whether scratch is enabled for this runner.
        scratch_dir (str): Path to scratch directory, if used.
        num_cores (int): Number of CPU cores allocated (from server).
        num_gpus (int): Number of GPUs allocated (from server).
        mem_gb (int): Memory allocation in gigabytes (from server).
    """
    
    JOBTYPES = ["thermochemistry", "boltzmann"]
    
    PROGRAM = "Thermochemistry"
    FAKE = False
    # Thermochemistry jobs are lightweight, scratch typically not needed
    SCRATCH = False

    def __init__(
        self, server, scratch=None, fake=False, scratch_dir=None, **kwargs
    ):
        """
        Initialize the thermochemistry job runner.
        
        Args:
            server: Server configuration for job execution
            scratch (bool, optional): Whether to use scratch directory
            fake (bool): Whether to run in fake mode for testing
            scratch_dir (str, optional): Path to scratch directory
            **kwargs: Additional keyword arguments for parent class
        """
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
        """
        Prepare the job environment before execution.
        
        Sets up file paths and directory structure required for
        the thermochemistry calculation.
        
        Args:
            job: Thermochemistry job instance to prepare
        """
        self._assign_variables(job)

    def _assign_variables(self, job):
        """
        Set up file paths for input, output, and error files.
        
        Configures the necessary file paths for job execution,
        choosing between scratch directory or job directory
        based on runner configuration.
        
        Args:
            job: Thermochemistry job instance to configure
        """
        self.job_outputfile = job.outputfile
        self.job_errfile = job.errfile

        if self.scratch and self.scratch_dir:
            self._set_up_variables_in_scratch(job)
        else:
            self._set_up_variables_in_job_directory(job)

    def _set_up_variables_in_scratch(self, job):
        """
        Set up file paths in a scratch directory.
        
        Creates a dedicated scratch directory for the job and
        configures all file paths to use the scratch location.
        
        Args:
            job: Thermochemistry job instance to configure
        """
        # Create scratch job directory
        scratch_job_dir = os.path.join(self.scratch_dir, job.label)
        if not os.path.exists(scratch_job_dir):
            with suppress(FileExistsError):
                os.makedirs(scratch_job_dir)
        self.running_directory = scratch_job_dir
        logger.debug(f"Running directory: {self.running_directory}")

        # Set up input file path in scratch
        job_inputfile = os.path.basename(job.inputfile)
        scratch_job_inputfile = os.path.join(scratch_job_dir, job_inputfile)
        self.job_inputfile = os.path.abspath(scratch_job_inputfile)

        # Set up output file path in scratch
        job_outfile = job.label + ".dat"
        scratch_job_outfile = os.path.join(scratch_job_dir, job_outfile)
        self.job_outputfile = os.path.abspath(scratch_job_outfile)

        # Set up error file path in scratch
        job_errfile = job.label + ".err"
        scratch_job_errfile = os.path.join(scratch_job_dir, job_errfile)
        self.job_errfile = os.path.abspath(scratch_job_errfile)

    def _set_up_variables_in_job_directory(self, job):
        """
        Set up file paths in the job's directory.
        
        Configures file paths to use the job's working directory
        for input, output, and error files.
        
        Args:
            job: Thermochemistry job instance to configure
        """
        self.running_directory = job.folder
        logger.debug(f"Running directory: {self.running_directory}")
        self.job_inputfile = os.path.abspath(job.inputfile)
        self.job_outputfile = os.path.abspath(job.outputfile)
        self.job_errfile = os.path.abspath(job.errfile)

    def _write_input(self, job):
        """
        Copy the input file to the running directory if using scratch.
        
        Ensures the input file is available in the working directory
        for thermochemistry analysis, copying from the original location
        if running in scratch mode.
        
        Args:
            job: Thermochemistry job instance with input file
        """
        if self.scratch and self.job_inputfile != job.inputfile:
            logger.info(
                f"Copying input file {job.inputfile} to {self.job_inputfile}"
            )
            copy(job.inputfile, self.job_inputfile)

    def _get_command(self, job):
        """
        Get the command to execute for thermochemistry jobs.
        
        Thermochemistry jobs are executed directly in Python,
        so no external command is needed.
        
        Args:
            job: Thermochemistry job instance
            
        Returns:
            None: No external command required
        """
        return None

    def _create_process(self, job, command, env):
        """
        Run the thermochemistry calculation directly.
        
        Executes the thermochemistry calculation using the job's
        compute method and handles any errors that occur during
        the calculation process.
        
        Args:
            job: Thermochemistry job instance to execute
            command: Command to execute (unused for thermochemistry)
            env: Environment variables (unused for thermochemistry)
            
        Returns:
            int: Exit code (0 for success, 1 for failure)
        """
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
        """
        Run the thermochemistry job.
        
        For thermochemistry jobs, the actual calculation is performed
        in _create_process, so this method is a no-op.
        
        Args:
            process: Process result from _create_process
            **kwargs: Additional keyword arguments (unused)
        """
        pass

    def _get_executable(self):
        """
        Get the executable for thermochemistry jobs.
        
        Thermochemistry jobs run directly in Python without
        external executables.
        
        Returns:
            None: No external executable required
        """
        return None

    def _postrun(self, job):
        """
        Handle post-run tasks, including file management and cleanup.
        
        Manages file copying from scratch directories and cleanup
        operations after job completion. Removes error files for
        successful jobs.
        
        Args:
            job: Thermochemistry job instance that was executed
        """
        
        if self.scratch:
            # Copy output and error files to job folder
            for file in [self.job_outputfile, self.job_errfile]:
                if not os.path.exists(file):
                    logger.info(f"Copying file {file} to {job.folder}")
                    copy(file, job.folder)

            # Clean up scratch directory
            # if os.path.exists(self.running_directory):
            #     logger.info(
            #         f"Removing scratch directory: {self.running_directory}"
            #     )
            #     rmtree(self.running_directory)

        if job.is_complete():
            self._remove_err_files(job)
