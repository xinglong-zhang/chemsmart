"""
xTB job runner implementation.

This module contains the XTBJobRunner class for executing xTB jobs.
It handles the setup of the execution environment, command generation,
and process management for xTB calculations.
"""

import logging
import os
import shlex
import subprocess
from contextlib import suppress
from functools import lru_cache
from glob import glob
from shutil import copy, rmtree

from chemsmart.io.xyz.xyzfile import XYZFile
from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.executable import XTBExecutable
from chemsmart.utils.periodictable import PeriodicTable

pt = PeriodicTable()

logger = logging.getLogger(__name__)


class XTBJobRunner(JobRunner):
    """
    Job runner for xTB calculations.

    Handles the execution of xTB jobs, including setting up the environment,
    writing input files, running the executable, and managing output files.

    Attributes:
        JOBTYPES (list): List of supported job types (e.g., 'xtbopt', 'xtbsp').
        PROGRAM (str): The program name ('xtb').
        FAKE (bool): Whether this is a fake runner (for testing).
    """

    # creates job runner process
    # combines information about server and program
    JOBTYPES = [
        "xtbopt",
        "xtbsp",
    ]

    PROGRAM = "xtb"

    FAKE = False

    def __init__(self, server, scratch=None, fake=False, **kwargs):
        """
        Initialize the XTBJobRunner.

        Args:
            server (Server): The server configuration.
            scratch (str, optional): Path to scratch directory.
            fake (bool, optional): Whether to run in fake mode.
            **kwargs: Additional arguments.
        """
        super().__init__(server=server, scratch=scratch, fake=fake, **kwargs)
        logger.debug(f"Jobrunner num cores: {self.num_cores}")
        logger.debug(f"Jobrunner num hours: {self.num_hours}")
        logger.debug(f"Jobrunner num gpus: {self.num_gpus}")
        logger.debug(f"Jobrunner mem gb: {self.mem_gb}")
        logger.debug(f"Jobrunner num threads: {self.num_threads}")

    @property
    @lru_cache(maxsize=12)
    def executable(self):
        """
        Get the xTB executable object.

        Returns:
            XTBExecutable: The executable object.

        Raises:
            TypeError: If the executable cannot be found.
        """
        try:
            executable_path = XTBExecutable().get_executable()
            if executable_path is None:
                import shutil

                executable_path = shutil.which("xtb")
                if executable_path is None:
                    raise TypeError(
                        "xtb executable not found in path or settings."
                    )

            executable_path = os.path.dirname(executable_path)
            executable = XTBExecutable(executable_path)
            logger.debug(f"Obtained xtb executable: {executable}")
            return executable
        except TypeError as e:
            logger.error(f"No xtb executable is found: {e}\n")
            raise

    def _prerun(self, job):
        """
        Prepare for job execution.

        Args:
            job (XTBJob): The job to run.
        """
        self._assign_variables(job)

    def _assign_variables(self, job):
        """
        Sets proper file paths for job input, output,
        and error files in scratch or not in scratch.

        Args:
            job (XTBJob): The job to run.
        """
        # keep job output file in job folder regardless of running in scratch or not
        self.job_outputfile = job.outputfile

        if self.scratch and self.scratch_dir:
            self._set_up_variables_in_scratch(job)
        else:
            self._set_up_variables_in_job_directory(job)

        if self.executable and self.executable.local_run is not None:
            logger.info(f"Local run is {self.executable.local_run}.")
            job.local = self.executable.local_run

    def _set_up_variables_in_scratch(self, job):
        """
        Set up variables when running in scratch directory.

        Args:
            job (XTBJob): The job to run.
        """
        scratch_job_dir = os.path.join(self.scratch_dir, str(job.label))
        if not os.path.exists(scratch_job_dir):
            with suppress(FileExistsError):
                os.makedirs(scratch_job_dir)
        self.running_directory = scratch_job_dir
        logger.debug(f"Running directory: {self.running_directory}")

        job_xyzfile = str(job.label) + ".xyz"
        scratch_job_xyzfile = (
            f"{os.path.join(str(scratch_job_dir), job_xyzfile)}"
        )
        self.job_xyzfile = os.path.abspath(scratch_job_xyzfile)  # type: ignore

        job_errfile = str(job.label) + ".err"
        scratch_job_errfile = (
            f"{os.path.join(str(scratch_job_dir), job_errfile)}"
        )
        self.job_errfile = os.path.abspath(scratch_job_errfile)  # type: ignore

    def _set_up_variables_in_job_directory(self, job):
        """
        Set up variables when running in the job directory (no scratch).

        Args:
            job (XTBJob): The job to run.
        """
        self.running_directory = job.folder
        logger.debug(f"Running directory: {self.running_directory}")
        self.job_xyzfile = os.path.abspath(job.xyzfile)
        self.job_errfile = os.path.abspath(job.errfile)

    def _write_input(self, job):
        """
        Write the geometry input file as .xyz for xtb jobs..

        Args:
            job (XTBJob): The job to run.
        """
        logger.info(f"Writting geometry input file to: {self.job_xyzfile}")
        job.molecule.write_xyz(self.job_xyzfile, mode="w")

    def _get_command(self, job=None, **kwargs):
        """
        Generate the command to run the job.

        Args:
            job (XTBJob, optional): The job object.
            **kwargs: Additional arguments.

        Returns:
            str: The command string.
        """
        exe = self._get_executable()
        command = f"{exe} {self.job_xyzfile} "
        return command

    def _create_process(self, job, command, env):
        """
        Create the subprocess to run the job.

        Args:
            job (XTBJob): The job object.
            command (str): The command to execute.
            env (dict): The environment variables.

        Returns:
            subprocess.Popen: The process object.
        """
        settings = job.settings
        # Convert multiplicity to uhf for xTB command line
        uhf = settings.multiplicity - 1
        if (
            settings.solvent_model is not None
            and settings.solvent_id is not None
        ):
            command += (
                f"--{settings.gfn_version} "
                f"--{settings.job_type} {settings.optimization_level} --chrg {settings.charge} "
                f"--uhf {uhf} --{settings.solvent_model} {settings.solvent_id}"
            )
        else:
            command += (
                f"--{settings.gfn_version} "
                f"--{settings.job_type} {settings.optimization_level} --chrg {settings.charge} "
                f"--uhf {uhf}"
            )

        # Add --grad flag if requested
        if settings.grad:
            command += " --grad"

        with (
            open(self.job_outputfile, "w") as out,
            open(self.job_errfile, "w") as err,
        ):
            logger.info(
                f"Command executed: {command}\n"
                f"Writing output file to: {self.job_outputfile}\n"
                f"And err file to: {self.job_errfile}"
            )
            logger.debug(f"Environments for running: {self.executable.env}")
            return subprocess.Popen(
                shlex.split(command),
                stdout=out,
                stderr=err,
                env=env,
                cwd=self.running_directory,
            )

    def _get_executable(self):
        """
        Get the path to the xTB executable.

        Returns:
            str: The executable path.
        """
        exe = self.executable.get_executable()
        logger.info(f"XTB executable: {exe}")
        return exe

    def _postrun(self, job, **kwargs):
        """
        Clean up after job execution.

        Args:
            job (XTBJob): The job that was run.
            **kwargs: Additional arguments.
        """
        if self.scratch:
            # if job was run in scratch, copy files to job folder except files starting with Gau-
            for file in glob(f"{self.running_directory}/{job.label}*"):
                if not file.startswith("Gau-"):
                    logger.info(
                        f"Copying file {file} from {self.running_directory} to {job.folder}"
                    )
                    copy(file, job.folder)

        if job.is_complete():
            # if job is completed, remove scratch directory and submit_script
            # and log.info and log.err files
            if self.scratch:
                logger.info(
                    f"Removing scratch directory: {self.running_directory}."
                )
                rmtree(self.running_directory)


class FakeXTBJobRunner(XTBJobRunner):
    """
    Fake job runner for testing purposes.

    Simulates the execution of xTB jobs without actually running the program.
    """

    # creates job runner process
    # combines information about server and program
    FAKE = True

    def __init__(self, server, scratch=None, fake=True, **kwargs):
        """
        Initialize the FakeXTBJobRunner.

        Args:
            server (Server): The server configuration.
            scratch (str, optional): Path to scratch directory.
            fake (bool, optional): Whether to run in fake mode.
            **kwargs: Additional arguments.
        """
        super().__init__(server=server, scratch=scratch, fake=fake, **kwargs)

    def run(self, job, **kwargs):
        """
        Simulate running the job.

        Args:
            job (XTBJob): The job to run.
            **kwargs: Additional arguments.
        """
        self._prerun(job=job)
        self._write_input(job=job)
        FakeXTB(self.job_xyzfile).run()
        self._postrun(job=job)

    def _set_up_variables_in_scratch(self, job):
        """
        Set up variables for fake run in scratch.

        Args:
            job (XTBJob): The job to run.
        """
        scratch_job_dir = os.path.join(self.scratch_dir, job.label)
        if not os.path.exists(scratch_job_dir):
            with suppress(FileExistsError):
                os.makedirs(scratch_job_dir)
        self.running_directory = scratch_job_dir
        logger.debug(f"Running directory: {self.running_directory}")

        # assign label with fake to differentiate from real job
        job.label = f"{job.label}_fake"
        logger.debug(f"Job label for fake job run: {job.label}")

        job_xyzfile = job.label + ".xyz"
        scratch_job_xyzfile = os.path.join(scratch_job_dir, job_xyzfile)
        self.job_xyzfile = os.path.abspath(scratch_job_xyzfile)

        job_errfile = job.label + ".err"
        scratch_job_errfile = os.path.join(scratch_job_dir, job_errfile)
        self.job_errfile = os.path.abspath(scratch_job_errfile)

    def _set_up_variables_in_job_directory(self, job):
        """
        Set up variables for fake run in job directory.

        Args:
            job (XTBJob): The job to run.
        """
        self.running_directory = job.folder
        logger.debug(f"Running directory: {self.running_directory}")
        job.label = f"{job.label}_fake"
        logger.debug(f"Job label for fake job run: {job.label}")
        self.job_xyzfile = os.path.abspath(job.xyzfile)
        self.job_errfile = os.path.abspath(job.errfile)


class FakeXTB:
    """
    Simulates the xTB program execution.
    """

    def __init__(self, file_to_run):
        """
        Initialize FakeXTB.

        Args:
            file_to_run (str): Path to the input file.
        """
        if not os.path.exists(file_to_run):
            raise FileNotFoundError(f"File {file_to_run} not found.")
        self.file_to_run = file_to_run
        self.input_object = XYZFile(filename=self.file_to_run)

    @property
    def file_folder(self):
        """Get the folder containing the file."""
        return os.path.dirname(self.file_to_run)

    @property
    def filename(self):
        """Get the filename."""
        return os.path.basename(self.file_to_run)

    @property
    def input_filepath(self):
        """Get the full input file path."""
        return self.file_to_run

    @property
    def output_filepath(self):
        """Get the output file path."""
        output_file = self.filename.split(".")[0] + ".log"
        return os.path.join(self.file_folder, output_file)

    @property
    def input_contents(self):
        """Get the contents of the input file."""
        return self.input_object.contents

    @property
    def molecule(self):
        """Get the molecule object."""
        return self.input_object.molecule

    @property
    def charge(self):
        """Get the charge."""
        return 0  # Default for XYZ

    @property
    def multiplicity(self):
        """Get the multiplicity."""
        return 1  # Default for XYZ

    @property
    def spin(self):
        """Get the spin state ('R' or 'U')."""
        if self.multiplicity == 1:
            return "R"
        return "U"

    @property
    def num_atoms(self):
        """Get the number of atoms."""
        return len(self.molecule.chemical_symbols)

    @property
    def atomic_symbols(self):
        """Get the atomic symbols."""
        return list(self.molecule.chemical_symbols)

    @property
    def atomic_numbers(self):
        """Get the atomic numbers."""
        return [pt.to_atomic_number(s) for s in self.atomic_symbols]

    @property
    def atomic_coordinates(self):
        """Get the atomic coordinates."""
        return self.molecule.positions

    @property
    def empirical_formula(self):
        """Get the empirical formula."""
        return self.molecule.get_chemical_formula(empirical=True)

    def run(self):
        """Simulate the execution and write output."""
        with open(self.output_filepath, "w") as g:
            g.write("-" * 100 + "\n")
            g.write("* xtb version 0.0.0 (Fake) \n")
            g.write("-" * 100 + "\n")
            g.write("\n")

            for line in self.input_contents:
                g.write(f" {line}\n")

            g.write("\n")
            g.write("* finished run (fake xtb) \n")
