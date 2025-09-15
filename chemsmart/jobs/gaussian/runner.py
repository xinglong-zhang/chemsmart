"""
Gaussian job runner for executing computational chemistry calculations.

This module provides the GaussianJobRunner class for managing the
execution of Gaussian calculations on various computing platforms.
It handles job submission, file management, scratch directory setup,
and result collection for Gaussian computational chemistry jobs.

The runner supports both local and remote execution with automatic
handling of checkpoint files, error recovery, and job monitoring.
"""

import logging
import os
import shlex
import subprocess
from contextlib import suppress
from datetime import datetime
from functools import lru_cache
from glob import glob
from random import random
from shutil import copy

from chemsmart.io.gaussian.input import Gaussian16Input
from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.executable import GaussianExecutable
from chemsmart.utils.periodictable import PeriodicTable

pt = PeriodicTable()

logger = logging.getLogger(__name__)


class GaussianJobRunner(JobRunner):
    """
    Job runner for executing Gaussian computational chemistry calculations.

    Manages the execution of Gaussian jobs on various computing platforms
    including local machines and remote servers. Handles file management,
    scratch directory setup, checkpoint file handling, and result collection.

    Supports multiple job types and provides automatic error recovery
    and job monitoring capabilities.

    Attributes:
        JOBTYPES (list): Supported Gaussian job type identifiers.
        PROGRAM (str): Program identifier ('gaussian').
        FAKE (bool): Whether to run in fake/test mode.
        SCRATCH (bool): Whether to use scratch directories by default.
        server: Server configuration used for execution.
        scratch (bool): Whether to use a scratch directory for I/O.
        scratch_dir (str): Path to the scratch directory when enabled.
        num_cores (int): Number of CPU cores allocated.
        num_gpus (int): Number of GPUs allocated.
        mem_gb (int): Memory allocation in gigabytes.
        num_hours (int): Maximum walltime from the server configuration.
        num_threads (int): Thread count from the server configuration.
        executable (GaussianExecutable): Executable configuration for the server.
        running_directory (str): Directory where the job is executed
            (scratch or job folder).
        job_inputfile (str): Full path to the input `.com` file for execution.
        job_outputfile (str): Full path to the output `.log` file.
        job_chkfile (str): Full path to the checkpoint `.chk` file.
        job_errfile (str): Full path to the error `.err` file.
    """

    # List of supported Gaussian job types
    JOBTYPES = [
        "g16crest",
        "g16job",
        "g16dias",
        "g16opt",
        "g16irc",
        "g16modred",
        "g16nci",
        "g16resp",
        "g16traj",
        "g16scan",
        "g16sp",
        "g16td",
        "g16ts",
        "g16uvvis",
        "g16wbi",
        "g16",
        "g16com",
        "g16link",
    ]

    PROGRAM = "gaussian"

    FAKE = False
    SCRATCH = True
    # Default to use scratch for Gaussian jobs - class attribute instead
    # of instance attribute so it needs not be set at
    # instance level - set during initialization (__init__).

    def __init__(
        self, server, scratch=None, fake=False, scratch_dir=None, **kwargs
    ):
        """
        Initialize Gaussian job runner with server and execution settings.

        Sets up the job runner with computational resources, scratch directory
        configuration, and execution mode. Inherits from base JobRunner class
        and applies Gaussian-specific defaults.

        Args:
            server: Server configuration object for job execution.
            scratch (bool, optional): Whether to use scratch directories.
                Defaults to class SCRATCH setting.
            fake (bool): Whether to run in fake/test mode.
            scratch_dir (str, optional): Custom scratch directory path.
            **kwargs: Additional arguments passed to parent JobRunner.
        """
        # Use default SCRATCH if scratch is not explicitly set
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
        logger.debug(f"Jobrunner num cores: {self.num_cores}")
        logger.debug(f"Jobrunner num hours: {self.num_hours}")
        logger.debug(f"Jobrunner num gpus: {self.num_gpus}")
        logger.debug(f"Jobrunner mem gb: {self.mem_gb}")
        logger.debug(f"Jobrunner num threads: {self.num_threads}")
        logger.debug(f"Jobrunner scratch: {self.scratch}")

    @property
    @lru_cache(maxsize=12)
    def executable(self):
        """
        Get the Gaussian executable configuration for the current server.

        Retrieves the appropriate Gaussian executable settings based on
        the server configuration. Uses LRU cache for performance.

        Returns:
            GaussianExecutable: Configured executable object for the server.

        Raises:
            FileNotFoundError: If server configuration file is not found.
        """
        try:
            logger.info(
                f"Obtaining executable from server: {self.server.name}"
            )
            executable = GaussianExecutable.from_servername(
                servername=self.server.name
            )
            return executable
        except FileNotFoundError as e:
            logger.error(
                f"No server file {self.server} is found: {e}\n"
                f"Available servers are: {GaussianExecutable.available_servers}"
            )
            raise

    def _prerun(self, job):
        """
        Perform pre-execution setup for the Gaussian job.

        Handles all necessary preparation before job execution including
        file path assignment and directory setup.

        Args:
            job: Job object containing calculation parameters and settings.
        """
        self._assign_variables(job)

    def _assign_variables(self, job):
        """
        Set up file paths for job execution in scratch or job directory.

        Configures appropriate file paths for input, output, checkpoint,
        and error files based on whether scratch directories are used.
        Also handles local execution settings.

        Args:
            job: Job object to configure file paths for.
        """
        if self.scratch and self.scratch_dir:
            logger.debug("Setting up job variables in scratch directory")
            self._set_up_variables_in_scratch(job)
        else:
            logger.debug("Setting up job variables in job directory")
            self._set_up_variables_in_job_directory(job)

        if self.executable and self.executable.local_run is not None:
            logger.info(f"Local run is {self.executable.local_run}.")
            job.local = self.executable.local_run

    def _set_up_variables_in_scratch(self, job):
        """
        Configure file paths for execution in a scratch directory.

        Sets up a dedicated scratch directory for the job and configures
        all file paths (input, output, checkpoint, error) to point to
        the scratch location for improved I/O performance.

        Args:
            job: Job object to configure scratch paths for.
        """
        scratch_job_dir = os.path.join(self.scratch_dir, job.label)
        if not os.path.exists(scratch_job_dir):
            with suppress(FileExistsError):
                os.makedirs(scratch_job_dir)
                logger.debug(f"Created scratch directory: {scratch_job_dir}")
        self.running_directory = scratch_job_dir
        logger.debug(f"Running directory: {self.running_directory}")

        job_inputfile = job.label + ".com"
        scratch_job_inputfile = os.path.join(scratch_job_dir, job_inputfile)
        self.job_inputfile = os.path.abspath(scratch_job_inputfile)

        job_chkfile = job.label + ".chk"
        scratch_job_chkfile = os.path.join(scratch_job_dir, job_chkfile)
        self.job_chkfile = os.path.abspath(scratch_job_chkfile)

        job_errfile = job.label + ".err"
        scratch_job_errfile = os.path.join(scratch_job_dir, job_errfile)
        self.job_errfile = os.path.abspath(scratch_job_errfile)

        job_outputfile = job.label + ".log"
        scratch_job_outputfile = os.path.join(scratch_job_dir, job_outputfile)
        self.job_outputfile = os.path.abspath(scratch_job_outputfile)

    def _set_up_variables_in_job_directory(self, job):
        """
        Configure file paths for execution in the job directory.

        Sets up all file paths to point to the job's directory instead
        of using a separate scratch directory for execution.

        Args:
            job: Job object to configure directory paths for.
        """
        self.running_directory = job.folder
        logger.debug(f"Running directory: {self.running_directory}")
        self.job_inputfile = os.path.abspath(job.inputfile)
        self.job_chkfile = os.path.abspath(job.chkfile)
        self.job_errfile = os.path.abspath(job.errfile)
        self.job_outputfile = os.path.abspath(job.outputfile)

    def _write_input(self, job):
        """
        Write the Gaussian input file to the execution directory.

        Creates a GaussianInputWriter instance and uses it to generate
        the input file (.com) in the appropriate directory for execution.

        Args:
            job: Job object containing calculation parameters to write.
        """
        from chemsmart.jobs.gaussian.writer import GaussianInputWriter

        input_writer = GaussianInputWriter(job=job)
        input_writer.write(target_directory=self.running_directory)

    def _get_command(self, job):
        """
        Generate the command string for executing the Gaussian job.

        Constructs the complete command line including the executable
        path and input file path for running the calculation.

        Args:
            job: Job object (used implicitly via self.job_inputfile).

        Returns:
            str: Complete command string for job execution.
        """
        exe = self._get_executable()
        command = f"{exe} {self.job_inputfile}"
        return command

    def _create_process(self, job, command, env):
        """
        Create and start a subprocess for executing the Gaussian job.

        Opens output and error files, then starts the Gaussian process
        with appropriate environment variables and working directory.

        Args:
            job: Job object containing execution parameters.
            command: Command string to execute.
            env: Environment variables for the process.

        Returns:
            subprocess.Popen: The created process object.
        """
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
        Get the executable path for Gaussian.

        Retrieves the configured Gaussian executable path and logs
        the information for debugging purposes.

        Returns:
            str: Path to the Gaussian executable.
        """
        exe = self.executable.get_executable()
        logger.info(f"Gaussian executable: {exe}")
        return exe

    def _postrun(self, job):
        """
        Perform post-execution cleanup and file management.

        If using scratch directories, copies output files from scratch
        back to the job directory, excluding temporary Gaussian files.

        Args:
            job: Job object containing folder and label information.
        """
        if self.scratch:
            # if job was run in scratch, copy files to job folder except
            # files starting with Gau-
            for file in glob(f"{self.running_directory}/{job.label}*"):
                if not file.startswith("Gau-"):
                    logger.info(
                        f"Copying file {file} from {self.running_directory} "
                        f"to {job.folder}"
                    )
                    try:
                        copy(file, job.folder)
                    except Exception as e:
                        logger.error(
                            f"File {file} cannot be copied to job folder "
                            f"{job.folder}: {e}"
                        )

        if job.is_complete():
            # if job is completed, remove scratch directory and submit_script
            # and log.info and log.err files
            # if self.scratch:
            #     logger.info(
            #         f"Removing scratch directory: {self.running_directory}."
            #     )
            #     rmtree(self.running_directory)

            self._remove_err_files(job)


class FakeGaussianJobRunner(GaussianJobRunner):
    """
    Fake Gaussian job runner for testing and simulation purposes.

    Provides a simulation environment for Gaussian calculations without
    actually running the Gaussian software, useful for testing workflows
    and development purposes.

    Attributes:
        FAKE (bool): Flag indicating fake/test mode (True for this runner).
        JOBTYPES (list): Supported Gaussian job types (inherited).
        PROGRAM (str): Program identifier ('gaussian').
        SCRATCH (bool): Whether to use scratch directories by default (inherited).
        server: Server configuration used for execution (inherited).
        scratch (bool): Whether to use a scratch directory (inherited).
        scratch_dir (str): Path to the scratch directory (inherited/derived).
        running_directory (str): Directory where the fake job is executed.
        job_inputfile (str): Full path to the input `.com` file.
        job_outputfile (str): Full path to the output `.log` file (created by fake run).
        job_chkfile (str): Full path to the checkpoint `.chk` file (path assigned).
        job_errfile (str): Full path to the error `.err` file (path assigned).
    """

    # creates job runner process
    # combines information about server and program
    FAKE = True

    def __init__(
        self, server, scratch=None, fake=True, scratch_dir=None, **kwargs
    ):
        """
        Initialize the fake Gaussian job runner.

        Sets up a simulation environment that mimics the behavior of
        the real GaussianJobRunner without executing actual calculations.

        Args:
            server: Server configuration for the fake runner.
            scratch: Whether to use scratch directories (default: None).
            fake: Flag indicating this is a fake runner (default: True).
            scratch_dir: Path to scratch directory (default: None).
            **kwargs: Additional arguments passed to parent class.
        """
        super().__init__(
            server=server,
            scratch=scratch,
            scratch_dir=scratch_dir,
            fake=fake,
            **kwargs,
        )

    def run(self, job, **kwargs):
        """
        Execute a fake Gaussian job simulation.

        Performs all the steps of a real job execution but uses a fake
        Gaussian process instead of the actual software for testing.

        Args:
            job: Job object to simulate execution for.
            **kwargs: Additional arguments (unused in fake execution).

        Returns:
            int: Return code from the fake Gaussian execution.
        """
        self._prerun(job=job)
        self._write_input(job=job)
        returncode = FakeGaussian(self.job_inputfile).run()
        self._postrun(job=job)
        return returncode

    def _set_up_variables_in_scratch(self, job):
        """
        Configure fake job file paths in scratch directory.

        Sets up scratch directory paths for the fake job execution,
        including special labeling to distinguish from real jobs.

        Args:
            job: Job object to configure fake scratch paths for.
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

        job_inputfile = job.label + ".com"
        scratch_job_inputfile = os.path.join(scratch_job_dir, job_inputfile)
        self.job_inputfile = os.path.abspath(scratch_job_inputfile)

        job_chkfile = job.label + ".chk"
        scratch_job_chkfile = os.path.join(scratch_job_dir, job_chkfile)
        self.job_chkfile = os.path.abspath(scratch_job_chkfile)

        job_errfile = job.label + ".err"
        scratch_job_errfile = os.path.join(scratch_job_dir, job_errfile)
        self.job_errfile = os.path.abspath(scratch_job_errfile)

    def _set_up_variables_in_job_directory(self, job):
        """
        Configure fake job file paths in job directory.

        Sets up file paths for fake job execution in the job's directory,
        with special labeling to distinguish from real jobs.

        Args:
            job: Job object to configure fake directory paths for.
        """
        self.running_directory = job.folder
        logger.debug(f"Running directory: {self.running_directory}")
        job.label = f"{job.label}_fake"
        logger.debug(f"Job label for fake job run: {job.label}")
        self.job_inputfile = os.path.abspath(job.inputfile)
        self.job_chkfile = os.path.abspath(job.chkfile)
        self.job_errfile = os.path.abspath(job.errfile)


class FakeGaussian:
    """
    Fake Gaussian simulation class for testing purposes.

    Simulates Gaussian software behavior by parsing input files and
    generating corresponding output files without actual calculations.
    Useful for testing workflows and development.

    Attributes:
        file_to_run (str): Path to the Gaussian input file to simulate.
        input_object (Gaussian16Input): Parsed representation of the
            input file, providing access to parsed blocks and metadata.

    Properties:
        file_folder (str): Directory path of the input file.
        filename (str): Basename of the input file.
        input_filepath (str): Full path to the input file.
        output_filepath (str): Path to the output log file produced by the
            fake run (same folder, `.log` extension).
        input_contents (str): Entire contents of the input file.
        input_blocks (dict): Parsed input content groups.
        molecule (Molecule): Molecular structure parsed from input.
        charge (int): Molecular charge.
        multiplicity (int): Spin multiplicity.
        spin (str): 'R' for singlet (restricted), 'U' otherwise.
        num_atoms (int): Number of atoms in the molecule.
        atomic_symbols (list[str]): Atomic symbols for all atoms.
        atomic_numbers (list[int]): Atomic numbers for all atoms.
        atomic_coordinates (array): Cartesian coordinates of atoms.
        empirical_formula (str): Empirical chemical formula string.
    """

    def __init__(self, file_to_run):
        """
        Initialize the fake Gaussian simulator.

        Sets up the simulation by parsing the input file and preparing
        for fake execution without running actual calculations.

        Args:
            file_to_run: Path to the Gaussian input file to simulate.

        Raises:
            FileNotFoundError: If the input file doesn't exist.
        """
        if not os.path.exists(file_to_run):
            raise FileNotFoundError(f"File {file_to_run} not found.")
        self.file_to_run = file_to_run
        self.input_object = Gaussian16Input(filename=self.file_to_run)

    @property
    def file_folder(self):
        """
        Get the directory containing the input file.

        Returns:
            str: Directory path of the input file.
        """
        return os.path.dirname(self.file_to_run)

    @property
    def filename(self):
        """
        Get the basename of the input file.

        Returns:
            str: Filename without directory path.
        """
        return os.path.basename(self.file_to_run)

    @property
    def input_filepath(self):
        """
        Get the full path to the input file.

        Returns:
            str: Complete path to the input file.
        """
        return self.file_to_run

    @property
    def output_filepath(self):
        """
        Get the path for the output log file.

        Constructs the output filename by taking the basename up to the
        first '.' and appending '.log' in the same directory.

        Returns:
            str: Path to the output log file.
        """
        output_file = self.filename.split(".")[0] + ".log"
        return os.path.join(self.file_folder, output_file)

    @property
    def input_contents(self):
        """
        Get the lines of the input file.

        Returns:
            list[str]: Lines of the input file.
        """
        return self.input_object.contents

    @property
    def input_blocks(self):
        """
        Get the parsed content blocks from the input file.

        Returns:
            list: Parsed content blocks; [0]=route, [1]=title,
                [2]=charge/coords, [3:]=additional sections.
        """
        return self.input_object.content_groups

    @property
    def molecule(self):
        """
        Get the molecular structure from the input.

        Returns:
            Molecule: Molecular structure object.
        """
        return self.input_object.molecule

    @property
    def charge(self):
        """
        Get the molecular charge from the input.

        Returns:
            int: Molecular charge value.
        """
        return self.input_object.charge

    @property
    def multiplicity(self):
        """
        Get the spin multiplicity from the input.

        Returns:
            int: Spin multiplicity value.
        """
        return self.input_object.multiplicity

    @property
    def spin(self):
        """
        Determine spin type based on multiplicity.

        Returns restricted (R) for singlets, unrestricted (U) for others.

        Returns:
            str: 'R' for restricted, 'U' for unrestricted calculations.
        """
        if self.multiplicity == 1:
            return "R"
        return "U"

    @property
    def num_atoms(self):
        """
        Get the number of atoms in the molecule.

        Returns:
            int: Number of atoms in the molecular structure.
        """
        return len(self.molecule.chemical_symbols)

    @property
    def atomic_symbols(self):
        """
        Get the atomic symbols of all atoms in the molecule.

        Returns:
            list: List of atomic symbols (e.g., ['C', 'H', 'O']).
        """
        return list(self.molecule.chemical_symbols)

    @property
    def atomic_numbers(self):
        """
        Get the atomic numbers of all atoms in the molecule.

        Returns:
            list: List of atomic numbers corresponding to the atoms.
        """
        return [pt.to_atomic_number(s) for s in self.atomic_symbols]

    @property
    def atomic_coordinates(self):
        """
        Get the Cartesian coordinates of all atoms.

        Returns:
            array: Array of atomic coordinates in Cartesian space.
        """
        return self.molecule.positions

    @property
    def empirical_formula(self):
        """
        Get the empirical formula of the molecule.

        Returns:
            str: Empirical formula string (e.g., 'H2O', 'CH4').
        """
        return self.molecule.get_chemical_formula(empirical=True)

    def run(self):
        """
        Execute the fake Gaussian simulation.

        Creates a fake output file that mimics the structure and content
        of a real Gaussian output log, including job information and
        basic molecular data for testing purposes.

        Returns:
            None
        """
        with open(self.output_filepath, "w") as g:
            g.write(" Entering Gaussian System, FakeGaussianRunner\n")
            g.write(f" Input={self.input_filepath}\n")
            g.write(f" Output={self.output_filepath}\n")
            g.write(" ******************************************\n")
            g.write(" Fake Gaussian Executable\n")
            g.write(" ******************************************\n")
            # write mem/nproc/chk information (%...)
            for line in self.input_contents:
                if line.startswith("%"):
                    g.write(f" {line}\n")
            # write route information
            for line in self.input_contents:
                if line.startswith("#"):
                    line_len = len(line)
                    g.write(" " + "-" * line_len + "\n")
                    g.write(f" {line}\n")
                    g.write(" " + "-" * line_len + "\n")

            # missing Z-matrix
            # write title
            g.write(" " + "-" * line_len + "\n")
            g.write(f" {self.input_blocks[1][0]}\n")
            g.write(" " + "-" * line_len + "\n")

            g.write(" Symbolic Z-matrix:\n")

            # write charge and multiplicity
            g.write(
                f" Charge =  {self.charge} Multiplicity = {self.multiplicity}\n"
            )
            for line in self.input_blocks[2][1:]:
                g.write(f" {line}\n")
            g.write("\n")

            if self.input_blocks[3:]:
                for line in self.input_blocks[3:]:
                    g.write(f"{line}\n")

            # missing Z-matrix

            # write input orientation
            g.write(
                "                          Input orientation:                             \n"
            )
            g.write(
                " ------------------------------------------------------------------------\n"
            )
            g.write(
                " Center     Atomic      Atomic             Coordinates (Angstroms)       \n"
            )
            g.write(
                " Number     Number       Type             X           Y           Z      \n"
            )
            g.write(
                " ------------------------------------------------------------------------\n"
            )

            # write coordinates
            for i in range(self.num_atoms):
                g.write(
                    f"{i + 1:>4} {self.atomic_numbers[i]:>10} {0!s:>10} "
                    f"{self.atomic_coordinates[i][0]:>20.6}"
                    f"{self.atomic_coordinates[i][1]:>14.6}"
                    f"{self.atomic_coordinates[i][2]:>14.6}\n"
                )
            g.write(
                " ------------------------------------------------------------------------\n"
            )
            g.write(f" Stoichiometry    {self.empirical_formula}")

            # write standard orientation
            g.write(
                "                       Standard orientation:                             \n"
            )
            g.write(
                " ------------------------------------------------------------------------\n"
            )
            g.write(
                " Center     Atomic      Atomic             Coordinates (Angstroms)       \n"
            )
            g.write(
                " Number     Number       Type             X           Y           Z      \n"
            )
            g.write(
                " ------------------------------------------------------------------------\n"
            )
            # write coordinates
            for i in range(self.num_atoms):
                g.write(
                    f"{i + 1:>4} {self.atomic_numbers[i]:>10} {0!s:>10} "
                    f"{self.atomic_coordinates[i][0]:>20.6}"
                    f"{self.atomic_coordinates[i][1]:>14.6}"
                    f"{self.atomic_coordinates[i][2]:>14.6}\n"
                )
            g.write(
                " ------------------------------------------------------------------------\n"
            )
            g.write(
                " ! Dummy Rotational constant values and dummy SCF energy value...\n"
            )
            g.write(
                " Rotational constants (GHZ):          11.4493930           9.4805599           5.3596246\n"
            )  # not real values
            g.write(f" Standard basis: {self.input_object.basis} (5D, 7F)\n")
            g.write(f" NAtoms=    {self.num_atoms}\n")
            g.write(
                f" SCF Done:  E({self.spin}{self.input_object.functional.upper()}) =  "
                f"-{1000 * random()}  A.U. after   14 cycles\n"
            )  # dummy energy
            g.write(" Mulliken charges:\n")
            g.write("               1\n")
            for i in range(self.num_atoms):
                g.write(
                    f"{i + 1:>7} {self.atomic_symbols[i]:>3} {random():>12.6}\n"
                )  # not real values
            g.write(" Elapsed time: xx\n")
            g.write(
                f" Normal termination of Gaussian 16 (fake executable) at {datetime.now()}."
            )
