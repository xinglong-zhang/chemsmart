"""
NCIPLOT job runner implementation.

This module contains job runners for executing NCIPLOT calculations,
including both real and fake execution modes for testing purposes.
"""

import logging
import os
import shlex
import subprocess
from contextlib import suppress
from datetime import datetime
from functools import lru_cache
from glob import glob
from shutil import SameFileError, copy

from chemsmart.io.converter import FileConverter
from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.executable import NCIPLOTExecutable

logger = logging.getLogger(__name__)


class NCIPLOTJobRunner(JobRunner):
    """
    Job runner for NCIPLOT jobs.

    This class handles the execution of NCIPLOT calculations, including
    file preparation, job execution, and post-processing.

    Attributes:
        PROGRAM (str): Program identifier ('NCIPLOT').
        JOBTYPES (list): Supported job types (['nciplot']).
        FAKE (bool): Whether this runner operates in fake/test mode.
        SCRATCH (bool): Whether to use scratch directories by default.
        server: Server configuration used for execution.
        scratch (bool): Whether scratch is enabled for this runner.
        scratch_dir (str): Path to scratch directory, if used.
        num_cores (int): Number of CPU cores allocated (from server).
        num_gpus (int): Number of GPUs allocated (from server).
        mem_gb (int): Memory allocation in gigabytes (from server).
    """

    PROGRAM = "NCIPLOT"
    JOBTYPES = [
        "nciplot",
    ]

    FAKE = False
    SCRATCH = True  # default to using scratch for NCIPLOT Jobs

    def __init__(
        self, server, scratch=None, fake=False, scratch_dir=None, **kwargs
    ):
        """
        Initialize the NCIPLOTJobRunner.

        Args:
            server: Server configuration for job execution
            scratch: Whether to use scratch directory (defaults to class SCRATCH)
            fake: Whether to run in fake mode for testing
            scratch_dir: Custom scratch directory path
            **kwargs: Additional keyword arguments
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
        logger.debug(f"Jobrunner scratch: {self.scratch}")
        logger.debug(f"Jobrunner fake mode: {self.fake}")
        logger.debug(f"Jobrunner delete_scratch: {self.delete_scratch}")

    @property
    @lru_cache(maxsize=12)
    def executable(self):
        """
        Get executable class object for NCIPLOT.

        Returns:
            NCIPLOTExecutable: Configured executable for NCIPLOT

        Raises:
            FileNotFoundError: If server configuration is not found
        """
        try:
            logger.info(
                f"Obtaining executable from server: {self.server.name}"
            )
            executable = NCIPLOTExecutable.from_servername(
                servername=self.server.name
            )
            logger.debug(f"Executable obtained: {executable}")
            return executable
        except FileNotFoundError as e:
            logger.error(
                f"No server file {self.server} is found: {e}\n"
                f"Available servers are: {NCIPLOTExecutable.available_servers}"
            )
            raise

    def _prerun(self, job):
        """
        Prepare the job environment before running.

        Args:
            job: NCIPLOTJob instance to prepare
        """
        logger.debug(f"Starting prerun preparation for job: {job.label}")
        self._assign_variables(job)
        self._prepare_files(job)
        logger.debug(f"Completed prerun preparation for job: {job.label}")

    def _assign_variables(self, job):
        """
        Set up file paths for input, output, and error files.

        Args:
            job: NCIPLOTJob instance to configure
        """
        if self.scratch and self.scratch_dir:
            logger.debug(
                f"Setting up job in scratch directory: {self.scratch_dir}"
            )
            self._set_up_variables_in_scratch(job)
        else:
            logger.debug("Setting up job in job directory")
            self._set_up_variables_in_job_directory(job)

    def _set_up_variables_in_scratch(self, job):
        """
        Set up file paths in a scratch directory.

        Args:
            job: NCIPLOTJob instance to configure
        """
        scratch_job_dir = os.path.join(self.scratch_dir, job.label)
        if not os.path.exists(scratch_job_dir):
            with suppress(FileExistsError):
                os.makedirs(scratch_job_dir)
        self.running_directory = scratch_job_dir
        logger.debug(f"Running directory: {self.running_directory}")

        # Set up input file path in scratch
        job_inputfile = job.label + ".nci"
        scratch_job_inputfile = os.path.join(scratch_job_dir, job_inputfile)
        self.job_inputfile = os.path.abspath(scratch_job_inputfile)
        logger.debug(f"Job input file in scratch: {self.job_inputfile}")

        # Set up output file path in scratch
        job_outputfile = job.label + ".nciout"
        scratch_job_outfile = os.path.join(scratch_job_dir, job_outputfile)
        self.job_outputfile = os.path.abspath(scratch_job_outfile)
        logger.debug(f"Job output file in scratch: {self.job_outputfile}")

        # Set up error file path in scratch
        job_errfile = job.label + ".ncierr"
        scratch_job_errfile = os.path.join(scratch_job_dir, job_errfile)
        self.job_errfile = os.path.abspath(scratch_job_errfile)
        logger.debug(f"Job error file in scratch: {self.job_errfile}")

    def _set_up_variables_in_job_directory(self, job):
        """
        Set up file paths in the job's directory.

        Args:
            job: NCIPLOTJob instance to configure
        """
        self.running_directory = job.folder
        logger.debug(f"Running directory: {self.running_directory}")
        self.job_inputfile = os.path.abspath(job.inputfile)
        logger.debug(f"Job input file in folder: {self.job_inputfile}")
        self.job_outputfile = os.path.abspath(job.outputfile)
        logger.debug(f"Job output file in folder: {self.job_outputfile}")
        self.job_errfile = os.path.abspath(job.errfile)
        logger.debug(f"Job error file in folder: {self.job_errfile}")

    def _prepare_files(self, job):
        """
        Prepare input files and write the input for the NCIPLOT job.

        Args:
            job: NCIPLOTJob instance to prepare files for
        """

        # Ensure running directory exists
        if job.molecule is not None:
            logger.debug("Writing XYZ from PubChem molecule")
            self._write_xyz_from_pubchem(job)
        else:
            # Validate filenames are provided
            assert (
                job.filenames is not None
            ), "No molecule provided and no filenames specified for NCIPLOT job."

            if not isinstance(job.filenames, (list, tuple)):
                raise TypeError(
                    f"Expected filenames to be a list or tuple, got "
                    f"{type(job.filenames).__name__}"
                )
            if len(job.filenames) == 0:
                raise ValueError(
                    "No filenames provided for NCIPLOT job. Please provide "
                    "at least one file."
                )
            else:
                # Check if all files are in supported formats
                if all(
                    f.endswith((".wfn", ".wfx", ".xyz")) for f in job.filenames
                ):
                    logger.debug("Copying supported input files")
                    self._copy_input_files(job)
                else:
                    logger.debug("Converting input files to XYZ format")
                    self._write_xyz_from_input_files(job)

    def _copy_input_files(self, job):
        """
        Copy input files to the running directory.

        Args:
            job: NCIPLOTJob instance with files to copy

        Raises:
            FileNotFoundError: If any input file doesn't exist
        """
        for filename in job.filenames:
            if not os.path.exists(filename):
                raise FileNotFoundError(
                    f"File {filename} does not exist for NCIPLOT job."
                )
            with suppress(SameFileError):
                logger.info(
                    f"Copying file {filename} to {self.running_directory}"
                )
                copy(filename, self.running_directory)

    def _write_xyz_from_input_files(self, job):
        """
        Convert input files to XYZ format and copy to running directory.

        Args:
            job: NCIPLOTJob instance with files to convert

        Raises:
            FileNotFoundError: If any input file doesn't exist
            ValueError: If file conversion fails
        """
        for filename in job.filenames:
            if not os.path.exists(filename):
                raise FileNotFoundError(
                    f"File {filename} does not exist for NCIPLOT job."
                )
            try:
                # Convert file types to .xyz
                logger.debug(f"Converting {filename} to XYZ format")
                converter = FileConverter(
                    filename=filename,
                    output_filetype="xyz",
                )
                converter.convert_files()

                with suppress(SameFileError):
                    # Copy the converted .xyz file to running directory
                    xyz_file = filename.rsplit(".", 1)[0] + ".xyz"
                    xyz_promolecular_file = (
                        filename.rsplit(".", 1)[0] + "_promolecular.xyz"
                    )
                    xyz_promolecular_filepath = os.path.join(
                        self.running_directory,
                        os.path.basename(xyz_promolecular_file),
                    )
                    copy(xyz_file, xyz_promolecular_filepath)
                    logger.info(
                        f"Copied file {xyz_file} to {xyz_promolecular_filepath}"
                    )
            except Exception as e:
                raise ValueError(
                    f"Could not convert file {filename} to .xyz format. "
                    f"Error: {e}\n"
                    f"Unsupported file format for NCIPLOT: {filename}.\n"
                    f"Supported formats are .xyz, .wfn, .wfx"
                )

    def _write_xyz_from_pubchem(self, job):
        """
        Write the molecule to an XYZ file if it is provided.

        Args:
            job: NCIPLOTJob instance with molecule to write
        """
        xyz_filepath = os.path.join(self.running_directory, f"{job.label}.xyz")
        job.molecule.write_xyz(filename=xyz_filepath, mode="w")
        logger.info(f"Wrote molecule to {xyz_filepath}")

    def _write_input(self, job):
        """
        Write the input file for NCIPLOT job.

        Args:
            job: NCIPLOTJob instance to write input for
        """
        from chemsmart.jobs.nciplot.writer import NCIPLOTInputWriter

        logger.debug(f"Writing input file for job: {job.label}")
        input_writer = NCIPLOTInputWriter(job=job)
        input_writer.write(target_directory=self.running_directory)
        logger.debug(f"Completed writing input file for job: {job.label}")

    def _get_command(self, job):
        """
        Get execution command for NCIPLOT jobs.

        Args:
            job: NCIPLOTJob instance to get command for

        Returns:
            str: Command string for execution
        """
        exe = self._get_executable()
        # Direct output to stdout without explicit output redirection
        command = f"{exe} {self.job_inputfile}"
        logger.debug(f"Generated command: {command}")
        return command

    def _create_process(self, job, command, env):
        """
        Run the NCIPLOT calculation directly.

        Args:
            job: NCIPLOTJob instance being executed
            command: Command string to execute
            env: Environment variables for execution

        Returns:
            subprocess.Popen: Process object for the running command
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
        Get executable for NCIPLOT.

        Returns:
            str: Path to NCIPLOT executable
        """
        exe = self.executable.get_executable()
        logger.info(f"NCIPLOT executable: {exe}")
        return exe

    def _postrun(self, job):
        """
        Handle post-run tasks, such as copying files from scratch and cleanup.

        Args:
            job: NCIPLOTJob instance that was executed
        """

        if self.scratch:
            # If job was run in scratch, copy files to job folder except temp files
            for file in glob(f"{self.running_directory}/*"):
                if not file.endswith(".tmp"):
                    logger.info(
                        f"Copying file {file} from {self.running_directory} "
                        f"to {job.folder}"
                    )
                    with suppress(SameFileError):
                        # Copy file to job folder
                        copy(file, job.folder)


class FakeNCIPLOTJobRunner(NCIPLOTJobRunner):
    """
    Fake NCIPLOT job runner for testing purposes.

    This class creates a mock job runner that simulates NCIPLOT execution
    without actually running the program, useful for testing and development.

    Attributes:
        FAKE (bool): True for this runner to indicate fake mode.
        SCRATCH (bool): Whether to use scratch directories (inherited default).
        server: Server configuration used for execution (retained/used normally).
        scratch (bool): Whether scratch is enabled for this runner.
        scratch_dir (str): Path to scratch directory, if used.
        num_cores (int): Number of CPU cores allocated (from server).
        num_gpus (int): Number of GPUs allocated (from server).
        mem_gb (int): Memory allocation in gigabytes (from server).
    """

    FAKE = True

    def __init__(
        self, server, scratch=None, fake=True, scratch_dir=None, **kwargs
    ):
        """
        Initialize FakeNCIPLOTJobRunner.

        Args:
            server: Server configuration (ignored in fake mode)
            scratch: Whether to use scratch directory
            fake: Always True for this class
            scratch_dir: Custom scratch directory path
            **kwargs: Additional keyword arguments
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
        Run a fake NCIPLOT job for testing.

        Args:
            job: NCIPLOTJob instance to run in fake mode
            **kwargs: Additional keyword arguments

        Returns:
            int: Return code (always 0 for successful fake runs)
        """
        # Append fake to label to distinguish from real jobs
        job.label = job.label + "_fake"

        self._prerun(job=job)
        self._write_input(job=job)
        returncode = FakeNCIPLOT(self.job_inputfile).run()
        self._postrun(job=job)
        self._postrun_cleanup(job=job)
        return returncode


class FakeNCIPLOT:
    """
    Fake NCIPLOT executable for testing purposes.

    This class simulates NCIPLOT execution by creating mock output files
    without running the actual program.

    Attributes:
        file_to_run (str): Path to the NCIPLOT input file to simulate.
        file_folder (str): Directory path containing the input file.
        filename (str): Basename of the input file.
        input_filepath (str): Full path to the input file.
        output_filepath (str): Full path to the simulated output file (.nciout).
    """

    def __init__(self, file_to_run):
        """
        Initialize FakeNCIPLOT.

        Args:
            file_to_run: Path to input file to simulate processing

        Raises:
            FileNotFoundError: If input file doesn't exist
        """
        if not os.path.exists(file_to_run):
            raise FileNotFoundError(f"File {file_to_run} not found.")
        self.file_to_run = file_to_run

    @property
    def file_folder(self):
        """
        Get the directory containing the input file.

        Returns:
            str: Directory path of the input file
        """
        return os.path.dirname(self.file_to_run)

    @property
    def filename(self):
        """
        Get the basename of the input file.

        Returns:
            str: Filename without directory path
        """
        return os.path.basename(self.file_to_run)

    @property
    def input_filepath(self):
        """
        Get the full path to the input file.

        Returns:
            str: Full path to input file
        """
        return self.file_to_run

    @property
    def output_filepath(self):
        """
        Get the expected output file path.

        Returns:
            str: Full path to output file
        """
        output_file = self.filename.split(".")[0] + ".nciout"
        return os.path.join(self.file_folder, output_file)

    def run(self):
        """
        Run the fake NCIPLOT calculation.

        Creates a mock output file with typical NCIPLOT header and footer
        information to simulate successful program execution.

        Returns:
            int: Return code (always 0 for success)
        """
        with open(self.output_filepath, "w") as g:
            # Write standard NCIPLOT header
            g.write(
                """# ----------------- NCIPLOT ------------------------
 # --- PLOTTING NON COVALENT INTERACTION REGIONS ----
 # ---             E.R. Johnson                  ----
 # ---          J. Contreras-Garcia              ----
 # ----------    Duke University         ------------
 #                                                   
 # ---             A. de la Roza                  ---
 # --------- University of California Merced --------
 #                                                   
 # ---               R. A. Boto                   ---
 # ---                 C. Quan                     --
 # --------  Université Pierre et Marie Curie -------
 # --------------------------------------------------
 # ---              Please cite                  ----
 # --J. Am. Chem. Soc., 2010, 132 (18), pp 6498–6506-
 # --------------------------------------------------
 # --------------------------------------------------
 # ---     Contributions for the wfn properties  ----
 # ---      from H. L. Schmider are acknowledged  ---
 # --------------------------------------------------
 # --------------------------------------------------
 # ---     Contributions for the wfx reader      ----
 # ---      from Dave Arias are acknowledged      ---
 # --------------------------------------------------
 # --------------------------------------------------
 # ---     Contributions for the integration --------
 # ---      algorithms from Erna Wieduwilt    --------
 # ---             are acknowledged           --------
 # ---------------------------------------------------
 #"""
            )
            g.write(f" # Start -- {datetime.now()}\n")
            g.write(
                """-----------------------------------------------------
      INPUT INFORMATION:
-----------------------------------------------------"""
            )
            # Include input file content in output
            for line in open(self.input_filepath, "r"):
                g.write(line)
            g.write("\n")
            # Write completion timestamp
            g.write(f"End -- {datetime.now()}\n")
