"""
ORCA job runner implementation.

This module contains the job runner classes for executing ORCA quantum chemistry
calculations on different computing environments, including both real and fake
runners for testing purposes.
"""

import logging
import os
import re
import shlex
import subprocess
from contextlib import suppress
from functools import lru_cache
from glob import glob
from shutil import copy

from chemsmart.io.orca.input import ORCAInput
from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.executable import ORCAExecutable
from chemsmart.utils.periodictable import PeriodicTable

pt = PeriodicTable()

logger = logging.getLogger(__name__)


class ORCAJobRunner(JobRunner):
    """
    ORCA-specific job runner.

    This class handles the execution of ORCA quantum chemistry calculations
    with support for scratch directories, file management, and various job types.

    Attributes:
        JOBTYPES (list): Supported job types handled by this runner.
        PROGRAM (str): Program identifier ('orca').
        FAKE (bool): Whether this runner operates in fake/test mode.
        SCRATCH (bool): Whether to use scratch directories by default.
        server: Server configuration used for execution.
        scratch (bool): Whether scratch is enabled for this runner.
        scratch_dir (str): Path to scratch directory, if used.
        num_cores (int): Number of CPU cores allocated (from server).
        num_gpus (int): Number of GPUs allocated (from server).
        mem_gb (int): Memory allocation in gigabytes (from server).
    """

    # creates job runner process
    # combines information about server and program

    JOBTYPES = [
        "orcajob",
        "orcainp",
        "orcaopt",
        "orcamodred",
        "orcaqrc",
        "orcascan",
        "orcats",
        "orcasp",
        "orcairc",
        "orcaqmmm",
        "orcaneb",
    ]

    PROGRAM = "orca"

    FAKE = False
    SCRATCH = True

    def __init__(
        self, server, scratch=None, fake=False, scratch_dir=None, **kwargs
    ):
        """
        Initialize the ORCA job runner.

        Args:
            server: Server configuration object
            scratch: Whether to use scratch directory (default: True)
            fake: Whether to use fake runner for testing
            scratch_dir: Path to scratch directory
            **kwargs: Additional keyword arguments
        """
        # Use default SCRATCH if scratch is not explicitly set
        if scratch is None:
            scratch = self.SCRATCH  # default to True for ORCA jobs
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
        logger.debug(f"Jobrunner delete_scratch: {self.delete_scratch}")

    @property
    @lru_cache(maxsize=12)
    def executable(self):
        """
        Get the executable class object for ORCA.

        Returns:
            ORCAExecutable: ORCA executable configuration object

        Raises:
            FileNotFoundError: If server configuration file is not found
        """
        try:
            logger.info(
                f"Obtaining executable from server: {self.server.name}"
            )
            executable = ORCAExecutable.from_servername(
                servername=self.server.name
            )
            return executable
        except FileNotFoundError as e:
            logger.error(
                f"No server file {self.server} is found: {e}\n"
                f"Available servers are: {ORCAExecutable.available_servers}"
            )
            raise

    def _prerun(self, job):
        """
        Perform pre-run setup for the job.

        Args:
            job: The job object to prepare for execution
        """
        self._assign_variables(job)

    def _assign_variables(self, job):
        """
        Set proper file paths for job input, output, and error files.

        Configures paths for execution in scratch directory or job directory.

        Args:
            job: The job object to configure paths for
        """
        if self.scratch and self.scratch_dir:
            self._set_up_variables_in_scratch(job)
        else:
            self._set_up_variables_in_job_directory(job)

        if self.executable and self.executable.local_run is not None:
            logger.info(f"Local run is {self.executable.local_run}.")
            job.local = self.executable.local_run

    def _set_up_variables_in_scratch(self, job):
        """
        Set up file paths for execution in scratch directory.

        Args:
            job: The job object to configure for scratch execution
        """
        scratch_job_dir = os.path.join(self.scratch_dir, job.label)
        if not os.path.exists(scratch_job_dir):
            with suppress(FileExistsError):
                os.makedirs(scratch_job_dir)
        self.running_directory = scratch_job_dir
        logger.debug(f"Running directory: {self.running_directory}")

        job_inputfile = job.label + ".inp"
        scratch_job_inputfile = os.path.join(scratch_job_dir, job_inputfile)
        self.job_inputfile = os.path.abspath(scratch_job_inputfile)

        job_gbwfile = job.label + ".gbw"
        scratch_job_gbwfile = os.path.join(scratch_job_dir, job_gbwfile)
        self.job_gbwfile = os.path.abspath(scratch_job_gbwfile)

        job_errfile = job.label + ".err"
        scratch_job_errfile = os.path.join(scratch_job_dir, job_errfile)
        self.job_errfile = os.path.abspath(scratch_job_errfile)

        job_outputfile = job.label + ".out"
        scratch_job_outputfile = os.path.join(scratch_job_dir, job_outputfile)
        self.job_outputfile = os.path.abspath(scratch_job_outputfile)

    def _set_up_variables_in_job_directory(self, job):
        """
        Set up file paths for execution in job directory.

        Args:
            job: The job object to configure for job directory execution
        """
        self.running_directory = job.folder
        logger.debug(f"Running directory: {self.running_directory}")
        self.job_inputfile = os.path.abspath(job.inputfile)
        self.job_gbwfile = os.path.abspath(job.gbwfile)
        self.job_errfile = os.path.abspath(job.errfile)
        self.job_outputfile = os.path.abspath(job.outputfile)

    def _copy_over_xyz_files(self, job):
        """
        Copy xyz files from job directory to scratch directory.

        This method searches for xyz file references in the input file
        and copies them to the scratch directory if running in scratch mode.
        Particularly important for NEB calculations which require multiple
        geometry files (ending_xyzfile, intermediate_xyzfile, restarting_xyzfile).

        Handles both relative and absolute file paths by:
        1. Resolving relative paths against the job folder
        2. Falling back to current working directory if not found
        3. Copying files to scratch using only basenames

        Note: Must be called AFTER _write_input() so the input file exists
        and can be parsed for XYZ file references.

        Args:
            job: The job object containing input file information

        Raises:
            FileNotFoundError: If referenced xyz file does not exist at resolved path
        """
        from chemsmart.utils.repattern import xyz_filename_pattern

        # Read from the scratch input file location
        input_file_to_read = (
            self.job_inputfile if self.scratch else job.inputfile
        )

        with open(input_file_to_read, "r") as f:
            for line in f:
                match = re.search(xyz_filename_pattern, line)
                if match:
                    xyz_file = match.group(1)

                    # Handle absolute and relative paths
                    if not os.path.isabs(xyz_file):
                        # Try relative to job folder first
                        xyz_file_path = os.path.join(job.folder, xyz_file)
                        if not os.path.exists(xyz_file_path):
                            # Try relative to current working directory
                            xyz_file_path = os.path.abspath(xyz_file)
                    else:
                        xyz_file_path = xyz_file

                    if not os.path.exists(xyz_file_path):
                        raise FileNotFoundError(
                            f"XYZ file {xyz_file} does not exist at {xyz_file_path}."
                        )

                    # copy to scratch if running in scratch
                    if self.scratch and self.scratch_dir:
                        # Use basename for scratch location
                        xyz_basename = os.path.basename(xyz_file_path)
                        xyz_file_scratch = os.path.join(
                            self.running_directory, xyz_basename
                        )
                        copy(xyz_file_path, xyz_file_scratch)
                        logger.info(
                            f"Copied {xyz_file_path} to {xyz_file_scratch}."
                        )

    def _write_input(self, job):
        """
        Write the input file for the job.

        Creates the ORCA input file in the running directory (scratch or job folder).
        After writing the input, automatically copies any referenced XYZ files to
        scratch directory if scratch is enabled. This is essential for NEB jobs
        that reference multiple geometry files.

        Args:
            job: The job object to write input for

        Note:
            XYZ file copying happens after input writing so the input file can
            be parsed to discover file references.
        """
        from chemsmart.jobs.orca.writer import ORCAInputWriter

        input_writer = ORCAInputWriter(job=job)
        input_writer.write(target_directory=self.running_directory)

        # Copy XYZ files to scratch after writing input file
        # Use the actual input file path (self.job_inputfile in scratch)
        if self.scratch and os.path.exists(self.job_inputfile):
            self._copy_over_xyz_files(job)

    def _get_command(self, job):
        """
        Get the command string to execute the ORCA job.

        Args:
            job: The job object to get command for

        Returns:
            str: Command string for job execution
        """
        exe = self._get_executable()
        command = f"{exe} {self.job_inputfile}"
        return command

    def _create_process(self, job, command, env):
        """
        Create subprocess for job execution.

        Args:
            job: The job object to execute
            command: Command string to execute
            env: Environment variables for execution

        Returns:
            subprocess.Popen: Process object for the running job
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
        Get executable path for ORCA.

        Returns:
            str: Path to ORCA executable
        """
        exe = self.executable.get_executable()
        logger.info(f"ORCA executable: {exe}")
        return exe

    def _postrun(self, job):
        """
        Perform post-run cleanup and file management.

        Args:
            job: The completed job object
        """
        logger.debug(f"Scratch: {self.scratch}")

        if self.scratch:
            logger.debug(f"Running directory: {self.running_directory}")
            # if job was run in scratch, copy files to job folder except files containing .tmp
            for file in glob(f"{self.running_directory}/{job.label}*"):
                if not file.endswith((".tmp", ".tmp.*")):
                    logger.info(
                        f"Copying file {file} from {self.running_directory} "
                        f"to {job.folder}"
                    )
                    try:
                        copy(file, job.folder)
                    except Exception as e:
                        logger.error(
                            f"Failed to copy file {file} to {job.folder}: {e}"
                        )


class FakeORCAJobRunner(ORCAJobRunner):
    """
    Fake ORCA job runner for testing purposes.

    This class simulates ORCA job execution without actually running
    calculations, useful for testing and development.

    Attributes:
        PROGRAM (str): Program identifier ('orca').
        JOBTYPES (list): Supported job types handled (inherits from ORCAJobRunner).
        FAKE (bool): True for this runner to indicate fake mode.
        SCRATCH (bool): Whether to use scratch directories (inherits default).
        server: Server configuration used for execution.
        scratch (bool): Whether scratch is enabled for this runner.
        scratch_dir (str): Path to scratch directory, if used.
        num_cores (int): Number of CPU cores allocated (from server).
        num_gpus (int): Number of GPUs allocated (from server).
        mem_gb (int): Memory allocation in gigabytes (from server).
    """

    # creates job runner process
    # combines information about server and program
    FAKE = True

    def __init__(self, server, scratch=None, fake=True, **kwargs):
        """
        Initialize the fake ORCA job runner.

        Args:
            server: Server configuration object
            scratch: Whether to use scratch directory
            fake: Always True for fake runner
            **kwargs: Additional keyword arguments
        """
        super().__init__(server=server, scratch=scratch, fake=fake, **kwargs)

    def run(self, job):
        """
        Run a fake ORCA job.

        Args:
            job: The job object to run

        Returns:
            int: Return code from fake execution
        """
        self._prerun(job=job)
        self._write_input(job=job)
        returncode = FakeORCA(self.job_inputfile).run()
        self._postrun(job=job)
        self._postrun_cleanup(job=job)
        return returncode


class FakeORCA:
    """
    Fake ORCA execution simulator.

    This class simulates ORCA program execution by generating
    fake output files for testing purposes.

    Attributes:
        file_to_run (str): Absolute path to the ORCA input file to simulate.
        input_object (ORCAInput): Parsed representation of the input file.

    Properties:
        file_folder (str): Directory containing the input file.
        filename (str): Basename of the input file.
        input_filepath (str): Absolute path to the input file.
        output_filepath (str): Absolute path to the simulated output (.out).
        input_contents (list[str]): Lines of the input file.
        molecule (Molecule): Parsed molecular structure from the input.
        charge (int): Molecular charge.
        multiplicity (int): Spin multiplicity.
        spin (str): 'R' for restricted (singlet), 'U' otherwise.
        num_atoms (int): Number of atoms in the molecule.
        atomic_symbols (list[str]): Chemical symbols of atoms.
        atomic_numbers (list[int]): Atomic numbers of atoms.
        atomic_coordinates (array): Cartesian coordinates of atoms.
        empirical_formula (str): Empirical chemical formula string.
    """

    def __init__(self, file_to_run):
        """
        Initialize the fake ORCA simulator.

        Args:
            file_to_run: Path to the input file to simulate

        Raises:
            FileNotFoundError: If input file does not exist
        """
        if not os.path.exists(file_to_run):
            raise FileNotFoundError(f"File {file_to_run} not found.")
        self.file_to_run = os.path.abspath(file_to_run)
        self.input_object = ORCAInput(filename=self.file_to_run)

    @property
    def file_folder(self):
        """
        Get the folder containing the input file.

        Returns:
            str: Directory path of the input file
        """
        return os.path.dirname(self.file_to_run)

    @property
    def filename(self):
        """
        Get the filename of the input file.

        Returns:
            str: Basename of the input file
        """
        return os.path.basename(self.file_to_run)

    @property
    def input_filepath(self):
        """
        Get the full path to the input file.

        Returns:
            str: Absolute path to the input file
        """
        return self.file_to_run

    @property
    def output_filepath(self):
        """
        Get the path to the output file.

        Returns:
            str: Path to the output file (.out extension)
        """
        output_file = self.filename.split(".")[0] + ".out"
        return os.path.join(self.file_folder, output_file)

    @property
    def input_contents(self):
        """
        Get the contents of the input file.

        Returns:
            str: Contents of the input file
        """
        return self.input_object.contents

    @property
    def molecule(self):
        """
        Get the molecule object from the input.

        Returns:
            Molecule: Molecule object from the input file
        """
        return self.input_object.molecule

    @property
    def charge(self):
        """
        Get the charge from the input.

        Returns:
            int: Molecular charge
        """
        return self.input_object.charge

    @property
    def multiplicity(self):
        """
        Get the multiplicity from the input.

        Returns:
            int: Spin multiplicity
        """
        return self.input_object.multiplicity

    @property
    def spin(self):
        """
        Get the spin type based on multiplicity.

        Returns:
            str: 'R' for restricted (singlet) or 'U' for unrestricted
        """
        if self.multiplicity == 1:
            return "R"
        return "U"

    @property
    def num_atoms(self):
        """
        Get the number of atoms in the molecule.

        Returns:
            int: Number of atoms
        """
        return len(self.molecule.chemical_symbols)

    @property
    def atomic_symbols(self):
        """
        Get the atomic symbols of the molecule.

        Returns:
            list: List of atomic symbols
        """
        return list(self.molecule.chemical_symbols)

    @property
    def atomic_numbers(self):
        """
        Get the atomic numbers of the molecule.

        Returns:
            list: List of atomic numbers
        """
        return [pt.to_atomic_number(s) for s in self.atomic_symbols]

    @property
    def atomic_coordinates(self):
        """
        Get the atomic coordinates of the molecule.

        Returns:
            array: Atomic coordinates
        """
        return self.molecule.positions

    @property
    def empirical_formula(self):
        """
        Get the empirical formula of the molecule.

        Returns:
            str: Empirical formula
        """
        return self.molecule.get_chemical_formula(empirical=True)

    def run(self):
        """
        Run the fake ORCA calculation.

        Generates a fake output file with standard ORCA output format
        for testing purposes.

        Returns:
            int: Return code (always 0 for successful fake run)
        """

        # get input lines
        with open(self.input_filepath) as f:
            input_lines = f.readlines()

        with open(self.output_filepath, "w") as g:
            g.write("\n")
            g.write("""                                 *****************
                                 * O   R   C   A *
                                 *****************\n""")
            g.write(
                "                         Fake Orca version 0.0.0 -  RELEASE  -                \n"
            )
            g.write(
                "================================================================================\n"
            )
            g.write(
                "                                       INPUT FILE                               \n"
            )
            g.write(
                "================================================================================\n"
            )
            g.write(f"NAME = {self.input_filepath}\n")
            for i, line in enumerate(input_lines):
                g.write(f"|  {i+1}> {line}")
            g.write(f"|  {len(input_lines)+1}> ")
            g.write(
                f"|  {len(input_lines) + 2}>                          ****END OF INPUT****\n"
            )
            g.write("                       *****************************\n")
            g.write("                       *        Fake Orca Run      *\n")
            g.write("                       *****************************\n")
            g.write(
                f"Number of atoms                         .... {self.num_atoms}\n"
            )
            g.write(
                f"  Total Charge           Charge          ....    {self.input_object.charge}\n"
            )
            g.write(
                f"  Multiplicity           Mult            ....    {self.input_object.multiplicity}\n"
            )
            g.write(
                "                     ***********************HURRAY********************\n"
            )
            g.write(
                "                     ***        THE OPTIMIZATION HAS CONVERGED     ***\n"
            )
            g.write(
                "                     ***********************HURRAY********************\n"
            )
            g.write(
                "                  *******************************************************\n"
            )
            g.write(
                "                  *** FINAL ENERGY EVALUATION AT THE STATIONARY POINT ***\n"
            )
            g.write(
                "                  *******************************************************\n"
            )
            g.write(" ---------------------------------\n")
            g.write(" CARTESIAN COORDINATES (ANGSTROEM)\n")
            g.write(" ---------------------------------\n")
            for line in input_lines[2:-2]:
                g.write(f"{line}")
            g.write(" ----------------\n")
            g.write(" TOTAL SCF ENERGY\n")
            g.write(" ----------------\n")
            g.write(
                " Total Energy       :          -76.32331101 Eh           -2076.86288 eV\n"
            )
            g.write(
                "                              ****ORCA TERMINATED NORMALLY****\n"
            )
            g.write(
                "TOTAL RUN TIME: 0 days 0 hours 0 minutes 0 seconds xx msec\n"
            )
