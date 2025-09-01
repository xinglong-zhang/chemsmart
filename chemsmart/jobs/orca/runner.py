import logging
import os
import re
import shlex
import subprocess
from contextlib import suppress
from functools import lru_cache
from glob import glob
from shutil import copy, rmtree

from chemsmart.io.orca.input import ORCAInput
from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.executable import ORCAExecutable
from chemsmart.utils.periodictable import PeriodicTable

pt = PeriodicTable()

logger = logging.getLogger(__name__)


class ORCAJobRunner(JobRunner):
    # creates job runner process
    # combines information about server and program

    JOBTYPES = [
        "orcajob",
        "orcainp",
        "orcaopt",
        "orcamodred",
        "orcascan",
        "orcats",
        "orcasp",
        "orcairc",
    ]

    PROGRAM = "orca"

    FAKE = False
    SCRATCH = True

    def __init__(
        self, server, scratch=None, fake=False, scratch_dir=None, **kwargs
    ):
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

    @property
    @lru_cache(maxsize=12)
    def executable(self):
        """Executable class object for ORCA."""
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
        self._assign_variables(job)
        if self.scratch and os.path.exists(job.inputfile):
            # copy input file to scratch directory
            self._copy_over_xyz_files(job)

    def _assign_variables(self, job):
        """Sets proper file paths for job input, output, and error files in scratch or not in scratch."""
        if self.scratch and self.scratch_dir:
            self._set_up_variables_in_scratch(job)
        else:
            self._set_up_variables_in_job_directory(job)

        if self.executable and self.executable.local_run is not None:
            logger.info(f"Local run is {self.executable.local_run}.")
            job.local = self.executable.local_run

    def _set_up_variables_in_scratch(self, job):
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
        self.running_directory = job.folder
        logger.debug(f"Running directory: {self.running_directory}")
        self.job_inputfile = os.path.abspath(job.inputfile)
        self.job_gbwfile = os.path.abspath(job.gbwfile)
        self.job_errfile = os.path.abspath(job.errfile)
        self.job_outputfile = os.path.abspath(job.outputfile)

    def _copy_over_xyz_files(self, job):
        """Copy over xyz files from run directory to scratch directory,
        if running in scratch."""
        from chemsmart.utils.repattern import xyz_filename_pattern

        with open(job.inputfile, "r") as f:
            for line in f:
                match = re.search(xyz_filename_pattern, line)
                if match:
                    xyz_file = match.group(1)
                    if not os.path.exists(xyz_file):
                        raise FileNotFoundError(
                            f"XYZ file {xyz_file} does not exist."
                        )

                    # copy to scratch if running in scratch
                    if self.scratch and self.scratch_dir:
                        xyz_file_scratch = os.path.join(
                            self.running_directory, xyz_file
                        )
                        copy(xyz_file, xyz_file_scratch)
                        logger.info(
                            f"Copied {xyz_file} to {self.running_directory}."
                        )

    def _write_input(self, job):
        from chemsmart.jobs.orca.writer import ORCAInputWriter

        input_writer = ORCAInputWriter(job=job)
        input_writer.write(target_directory=self.running_directory)

    def _get_command(self, job):
        exe = self._get_executable()
        command = f"{exe} {self.job_inputfile}"
        return command

    def _create_process(self, job, command, env):
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
        """Get executable for Gaussian."""
        exe = self.executable.get_executable()
        logger.info(f"ORCA executable: {exe}")
        return exe

    def _postrun(self, job):
        logger.debug(f"Scratch: {self.scratch}")
        if self.scratch:
            # if job was run in scratch, copy files to job folder except files containing .tmp
            for file in glob(f"{self.running_directory}/{job.label}*"):
                if not file.endswith((".tmp", ".tmp.*")):
                    if not os.path.exists(file):
                        logger.info(
                            f"Copying file {file} from {self.running_directory} to {job.folder}"
                        )
                        try:
                            copy(file, job.folder)
                        except Exception as e:
                            logger.error(
                                f"Failed to copy file {file} to {job.folder}: {e}"
                            )

        if job.is_complete():
            # if job is completed, remove scratch directory
            # and err files
            if self.scratch:
                logger.info(
                    f"Removing scratch directory: {self.running_directory}."
                )
                rmtree(self.running_directory)

            self._remove_err_files(job)


class FakeORCAJobRunner(ORCAJobRunner):
    # creates job runner process
    # combines information about server and program
    FAKE = True

    def __init__(self, server, scratch=None, fake=True, **kwargs):
        super().__init__(server=server, scratch=scratch, fake=fake, **kwargs)

    def run(self, job):
        self._prerun(job=job)
        self._write_input(job=job)
        returncode = FakeORCA(self.job_inputfile).run()
        self._postrun(job=job)
        return returncode


class FakeORCA:
    def __init__(self, file_to_run):
        if not os.path.exists(file_to_run):
            raise FileNotFoundError(f"File {file_to_run} not found.")
        self.file_to_run = os.path.abspath(file_to_run)
        self.input_object = ORCAInput(filename=self.file_to_run)

    @property
    def file_folder(self):
        return os.path.dirname(self.file_to_run)

    @property
    def filename(self):
        return os.path.basename(self.file_to_run)

    @property
    def input_filepath(self):
        return self.file_to_run

    @property
    def output_filepath(self):
        output_file = self.filename.split(".")[0] + ".out"
        return os.path.join(self.file_folder, output_file)

    @property
    def input_contents(self):
        return self.input_object.contents

    @property
    def molecule(self):
        return self.input_object.molecule

    @property
    def charge(self):
        return self.input_object.charge

    @property
    def multiplicity(self):
        return self.input_object.multiplicity

    @property
    def spin(self):
        if self.multiplicity == 1:
            return "R"
        return "U"

    @property
    def num_atoms(self):
        return len(self.molecule.chemical_symbols)

    @property
    def atomic_symbols(self):
        return list(self.molecule.chemical_symbols)

    @property
    def atomic_numbers(self):
        return [pt.to_atomic_number(s) for s in self.atomic_symbols]

    @property
    def atomic_coordinates(self):
        return self.molecule.positions

    @property
    def empirical_formula(self):
        return self.molecule.get_chemical_formula(empirical=True)

    def run(self):
        # get input lines
        with open(self.input_filepath) as f:
            input_lines = f.readlines()

        with open(self.output_filepath, "w") as g:
            g.write("\n")
            g.write(
                """                                 *****************
                                 * O   R   C   A *
                                 *****************\n"""
            )
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
