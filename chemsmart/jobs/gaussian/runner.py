import logging
import os
import sys
import shlex
import subprocess
from contextlib import suppress
from functools import lru_cache
from datetime import datetime
from glob import glob
from random import random
from shutil import copy, rmtree

from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.executable import GaussianExecutable
from chemsmart.io.gaussian.input import Gaussian16Input
from chemsmart.utils.periodictable import PeriodicTable

pt = PeriodicTable()

if sys.version_info >= (3, 10):
    from shutil import _USE_CP_SENDFILE  # noqa F811

    _USE_CP_SENDFILE = False  # noqa F811
    # to avoid "BlockingIOError: [Errno 11] Resource temporarily unavailable:" Error when copying
    # only works in Python 3.10

logger = logging.getLogger(__name__)


class GaussianJobRunner(JobRunner):
    # creates job runner process
    # combines information about server and program
    JOBTYPES = [
        "g16crest",
        "g16crestopt",
        "g16crestts",
        "g16job",
        "g16dias",
        "g16opt",
        "g16irc",
        "g16modred",
        "g16nci",
        "g16resp",
        "g16saopt",
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

    def __init__(self, server, scratch=True, fake=False, **kwargs):
        super().__init__(server=server, scratch=scratch, fake=fake, **kwargs)
        logger.debug(f"Jobrunner num cores: {self.num_cores}")
        logger.debug(f"Jobrunner num hours: {self.num_hours}")
        logger.debug(f"Jobrunner num gpus: {self.num_gpus}")
        logger.debug(f"Jobrunner mem gb: {self.mem_gb}")
        logger.debug(f"Jobrunner num threads: {self.num_threads}")

    @property
    @lru_cache(maxsize=12)
    def executable(self):
        """Executable class object for Gaussian."""
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
        self._assign_variables(job)

    def _assign_variables(self, job):
        """Sets proper file paths for job input, output, and error files in scratch or not in scratch."""
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
        scratch_job_dir = os.path.join(self.scratch_dir, job.label)
        if not os.path.exists(scratch_job_dir):
            with suppress(FileExistsError):
                os.makedirs(scratch_job_dir)
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

    def _set_up_variables_in_job_directory(self, job):
        self.running_directory = job.folder
        logger.debug(f"Running directory: {self.running_directory}")
        self.job_inputfile = os.path.abspath(job.inputfile)
        self.job_chkfile = os.path.abspath(job.chkfile)
        self.job_errfile = os.path.abspath(job.errfile)

    def _write_input(self, job):
        from chemsmart.jobs.gaussian.writer import GaussianInputWriter

        input_writer = GaussianInputWriter(job=job, jobrunner=self)
        input_writer.write(target_directory=self.running_directory)

    def _get_command(self):
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
        logger.info(f"Gaussian executable: {exe}")
        return exe

    def _postrun(self, job):
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

            # writer = SubmitscriptWriter(job)
            # submit_script = writer.job_submit_script
            # run_script = writer.job_run_script
            # err_filepath = os.path.join(job.folder, f"{job.errfile}")
            # joblogerr_filepath = os.path.join(job.folder, "log.err")
            # jobloginfo_filepath = os.path.join(job.folder, "log.info")
            # pbs_errfile = os.path.join(job.folder, "pbs.err")
            # pbs_infofile = os.path.join(job.folder, "pbs.info")
            #
            # files_to_remove = [
            #     submit_script,
            #     run_script,
            #     err_filepath,
            #     joblogerr_filepath,
            #     jobloginfo_filepath,
            #     pbs_errfile,
            #     pbs_infofile,
            # ]
            #
            # for f in files_to_remove:
            #     with suppress(FileNotFoundError):
            #         os.remove(f)


class FakeGaussianJobRunner(GaussianJobRunner):
    # creates job runner process
    # combines information about server and program
    FAKE = True

    def __init__(self, server, scratch=True, fake=True, **kwargs):
        super().__init__(server=server, scratch=scratch, fake=fake, **kwargs)

    def run(self, job, **kwargs):
        self._prerun(job=job)
        self._write_input(job=job)
        returncode = FakeGaussian(self.job_inputfile).run()
        self._postrun(job=job)
        return returncode

    def _set_up_variables_in_scratch(self, job):
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
        self.running_directory = job.folder
        logger.debug(f"Running directory: {self.running_directory}")
        job.label = f"{job.label}_fake"
        logger.debug(f"Job label for fake job run: {job.label}")
        self.job_inputfile = os.path.abspath(job.inputfile)
        self.job_chkfile = os.path.abspath(job.chkfile)
        self.job_errfile = os.path.abspath(job.errfile)


class FakeGaussian:
    def __init__(self, file_to_run):
        if not os.path.exists(file_to_run):
            raise FileNotFoundError(f"File {file_to_run} not found.")
        self.file_to_run = file_to_run
        self.input_object = Gaussian16Input(filename=self.file_to_run)

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
        output_file = self.filename.split(".")[0] + ".log"
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
            lines = f.readlines()

        with open(self.output_filepath, "w") as g:
            g.write(" Entering Gaussian System, FakeGaussianRunner\n")
            g.write(f" Input={self.input_filepath}\n")
            g.write(f" Output={self.output_filepath}\n")
            g.write(" ******************************************\n")
            g.write(" Fake Gaussian Executable\n")
            g.write(" ******************************************\n")
            # write mem/nproc/chk information (%...)
            for line in lines:
                if line.startswith("%"):
                    g.write(f" {line}")
            # write route information
            for line in lines:
                if line.startswith("#"):
                    line_len = len(line)
                    g.write(" " + "-" * line_len + "\n")
                    g.write(f" {line}")
                    g.write(" " + "-" * line_len + "\n")
            # write charge and multiplicity
            g.write(
                f" Charge =  {self.charge} Multiplicity = {self.multiplicity}\n"
            )
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
