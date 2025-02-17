import logging
import os
import shlex
import subprocess
from contextlib import suppress
from functools import lru_cache
from glob import glob
from shutil import copy, rmtree

from chemsmart.io.gaussian.input import Gaussian16Input
from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.executable import XTBExecutable
from chemsmart.utils.periodictable import PeriodicTable

pt = PeriodicTable()

logger = logging.getLogger(__name__)


class XTBJobRunner(JobRunner):
    # creates job runner process
    # combines information about server and program
    JOBTYPES = [
        "xtbopt",
    ]

    PROGRAM = "XTB"

    FAKE = False

    def __init__(self, server, scratch=None, fake=False, **kwargs):
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
            executable = XTBExecutable.from_servername(
                servername=self.server.name
            )
            return executable
        except FileNotFoundError as e:
            logger.error(
                f"No server file {self.server} is found: {e}\n"
                f"Available servers are: {XTBExecutable.available_servers}"
            )
            raise

    def _prerun(self, job):
        self._assign_variables(job)

    def _assign_variables(self, job):
        """Sets proper file paths for job input, output,
        and error files in scratch or not in scratch."""
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

        job_inputfile = job.label + ".xyz"
        scratch_job_inputfile = os.path.join(scratch_job_dir, job_inputfile)
        self.job_inputfile = os.path.abspath(scratch_job_inputfile)

        job_errfile = job.label + ".err"
        scratch_job_errfile = os.path.join(scratch_job_dir, job_errfile)
        self.job_errfile = os.path.abspath(scratch_job_errfile)

    def _set_up_variables_in_job_directory(self, job):
        self.running_directory = job.folder
        logger.debug(f"Running directory: {self.running_directory}")
        self.job_inputfile = os.path.abspath(job.inputfile)
        self.job_errfile = os.path.abspath(job.errfile)

    def _write_input(self, job):
        pass

    def _get_command(self, settings):
        exe = self._get_executable()
        if (
            settings.solvent_model is not None
            and settings.solvent_id is not None
        ):
            command = (
                f"{exe} {self.job_inputfile} --{settings.xtb_version} "
                f"--{settings.jobtype} {settings.optimization_level}  --chrg {settings.charge} "
                f"--uhf {settings.unpair_electrons} --{settings.solvent_model} {settings.solvent_id}"
            )
        else:
            command = (
                f"{exe} {self.job_inputfile} --{settings.xtb_version} "
                f"--{settings.jobtype} {settings.optimization_level}  --chrg {settings.charge} "
                f"--uhf {settings.unpair_electrons}"
            )
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


class FakeXTBJobRunner(XTBJobRunner):
    # creates job runner process
    # combines information about server and program
    FAKE = True

    def __init__(self, server, scratch=None, fake=True, **kwargs):
        super().__init__(server=server, scratch=scratch, fake=fake, **kwargs)

    def run(self, job, **kwargs):
        self._prerun(job=job)
        self._write_input(job=job)
        returncode = FakeXTB(self.job_inputfile).run()
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


class FakeXTB:
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
        with open(self.output_filepath, "w") as g:
            g.write("-" * 100 + "\n")
            g.write("* xtb version 0.0.0 (Fake) \n")
            g.write("-" * 100 + "\n")
            g.write("\n")

            for line in self.input_contents:
                g.write(f" {line}\n")

            g.write("\n")
            g.write("* finished run (fake xtb) \n")
