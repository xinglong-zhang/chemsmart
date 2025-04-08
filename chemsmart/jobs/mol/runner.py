import logging
import os
import shlex
import shutil
import subprocess
from glob import glob
from pathlib import Path
from shutil import copy, rmtree

from chemsmart.jobs.runner import JobRunner
from chemsmart.utils.periodictable import PeriodicTable
from chemsmart.utils.utils import run_command

pt = PeriodicTable()

logger = logging.getLogger(__name__)


class PyMOLJobRunner(JobRunner):
    # creates job runner process
    # combines information about server and program
    JOBTYPES = [
        "pymol_visualization",
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
    SCRATCH = False
    # default to not use scratch for pymol jobs

    def __init__(self, server, scratch=None, fake=False, **kwargs):
        # Use default SCRATCH if scratch is not explicitly set
        if scratch is None:
            scratch = self.SCRATCH
        super().__init__(server=server, scratch=scratch, fake=fake, **kwargs)
        logger.debug(f"Jobrunner server: {self.server}")
        logger.debug(f"Jobrunner num cores: {self.num_cores}")
        logger.debug(f"Jobrunner num hours: {self.num_hours}")
        logger.debug(f"Jobrunner num gpus: {self.num_gpus}")
        logger.debug(f"Jobrunner mem gb: {self.mem_gb}")
        logger.debug(f"Jobrunner num threads: {self.num_threads}")
        logger.debug(f"Jobrunner scratch: {self.scratch}")

    @property
    def executable(self):
        """Define the path for the PyMOL executable."""
        pymol_path = run_command(["which", "pymol"])
        if not os.path.exists(pymol_path):
            raise FileNotFoundError(
                f"PyMOL not found: {pymol_path}. Please install PyMOL!"
            )
        return pymol_path

    @property
    def pymol_templates_path(self):
        """Define the path for PyMOL templates directory."""
        return Path(__file__).resolve().parent / "templates"

    def _prerun(self, job):
        self._generate_visualization_style_script()

    def _generate_visualization_style_script(self):
        # Define the source and destination file paths
        source_style_file = (
            self.pymol_templates_path / "zhang_group_pymol_style.py"
        )
        cwd = os.getcwd()
        dest_style_file = os.path.join(cwd, "zhang_group_pymol_style.py")

        # Check if the style file already exists in the current working directory
        if not os.path.exists(dest_style_file):
            # Copy the file from templates to the current working directory
            logger.debug(
                f"Copying file from {source_style_file} to {dest_style_file}."
            )
            shutil.copy(source_style_file, dest_style_file)

    def _write_input(self, job):
        # write to .xyz file if the supplied file is not .xyz
        if not os.path.exists(job.inputfile):
            mol = job.molecule
            logger.info(f"Writing Molecule to {job.inputfile}.")
            mol.write(job.inputfile, format="xyz")

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

            self._remove_err_files(job)
