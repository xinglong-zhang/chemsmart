import logging
import os
import shlex
import shutil
import subprocess
from pathlib import Path
from shutil import rmtree

from chemsmart.jobs.runner import JobRunner
from chemsmart.utils.periodictable import PeriodicTable
from chemsmart.utils.utils import run_command

pt = PeriodicTable()

logger = logging.getLogger(__name__)


class PyMOLJobRunner(JobRunner):
    # creates job runner process
    # combines information about server and program
    JOBTYPES = []

    PROGRAM = "pymol"

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
        self._assign_variables(job)

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

    def _assign_variables(self, job):
        """Sets proper file paths for job input, output, and error files."""
        self.running_directory = job.folder
        logger.debug(f"Running directory: {self.running_directory}")
        self.job_inputfile = os.path.abspath(job.inputfile)
        self.job_logfile = os.path.abspath(job.logfile)
        self.job_outputfile = os.path.abspath(job.outputfile)
        self.job_errfile = os.path.abspath(job.errfile)
        self.job_pymol_script = job.pymol_script
        self.job_quite_mode = job.quite_mode
        self.job_command_line_only = job.command_line_only

    def _write_input(self, job):
        # write to .xyz file if the supplied file is not .xyz
        if not os.path.exists(job.inputfile):
            mol = job.molecule
            logger.info(f"Writing Molecule to {job.inputfile}.")
            mol.write(job.inputfile, format="xyz")

    def _get_command(self):
        # subclass to implement
        pass

    def _update_os_environ(self, job):
        # no envs to update for pymol
        pass

    def _create_process(self, job, command, env):
        with (
            open(self.job_errfile, "w") as err,
            open(self.job_outputfile, "w") as out,
        ):
            logger.info(
                f"Command executed: {command}\n"
                f"Writing output file to: {self.job_logfile}\n"
                f"And err file to: {self.job_errfile}"
            )
            return subprocess.Popen(
                shlex.split(command),
                stdout=out,
                stderr=err,
                env=env,
                cwd=self.running_directory,
            )

    def _postrun(self, job):
        if job.is_complete():
            # if job is completed, remove scratch directory and submit_script
            # and log.info and log.err files
            if self.scratch:
                logger.info(
                    f"Removing scratch directory: {self.running_directory}."
                )
                rmtree(self.running_directory)

            self._remove_err_files(job)


class PyMOLVisualizationJobRunner(PyMOLJobRunner):
    JOBTYPES = [
        "pymol_visualization",
    ]

    def _get_command(self):
        exe = self.executable
        command = f"{exe} {self.job_inputfile}"
        # get style file
        if self.job_pymol_script is None:
            if os.path.exists("zhang_group_pymol_style.py"):
                job_pymol_script = os.path.abspath(
                    "zhang_group_pymol_style.py"
                )
                if os.path.exists(job_pymol_script):
                    command += f" -r {job_pymol_script}"
            #     pass
            else:
                raise ValueError(
                    "No PyMOL style file can be found!\n Please specify via -s flag."
                )

        if self.job_quite_mode:
            command += " -q"
        if self.job_command_line_only:
            command += " -c"
        # save pse file by default
        command += f' -d "save {self.job_outputfile}; quit"'

        return command
