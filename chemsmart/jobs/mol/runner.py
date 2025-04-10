import logging
import os
import shlex
import shutil
import subprocess
import sys  # Add this import for sys.platform
from pathlib import Path
from shutil import rmtree

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.runner import JobRunner
from chemsmart.utils.periodictable import PeriodicTable

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
        """Define the path for the PyMOL executable, handling Windows-specific naming."""
        # Use pymol.exe on Windows, pymol on Unix-like systems
        pymol_cmd = "pymol.exe" if sys.platform == "win32" else "pymol"
        pymol_path = shutil.which(pymol_cmd)
        if pymol_path is None or not os.path.exists(pymol_path):
            raise FileNotFoundError(
                f"PyMOL executable '{pymol_cmd}' not found in PATH. Please install PyMOL!"
            )
        return pymol_path

    @property
    def pymol_templates_path(self):
        """Define the path for PyMOL templates directory."""
        return Path(__file__).resolve().parent / "templates"

    def _prerun(self, job):
        self._assign_variables(job)

    def _assign_variables(self, job):
        """Sets proper file paths for job input, output, and error files."""
        self.running_directory = job.folder
        logger.debug(f"Running directory: {self.running_directory}")
        self.job_basename = job.label
        self.job_inputfile = os.path.abspath(job.inputfile)
        self.job_logfile = os.path.abspath(job.logfile)
        self.job_outputfile = os.path.abspath(job.outputfile)
        self.job_errfile = os.path.abspath(job.errfile)

    def _generate_visualization_style_script(self, job):
        # Define the source and destination file paths
        source_style_file = (
            self.pymol_templates_path / "zhang_group_pymol_style.py"
        )
        dest_style_file = os.path.join(
            job.folder, "zhang_group_pymol_style.py"
        )

        # Check if the style file already exists in the current working directory
        if not os.path.exists(dest_style_file):
            # Copy the file from templates to the current working directory
            logger.debug(
                f"Copying file from {source_style_file} to {dest_style_file}."
            )
            shutil.copy(source_style_file, dest_style_file)
        return dest_style_file

    def _write_input(self, job):
        # write to .xyz file if the supplied file is not .xyz
        if not os.path.exists(job.inputfile):
            mol = job.molecule
            # if mol is a list of molecules, then write to .xyz for all molecules
            if isinstance(mol, list):
                if isinstance(mol[0], Molecule):
                    logger.info(
                        f"Writing list of molecules: {mol} to {job.inputfile}"
                    )
                    for m in mol:
                        m.write(job.inputfile, format="xyz", mode="a")
                else:
                    raise ValueError(
                        f"Object {mol[0]} is not of Molecule type!"
                    )
            elif isinstance(mol, Molecule):
                logger.info(f"Writing Molecule to {job.inputfile}.")
                mol.write(job.inputfile, format="xyz", mode="w")
            else:
                raise ValueError(f"Object {mol} is not of Molecule type!")
        else:
            logger.warning(
                f"File {job.inputfile} already exists!\n"
                f"Will proceed to visualize this file instead!"
            )

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

    def _get_command(self, job):
        exe = self.executable
        command = f"{exe} {job.inputfile}"
        # get style file
        if job.pymol_script is None:
            if os.path.exists(
                os.path.join(job.folder, "zhang_group_pymol_style.py")
            ):
                job_pymol_script = os.path.join(
                    job.folder, "zhang_group_pymol_style.py"
                )
            else:
                logger.info(
                    "Using default zhang_group_pymol_style for rendering."
                )
                job_pymol_script = self._generate_visualization_style_script(
                    job
                )
            if os.path.exists(job_pymol_script):
                command += f" -r {job_pymol_script}"
        else:
            # using user-defined style file
            assert os.path.exists(
                job.pymol_script
            ), f"Supplied PyMOL Style file {job.pymol_script} does not exist!"
            command += f" -r {job.pymol_script}"

        if job.quiet_mode:
            command += " -q"
        if job.command_line_only:
            command += " -c"
        if job.render_style is None:
            if os.path.exists("zhang_group_pymol_style.py"):
                # defaults to using pymol_style if not specified
                command += f' -d "pymol_style {self.job_basename}'
            else:
                # no render style and no style file present
                command += ' -d "'
        else:
            if job.render_style.lower() == "pymol":
                command += f' -d "pymol_style {self.job_basename}'
            elif job.render_style.lower() == "cylview":
                command += f' -d "cylview_style {self.job_basename}'
            else:
                raise ValueError(
                    f"The style {job.render_style} is not available!"
                )

        if job.vdw:
            command += f"; add_vdw {self.job_basename}"

        # save pse file by default
        command += f'; zoom; save {job.outputfile}; quit"'

        return command
