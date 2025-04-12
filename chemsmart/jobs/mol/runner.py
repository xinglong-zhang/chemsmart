import glob
import logging
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.runner import JobRunner
from chemsmart.utils.periodictable import PeriodicTable
from chemsmart.utils.utils import quote_path

pt = PeriodicTable()

logger = logging.getLogger(__name__)


class PyMOLJobRunner(JobRunner):
    """Class for creating pymol job runner processs.
    Default to not use scratch for pymol jobs."""

    JOBTYPES = []
    PROGRAM = "pymol"
    FAKE = False
    SCRATCH = False

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
        """Define the path for the PyMOL executable, handling Windows-specific naming.
        Use pymol.exe on Windows, pymol on Unix-like systems."""
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

    def _get_visualization_command(self, job):
        exe = quote_path(self.executable)
        input_file = quote_path(job.inputfile)
        command = f"{exe} {input_file}"

        # Get style file
        if job.pymol_script is None:
            style_file_path = os.path.join(
                job.folder, "zhang_group_pymol_style.py"
            )
            if os.path.exists(style_file_path):
                job_pymol_script = style_file_path
            else:
                logger.info(
                    "Using default zhang_group_pymol_style for rendering."
                )
                job_pymol_script = self._generate_visualization_style_script(
                    job
                )
            if os.path.exists(job_pymol_script):
                command += f" -r {quote_path(job_pymol_script)}"
        else:
            # using user-defined style file
            assert os.path.exists(
                job.pymol_script
            ), f"Supplied PyMOL Style file {job.pymol_script} does not exist!"
            command += f" -r {quote_path(job.pymol_script)}"

        if job.quiet_mode:
            command += " -q"
        if job.command_line_only:
            command += " -c"

        # Handle the -d argument (PyMOL commands)
        if job.render is None:
            # defaults to using zhang_group_pymol_style if not specified
            if os.path.exists("zhang_group_pymol_style.py"):
                command += f' -d "pymol_style {self.job_basename}'
            else:
                # no render style and no style file present
                command += ' -d "'
        else:
            if job.render.lower() == "pymol":
                command += f' -d "pymol_style {self.job_basename}'
            elif job.render.lower() == "cylview":
                command += f' -d "cylview_style {self.job_basename}'
            else:
                raise ValueError(f"The style {job.render} is not available!")

        if job.vdw:
            command += f"; add_vdw {self.job_basename}"

        # zoom
        command += "; zoom"

        if job.trace:
            command += "; ray"

        return command

    def _save_pse_command(self, job, command):
        # Append the final PyMOL commands, quoting the output file path
        command += f"; save {quote_path(job.outputfile)}"

        return command

    def _quit_command(self, job, command):
        # Append the quit command to the PyMOL command
        command += '; quit"'
        return command

    def _create_process(self, job, command, env):
        # Open files for stdout/stderr
        with (
            open(self.job_errfile, "w") as err,
            open(self.job_outputfile, "w") as out,
        ):
            logger.info(
                f"Command executed: {command}\n"
                f"Writing output file to: {self.job_logfile}\n"
                f"And err file to: {self.job_errfile}"
            )
            # Start PyMOL process
            process = subprocess.Popen(
                shlex.split(command),
                stdin=subprocess.DEVNULL,  # prevent hanging
                stdout=out,
                stderr=err,
                env=env,
                cwd=self.running_directory,
            )
            # Wait for process to complete
            # process.wait()
            returncode = process.wait()
            logger.debug(f"PyMOL process exited with code {returncode}")

            if process.returncode != 0:
                logger.error(
                    f"PyMOL process failed with return code {process.returncode}"
                )
                raise subprocess.CalledProcessError(
                    process.returncode, command
                )
        return process


class PyMOLVisualizationJobRunner(PyMOLJobRunner):
    JOBTYPES = ["pymol_visualization"]

    def _get_command(self, job):
        command = self._get_visualization_command(job)
        command = self._save_pse_command(job, command)
        command = self._quit_command(job, command)
        return command


class PyMOLMovieJobRunner(PyMOLJobRunner):
    JOBTYPES = ["pymol_movie"]

    def _get_rotation_command(self, job, command):
        # rotation commands
        command += "; mset 1, 180"  # 360-degree rotation over 180 frames
        command += "; util.mroll 1, 180, 1"  # rotate about the z-axis
        if job.trace:
            command += "; set ray_trace_frames, 1"
        command += "; set cache_frames, 0"
        return command

    def _export_movie_command(self, job, command):
        """Export movie frames (use a temporary prefix to avoid conflicts)"""
        frame_prefix = os.path.join(job.folder, f"{self.job_basename}_frame_")
        command += f"; mpng {quote_path(frame_prefix)}"
        return command

    def _get_command(self, job):
        command = self._get_visualization_command(job)
        command = self._get_rotation_command(job, command)
        command = self._export_movie_command(job, command)
        command = self._save_pse_command(job, command)
        command = self._quit_command(job, command)
        return command

    def _postrun(self, job, framerate=30):
        """Convert to MP4 using ffmpeg."""
        frame_prefix = os.path.join(job.folder, f"{self.job_basename}_frame_")
        frame_pattern = f"{frame_prefix}%04d.png"
        output_mp4 = (
            os.path.splitext(job.outputfile)[0] + ".mp4"
        )  # Avoid .pse conflict

        # Check for PNG files
        png_files = glob.glob(f"{frame_prefix}*.png")
        if not png_files:
            raise FileNotFoundError(
                f"No PNG frames found matching '{frame_prefix}*.png'."
            )

        # Run ffmpeg command
        ffmpeg_cmd = [
            "ffmpeg",
            "-framerate",
            str(framerate),
            "-i",
            os.path.normpath(frame_pattern),  # Use %04d for ffmpeg
            "-c:v",
            "libx264",  # Uses H.264 codec for compatibility.
            "-pix_fmt",
            "yuv420p",  # ensures compatibility with most video players
            os.path.normpath(output_mp4),
        ]
        logger.info(f"Executing FFmpeg command: {' '.join(ffmpeg_cmd)}")
        try:
            subprocess.run(
                ffmpeg_cmd,
                check=True,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            logger.info(f"Movie saved as {output_mp4}")
        except subprocess.CalledProcessError as e:
            logger.error(f"FFmpeg error: {e.stderr}")
            raise
        except FileNotFoundError:
            raise FileNotFoundError(
                "FFmpeg not found. Ensure it’s installed in the environment."
            )

        # Clean up PNG files
        for png_file in png_files:
            os.remove(png_file)
        logger.info("Cleaned up PNG frames.")
