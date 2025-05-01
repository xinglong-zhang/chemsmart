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
from chemsmart.settings.executable import GaussianExecutable
from chemsmart.utils.periodictable import PeriodicTable
from chemsmart.utils.utils import (
    get_prepend_string_list_from_modred_free_format,
    quote_path,
    run_command,
)

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
        return command

    def _setup_style(self, job, command):
        # Handle the -d argument (PyMOL commands)
        if job.style is None:
            # defaults to using zhang_group_pymol_style if not specified
            if os.path.exists("zhang_group_pymol_style.py"):
                command += f' -d "pymol_style {self.job_basename}'
            else:
                # no render style and no style file present
                command += ' -d "'
        else:
            if job.style.lower() == "pymol":
                command += f' -d "pymol_style {self.job_basename}'
            elif job.style.lower() == "cylview":
                command += f' -d "cylview_style {self.job_basename}'
            else:
                raise ValueError(f"The style {job.style} is not available!")

        return command

    def _setup_viewport(self, command):
        command += "; viewport 800,600"
        return command

    def _add_vdw(self, job, command):
        if job.vdw:
            command += f"; add_vdw {self.job_basename}"

        return command

    def _add_coordinates_labels(self, job, command):
        distances = []
        angles = []
        dihedrals = []
        if job.coordinates:
            prepend_string_list = (
                get_prepend_string_list_from_modred_free_format(
                    input_modred=job.coordinates, program="pymol"
                )
            )
            for prepend_string in prepend_string_list:
                if prepend_string.startswith("B"):
                    distances.append(
                        [int(i) for i in prepend_string.split()[1:]]
                    )
                elif prepend_string.startswith("A"):
                    angles.append([int(i) for i in prepend_string.split()[1:]])
                elif prepend_string.startswith("D"):
                    dihedrals.append(
                        [int(i) for i in prepend_string.split()[1:]]
                    )

            for i, distance in enumerate(distances):
                command += (
                    f"; distance d{i+1}, id {distance[0]}, id {distance[1]}"
                )
            for i, angle in enumerate(angles):
                command += f"; angle a{i+1}, id {angle[0]}, id {angle[1]}, id {angle[2]}"
            for i, dihedral in enumerate(dihedrals):
                command += f"; dihedral di{i+1}, id {dihedral[0]}, id {dihedral[1]}, id {dihedral[2]}, id {dihedral[3]}"
        return command

    def _offset_labels(self, job, command, x=-1.2):
        command += f"; set label_position, ({x},0,0)"
        return command

    def _add_zoom_command(self, job, command):
        # zoom
        command += "; zoom"
        return command

    def _add_ray_command(self, job, command):
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
        command = self._setup_style(job, command)
        command = self._setup_viewport(command)
        command = self._add_coordinates_labels(job, command)
        command = self._offset_labels(job, command)
        command = self._add_vdw(job, command)
        command = self._add_zoom_command(job, command)
        command = self._job_specific_commands(job, command)
        command = self._save_pse_command(job, command)
        command = self._quit_command(job, command)
        return command

    def _hide_labels(self, job, command):
        """Hide labels for NCI analysis."""
        command += "; hide labels"
        return command

    def _job_specific_commands(self, job, command):
        """Job specific commands."""
        command = self._add_ray_command(job, command)
        return command


class PyMOLMovieJobRunner(PyMOLVisualizationJobRunner):
    JOBTYPES = ["pymol_movie"]

    def _setup_style(self, job, command):
        if os.path.exists("zhang_group_pymol_style.py"):
            command += f' -d "movie_style {self.job_basename}'
        else:
            # no render style and no style file present
            command += ' -d "'
        return command

    def _job_specific_commands(self, job, command):
        """Job specific commands."""
        command = self._get_rotation_command(job, command)
        command = self._set_ray_trace_frames(job, command)
        command = self._export_movie_command(job, command)
        return command

    def _get_rotation_command(self, job, command):
        # rotation commands
        command += "; mset 1, 180"  # 360-degree rotation over 180 frames
        command += "; util.mroll 1, 180, 1"  # rotate about the z-axis
        command += "; set cache_frames, 0"
        return command

    def _set_ray_trace_frames(self, job, command):
        if job.trace:
            command += "; set ray_trace_frames, 1; set ray_trace_mode, 1"
        return command

    def _export_movie_command(self, job, command):
        """Export movie frames (use a temporary prefix to avoid conflicts)"""
        frame_prefix = os.path.join(job.folder, f"{self.job_basename}_frame_")
        command += f"; mpng {frame_prefix}"
        return command

    def _postrun(self, job, framerate=30):
        """Convert to MP4 using ffmpeg."""
        self._create_movie(job=job, framerate=framerate)

    def _create_movie(self, job, framerate=30, constant_rate_factor=18):
        frame_prefix = os.path.join(job.folder, f"{self.job_basename}_frame_")
        frame_pattern = f"{frame_prefix}%04d.png"
        output_mp4 = (
            os.path.splitext(job.outputfile)[0] + ".mp4"
        )  # Avoid .pse conflict

        # Check if mp4 already exists and include options for overwriting
        if os.path.exists(output_mp4):
            if not job.overwrite:
                logger.warning(
                    f"{output_mp4} already exists. Skipping movie creation."
                )
                return
            else:
                logger.info(f"{output_mp4} exists, but will be overwritten.")

        # Check for PNG files
        png_files = glob.glob(f"{frame_prefix}*.png")
        if not png_files:
            raise FileNotFoundError(
                f"No PNG frames found matching '{frame_prefix}*.png'."
            )

        # Run ffmpeg command
        ffmpeg_cmd = [
            "ffmpeg",
            "-y",  # Overwrite output file if it exists
            "-framerate",
            str(framerate),
            "-i",
            os.path.normpath(frame_pattern),  # Use %04d for ffmpeg
            "-c:v",
            "libx264",  # Uses H.264 codec for compatibility.
            "-pix_fmt",
            "yuv420p",  # ensures compatibility with most video players
            "-crf",
            str(
                constant_rate_factor
            ),  # High quality, adjust as needed (0-51, lower = better)
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


class PyMOLIRCMovieJobRunner(PyMOLMovieJobRunner):
    JOBTYPES = ["pymol_ircmovie"]

    def _get_rotation_command(self, job, command):
        # no rotation commands for ircmovie
        return command


class PyMOLNCIJobRunner(PyMOLVisualizationJobRunner):
    JOBTYPES = ["pymol_nci"]

    def _job_specific_commands(self, job, command):
        """Job specific commands."""
        command = self._hide_labels(job, command)
        command = self._load_cube_files(job, command)
        command = self._run_nci_command(job, command)
        command = self._add_ray_command(job, command)
        return command

    def _load_cube_files(self, job, command):
        """Load cube files for NCI analysis."""
        dens_file = os.path.join(job.folder, f"{self.job_basename}-dens.cube")
        grad_file = os.path.join(job.folder, f"{self.job_basename}-grad.cube")
        assert os.path.exists(
            dens_file
        ), f"Density cube file {dens_file} not found!"
        assert os.path.exists(
            grad_file
        ), f"Gradient cube file {grad_file} not found!"

        command += f"; load {quote_path(dens_file)}"
        command += f"; load {quote_path(grad_file)}"

        return command

    def _run_nci_command(self, job, command):
        """Run the nci command"""
        if job.binary:
            command += f"; nci_binary {self.job_basename}"
        else:
            command += f"; nci {self.job_basename}"
        return command


class PyMOLMOJobRunner(PyMOLVisualizationJobRunner):
    JOBTYPES = ["pymol_mo"]

    def _get_gaussian_executable(self, job):
        """Get the Gaussian executable for the job."""
        logger.info(
            f"Obtaining Gaussian executable from server: {self.server.name}"
        )
        gaussian_exe = GaussianExecutable.from_servername(self.server.name)
        gaussian_exe_path = gaussian_exe.executable_folder
        return gaussian_exe_path

    def _prerun(self, job):
        self._assign_variables(job)
        self._generate_fchk_file(job)
        self._generate_mo_cube_file(job)
        self._write_molecular_orbital_pml(job)

    def _generate_fchk_file(self, job):
        """Generate the fchk file from the chk file."""
        assert os.path.exists(
            os.path.join(job.folder, f"{self.job_basename}.chk")
        ), ".chk file is required but not found!"

        gaussian_exe = self._get_gaussian_executable(job)
        if os.path.exists(
            os.path.join(job.folder, f"{self.job_basename}.fchk")
        ):
            pass
        else:
            # generate .fchk file from .chk file
            logger.info(f"Generating .fchk file from {self.job_basename}.chk")
            fchk_command = f"{gaussian_exe}/formchk {self.job_basename}.chk"
            run_command(fchk_command)

    def _generate_mo_cube_file(self, job):
        """Generate the MO cube file."""
        gaussian_exe = self._get_gaussian_executable(job)
        if [job.number, job.homo, job.lumo].count(None) < 1:
            raise ValueError(
                "Cannot specify more than two of MO number, HOMO, and LUMO together."
            )

        if job.number:
            if os.path.exists(f"{self.job_basename}_MO{job.number}.cube"):
                logger.info(
                    f"cube file {self.job_basename}_MO{job.number}.cube already exists."
                )
                pass
            else:
                cubegen_command = f"{gaussian_exe}/cubegen 0 MO={job.number} {self.job_basename}.fchk {self.job_basename}_MO{job.number}.cube 0 h"
                run_command(cubegen_command)

        if job.homo:
            if os.path.exists(f"{self.job_basename}_HOMO.cube"):
                logger.info(
                    f"cube file {self.job_basename}_HOMO.cube already exists."
                )
                pass
            else:
                cubegen_command = f"{gaussian_exe}/cubegen 0 MO=HOMO {self.job_basename}.fchk {self.job_basename}_HOMO.cube 0 h"
                run_command(cubegen_command)
        if job.lumo:
            if os.path.exists(f"{self.job_basename}_LUMO.cube"):
                logger.info(
                    f"cube file {self.job_basename}_LUMO.cube already exists."
                )
                pass
            else:
                cubegen_command = f"{gaussian_exe}/cubegen 0 MO=LUMO {self.job_basename}.fchk {self.job_basename}_LUMO.cube 0 h"
                run_command(cubegen_command)

    def _write_molecular_orbital_pml(
        self, job, isosurface=0.05, transparency=0.2
    ):
        pml_file = os.path.join(job.folder, f"{self.job_basename}.pml")
        if not os.path.exists(pml_file):
            with open(pml_file, "w") as f:
                f.write(f"load {self.job_basename}.cube\n")
                f.write(
                    f"isosurface pos_iso, {self.job_basename}, {isosurface}\n"
                )
                f.write(
                    f"isosurface neg_iso, {self.job_basename}, {-isosurface}\n"
                )
                f.write("print(pos_iso)\n")
                f.write("print(neg_iso)\n")
                f.write("set surface_color, blue, pos_iso\n")
                f.write("set surface_color, red, neg_iso\n")
                f.write(f"set transparency, {transparency}\n")
            logger.info(f"Wrote PML file: {pml_file}")

    def _job_specific_commands(self, job, command):
        """Job specific commands."""
        command = self._hide_labels(job, command)
        command = self._call_pml(job, command)
        command = self._add_ray_command(job, command)
        return command

    def _call_pml(self, job, command):
        """Call the PML file for visualization."""
        pml_file = os.path.join(job.folder, f"{self.job_basename}.pml")
        command += f"; load {quote_path(pml_file)}"
        return command
