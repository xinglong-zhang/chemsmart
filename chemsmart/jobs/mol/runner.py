"""
PyMOL job runners for molecular visualization and analysis.

This module implements runners that execute PyMOL-based jobs, including
static visualizations, rotating movies, IRC trajectory animations, NCI
isosurface visualizations, MO plots, and spin density renderings. The
base PyMOLJobRunner provides command construction, environment setup,
and process management; specialized runners extend it with job-specific
commands and post-processing (e.g., FFmpeg movie creation).
"""

import glob
import logging
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.executable import GaussianExecutable
from chemsmart.utils.periodictable import PeriodicTable
from chemsmart.utils.repattern import (
    pymol_color_range_pattern,
    pymol_isosurface_pattern,
)
from chemsmart.utils.utils import (
    get_prepend_string_list_from_modred_free_format,
    quote_path,
    run_command,
)

pt = PeriodicTable()

logger = logging.getLogger(__name__)


class PyMOLJobRunner(JobRunner):
    """
    Job runner for executing PyMOL molecular visualization tasks.

    Provides comprehensive functionality for running PyMOL jobs including
    molecular visualization, animation creation, and specialized analysis
    visualizations. Handles executable detection, command generation,
    and process management for various PyMOL-based tasks.

    Attributes:
        JOBTYPES (list): Supported job type identifiers (set by subclasses).
        PROGRAM (str): Program identifier ('pymol').
        FAKE (bool): Whether this runner operates in fake/test mode.
        SCRATCH (bool): Whether to use scratch directories by default.
        server: Server configuration used for execution.
        scratch (bool): Whether scratch is enabled for this runner.
        scratch_dir (str): Path to scratch directory, if used.
        num_cores (int): Number of CPU cores allocated (from server).
        num_gpus (int): Number of GPUs allocated (from server).
        mem_gb (int): Memory allocation in gigabytes (from server).
    """

    JOBTYPES = []
    PROGRAM = "pymol"
    FAKE = False
    SCRATCH = False

    def __init__(
        self,
        server,
        scratch=None,
        fake=False,
        **kwargs,
    ):
        """
        Initialize the PyMOL job runner.

        Sets up the runner with server configuration and execution
        parameters, defaulting to no scratch directory usage for
        PyMOL jobs unless explicitly specified.

        Args:
            server: Server configuration for job execution.
            scratch: Whether to use scratch directories (default: None).
            fake: Whether this is a fake runner for testing (default: False).
            **kwargs: Additional arguments passed to parent JobRunner.
        """
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
        """
        Get the path to the PyMOL executable.

        Automatically detects PyMOL installation on the system,
        handling platform-specific executable names (pymol.exe on
        Windows, pymol on Unix-like systems).

        Returns:
            str: Path to the PyMOL executable.

        Raises:
            FileNotFoundError: If PyMOL is not found in system PATH.
        """
        pymol_cmd = "pymol.exe" if sys.platform == "win32" else "pymol"
        pymol_path = shutil.which(pymol_cmd)
        if pymol_path is None or not os.path.exists(pymol_path):
            raise FileNotFoundError(
                f"PyMOL executable '{pymol_cmd}' not found in PATH. "
                f"Please install PyMOL!"
            )
        return pymol_path

    @property
    def pymol_templates_path(self):
        """
        Get the path to the PyMOL templates directory.

        Returns:
            Path: Absolute path to the templates directory containing
                  PyMOL visualization style scripts.
        """
        return Path(__file__).resolve().parent / "templates"

    def _prerun(self, job):
        """
        Perform pre-execution setup for the PyMOL job.

        Configures necessary file paths and variables before job
        execution begins.

        Args:
            job: PyMOL job object to configure.
        """
        self._assign_variables(job)

    def _assign_variables(self, job):
        """
        Set up file paths for PyMOL job execution.

        Configures absolute paths for input, output, log, and error
        files based on the job configuration and folder structure.

        Args:
            job: PyMOL job object containing file path information.
        """
        self.running_directory = job.folder
        logger.debug(f"Running directory: {self.running_directory}")

    def _generate_visualization_style_script(self, job):
        """
        Generate or copy the PyMOL visualization style script.

        Copies the default Zhang group PyMOL style script to the job
        directory. Modify the PyMOL style script for job (calls
        self._modify_job_pymol_script() if job.isosurface_value
        or job.range is not None.
        Works for using default zhang_group_pymol_style.py.
        Overwrites existing zhang_group_pymol_style.py style.

        Args:
            job: PyMOL job object containing folder information.

        Returns:
            str: Path to the style script file in the job directory.
        """
        # Define the source and destination file paths
        source_style_file = (
            self.pymol_templates_path / "zhang_group_pymol_style.py"
        )
        dest_style_file = os.path.join(
            job.folder, "zhang_group_pymol_style.py"
        )

        # Check if the style file already exists in the current working
        # directory
        if os.path.exists(dest_style_file):
            logger.warning(
                f"Style file {dest_style_file} already exists! Overwriting..."
            )
        logger.debug(
            f"Copying file from {source_style_file} to {dest_style_file}."
        )
        shutil.copy(source_style_file, dest_style_file)

        # check if job has attribute of isosurface_value or color_range
        if not hasattr(job, "isosurface_value"):
            job.isosurface_value = None
        if not hasattr(job, "color_range"):
            job.color_range = None

        logger.debug(f"Job isosurface_value: {job.isosurface_value}")
        logger.debug(f"Job color_range: {job.color_range}")

        if job.isosurface_value is None and job.color_range is None:
            return dest_style_file
        else:
            return self._modify_job_pymol_script(job, dest_style_file)

    def _modify_job_pymol_script(self, job, style_file_path=None):
        """Modify the PyMOL style script for the job.
        Modifies the PyMOL style script to set the isosurface value
        and color range if they are provided in the job.
        Args:
            job: PyMOL job object containing isosurface and color range.
            style_file_path: Path to the style script file to modify.
            If None, defaults to zhang_group_pymol_style.py in the job folder.
        Returns:
            str: Path to the modified style script file."""

        if style_file_path is None:
            style_file_path = os.path.join(
                job.folder, "zhang_group_pymol_style.py"
            )
        if not os.path.exists(style_file_path):
            raise FileNotFoundError(
                f"Style file {style_file_path} does not exist!"
            )

        with open(style_file_path, "r") as file:
            style_script = file.read()

        new_script = style_script

        if job.isosurface_value is not None:
            new_script = re.sub(
                pymol_isosurface_pattern,
                f"isosurface={job.isosurface_value}",
                new_script,
            )

        if job.color_range is not None:
            new_script = re.sub(
                pymol_color_range_pattern,
                f"range={job.color_range}",
                new_script,
            )

        if new_script != style_script:
            logger.warning(
                f"Modifying style file {style_file_path} "
                f"for job {job.label}."
            )
            with open(style_file_path, "w") as file:
                file.write(new_script)

        return style_file_path

    def _get_gaussian_executable(self, job):
        """
        Get the Gaussian executable path for cube file generation.

        Retrieves the Gaussian executable configuration needed for
        generating cube files from formatted checkpoint files.

        Args:
            job: PyMOL MO job object.

        Returns:
            str: Path to the Gaussian executable directory.
        """
        logger.info(
            f"Obtaining Gaussian executable from server: {self.server.name}"
        )
        gaussian_exe = GaussianExecutable.from_servername(self.server.name)
        gaussian_exe_path = gaussian_exe.executable_folder
        return gaussian_exe_path

    def _generate_fchk_file(self, job):
        """
        Generate the formatted checkpoint file from Gaussian checkpoint.

        Uses Gaussian's formchk utility to convert binary checkpoint
        files to formatted checkpoint files required for cube generation.

        Args:
            job: Job object containing checkpoint file information.

        Raises:
            FileNotFoundError: If the required .chk file is not found.
        """
        chk_file_path = os.path.join(job.folder, f"{job.label}.chk")
        if not os.path.exists(chk_file_path):
            raise FileNotFoundError(
                f".chk file is required but not found at {chk_file_path}!"
            )

        gaussian_exe = self._get_gaussian_executable(job)
        if os.path.exists(os.path.join(job.folder, f"{job.label}.fchk")):
            logger.info(
                f".fchk file {job.label}.fchk already exists.\n"
                f"Skipping generation of .fchk file."
            )
            pass
        else:
            # generate .fchk file from .chk file
            logger.info(f"Generating .fchk file from {job.label}.chk")
            fchk_command = f"{gaussian_exe}/formchk {job.label}.chk"
            run_command(fchk_command)

    def _write_input(self, job):
        """
        Write molecular structure data to XYZ input file.

        Creates an XYZ coordinate file from the job's molecular data,
        handling both single molecules and lists of molecules with
        appropriate formatting and logging.

        Args:
            job: PyMOL job object containing molecular structure data.

        Raises:
            ValueError: If the molecule object is not of correct type.
        """
        # write to .xyz file if the supplied file is not .xyz
        if not os.path.exists(job.inputfile):
            mol = job.molecule
            # if mol is a list of molecules, then write to .xyz for all
            # molecules
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
        """
        Update operating system environment variables.

        PyMOL jobs typically don't require special environment
        variable configuration, so this method is a no-op.

        Args:
            job: PyMOL job object (unused for PyMOL).
        """
        # no envs to update for pymol
        pass

    def _get_visualization_command(self, job):
        """
        Generate the base PyMOL visualization command.

        Constructs the initial PyMOL command with input file and
        style script parameters, handling file path quoting and
        command-line options.

        Args:
            job: PyMOL job object with visualization parameters.

        Returns:
            str: Base PyMOL command string for execution.
        """
        exe = quote_path(self.executable)
        input_file = quote_path(job.inputfile)
        command = f"{exe} {input_file}"

        # Get style file
        if job.pymol_script is None:
            logger.info(
                f"No style file supplied for job {job.label}."
                f"Using default zhang_group_pymol_style.py."
            )
            job_pymol_script = self._generate_visualization_style_script(job)
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
        """
        Configure PyMOL visualization style commands.

        Adds style-specific PyMOL commands to the command string,
        supporting different visualization styles like pymol and cylview.

        Args:
            job: PyMOL job object with style configuration.
            command: Base command string to extend.

        Returns:
            str: Command string with style configuration added.

        Raises:
            ValueError: If an unsupported style is specified.
        """
        # Handle the -d argument (PyMOL commands)

        if job.style is None:
            # defaults to using zhang_group_pymol_style if not specified
            if os.path.exists("zhang_group_pymol_style.py"):
                command += f' -d "pymol_style {job.label}'
            else:
                # no render style and no style file present
                command += ' -d "'
        else:
            if job.style.lower() == "pymol":
                command += f' -d "pymol_style {job.label}'
            elif job.style.lower() == "cylview":
                command += f' -d "cylview_style {job.label}'
            else:
                raise ValueError(f"The style {job.style} is not available!")

        return command

    def _setup_viewport(self, command):
        """
        Configure PyMOL viewport dimensions.

        Sets the PyMOL viewport to a standard size for consistent
        visualization output across different systems.

        Args:
            command: Command string to extend with viewport settings.

        Returns:
            str: Command string with viewport configuration added.
        """
        command += "; viewport 800,600"
        return command

    def _add_vdw(self, job, command):
        """
        Add van der Waals surface visualization to PyMOL command.

        Conditionally adds VDW surface representation to the
        molecular visualization if requested in the job configuration.

        Args:
            job: PyMOL job object with VDW settings.
            command: Command string to extend.

        Returns:
            str: Command string with VDW surface added if requested.
        """
        if job.vdw:
            command += f"; add_vdw {job.label}"

        return command

    def _add_coordinates_labels(self, job, command):
        """
        Add coordinate labels for distances, angles, and dihedrals.

        Parses coordinate specifications from the job and adds
        appropriate PyMOL commands for labeling geometric measurements
        like bond distances, angles, and dihedral angles.

        Args:
            job: PyMOL job object with coordinate specifications.
            command: Command string to extend.

        Returns:
            str: Command string with coordinate labels added.
        """
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
        """
        Set label position offset in PyMOL visualization.

        Adjusts the position of text labels relative to atoms or
        measurements for better visual clarity in the output.

        Args:
            job: PyMOL job object (used for consistency).
            command: Command string to extend.
            x: X-axis offset for label positioning (default: -1.2).

        Returns:
            str: Command string with label positioning added.
        """
        command += f"; set label_position, ({x},0,0)"
        return command

    def _add_zoom_command(self, job, command):
        """
        Add zoom command to fit molecule in PyMOL viewport.

        Ensures the molecular structure is properly centered and
        scaled to fit within the visualization viewport.

        Args:
            job: PyMOL job object (used for consistency).
            command: Command string to extend.

        Returns:
            str: Command string with zoom command added.
        """
        # zoom
        command += "; zoom"
        return command

    def _add_ray_command(self, job, command):
        """
        Add high-resolution ray tracing command if requested.

        Conditionally adds ray tracing for high-quality rendering
        when the trace option is enabled in the job configuration.

        Args:
            job: PyMOL job object with trace settings.
            command: Command string to extend.

        Returns:
            str: Command string with ray tracing added if requested.
        """
        if job.trace:
            command += "; ray 2400,1800"

        return command

    def _add_refresh_command(self, job, command):
        """
        Add refresh command to ensure command completion.

        Forces PyMOL to complete all pending commands before
        proceeding to the next operation.

        Args:
            job: PyMOL job object (used for consistency).
            command: Command string to extend.

        Returns:
            str: Command string with refresh command added.
        """
        command += "; refresh"
        return command

    def _save_pse_command(self, job, command):
        """
        Add save command to write PyMOL session file.

        Appends the command to save the current PyMOL session
        to the specified output file with proper path quoting.

        Args:
            job: PyMOL job object with output file path.
            command: Command string to extend.

        Returns:
            str: Command string with save command added.
        """
        # Append the final PyMOL commands, quoting the output file path
        command += f"; save {quote_path(job.outputfile)}"

        return command

    def _quit_command(self, job, command):
        """
        Add quit command to terminate PyMOL session.

        Appends the final quit command to cleanly exit PyMOL
        after all visualization tasks are completed.

        Args:
            job: PyMOL job object (used for consistency).
            command: Command string to extend.

        Returns:
            str: Command string with quit command added.
        """
        # Append the quit command to the PyMOL command
        command += '; quit"'
        return command

    def _create_process(self, job, command, env):
        """
        Create and execute the PyMOL subprocess.

        Opens output and error files, starts the PyMOL process with
        the constructed command, and waits for completion while
        handling error conditions and logging.

        Args:
            job: PyMOL job object with file paths.
            command: Complete PyMOL command string to execute.
            env: Environment variables for the process.

        Returns:
            subprocess.Popen: The completed process object.

        Raises:
            subprocess.CalledProcessError: If PyMOL process fails.
        """
        # Open files for stdout/stderr
        job_errfile = os.path.abspath(job.errfile)
        job_outputfile = os.path.abspath(job.outputfile)
        with (
            open(job_errfile, "w") as err,
            open(job_outputfile, "w") as out,
        ):
            logger.info(
                f"Command executed: {command}\n"
                f"Writing output file to: {os.path.abspath(job.logfile)}\n"
                f"And err file to: {job_errfile}"
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

    def _get_base_filepath_to_remove(self, job):
        """Get the base filepath for the job to assist in file removal."""
        return Path(job.folder) / job.job_basename


class PyMOLVisualizationJobRunner(PyMOLJobRunner):
    """
    Specialized PyMOL job runner for basic molecular visualization.

    Extends the base PyMOL runner to provide standard molecular
    visualization capabilities with customizable styling, labeling,
    and rendering options for static molecular images.
    """

    JOBTYPES = ["pymol_visualization"]

    def _get_command(self, job):
        """
        Build the complete PyMOL command for visualization.

        Assembles all visualization components including style setup,
        viewport configuration, labels, and rendering options into
        a single PyMOL command string.

        Args:
            job: PyMOL visualization job object.

        Returns:
            str: Complete PyMOL command string for execution.
        """
        command = self._get_visualization_command(job)
        command = self._setup_style(job, command)
        command = self._setup_viewport(command)
        command = self._add_coordinates_labels(job, command)
        command = self._offset_labels(job, command, x=-1.2)
        command = self._add_vdw(job, command)
        command = self._add_zoom_command(job, command)
        command = self._job_specific_commands(job, command)
        command = self._save_pse_command(job, command)
        command = self._quit_command(job, command)
        return command

    def _hide_labels(self, job, command):
        """
        Hide atom and measurement labels in PyMOL visualization.

        Adds command to hide all labels for cleaner visualization,
        particularly useful for specialized analysis views.

        Args:
            job: PyMOL job object (used for consistency).
            command: Command string to extend.

        Returns:
            str: Command string with label hiding added.
        """
        command += "; hide labels"
        return command

    def _job_specific_commands(self, job, command):
        """
        Add job-specific commands for basic visualization.

        Incorporates visualization-specific commands like ray tracing
        that are unique to basic molecular visualization tasks.

        Args:
            job: PyMOL visualization job object.
            command: Command string to extend.

        Returns:
            str: Command string with visualization-specific commands.
        """
        command = self._add_ray_command(job, command)
        return command


class PyMOLMovieJobRunner(PyMOLVisualizationJobRunner):
    """
    Specialized PyMOL job runner for creating molecular animations.

    Extends the visualization runner to generate rotating molecular
    movies with optional ray tracing, handling frame generation
    and video conversion using FFmpeg.
    """

    JOBTYPES = ["pymol_movie"]

    def _setup_style(self, job, command):
        """
        Configure PyMOL style specifically for movie generation.

        Sets up movie-specific styling that works well for animated
        visualizations, using appropriate templates if available.

        Args:
            job: PyMOL movie job object.
            command: Command string to extend.

        Returns:
            str: Command string with movie style configuration.
        """
        if os.path.exists("zhang_group_pymol_style.py"):
            command += f' -d "movie_style {job.job_basename}'
        else:
            # no render style and no style file present
            command += ' -d "'
        return command

    def _job_specific_commands(self, job, command):
        """
        Add movie-specific commands for animation generation.

        Includes rotation setup, frame generation, ray tracing
        configuration, and movie export commands specific to
        creating molecular animations.

        Args:
            job: PyMOL movie job object.
            command: Command string to extend.

        Returns:
            str: Command string with movie-specific commands.
        """
        command = self._get_rotation_command(job, command)
        command = self._set_ray_trace_frames(job, command)
        command = self._add_refresh_command(job, command)
        command = self._export_movie_command(job, command)
        return command

    def _get_rotation_command(self, job, command):
        """
        Add rotation commands for creating molecular movies.

        Sets up PyMOL animation commands to create a 360-degree
        rotation over 180 frames for smooth molecular visualization.

        Args:
            job: PyMOL movie job object.
            command: Command string to extend.

        Returns:
            str: Command string with rotation commands added.
        """
        # rotation commands
        command += "; mset 1, 180"  # 360-degree rotation over 180 frames
        command += "; util.mroll 1, 180, 1"  # rotate about the z-axis
        command += "; set cache_frames, 0"
        return command

    def _set_ray_trace_frames(self, job, command):
        """
        Configure ray tracing for movie frames if requested.

        Enables high-quality ray tracing for each frame when the
        trace option is enabled, producing higher quality animations.

        Args:
            job: PyMOL movie job object with trace settings.
            command: Command string to extend.

        Returns:
            str: Command string with ray tracing configuration.
        """
        if job.trace:
            command += "; set ray_trace_frames, 1; set ray_trace_mode, 1"
        return command

    def _export_movie_command(self, job, command):
        """
        Add command to export individual movie frames as PNG files.

        Configures PyMOL to export each frame of the animation as
        separate PNG files with a consistent naming convention.

        Args:
            job: PyMOL movie job object.
            command: Command string to extend.

        Returns:
            str: Command string with frame export command.
        """
        frame_prefix = os.path.join(job.folder, f"{job.job_basename}_frame_")
        command += f"; mpng {frame_prefix}"
        return command

    def _postrun(self, job, framerate=30):
        """
        Perform post-processing to convert frames to MP4 video.

        Executes FFmpeg conversion of individual PNG frames into
        a compressed MP4 video file with specified frame rate.

        Args:
            job: PyMOL movie job object.
            framerate: Video frame rate in FPS (default: 30).
        """
        self._create_movie(job=job, framerate=framerate)

    def _create_movie(self, job, framerate=30, constant_rate_factor=18):
        """
        Create an MP4 movie from exported PNG frames using FFmpeg.

        Builds and executes an FFmpeg command to convert the per-frame
        PNG images written by PyMOL (via `mpng`) into a compressed MP4
        video next to the session file. Optionally overwrites existing
        MP4 files and removes PNG frames after successful conversion.

        Args:
            job: PyMOL movie job holding folder, label, and overwrite flag.
            framerate (int): Output video framerate in frames per second.
            constant_rate_factor (int): FFmpeg CRF value (0–51, lower = better).

        Raises:
            FileNotFoundError: If no PNG frames are found or FFmpeg is missing.
            subprocess.CalledProcessError: If FFmpeg fails to encode the video.
        """
        frame_prefix = os.path.join(job.folder, f"{job.job_basename}_frame_")
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
    """
    Specialized PyMOL job runner for IRC trajectory animations.

    Extends the movie runner for visualizing Intrinsic Reaction
    Coordinate (IRC) trajectories without rotation, focusing on
    the molecular changes along the reaction path.
    """

    JOBTYPES = ["pymol_ircmovie"]

    def _get_rotation_command(self, job, command):
        """
        Override rotation for IRC movies.

        IRC movies display reaction trajectories and don't need
        rotation effects, so this method returns the command unchanged.

        Args:
            job: PyMOL IRC movie job object.
            command: Command string to return unchanged.

        Returns:
            str: Unmodified command string (no rotation for IRC).
        """
        # no rotation commands for ircmovie
        return command


class PyMOLNCIJobRunner(PyMOLVisualizationJobRunner):
    """
    Specialized PyMOL job runner for NCI (Non-Covalent Interactions) analysis.

    Provides visualization capabilities for non-covalent interaction
    analysis using density and gradient cube files with specialized
    PyMOL commands for isosurface generation and coloring.
    """

    JOBTYPES = ["pymol_nci"]

    def _job_specific_commands(self, job, command):
        """
        Add NCI-specific visualization commands.

        Incorporates commands for loading cube files, hiding labels,
        and generating NCI isosurfaces with appropriate coloring
        for interaction analysis.

        Args:
            job: PyMOL NCI job object.
            command: Command string to extend.

        Returns:
            str: Command string with NCI-specific commands.
        """
        command = self._hide_labels(job, command)
        command = self._load_cube_files(job, command)
        command = self._add_refresh_command(job, command)  # adding sync fails
        command = self._run_nci_command(job, command)
        command = self._add_ray_command(job, command)
        return command

    def _load_cube_files(self, job, command):
        """
        Load density and gradient cube files for NCI analysis.

        Loads the required density and gradient cube files that
        contain the electronic structure data needed for NCI
        visualization and analysis.

        Args:
            job: PyMOL NCI job object.
            command: Command string to extend.

        Returns:
            str: Command string with cube file loading commands.

        Raises:
            AssertionError: If required cube files are not found.
        """
        dens_file = os.path.join(job.folder, f"{job.label}-dens.cube")
        grad_file = os.path.join(job.folder, f"{job.label}-grad.cube")
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
        """
        Execute the appropriate NCI visualization command.

        Runs the NCI analysis command with the appropriate mode
        (binary, intermediate, or standard) based on job configuration
        to generate non-covalent interaction visualizations.

        Args:
            job: PyMOL NCI job object with analysis mode settings.
            command: Command string to extend.

        Returns:
            str: Command string with NCI visualization command.
        """
        if job.binary:
            command += f"; nci_binary {job.label}"
        elif job.intermediate:
            command += f"; nci_intermediate {job.label}"
        else:
            command += f"; nci {job.label}"
        return command


class PyMOLMOJobRunner(PyMOLVisualizationJobRunner):
    """
    Specialized PyMOL job runner for molecular orbital visualization.

    Provides visualization capabilities for molecular orbitals by
    integrating with Gaussian calculations to generate cube files
    and create isosurface visualizations of HOMO, LUMO, or specific
    molecular orbitals.
    """

    JOBTYPES = ["pymol_mo"]

    def _prerun(self, job):
        """
        Perform pre-execution setup for molecular orbital visualization.

        Prepares all required artifacts before launching PyMOL by
        assigning file paths, converting the binary checkpoint to
        a formatted checkpoint (.fchk), generating the MO cube file,
        and writing a PML script that defines the isosurfaces.

        Args:
            job: PyMOL MO job object to prepare.
        """
        self._assign_variables(job)
        self._generate_fchk_file(job)
        self._generate_mo_cube_file(job)
        self._write_molecular_orbital_pml(job)

    def _generate_mo_cube_file(self, job):
        """
        Generate a molecular orbital cube file via Gaussian's cubegen.

        Generates a cube file for either a specific MO number, HOMO,
        or LUMO, based on the job settings. Exactly one of these
        options must be selected; otherwise a ValueError is raised.

        Args:
            job: PyMOL MO job object with orbital selection flags.

        Raises:
            ValueError: If zero or multiple orbital selections are provided.
        """
        gaussian_exe = self._get_gaussian_executable(job)

        selected = sum(
            [
                job.number is not None,
                job.homo,
                job.lumo,
            ]
        )

        if selected != 1:
            raise ValueError(
                "You must specify exactly one of --number, --homo, or --lumo."
            )

        if os.path.exists(f"{job.job_basename}.cube"):
            logger.info(f"cube file {job.job_basename}.cube already exists.")
            return

        mo_type = None
        if job.number:
            mo_type = f"{job.number}"

        if job.homo:
            mo_type = "HOMO"

        if job.lumo:
            mo_type = "LUMO"

        if mo_type is None:
            raise ValueError(
                f"No MO specified, MO type/number given is: {mo_type}"
            )

        cubegen_command = (
            f"{gaussian_exe}/cubegen 0 MO={mo_type} {job.label}.fchk "
            f"{job.job_basename}.cube 0 h"
        )

        run_command(cubegen_command)

    def _write_molecular_orbital_pml(
        self,
        job,
    ):
        """
        Write a PML script to visualize the MO isosurfaces.

        Creates a PyMOL script that loads the generated MO cube file
        and defines positive/negative isosurfaces with distinct colors
        and a configurable transparency level.

        Args:
            job: PyMOL MO job object.
        """

        pml_file = os.path.join(job.folder, f"{job.mo_basename}.pml")
        if os.path.exists(pml_file):
            logger.warning(f"PML file {pml_file} already exists! Overwriting.")
        with open(pml_file, "w") as f:
            f.write(f"load {job.mo_basename}.cube\n")
            f.write(
                f"isosurface pos_iso, {job.mo_basename}, {job.isosurface_value}\n"
            )
            f.write(
                f"isosurface neg_iso, {job.mo_basename}, {-job.isosurface_value}\n"
            )
            f.write("print(pos_iso)\n")
            f.write("print(neg_iso)\n")
            f.write("set surface_color, blue, pos_iso\n")
            f.write("set surface_color, red, neg_iso\n")
            f.write(f"set transparency, {job.transparency_value}\n")
            f.write(f"set surface_quality, {job.surface_quality}\n")
            f.write(f"set antialias, {job.antialias_value}\n")
            logger.info(f"Wrote PML file: {pml_file}")

    def _job_specific_commands(self, job, command):
        """
        Add MO-specific commands to the PyMOL command string.

        Hides labels, loads the generated PML script, and enables
        ray tracing for higher-quality rendering.

        Args:
            job: PyMOL MO job object.
            command: Existing command string to extend.

        Returns:
            str: Extended command string with MO-specific commands.
        """
        command = self._hide_labels(job, command)
        command = self._call_pml(job, command)
        command = self._add_ray_command(job, command)
        return command

    def _offset_labels(self, job, command, x=-1.2):
        """
        No label offsetting for MO visualization.

        Labels are hidden in MO visualizations, so this method
        returns the command unchanged.

        Args:
            job: PyMOL MO job object (unused).
            command: Current command string.

        Returns:
            str: Unmodified command string.
        """
        return command

    def _call_pml(self, job, command):
        """
        Load the generated PML script for MO visualization.

        Adds a PyMOL command to load the PML script that configures
        MO isosurfaces and styling.

        Args:
            job: PyMOL MO job object.
            command: Current command string to extend.

        Returns:
            str: Extended command string with PML load command.
        """
        pml_file = os.path.join(job.folder, f"{job.mo_basename}.pml")
        command += f"; load {quote_path(pml_file)}"
        return command


class PyMOLSpinJobRunner(PyMOLVisualizationJobRunner):
    """
    PyMOL job runner for spin density visualization.

    Specialized runner for creating spin density visualizations using
    PyMOL. Generates spin density cube files from Gaussian calculations
    and creates PyMOL sessions with isosurface representations.
    """

    JOBTYPES = ["pymol_spin"]

    def _get_gaussian_executable(self, job):
        """
        Get the Gaussian executable for the job.

        Retrieves the Gaussian executable path from the server
        configuration for generating spin density cube files.

        Args:
            job: PyMOL spin job instance.

        Returns:
            str: Path to Gaussian executable folder.
        """
        logger.info(
            f"Obtaining Gaussian executable from server: {self.server.name}"
        )
        gaussian_exe = GaussianExecutable.from_servername(self.server.name)
        gaussian_exe_path = gaussian_exe.executable_folder
        return gaussian_exe_path

    def _prerun(self, job):
        """
        Execute pre-run setup for spin density visualization.

        Performs all necessary setup steps including variable assignment,
        fchk file generation, spin cube file generation, and PML script
        creation before PyMOL execution.

        Args:
            job: PyMOL spin job instance to prepare.
        """
        self._assign_variables(job)
        self._generate_fchk_file(job)
        self._generate_spin_cube_file(job)
        self._write_spin_density_pml(job)

    def _generate_spin_cube_file(self, job):
        """
        Generate the spin density cube file.

        Uses Gaussian's cubegen utility to create a cube file containing
        spin density data from the formatted checkpoint file.

        Args:
            job: PyMOL spin job instance containing grid parameters.
        """
        gaussian_exe = self._get_gaussian_executable(job)

        cubegen_command = f"{gaussian_exe}/cubegen 0 spin {job.label}.fchk {job.job_basename}.cube {job.npts}"
        run_command(cubegen_command)

    def _write_spin_density_pml(self, job):
        """
        Write the .pml file based on the .cube file.

        Creates a PyMOL script file that loads the spin density cube
        file and sets up positive/negative isosurfaces with appropriate
        coloring and transparency settings.

        Args:
            job: PyMOL spin job instance.
        """

        pml_file = os.path.join(job.folder, f"{job.spin_basename}.pml")
        if os.path.exists(pml_file):
            logger.warning(f"PML file {pml_file} already exists. Overwriting.")
        with open(pml_file, "w") as f:
            f.write(f"load {job.spin_basename}.cube\n")
            f.write(
                f"isosurface pos_iso_spin, {job.spin_basename}, {job.isosurface_value}\n"
            )
            f.write(
                f"isosurface neg_iso_spin, {job.spin_basename}, {-job.isosurface_value}\n"
            )
            f.write(
                f"ramp_new ramp, {job.spin_basename}, [{-job.isosurface_value},{job.isosurface_value}], [red, blue]\n"
            )
            f.write("set surface_color, ramp, pos_iso_spin\n")
            f.write("set surface_color, ramp, neg_iso_spin\n")
            f.write(f"set transparency, {job.transparency_value}\n")
            f.write(f"set surface_quality, {job.surface_quality}\n")
            f.write(f"set antialias, {job.antialias_value}\n")
            f.write(f"set ray_trace_mode, {job.ray_trace_mode}\n")
            logger.info(f"Wrote PML file: {pml_file}")

    def _job_specific_commands(self, job, command):
        """
        Add job-specific commands for spin density visualization.

        Configures PyMOL commands specific to spin density visualization
        including label hiding, PML file loading, and ray tracing setup.

        Args:
            job: PyMOL spin job instance.
            command: Current command string to extend.

        Returns:
            str: Extended command string with spin-specific commands.
        """
        command = self._hide_labels(job, command)
        command = self._call_pml(job, command)
        command = self._add_ray_command(job, command)
        return command

    def _offset_labels(self, job, command, x=-1.2):
        """
        Handle label offset for spin density visualization.

        Label offsetting is not needed for spin density visualization
        as labels are typically hidden for this type of visualization.

        Args:
            job: PyMOL spin job instance.
            command: Current command string.

        Returns:
            str: Unchanged command string.
        """
        # not needed for spin density visualization
        return command

    def _call_pml(self, job, command):
        """
        Call the PML file for spin density visualization.

        Adds PyMOL command to load the generated PML script containing
        spin density isosurface setup and visualization settings.

        Args:
            job: PyMOL spin job instance.
            command: Current command string to extend.

        Returns:
            str: Extended command string with PML loading command.
        """
        pml_file = os.path.join(job.folder, f"{job.spin_basename}.pml")
        command += f"; load {quote_path(pml_file)}"
        return command
