import logging
import os
import shlex
import subprocess
from contextlib import suppress
from datetime import datetime
from functools import lru_cache
from glob import glob
from shutil import SameFileError, copy

from chemsmart.io.converter import FileConverter
from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.executable import NCIPLOTExecutable

logger = logging.getLogger(__name__)


class NCIPLOTJobRunner(JobRunner):
    """Job runner for NCIPLOT jobs."""

    PROGRAM = "NCIPLOT"
    JOBTYPES = [
        "nciplot",
    ]

    FAKE = False
    SCRATCH = True  # default to using scratch for NCIPLOT Jobs

    def __init__(
        self, server, scratch=None, fake=False, scratch_dir=None, **kwargs
    ):
        """Initialize the NCIPLOTJobRunner."""
        # Use default SCRATCH if scratch is not explicitly set
        if scratch is None:
            scratch = self.SCRATCH
        super().__init__(
            server=server,
            scratch=scratch,
            scratch_dir=scratch_dir,
            fake=fake,
            **kwargs,
        )
        logger.debug(f"Jobrunner server: {self.server}")
        logger.debug(f"Jobrunner scratch: {self.scratch}")
        logger.debug(f"Jobrunner fake mode: {self.fake}")

    @property
    @lru_cache(maxsize=12)
    def executable(self):
        """Executable class object for NCIPLOT."""
        try:
            logger.info(
                f"Obtaining executable from server: {self.server.name}"
            )
            executable = NCIPLOTExecutable.from_servername(
                servername=self.server.name
            )
            logger.debug(f"Executable obtained: {executable}")
            return executable
        except FileNotFoundError as e:
            logger.error(
                f"No server file {self.server} is found: {e}\n"
                f"Available servers are: {NCIPLOTExecutable.available_servers}"
            )
            raise

    def _prerun(self, job):
        """Prepare the job environment before running."""
        self._assign_variables(job)
        self._prepare_files(job)

    def _assign_variables(self, job):
        """Set up file paths for input, output, and error files."""
        if self.scratch and self.scratch_dir:
            logger.debug(
                f"Setting up job in scratch directory: {self.scratch_dir}"
            )
            self._set_up_variables_in_scratch(job)
        else:
            self._set_up_variables_in_job_directory(job)

    def _set_up_variables_in_scratch(self, job):
        """Set up file paths in a scratch directory."""
        scratch_job_dir = os.path.join(self.scratch_dir, job.label)
        if not os.path.exists(scratch_job_dir):
            with suppress(FileExistsError):
                os.makedirs(scratch_job_dir)
        self.running_directory = scratch_job_dir
        logger.debug(f"Running directory: {self.running_directory}")

        job_inputfile = job.label + ".nci"
        scratch_job_inputfile = os.path.join(scratch_job_dir, job_inputfile)
        self.job_inputfile = os.path.abspath(scratch_job_inputfile)
        logger.debug(f"Job input file in scratch: {self.job_inputfile}")

        job_outputfile = job.label + ".nciout"
        scratch_job_outfile = os.path.join(scratch_job_dir, job_outputfile)
        self.job_outputfile = os.path.abspath(scratch_job_outfile)
        logger.debug(f"Job output file in scratch: {self.job_outputfile}")

        job_errfile = job.label + ".ncierr"
        scratch_job_errfile = os.path.join(scratch_job_dir, job_errfile)
        self.job_errfile = os.path.abspath(scratch_job_errfile)
        logger.debug(f"Job error file in scratch: {self.job_errfile}")

    def _set_up_variables_in_job_directory(self, job):
        """Set up file paths in the job's directory."""
        self.running_directory = job.folder
        logger.debug(f"Running directory: {self.running_directory}")
        self.job_inputfile = os.path.abspath(job.inputfile)
        logger.debug(f"Job input file in folder: {self.job_inputfile}")
        self.job_outputfile = os.path.abspath(job.outputfile)
        logger.debug(f"Job output file in folder: {self.job_outputfile}")
        self.job_errfile = os.path.abspath(job.errfile)
        logger.debug(f"Job error file in folder: {self.job_errfile}")

    def _prepare_files(self, job):
        """Prepare input files and write the input for the NCIPLOT job."""
        # Ensure running directory exists
        if job.molecule is not None:
            self._write_xyz_from_pubchem(job)
        else:
            assert (
                job.filenames is not None
            ), "No molecule provided and no filenames specified for NCIPLOT job."
            if not isinstance(job.filenames, (list, tuple)):
                raise TypeError(
                    f"Expected filenames to be a list or tuple, got {type(job.filenames).__name__}"
                )
            if len(job.filenames) == 0:
                raise ValueError(
                    "No filenames provided for NCIPLOT job. Please provide at least one file."
                )
            else:
                if all(
                    f.endswith((".wfn", ".wfx", ".xyz")) for f in job.filenames
                ):
                    self._copy_input_files(job)
                else:
                    self._write_xyz_from_input_files(job)

    def _copy_input_files(self, job):
        """Copy input files to the running directory."""
        for filename in job.filenames:
            if not os.path.exists(filename):
                raise FileNotFoundError(
                    f"File {filename} does not exist for NCIPLOT job."
                )
            with suppress(SameFileError):
                logger.info(
                    f"Copying file {filename} to {self.running_directory}"
                )
                copy(filename, self.running_directory)

    def _write_xyz_from_input_files(self, job):
        for filename in job.filenames:
            if not os.path.exists(filename):
                raise FileNotFoundError(
                    f"File {filename} does not exist for NCIPLOT job."
                )
            try:
                # convert file types to .xyz
                converter = FileConverter(
                    filename=filename,
                    output_filetype="xyz",
                )
                converter.convert_files()
                with suppress(SameFileError):
                    # copy the converted .xyz file to running directory
                    xyz_file = filename.rsplit(".", 1)[0] + ".xyz"
                    xyz_promolecular_file = (
                        filename.rsplit(".", 1)[0] + "_promolecular.xyz"
                    )
                    xyz_promolecular_filepath = os.path.join(
                        self.running_directory,
                        os.path.basename(xyz_promolecular_file),
                    )
                    copy(xyz_file, xyz_promolecular_filepath)
                    logger.info(
                        f"Copied file {xyz_file} to {xyz_promolecular_filepath}"
                    )
            except Exception as e:
                raise ValueError(
                    f"Could not convert file {filename} to .xyz format. Error: {e}\n"
                    f"Unsupported file format for NCIPLOT: {filename}.\n"
                    f"Supported formats are .xyz, .wfn, .wfx"
                )

    def _write_xyz_from_pubchem(self, job):
        """Write the molecule to an XYZ file if it is provided."""
        xyz_filepath = os.path.join(self.running_directory, f"{job.label}.xyz")
        job.molecule.write_xyz(filename=xyz_filepath, mode="w")
        logger.info(f"Wrote molecule to {xyz_filepath}")

    def _write_input(self, job):
        """Write the input file for NCIPLOT job."""
        from chemsmart.jobs.nciplot.writer import NCIPLOTInputWriter

        input_writer = NCIPLOTInputWriter(job=job)
        input_writer.write(target_directory=self.running_directory)

    def _get_command(self, job):
        """Get execution command for NCIPLOT jobs."""
        exe = self._get_executable()
        command = f"{exe} {self.job_inputfile}"  # without output, so direct output to stdout
        return command

    def _create_process(self, job, command, env):
        """Run the NCIPLOT calculation directly."""
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
        """Get executable for NCIPLOT."""
        exe = self.executable.get_executable()
        logger.info(f"NCIPLOT executable: {exe}")
        return exe

    def _postrun(self, job):
        """Handle post-run tasks, such as copying files from scratch and cleanup."""
        if self.scratch:
            # if job was run in scratch, copy files to job folder except files starting with Gau-
            for file in glob(f"{self.running_directory}/*"):
                if not file.endswith(".tmp"):
                    logger.info(
                        f"Copying file {file} from {self.running_directory} to {job.folder}"
                    )
                    with suppress(SameFileError):
                        # Copy file to job folder
                        copy(file, job.folder)

        if job.is_complete():
            # # if job is completed, remove scratch directory and submit_script
            # # and log.info and log.err files
            # if self.scratch:
            #     logger.info(
            #         f"Removing scratch directory: {self.running_directory}."
            #     )
            #     rmtree(self.running_directory)

            self._remove_err_files(job)


class FakeNCIPLOTJobRunner(NCIPLOTJobRunner):
    # creates job runner process
    # combines information about server and program
    FAKE = True

    def __init__(
        self, server, scratch=None, fake=True, scratch_dir=None, **kwargs
    ):
        super().__init__(
            server=server,
            scratch=scratch,
            scratch_dir=scratch_dir,
            fake=fake,
            **kwargs,
        )

    def run(self, job, **kwargs):
        job.label = job.label + "_fake"  # append fake to label
        self._prerun(job=job)
        self._write_input(job=job)
        returncode = FakeNCIPLOT(self.job_inputfile).run()
        self._postrun(job=job)
        return returncode


class FakeNCIPLOT:
    def __init__(self, file_to_run):
        if not os.path.exists(file_to_run):
            raise FileNotFoundError(f"File {file_to_run} not found.")
        self.file_to_run = file_to_run

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
        output_file = self.filename.split(".")[0] + ".nciout"
        return os.path.join(self.file_folder, output_file)

    def run(self):
        with open(self.output_filepath, "w") as g:
            g.write(
                """# ----------------- NCIPLOT ------------------------
 # --- PLOTTING NON COVALENT INTERACTION REGIONS ----
 # ---             E.R. Johnson                  ----
 # ---          J. Contreras-Garcia              ----
 # ----------    Duke University         ------------
 #                                                   
 # ---             A. de la Roza                  ---
 # --------- University of California Merced --------
 #                                                   
 # ---               R. A. Boto                   ---
 # ---                 C. Quan                     --
 # --------  Université Pierre et Marie Curie -------
 # --------------------------------------------------
 # ---              Please cite                  ----
 # --J. Am. Chem. Soc., 2010, 132 (18), pp 6498–6506-
 # --------------------------------------------------
 # --------------------------------------------------
 # ---     Contributions for the wfn properties  ----
 # ---      from H. L. Schmider are acknowledged  ---
 # --------------------------------------------------
 # --------------------------------------------------
 # ---     Contributions for the wfx reader      ----
 # ---      from Dave Arias are acknowledged      ---
 # --------------------------------------------------
 # --------------------------------------------------
 # ---     Contributions for the integration --------
 # ---      algorithms from Erna Wieduwilt    --------
 # ---             are acknowledged           --------
 # ---------------------------------------------------
 #"""
            )
            g.write(f" # Start -- {datetime.now()}\n")
            g.write(
                """-----------------------------------------------------
      INPUT INFORMATION:
-----------------------------------------------------"""
            )
            for line in open(self.input_filepath, "r"):
                g.write(line)
            g.write("\n")
            g.write(f"End -- {datetime.now()}\n")
