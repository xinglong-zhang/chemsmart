import logging
import os
import subprocess
from contextlib import suppress
from functools import lru_cache
from glob import glob
from shutil import copy

from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.executable import XTBExecutable

logger = logging.getLogger(__name__)


class XTBJobRunner(JobRunner):
    """Job runner for xTB command-line calculations."""

    JOBTYPES = ["xtbopt", "xtbsp", "xtbhess"]
    PROGRAM = "xtb"
    FAKE = False
    SCRATCH = True

    def __init__(
        self, server, scratch=None, fake=False, scratch_dir=None, **kwargs
    ):
        if scratch is None:
            scratch = self.SCRATCH
        super().__init__(
            server=server,
            scratch=scratch,
            scratch_dir=scratch_dir,
            fake=fake,
            **kwargs,
        )
        logger.debug(f"xTB jobrunner server: {self.server}")
        logger.debug(f"xTB jobrunner num cores: {self.num_cores}")
        logger.debug(f"xTB jobrunner num hours: {self.num_hours}")
        logger.debug(f"xTB jobrunner num gpus: {self.num_gpus}")
        logger.debug(f"xTB jobrunner mem gb: {self.mem_gb}")
        logger.debug(f"xTB jobrunner num threads: {self.num_threads}")
        logger.debug(f"xTB jobrunner scratch: {self.scratch}")
        logger.debug(f"xTB jobrunner delete_scratch: {self.delete_scratch}")

    @property
    @lru_cache(maxsize=12)
    def executable(self):
        try:
            logger.info(
                f"Obtaining xTB executable from server: {self.server.name}"
            )
            return XTBExecutable.from_servername(servername=self.server.name)
        except FileNotFoundError as exc:
            logger.error(
                f"No server file {self.server} is found for xTB: {exc}\n"
                f"Available servers: {XTBExecutable.available_servers}"
            )
            raise

    def _prerun(self, job):
        self._assign_variables(job)

    def _assign_variables(self, job):
        os.makedirs(job.folder, exist_ok=True)
        self.job_outputfile = os.path.abspath(job.outputfile)

        if self.scratch and self.scratch_dir:
            self._set_up_variables_in_scratch(job)
        else:
            self._set_up_variables_in_job_directory(job)

        if self.executable and self.executable.local_run is not None:
            logger.info(f"xTB local run is {self.executable.local_run}.")
            job.local = self.executable.local_run

        logger.debug(f"xTB running directory: {self.running_directory}")
        logger.debug(f"xTB geometry input path: {self.job_xyzfile}")
        logger.debug(f"xTB output path: {self.job_outputfile}")
        logger.debug(f"xTB error path: {self.job_errfile}")

    def _set_up_variables_in_scratch(self, job):
        scratch_job_dir = os.path.join(self.scratch_dir, job.label)
        if not os.path.exists(scratch_job_dir):
            with suppress(FileExistsError):
                os.makedirs(scratch_job_dir)
        self.running_directory = scratch_job_dir
        self.job_xyzfile = os.path.abspath(
            os.path.join(scratch_job_dir, f"{job.label}.xyz")
        )
        self.job_errfile = os.path.abspath(job.errfile)

    def _set_up_variables_in_job_directory(self, job):
        self.running_directory = job.folder
        self.job_xyzfile = os.path.abspath(job.xyzfile)
        self.job_errfile = os.path.abspath(job.errfile)

    def _write_input(self, job):
        logger.info(f"Writing xTB geometry input file to: {self.job_xyzfile}")
        job.molecule.write_xyz(self.job_xyzfile, mode="w")

    def _get_command(self, job):
        executable = self.executable.get_executable()
        command = [executable, self.job_xyzfile]
        command.extend(self._settings_args(job.settings))
        logger.debug(f"Generated xTB command: {command}")
        return command

    def _settings_args(self, settings):
        # Keep command rendering centralized so run/sub/fake paths cannot drift
        # when new xTB flags or jobtypes are added.
        args = []
        args.extend(self._gfn_args(settings.gfn_version))

        if settings.jobtype == "opt":
            args.extend(["--opt", settings.optimization_level])
        elif settings.jobtype == "hess":
            args.append("--hess")
        elif settings.jobtype != "sp":
            raise ValueError(f"Unsupported xTB jobtype: {settings.jobtype}")

        args.extend(["--chrg", str(settings.charge)])
        args.extend(["--uhf", str(settings.multiplicity - 1)])

        if (
            settings.solvent_model is not None
            and settings.solvent_id is not None
        ):
            args.extend([f"--{settings.solvent_model}", settings.solvent_id])
        if settings.grad:
            args.append("--grad")
        return args

    @staticmethod
    def _gfn_args(gfn_version):
        if gfn_version is None:
            return []
        if gfn_version.startswith("gfn") and gfn_version[-1].isdigit():
            return ["--gfn", gfn_version[-1]]
        if gfn_version == "gfnff":
            return ["--gfnff"]
        return [f"--{gfn_version}"]

    def _create_process(self, job, command, env):
        logger.info(
            f"Executing xTB command: {' '.join(command)}\n"
            f"Writing output file to: {self.job_outputfile}\n"
            f"Writing err file to: {self.job_errfile}"
        )
        logger.debug(f"xTB run environment updates: {self.executable.env}")
        with (
            open(self.job_outputfile, "w") as out,
            open(self.job_errfile, "w") as err,
        ):
            return subprocess.Popen(
                command,
                stdout=out,
                stderr=err,
                env=env,
                cwd=self.running_directory,
            )

    def _postrun(self, job, **kwargs):
        if not self.scratch:
            return
        logger.debug(
            f"Copying xTB scratch files from: {self.running_directory}"
        )
        for filepath in glob(os.path.join(self.running_directory, "*")):
            destination = os.path.join(job.folder, os.path.basename(filepath))
            if os.path.abspath(filepath) == os.path.abspath(destination):
                continue
            logger.info(f"Copying xTB file {filepath} to {job.folder}")
            with suppress(IsADirectoryError):
                copy(filepath, job.folder)


class FakeXTBJobRunner(XTBJobRunner):
    """Fake xTB runner for CLI, submit, and parser smoke tests."""

    FAKE = True

    def __init__(self, server, scratch=None, fake=True, **kwargs):
        super().__init__(server=server, scratch=scratch, fake=fake, **kwargs)

    @property
    @lru_cache(maxsize=12)
    def executable(self):
        return XTBExecutable(executable_folder=None, local_run=True)

    def _set_up_variables_in_scratch(self, job):
        scratch_job_dir = os.path.join(self.scratch_dir, job.label)
        if not os.path.exists(scratch_job_dir):
            with suppress(FileExistsError):
                os.makedirs(scratch_job_dir)
        self.running_directory = scratch_job_dir
        self._append_suffix_to_job_label(job, "_fake")
        self.job_xyzfile = os.path.abspath(
            os.path.join(scratch_job_dir, f"{job.label}.xyz")
        )
        self.job_outputfile = os.path.abspath(job.outputfile)
        self.job_errfile = os.path.abspath(job.errfile)

    def _set_up_variables_in_job_directory(self, job):
        self.running_directory = job.folder
        self._append_suffix_to_job_label(job, "_fake")
        self.job_xyzfile = os.path.abspath(job.xyzfile)
        self.job_outputfile = os.path.abspath(job.outputfile)
        self.job_errfile = os.path.abspath(job.errfile)

    def run(self, job, **kwargs):
        self._prerun(job=job)
        self._write_input(job=job)
        command = self._get_command(job)
        # Exercise the same setup, command rendering, and parser-completion
        # path as real xTB without requiring an executable on CI/developer hosts.
        FakeXTB(
            xyzfile=self.job_xyzfile,
            outputfile=self.job_outputfile,
            errfile=self.job_errfile,
            command=command,
        ).run()
        self._postrun(job=job)
        self._postrun_cleanup(job=job)
        return 0


class FakeXTB:
    """Small fake xTB execution simulator."""

    def __init__(self, xyzfile, outputfile, errfile, command):
        self.xyzfile = xyzfile
        self.outputfile = outputfile
        self.errfile = errfile
        self.command = command

    def run(self):
        if not os.path.exists(self.xyzfile):
            raise FileNotFoundError(f"File {self.xyzfile} not found.")
        with open(self.outputfile, "w") as out:
            out.write(
                "------------------------------------------------------------\n"
            )
            out.write("* xtb version 0.0.0 (Fake)\n")
            out.write(
                f"program call               : {' '.join(self.command)}\n"
            )
            out.write(
                "------------------------------------------------------------\n"
            )
            out.write("* finished run (fake xtb)\n")
        with open(self.errfile, "w") as err:
            err.write("")
        if "--opt" in self.command:
            # Real xtb writes the relaxed structure to xtbopt.xyz in the run
            # directory; downstream geometry extraction depends on it.
            copy(
                self.xyzfile,
                os.path.join(os.path.dirname(self.outputfile), "xtbopt.xyz"),
            )
        return 0
