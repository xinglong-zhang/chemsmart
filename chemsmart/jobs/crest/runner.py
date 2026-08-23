import logging
import os
import subprocess
from contextlib import suppress
from functools import lru_cache
from glob import glob
from shutil import copy

from chemsmart.jobs.runner import JobRunner
from chemsmart.settings.executable import CRESTExecutable

logger = logging.getLogger(__name__)


class CRESTJobRunner(JobRunner):
    """Job runner for CREST conformational search calculations."""

    JOBTYPES = ["crestconformers"]
    PROGRAM = "crest"
    FAKE = False
    SCRATCH = False

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
        logger.debug(f"CREST Jobrunner server: {self.server}")
        logger.debug(f"CREST Jobrunner num cores: {self.num_cores}")
        logger.debug(f"CREST Jobrunner num hours: {self.num_hours}")
        logger.debug(f"CREST Jobrunner num gpus: {self.num_gpus}")
        logger.debug(f"CREST Jobrunner mem gb: {self.mem_gb}")
        logger.debug(f"CREST Jobrunner num threads: {self.num_threads}")
        logger.debug(f"CREST Jobrunner scratch: {self.scratch}")
        logger.debug(f"CREST Jobrunner delete_scratch: {self.delete_scratch}")

    @property
    @lru_cache(maxsize=12)
    def executable(self):
        try:
            logger.info(
                f"Obtaining CREST executable from server: {self.server.name}"
            )
            return CRESTExecutable.from_servername(servername=self.server.name)
        except FileNotFoundError as exc:
            logger.error(
                f"No server file {self.server} is found for CREST: {exc}\n"
                f"Available servers: {CRESTExecutable.available_servers}"
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
            logger.info(f"CREST local run is {self.executable.local_run}.")
            job.local = self.executable.local_run

        logger.debug(f"CREST running directory: {self.running_directory}")
        logger.debug(f"CREST geometry input path: {self.job_xyzfile}")
        logger.debug(f"CREST output path: {self.job_outputfile}")
        logger.debug(f"CREST error path: {self.job_errfile}")

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
        logger.info(
            f"Writing CREST geometry input file to: {self.job_xyzfile}"
        )
        job.molecule.write_xyz(self.job_xyzfile, mode="w")

        from chemsmart.jobs.crest.writer import CRESTInputWriter

        input_writer = CRESTInputWriter(job=job)
        input_writer.write(target_directory=self.running_directory)

    def _get_command(self, job):
        executable = self.executable.get_executable()
        command = [executable, self.job_xyzfile]
        command.extend(self.get_settings_args(job.settings))
        logger.debug(f"Generated CREST command: {command}")
        return command

    def get_settings_args(self, settings):
        args = []
        args.extend(self._gfn_args(settings.gfn_version))

        if settings.charge is not None:
            args.extend(["--chrg", str(settings.charge)])
        if settings.multiplicity is not None:
            args.extend(["--uhf", str(settings.multiplicity - 1)])

        if (
            settings.solvent_model is not None
            and settings.solvent_id is not None
        ):
            args.extend([f"--{settings.solvent_model}", settings.solvent_id])

        if settings.energy_window is not None:
            args.extend(["--ewin", str(settings.energy_window)])
        if settings.rmsd_threshold is not None:
            args.extend(["--rthr", str(settings.rmsd_threshold)])
        if settings.energy_threshold is not None:
            args.extend(["--ethr", str(settings.energy_threshold)])
        if settings.bconst_threshold is not None:
            args.extend(["--bthr", str(settings.bconst_threshold)])
        if settings.population_threshold is not None:
            args.extend(["--pthr", str(settings.population_threshold)])
        if settings.temperature is not None:
            args.extend(["--temp", str(settings.temperature)])
        if settings.md_timestep is not None:
            args.extend(["--tstep", str(settings.md_timestep)])
        if settings.md_length is not None:
            args.extend(["--mdlen", str(settings.md_length)])
        if settings.md_dump_step is not None:
            args.extend(["--mddump", str(settings.md_dump_step)])
        if settings.vbias_dump_interval is not None:
            args.extend(["--vbdump", str(settings.vbias_dump_interval)])
        if settings.additional_md_temperature is not None:
            args.extend(["--tnmd", str(settings.additional_md_temperature)])
        if settings.no_topology_check:
            args.append("--notopo")
        if settings.no_reference_topology_check:
            args.append("--noreftopo")

        if settings.optimization_level is not None:
            args.extend(["--optlev", settings.optimization_level])

        if settings.quick_mode is not None:
            args.append(f"--{settings.quick_mode}")

        if settings.nci:
            args.append("--nci")

        if settings.nprocs is not None:
            args.extend(["--T", str(settings.nprocs)])

        if settings.constraints:
            args.extend(["--cinp", "constraints.inp"])

        if settings.additional_flags is not None:
            args.extend(settings.additional_flags.split())

        return args

    @staticmethod
    def _gfn_args(gfn_version):
        if gfn_version is None:
            return []
        return [f"--{gfn_version}"]

    def _create_process(self, job, command, env):
        logger.info(
            f"Executing CREST command: {' '.join(command)}\n"
            f"Writing output file to: {self.job_outputfile}\n"
            f"Writing err file to: {self.job_errfile}"
        )
        logger.debug(f"CREST run environment updates: {self.executable.env}")
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
            f"Copying CREST scratch files from: {self.running_directory}"
        )
        for filepath in glob(os.path.join(self.running_directory, "*")):
            destination = os.path.join(job.folder, os.path.basename(filepath))
            if os.path.abspath(filepath) == os.path.abspath(destination):
                continue
            logger.info(f"Copying CREST file {filepath} to {job.folder}")
            with suppress(IsADirectoryError):
                copy(filepath, job.folder)


class FakeCRESTJobRunner(CRESTJobRunner):
    """Fake CREST runner for CLI, submit, and parser smoke tests."""

    FAKE = True

    def __init__(self, server, scratch=None, fake=True, **kwargs):
        super().__init__(server=server, scratch=scratch, fake=fake, **kwargs)

    @property
    @lru_cache(maxsize=12)
    def executable(self):
        return CRESTExecutable(executable_folder=None, local_run=True)

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
        FakeCREST(
            xyzfile=self.job_xyzfile,
            outputfile=self.job_outputfile,
            errfile=self.job_errfile,
            command=command,
        ).run()
        self._postrun(job=job)
        self._postrun_cleanup(job=job)
        return 0


class FakeCREST:
    """Small fake CREST execution simulator."""

    def __init__(self, xyzfile, outputfile, errfile, command):
        self.xyzfile = xyzfile
        self.outputfile = outputfile
        self.errfile = errfile
        self.command = command

    def run(self):
        if not os.path.exists(self.xyzfile):
            raise FileNotFoundError(f"File {self.xyzfile} not found.")
        with open(self.outputfile, "w") as out:
            out.write(" ╔════════════════════════════════════════════╗\n")
            out.write(" ║            ___ ___ ___ ___ _____           ║\n")
            out.write(" ║           / __| _ \ __/ __|_   _|          ║\n")
            out.write(" ║          | (__|   / _|\__ \ | |            ║\n")
            out.write(" ║           \___|_|_\___|___/ |_|            ║\n")
            out.write(" ║                                            ║\n")
            out.write(" ║  Conformer-Rotamer Ensemble Sampling Tool  ║\n")
            out.write(" ║          based on the xTB methods          ║\n")
            out.write(" ║                   (Fake)                   ║\n")
            out.write(" ╚════════════════════════════════════════════╝\n")
            out.write(f" command: {' '.join(self.command)}\n")
            out.write("\n")
            out.write(" CREST terminated normally.\n")
        with open(self.errfile, "w") as err:
            err.write("")
        return 0
