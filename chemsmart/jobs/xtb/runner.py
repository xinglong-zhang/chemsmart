"""Argv-only local and fake runners for ChemSmart xTB jobs."""

import logging
import os
import subprocess
import tempfile
from contextlib import suppress
from functools import lru_cache
from glob import glob
from pathlib import Path
from shutil import copy

from chemsmart.jobs.runner import JobRunner
from chemsmart.jobs.xtb.validation import (
    XTBEnvironmentError,
    XTBResultValidationError,
    bind_xtb_execution_input,
    build_preview_receipt,
    probe_xtb_environment,
    validate_xtb_result,
    verify_xtb_provenance_binding,
)

logger = logging.getLogger(__name__)


class XTBJobRunner(JobRunner):
    """Execute the strict xTB ``sp``/``opt``/``hess`` command contract."""

    JOBTYPES = ["xtbopt", "xtbsp", "xtbhess"]
    PROGRAM = "xtb"
    FAKE = False
    SCRATCH = True

    _STALE_EXACT_NAMES = frozenset(
        {
            ".xtboptok",
            "charges",
            "energy",
            "g98.out",
            "gradient",
            "hessian",
            "vibspectrum",
            "wbo",
            "xtbopt.log",
            "xtbopt.EIn",
            "xtbopt.coord",
            "xtbopt.gen",
            "xtbopt.pdb",
            "xtbopt.poscar",
            "xtbopt.sdf",
            "xtbopt.xyz",
            "xtbrestart",
            "xtbtopo.mol",
        }
    )
    _STALE_SUFFIXES = (".out", ".err", ".engrad")

    def __init__(
        self,
        server,
        scratch=None,
        fake=False,
        scratch_dir=None,
        num_gpus=None,
        **kwargs,
    ):
        explicit_gpu_request = num_gpus is not None
        self._validate_cpu_only_gpus(num_gpus, source="explicit request")
        if scratch is None:
            scratch = self.SCRATCH
        super().__init__(
            server=server,
            scratch=scratch,
            scratch_dir=scratch_dir,
            fake=fake,
            num_gpus=num_gpus,
            **kwargs,
        )
        self.detected_server_num_gpus = self.num_gpus
        self.xtb_gpu_request_explicit = explicit_gpu_request
        if not explicit_gpu_request:
            # A server may expose GPUs while this program node still requests
            # CPU execution.  Host inventory is not an implicit GPU request.
            self.num_gpus = 0

    @staticmethod
    def _validate_cpu_only_gpus(num_gpus, *, source):
        if num_gpus is None:
            return
        if isinstance(num_gpus, bool) or not isinstance(num_gpus, int):
            raise TypeError(
                f"xTB num_gpus from {source} must be the integer 0."
            )
        if num_gpus != 0:
            raise ValueError(
                "xTB is CPU-only in ChemSmart; "
                f"{source} resolved num_gpus={num_gpus}."
            )

    @property
    @lru_cache(maxsize=12)
    def executable(self):
        # Imported lazily so parser/settings-only use does not require the
        # coordinator-owned executable registry integration to be complete.
        from chemsmart.settings.executable import XTBExecutable

        try:
            return XTBExecutable.from_servername(servername=self.server.name)
        except FileNotFoundError as exc:
            logger.error(
                f"No server file {self.server} is available for xTB: {exc}\n"
                f"Available servers: {XTBExecutable.available_servers}"
            )
            raise

    @classmethod
    def _stale_artifacts(cls, directory):
        root = Path(directory)
        if not root.exists():
            return []
        stale = []
        for path in root.iterdir():
            if not path.is_file():
                continue
            if path.name in cls._STALE_EXACT_NAMES or path.name.endswith(
                cls._STALE_SUFFIXES
            ):
                stale.append(path)
        return sorted(stale)

    @classmethod
    def _reject_stale_artifacts(cls, directory):
        stale = cls._stale_artifacts(directory)
        if stale:
            names = ", ".join(path.name for path in stale)
            raise RuntimeError(
                "Refusing to run xTB in a directory containing stale output "
                f"artifacts: {names}. Archive or remove them explicitly first."
            )

    @staticmethod
    def _private_workspace(parent, *, namespace, prefix):
        """Create one non-reused workspace without following a container link."""

        root = Path(parent)
        root.mkdir(parents=True, exist_ok=True)
        if root.is_symlink():
            raise ValueError(f"Refusing symlink xTB workspace root: {root}")
        container = root / namespace
        container.mkdir(mode=0o700, parents=False, exist_ok=True)
        if container.is_symlink():
            raise ValueError(
                f"Refusing symlink xTB workspace container: {container}"
            )
        workspace = Path(tempfile.mkdtemp(prefix=prefix, dir=container))
        os.chmod(workspace, 0o700)
        return workspace

    def _prepare_job_workspace(self, job, *, preview=False):
        namespace = ".chemsmart-previews" if preview else ".chemsmart-runs"
        prefix = f"{job.label}-preview-" if preview else f"{job.label}-run-"
        workspace = self._private_workspace(
            job.execution_root,
            namespace=namespace,
            prefix=prefix,
        )
        job.bind_workspace(workspace)
        return workspace

    def _prerun(self, job):
        self._prepare_job_workspace(job)
        # This is now an invariant check on a newly allocated directory, not
        # a rejection of outputs preserved from an earlier run.
        self._reject_stale_artifacts(job.folder)
        self._assign_variables(job)
        if os.path.abspath(self.running_directory) != os.path.abspath(
            job.folder
        ):
            self._reject_stale_artifacts(self.running_directory)

    def _assign_variables(self, job):
        os.makedirs(job.folder, exist_ok=True)
        self.job_outputfile = os.path.abspath(job.outputfile)

        if self.scratch and self.scratch_dir:
            self._set_up_variables_in_scratch(job)
        else:
            self._set_up_variables_in_job_directory(job)

        if self.executable.local_run is not None:
            job.local = self.executable.local_run

    def _set_up_variables_in_scratch(self, job):
        scratch_job_dir = self._private_workspace(
            self.scratch_dir,
            namespace=".chemsmart-xtb-runs",
            prefix=f"{job.label}-",
        )
        self.running_directory = str(scratch_job_dir)
        self.job_xyzfile = os.path.abspath(
            os.path.join(scratch_job_dir, f"{job.label}.xyz")
        )
        self.job_errfile = os.path.abspath(job.errfile)

    def _set_up_variables_in_job_directory(self, job):
        self.running_directory = job.folder
        self.job_xyzfile = os.path.abspath(job.xyzfile)
        self.job_errfile = os.path.abspath(job.errfile)

    def _write_input(self, job):
        logger.info(f"Writing xTB geometry input to {self.job_xyzfile}")
        destination = Path(self.job_xyzfile)
        if destination.is_symlink():
            raise ValueError(
                f"Refusing symlink xTB input destination: {self.job_xyzfile}"
            )
        destination_resolved = destination.resolve(strict=False)
        for role in ("source_artifact", "project_artifact"):
            record = job.declared_provenance_binding.get(role)
            if record and destination_resolved == Path(
                record["resolved_path"]
            ):
                raise ValueError(
                    "Refusing to overwrite a provenance-bound xTB source "
                    f"artifact: {record['declared_path']}"
                )
        job.molecule.write_xyz(self.job_xyzfile, mode="w")

    def _get_command(self, job):
        job.settings.validate(
            expected_jobtype=job._TYPE_TO_JOBTYPE.get(job.TYPE)
        )
        executable = self.executable.get_executable()
        if not executable:
            raise FileNotFoundError(
                "No xTB executable was configured or discovered on PATH."
            )
        command = [str(executable), self.job_xyzfile]
        command.extend(self._settings_args(job.settings))
        return command

    def _execution_environment(self, job):
        """Bind xTB and its numerical libraries to the requested CPU count."""

        if (
            isinstance(self.num_cores, bool)
            or not isinstance(self.num_cores, int)
            or self.num_cores < 1
        ):
            raise ValueError("xTB num_cores must be a positive integer.")
        env = self._update_os_environ(job)
        for variable in (
            "OMP_NUM_THREADS",
            "MKL_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
        ):
            env[variable] = str(self.num_cores)
        return env

    @staticmethod
    def _settings_args(settings):
        settings.validate()
        args = XTBJobRunner._gfn_args(settings.gfn_version)

        if settings.jobtype == "opt":
            args.extend(["--opt", settings.optimization_level])
        elif settings.jobtype == "hess":
            args.append("--hess")
        elif settings.jobtype != "sp":
            raise ValueError(f"Unsupported xTB jobtype: {settings.jobtype!r}")

        args.extend(["--chrg", str(settings.charge)])
        args.extend(["--uhf", str(settings.multiplicity - 1)])
        if settings.solvent_model is not None:
            args.extend([f"--{settings.solvent_model}", settings.solvent_id])
        if settings.grad:
            args.append("--grad")
        return args

    @staticmethod
    def _gfn_args(gfn_version):
        if gfn_version in ("gfn0", "gfn1", "gfn2"):
            return ["--gfn", gfn_version[-1]]
        if gfn_version == "gfnff":
            return ["--gfnff"]
        raise ValueError(f"Unsupported xTB gfn_version: {gfn_version!r}")

    def _create_process(self, job, command, env):
        if not isinstance(command, list) or not all(
            isinstance(part, str) for part in command
        ):
            raise TypeError(
                "xTB commands must be a list of string argv items."
            )
        logger.info(
            f"Executing xTB argv {command!r}; output={self.job_outputfile}; "
            f"stderr={self.job_errfile}"
        )
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
        for filepath in glob(os.path.join(self.running_directory, "*")):
            destination = os.path.join(job.folder, os.path.basename(filepath))
            if os.path.abspath(filepath) == os.path.abspath(destination):
                continue
            with suppress(IsADirectoryError):
                copy(filepath, job.folder)

    def _postrun_cleanup(self, job):
        """Preserve stderr evidence; optionally remove validated scratch."""

        if job.is_complete() and self.scratch and self.delete_scratch:
            self._delete_scratch_directory()

    def run(self, job, **kwargs):
        """Execute xTB and require exit, environment, and result evidence."""

        provenance_binding = dict(job.declared_provenance_binding)
        provenance_findings = verify_xtb_provenance_binding(provenance_binding)
        if provenance_findings:
            raise XTBResultValidationError(provenance_findings)
        self._prerun(job=job)
        self._write_input(job=job)
        provenance_binding["execution_input_artifact"] = (
            bind_xtb_execution_input(self.job_xyzfile)
        )
        self.provenance_binding = provenance_binding
        command = self._get_command(job)
        env = self._execution_environment(job)
        environment_receipt = probe_xtb_environment(
            executable=command[0],
            num_cores=self.num_cores,
            num_gpus=self.num_gpus,
            mem_gb=self.mem_gb,
            env=env,
            receipt_path=job.environment_receiptfile,
            settings=job.settings,
        )
        self.environment_receipt = environment_receipt
        self.last_command = tuple(command)
        if environment_receipt.get("execution_ready") is not True:
            raise XTBEnvironmentError(
                environment_receipt.get("findings", ()),
                receipt_path=job.environment_receiptfile,
            )

        process = self._create_process(job, command=command, env=env)
        returncode = None
        try:
            returncode = self._run(process, **kwargs)
        finally:
            self._postrun(job=job)
        self.last_returncode = returncode
        try:
            validation_receipt = validate_xtb_result(
                job=job,
                command=command,
                environment_receipt=environment_receipt,
                provenance_binding=provenance_binding,
                returncode=returncode,
                receipt_path=job.result_receiptfile,
            )
        except Exception as exc:
            if returncode != 0:
                raise subprocess.CalledProcessError(
                    returncode, command
                ) from exc
            raise
        self.validation_receipt = validation_receipt
        job.validation_receipt = validation_receipt
        if returncode != 0:
            raise subprocess.CalledProcessError(returncode, command)
        if validation_receipt.get("ready") is not True:
            raise XTBResultValidationError(
                validation_receipt.get("findings", ()),
                receipt_path=job.result_receiptfile,
                returncode=returncode,
            )
        self._postrun_cleanup(job=job)
        return returncode


class _FakeXTBExecutable:
    """No-dependency executable facade used only by the fake runner."""

    local_run = True
    env = None
    scratch_dir = None

    @staticmethod
    def get_executable():
        return "xtb"


class FakeXTBJobRunner(XTBJobRunner):
    """Render an isolated xTB preview without creating completion artifacts."""

    FAKE = True

    def __init__(self, server, scratch=None, fake=True, **kwargs):
        super().__init__(server=server, scratch=False, fake=fake, **kwargs)

    @property
    @lru_cache(maxsize=12)
    def executable(self):
        return _FakeXTBExecutable()

    def run(self, job, **kwargs):
        preview_dir = self._prepare_job_workspace(job, preview=True)
        self.running_directory = str(preview_dir)
        self.job_xyzfile = os.path.join(preview_dir, f"{job.label}.xyz")
        self.job_outputfile = os.path.join(
            preview_dir, "preview-transcript.txt"
        )
        self.job_errfile = os.path.join(preview_dir, "preview-stderr.txt")
        self._write_input(job=job)
        command = ["xtb", self.job_xyzfile]
        command.extend(self._settings_args(job.settings))
        FakeXTB(
            xyzfile=self.job_xyzfile,
            outputfile=self.job_outputfile,
            errfile=self.job_errfile,
            command=command,
        ).run()
        receipt = build_preview_receipt(
            job=job,
            command=command,
            preview_dir=preview_dir,
            receipt_path=job.preview_receiptfile,
        )
        self.preview_receipt = receipt
        job.preview_receipt = receipt
        self.last_command = tuple(command)
        return 0


class FakeXTB:
    """Write a visibly non-executed preview transcript."""

    def __init__(self, xyzfile, outputfile, errfile, command):
        self.xyzfile = xyzfile
        self.outputfile = outputfile
        self.errfile = errfile
        self.command = command

    def run(self):
        if not os.path.exists(self.xyzfile):
            raise FileNotFoundError(f"File {self.xyzfile} not found.")
        with open(self.outputfile, "w") as out:
            out.write("CHEMSMART XTB PREVIEW ONLY\n")
            out.write("No xTB process or chemistry engine was invoked.\n")
            out.write(f"canonical argv: {self.command!r}\n")
        with open(self.errfile, "w") as err:
            err.write("")
        return 0
