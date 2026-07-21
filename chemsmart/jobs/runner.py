import copy
import logging
import os
from abc import abstractmethod
from contextlib import suppress
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from shutil import rmtree
from typing import Callable, Optional, Sequence

from chemsmart.jobs.job import Job
from chemsmart.settings.server import Server
from chemsmart.settings.user import ChemsmartUserSettings
from chemsmart.utils.mixins import RegistryMixin

user_settings = ChemsmartUserSettings()


logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class SerialMode:
    """Simple view of serial-mode flags derived from a jobrunner."""

    no_run_in_parallel: bool
    run_in_parallel: bool


@dataclass(frozen=True)
class PhaseTransitionDecision:
    """Decision payload for moving from one workflow phase to the next."""

    proceed: bool
    should_raise: bool
    message: Optional[str] = None


def run_phase_jobs(
    *,
    parent_runner,
    serial_mode: Optional[SerialMode] = None,
    jobs: Optional[Sequence] = None,
    jobs_factory: Optional[Callable[[], Optional[Sequence]]] = None,
    stop_on_incomplete: bool = False,
    before_run: Optional[Callable[[], None]] = None,
    logger_obj=None,
    phase_label: str = "phase",
) -> None:
    """Run a workflow phase of child jobs through ``Job._execute_phase_jobs``.

    Phase siblings are always executed sequentially. The derived ``serial_mode``
    (from ``parent_runner`` when omitted) only controls whether an incomplete
    child stops the phase early when ``stop_on_incomplete`` is true; it does not
    run HA/A (or other intra-phase sub-jobs) concurrently.

    For pKa, intra-molecule phases (gas opt, solvation SP, reference legs, etc.)
    never run in parallel by design. ``--run-in-parallel`` applies to separate
    pKa target jobs wrapped in ``BatchJob``, not to sub-jobs inside one pKa
    thermodynamic cycle.
    """
    if serial_mode is None:
        serial_mode = get_serial_mode(parent_runner)
    Job._execute_phase_jobs(
        parent_runner=parent_runner,
        jobs=jobs,
        jobs_factory=jobs_factory,
        no_run_in_parallel=serial_mode.no_run_in_parallel,
        stop_on_incomplete=stop_on_incomplete,
        before_run=before_run,
        logger_obj=logger_obj,
        phase_label=phase_label,
    )


def get_serial_mode(jobrunner) -> SerialMode:
    """Return serial-mode flags from a jobrunner."""
    if jobrunner is None:
        no_run_in_parallel = False
    else:
        no_run_in_parallel = bool(jobrunner.no_run_in_parallel)
    return SerialMode(
        no_run_in_parallel=no_run_in_parallel,
        run_in_parallel=not no_run_in_parallel,
    )


def _positive_int_or_none(value) -> Optional[int]:
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return None
    if parsed <= 0:
        return None
    return parsed


def get_configured_array_concurrency_limit(
    jobrunner=None,
) -> Optional[int]:
    """Return optional SLURM array throttle from env or server policy.

    Used for ``chemsmart sub --run-in-parallel`` when ``-N`` is not passed.
    Does not fall back to ``num_cores`` (cores per task and array concurrency
    are unrelated).

    Resolution order:
    1. ``CHEMSMART_MAX_SUBMITTERS`` environment variable
    2. ``jobrunner.max_submitters`` (if present)
    3. ``jobrunner.server.max_submitters`` (if present)

    Returns ``None`` when unset (caller may run all array tasks at once).
    """
    env_value = _positive_int_or_none(
        os.environ.get("CHEMSMART_MAX_SUBMITTERS")
    )
    if env_value is not None:
        return env_value

    if jobrunner is None:
        return None

    try:
        runner_value = _positive_int_or_none(jobrunner.max_submitters)
    except AttributeError:
        runner_value = None
    if runner_value is not None:
        return runner_value

    if isinstance(jobrunner, JobRunner):
        try:
            server_value = _positive_int_or_none(
                jobrunner.server.max_submitters
            )
        except AttributeError:
            server_value = None
        if server_value is not None:
            return server_value

    return None


def decide_phase_transition(
    *,
    phase_name: str,
    failures: Optional[Sequence[str]] = None,
    require_complete: bool = False,
    is_complete: Optional[bool] = None,
    stop_message: Optional[str] = None,
) -> PhaseTransitionDecision:
    """Return a shared decision for whether a workflow should enter next phase."""
    phase_failures = [f for f in (failures or []) if f]
    if phase_failures:
        summary = (
            f"{phase_name} phase failed in {len(phase_failures)} worker(s):\n"
            + "\n".join(f"  - {item}" for item in phase_failures)
        )
        return PhaseTransitionDecision(
            proceed=False,
            should_raise=True,
            message=summary,
        )

    if require_complete and is_complete is False:
        message = (
            stop_message or f"{phase_name} jobs incomplete, halting execution."
        )
        return PhaseTransitionDecision(
            proceed=False,
            should_raise=False,
            message=message,
        )

    return PhaseTransitionDecision(proceed=True, should_raise=False)


class JobRunner(RegistryMixin):
    """Abstract base class for job runner for running a job on a server.

    Args:
        server (Server): Server to run the job on.
        scratch (bool): Whether to use scratch directory.
        scratch_dir (str): Path to scratch directory.
        delete_scratch (bool): whether to delete scratch after
            job finishes normally.
        fake (bool): Whether to use fake job runner.
        **kwargs: Additional keyword arguments.
    """

    JOBTYPES: list = NotImplemented
    PROGRAM: str = NotImplemented
    FAKE: bool = False

    def __init__(
        self,
        server,
        scratch=None,
        scratch_dir=None,  # Explicit scratch directory
        delete_scratch=False,
        fake=False,
        no_run_in_parallel=False,
        num_cores=None,
        array_concurrency=None,
        num_nodes=None,
        num_gpus=None,
        mem_gb=None,
        **kwargs,
    ):
        if server is None:
            server = Server.current()

        if isinstance(server, str):
            server = Server.from_servername(server)

        if not isinstance(server, Server):
            raise ValueError(
                f"server must be instance of Server. Instead was: {server}"
            )

        self.server = server
        self.scratch = scratch
        self._scratch_dir = scratch_dir  # Store user-defined scratch_dir
        self.delete_scratch = delete_scratch
        self.no_run_in_parallel = no_run_in_parallel

        if self.scratch:
            self._set_scratch()

        self.fake = fake

        if num_cores is not None:
            self.num_cores = num_cores
        else:
            self.num_cores = self.server.num_cores

        # CLI -N/--array-concurrency only; never overrides server.num_nodes.
        self.cli_array_concurrency = array_concurrency
        if num_nodes is not None:
            self.num_nodes = num_nodes
        else:
            self.num_nodes = self.server.num_nodes

        if num_gpus is not None:
            self.num_gpus = num_gpus
        else:
            self.num_gpus = self.server.num_gpus

        if mem_gb is not None:
            self.mem_gb = mem_gb
        else:
            self.mem_gb = self.server.mem_gb

        self.kwargs = kwargs

    @property
    def scratch_dir(self):
        """Return the scratch directory, setting it if necessary."""
        if self._scratch_dir is None:
            self._scratch_dir = self._set_scratch()
        return self._scratch_dir

    @scratch_dir.setter
    def scratch_dir(self, value):
        """Allow explicit setting of scratch_dir."""
        if value is not None:
            value = os.path.expanduser(value)  # Expand '~' to absolute path
            if not os.path.exists(value):
                raise FileNotFoundError(
                    f"Specified scratch dir does not exist: {value}"
                )
        self._scratch_dir = value

    @lru_cache(maxsize=12)
    def _set_scratch(self):
        """Determine the scratch directory, considering multiple sources."""
        if self._scratch_dir is not None:
            return self._scratch_dir  # Use explicitly set directory

        scratch_dir = None
        if self.executable is not None:
            scratch_dir = self.executable.scratch_dir
            logger.debug(f"Setting scratch dir from executable: {scratch_dir}")
        # (2) then try to get from server specific environment variable
        if scratch_dir is None:
            scratch_dir = self.server.scratch_dir
            logger.debug(
                f"Setting scratch dir from server specific env: {scratch_dir}"
            )

        # (3) then try to get from user settings
        if scratch_dir is None:
            scratch_dir = user_settings.scratch
            logger.debug(
                f"Setting scratch dir from user settings: {scratch_dir}"
            )

        # (4) finally, if scratch_dir is still None, then disable scratch
        if scratch_dir is None:
            logger.warning(
                f"Could not determine scratch dir for {self}. Not using scratch."
            )
            self.scratch = False
        else:
            # check that the scratch folder exists
            scratch_dir = os.path.expanduser(scratch_dir)
            if not os.path.exists(scratch_dir):
                raise FileNotFoundError(
                    f"Specified scratch dir does not exist: {scratch_dir}"
                )
        return scratch_dir

    def __repr__(self):
        return f"{self.__class__.__qualname__}<server={self.server}>"

    @property
    def servername(self):
        return self.server.name

    @property
    def num_hours(self):
        return self.server.num_hours

    @property
    def num_threads(self):
        return self.server.num_threads

    @property
    @abstractmethod
    def executable(self):
        """Subclasses to implement."""
        pass

    def _prerun(self, job):
        # Subclasses can implement
        pass

    def _write_input(self, job):
        # Subclasses can implement
        pass

    def _run(self, process, **kwargs):
        process.communicate()
        return process.poll()

    def _postrun(self, job, **kwargs):
        # Subclasses can implement
        pass

    def _postrun_cleanup(self, job):
        """Perform cleanup tasks after job completion.
        This includes removing error files for successful jobs and
        deleting scratch directories if applicable."""
        if job.is_complete():
            logger.debug("Job completed successfully, deleting .err files")
            self._remove_err_files(job)

            # Delete scratch directory if requested and scratch was used
            if self.scratch and self.delete_scratch:
                logger.debug(
                    "Job completed successfully and delete_scratch is enabled"
                )
                self._delete_scratch_directory()

    @abstractmethod
    def _get_command(self, job):
        raise NotImplementedError

    @abstractmethod
    def _create_process(self, job, command, env):
        raise NotImplementedError

    def _update_os_environ(self, job):
        env = os.environ.copy()
        env_vars = self.executable.env if self.executable else None
        if not env_vars:
            return env
        logger.debug(f"Environment variables to update: \n{env_vars}")
        for k, v in env_vars.items():
            if isinstance(v, str):
                v = os.path.expanduser(v)
            env[k] = str(v)
        return env

    def run(self, job, **kwargs):
        """Main method to run a job. The run consists of
        several steps: prerun, write input, get command,
        create process, run process, postrun, and postrun cleanup.
        prerun and postrun are hooks for subclasses to implement.
        prerun consist of any setup needed before running the job,
        such as creating scratch directories or copying additional
        files into scratch (e.g., in ORCA copying .xyz files).
        postrun consist of e.g., copying files back from scratch to job
        folder (this may be different in different subclasses).
        Args:
            job: Job instance to run.
            **kwargs: Additional keyword arguments for the run method.
        """
        logger.debug(f"Running job {job} with runner {self}")
        logger.debug(f"Prerunning job: {job}")
        self._prerun(job)
        logger.debug(f"Writing input for job: {job}")
        self._write_input(job)
        logger.debug(f"Obtaining command for job: {job}")
        command = self._get_command(job)
        logger.debug(f"Command obtained for job {job}: {command}")
        logger.debug(f"Obtaining environment for job: {job}")
        env = self._update_os_environ(job)
        logger.debug(f"Environment obtained for job {job}: {env}")
        logger.debug(f"Creating process for job: {job}")
        process = self._create_process(job, command=command, env=env)
        logger.debug(f"Process created for job {job}: {process}")
        logger.debug(f"Running process for job: {job}")
        self._run(process, **kwargs)
        logger.debug(f"Postrunning job: {job}")
        self._postrun(job)
        logger.debug(f"Postrun cleanup for job: {job}")
        self._postrun_cleanup(job)

    def copy(self):
        return copy.copy(self)

    @classmethod
    def from_job(cls, job, server, scratch=None, fake=False, **kwargs):
        runners = cls.subclasses()
        logger.debug(f"Available runners: {runners}")
        jobtype = job.TYPE
        candidate_runners = []

        for runner in runners:
            logger.debug(f"Checking runner: {runner} for job: {job}")
            runner_jobtypes = runner.JOBTYPES
            logger.debug(f"Runner jobtypes: {runner_jobtypes}")

            if runner_jobtypes is NotImplemented:
                runner_jobtypes = []

            if jobtype in runner_jobtypes:
                candidate_runners.append(runner)

        if candidate_runners:
            selected_runner = None
            for runner in candidate_runners:
                if runner.FAKE == fake:
                    selected_runner = runner
                    break

            if selected_runner is None:
                selected_runner = candidate_runners[0]

            logger.info(f"Using job runner: {selected_runner} for job: {job}")

            # If scratch is None, use the runner's default scratch value
            if scratch is not None:
                scratch = scratch
            else:
                scratch = selected_runner.SCRATCH
            logger.info(
                f"Using scratch={scratch} for job runner: {selected_runner}"
            )

            return selected_runner(
                server=server, scratch=scratch, fake=fake, **kwargs
            )

        raise ValueError(
            f"Could not find any runners for job: {job}. \n"
            f"Runners in registry: {runners}. \n "
            f"Fake: {fake}"
        )

    def _remove_err_files(self, job):
        """Remove error files associated with the job.
        Files that end with .err, .pbserr, .slurmerr are removed."""

        basefilepath = self._get_base_filepath_to_remove(job)
        patterns = [".err", ".pbserr", ".slurmerr"]

        files_to_be_removed = [
            basefilepath.with_suffix(pattern) for pattern in patterns
        ]

        for file in files_to_be_removed:
            with suppress(FileNotFoundError):
                logger.info(f"Removing file {file}.")
                os.remove(file)

    def _get_base_filepath_to_remove(self, job):
        """Get the base filepath for the job to assist in file removal."""
        return Path(job.folder) / job.label

    def _append_suffix_to_job_label(self, job, suffix):
        """Append ``suffix`` to ``job.label`` once."""
        if suffix and not job.label.endswith(suffix):
            job.label = f"{job.label}{suffix}"
        logger.debug(f"Job label: {job.label}")

    def _delete_scratch_directory(self):
        """
        Delete the scratch directory if it exists.

        This method safely removes the scratch directory and all its contents
        after the job has completed successfully. Only deletes if the
        running_directory is actually within the scratch_dir.
        """
        if (
            hasattr(self, "running_directory")
            and hasattr(self, "scratch_dir")
            and self.scratch_dir
            and os.path.exists(self.running_directory)
        ):

            # Check if running_directory is actually within scratch_dir
            # to avoid accidentally deleting non-scratch directories
            # use resolve() to handle .. and symlinks
            rd = Path(self.running_directory).resolve()
            sd = Path(self.scratch_dir).resolve()

            # Basic sanity checks
            if not sd.exists() or not sd.is_dir():
                logger.error(
                    "scratch_dir %s doesn't exist or is not a directory; "
                    "refusing to proceed.",
                    sd,
                )
            elif rd == sd:
                logger.warning(
                    "Refusing to delete the scratch root itself: %s", sd
                )
            # Python 3.9+: Path.is_relative_to
            elif rd.is_relative_to(sd):
                try:
                    logger.info(
                        f"Deleting scratch directory: {self.running_directory}"
                    )
                    rmtree(self.running_directory)
                    logger.info(
                        f"Successfully deleted scratch directory: {self.running_directory}"
                    )
                except Exception as e:
                    logger.error(
                        f"Failed to delete scratch directory {self.running_directory}: {e}"
                    )
            else:
                logger.debug(
                    f"Running directory {self.running_directory} is not in scratch, "
                    f"skipping deletion."
                )
        else:
            logger.debug(
                "No scratch directory to delete or directory does not exist."
            )
