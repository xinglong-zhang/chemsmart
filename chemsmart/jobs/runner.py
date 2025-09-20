import copy
import logging
import os
from abc import abstractmethod
from contextlib import suppress
from functools import lru_cache
from shutil import rmtree

from chemsmart.settings.server import Server
from chemsmart.settings.user import ChemsmartUserSettings
from chemsmart.utils.mixins import RegistryMixin

user_settings = ChemsmartUserSettings()


logger = logging.getLogger(__name__)


class JobRunner(RegistryMixin):
    """Abstract base class for job runner for running a job on a server.

    Args:
        server (Server): Server to run the job on.
        scratch (bool): Whether to use scratch directory.
        scratch_dir (str): Path to scratch directory.
        delete_scratch: Boolean to delete scratch after
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
        num_cores=None,
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

        if self.scratch:
            self._set_scratch()

        self.fake = fake

        if num_cores is not None:
            self.num_cores = num_cores
        else:
            self.num_cores = self.server.num_cores

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

        for runner in runners:
            logger.debug(f"Checking runner: {runner} for job: {job}")
            runner_jobtypes = runner.JOBTYPES
            logger.debug(f"Runner jobtypes: {runner_jobtypes}")

            if runner_jobtypes is NotImplemented:
                runner_jobtypes = []

            if jobtype in runner_jobtypes:
                logger.info(f"Using job runner: {runner} for job: {job}")

                # If scratch is None, use the runner's default scratch value
                scratch = (
                    scratch
                    if scratch is not None
                    else getattr(runner, "SCRATCH", None)
                )
                logger.info(
                    f"Using scratch={scratch} for job runner: {runner}"
                )

                return runner(server=server, scratch=scratch, **kwargs)

        raise ValueError(
            f"Could not find any runners for job: {job}. \n"
            f"Runners in registry: {runners}. \n "
            f"Fake: {fake}"
        )

    def _remove_err_files(self, job):
        # also remove .err and .pbs* and .slurm* files if job is complete
        err_file = f"{job.folder}/{job.label}.err"
        pbs_err_file = f"{job.folder}/{job.label}.pbserr"
        slurm_err_file = f"{job.folder}/{job.label}.slurmerr"

        files_to_be_removed = [err_file, pbs_err_file, slurm_err_file]
        for file in files_to_be_removed:
            with suppress(FileNotFoundError):
                logger.info(f"Removing file {file}.")
                os.remove(file)

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
            rd = os.path.realpath(self.running_directory)
            sd = os.path.realpath(self.scratch_dir)
            if os.path.commonpath([rd, sd]) == sd:
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
                    f"Running directory {self.running_directory} is not in scratch, skipping deletion"
                )
        else:
            logger.debug(
                "No scratch directory to delete or directory does not exist"
            )
