import copy
import logging
import os
from abc import abstractmethod
from functools import lru_cache

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
    @lru_cache(maxsize=12)
    def scratch_dir(self):
        return self._set_scratch()

    @lru_cache(maxsize=12)
    def _set_scratch(self):
        scratch_dir = None
        if self.executable is not None:
            logger.info(f"Setting scratch dir for {self} from executable.")
            scratch_dir = self.executable.scratch_dir
            logger.info(f"Scratch dir set to: {scratch_dir}")
        # (2) then try to get from server specific environment variable
        if scratch_dir is None:
            scratch_dir = self.server.scratch_dir

        # (3) then try to get from user settings
        if scratch_dir is None:
            scratch_dir = user_settings.scratch

        if scratch_dir is not None:
            # check that the scratch folder exists
            scratch_dir = os.path.expanduser(scratch_dir)
            if not os.path.exists(scratch_dir):
                raise FileNotFoundError(
                    f"Specified scratch dir does not exist: {scratch_dir}"
                )
            return scratch_dir
        else:
            return None

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

    def _postrun(self, job):
        # Subclasses can implement
        pass

    @abstractmethod
    def _get_command(self, **kwargs):
        raise NotImplementedError

    @abstractmethod
    def _create_process(self, job, command, env):
        raise NotImplementedError

    def _update_os_environ(self, job):
        env = os.environ.copy()
        env_vars = self.executable.env
        if not env_vars:
            return env
        logger.debug(f"Environment variables to update: \n{env_vars}")
        for k, v in env_vars.items():
            if isinstance(v, str):
                v = os.path.expanduser(v)
            env[k] = str(v)
        return env

    def run(self, job, **kwargs):
        self._prerun(job)
        self._write_input(job)
        command = self._get_command(**kwargs)
        env = self._update_os_environ(job)
        process = self._create_process(job, command=command, env=env)
        self._run(process, **kwargs)
        self._postrun(job)

    def copy(self):
        return copy.copy(self)

    @classmethod
    def from_job(cls, job, server, scratch=None, fake=False, **kwargs):
        runners = cls.subclasses()
        logger.debug(f"All available runners: {runners}")
        jobtype = job.TYPE
        logger.debug(f"Running job type: {jobtype}")

        for runner in runners:
            logger.debug(f"Checking runner: {runner}")
            runner_jobtypes = runner.JOBTYPES
            logger.debug(f"Runner jobtypes: {runner_jobtypes}")

            if runner_jobtypes is NotImplemented:
                runner_jobtypes = []

            if jobtype in runner_jobtypes and fake == runner.FAKE:
                logger.info(f"Using job runner: {runner} for job: {job}")
                return runner(server=server, scratch=scratch, **kwargs)

        raise ValueError(
            f"Could not find any runners for job: {job}. \n"
            f"Runners in registry: {runners}. \n "
            f"Fake: {fake}"
        )
