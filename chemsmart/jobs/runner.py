import copy
import logging
import os
from abc import abstractmethod
from chemsmart.utils.mixins import RegistryMixin
from chemsmart.settings.server import Server
from chemsmart.settings.user import ChemsmartUserSettings
from debugpy.launcher.debuggee import process

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

    def __init__(
        self, server, scratch=False, scratch_dir=None, fake=False, **kwargs
    ):
        if server is None:
            server = Server.current()

        if isinstance(server, str):
            server = Server.from_server(server)

        if not isinstance(server, Server):
            raise ValueError(
                f"server must be instance of Server. Instead was: {server}"
            )

        self.server = server
        self.scratch = scratch
        self.scratch_dir = scratch_dir
        self.fake = fake
        self.kwargs = kwargs

        if self.scratch:
            self._set_scratch()

    def _set_scratch(self):
        scratch_dir = self.executable.scratch_dir
        # (2) then try to get from server specific environment variable
        if scratch_dir is None:
            scratch_dir = self.server.scratch

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
        else:
            raise ValueError("No valid scratch directory could be determined.")

        self.scratch_dir = scratch_dir

    def __repr__(self):
        return f"{self.__class__.__qualname__}<server={self.server}>"

    @property
    def servername(self):
        return self.server.name

    @property
    def num_hours(self):
        return self.server.num_hours

    @property
    def mem_gb(self):
        return self.server.mem_gb

    @property
    def num_cores(self):
        return self.server.num_cores

    @property
    def num_gpus(self):
        return self.server.num_gpus

    @property
    def num_threads(self):
        return self.server.num_threads

    @property
    @abstractmethod
    def executable(self):
        """Subclasses to implement. Return None if no executable."""
        raise NotImplementedError

    def _prerun(self, job):
        # Subclasses can implement
        pass

    @abstractmethod
    def _run(self, job, process, **kwargs):
        raise NotImplementedError

    def _postrun(self, job):
        # Subclasses can implement
        pass

    @abstractmethod
    def _get_command(self):
        command = self.server
        raise NotImplementedError

    @abstractmethod
    def _create_process(self, job, command, env):
        raise NotImplementedError

    def run(self, job, **kwargs):
        self._prerun(job)
        command = self._get_command()
        process = self._create_process(job, command=command, env=self.executable.env)
        # self._create_jobrunner(job, **kwargs)
        self._run(job, process, **kwargs)
        self._postrun(job)

    def copy(self):
        return copy.copy(self)

    @classmethod
    def from_jobtype(cls, job, server, scratch=False, scratch_dir=None, fake=False, **kwargs):
        runners = cls.subclasses()
        jobtype = job.TYPE

        for runner in runners:
            runner_jobtypes = runner.JOBTYPES

            if runner_jobtypes is NotImplemented:
                runner_jobtypes = []

            if jobtype in runner_jobtypes:
                return runner(server=server, scratch=scratch, **kwargs)

        raise ValueError(
            f'Could not find any runners for job: {job}. \n'
            f'Runners in registry: {cls.subclasses()}. \n '
            f'Fake: {fake}'
        )