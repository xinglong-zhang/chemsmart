import copy
import logging
import os
from abc import abstractmethod
from chemsmart.utils.mixins import RegistryMixin
from chemsmart.settings.server import Server
from chemsmart.settings.user import ChemsmartUserSettings
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

    def __init__(self, server, scratch=False, scratch_dir=None, fake=False, **kwargs):
        self.server = server
        self.scratch = scratch
        self.scratch_dir = scratch_dir
        self.fake = fake
        self.kwargs = kwargs

        if server is None:
            server = Server.current()
            logger.info(
                f"No server is specified, using current server {server} instead"
            )

        if isinstance(server, str):
            server = Server.from_name(server)

        if not isinstance(server, Server):
            raise ValueError(
                f"server must be instance of Server. Instead was: {server}"
            )

        if self.scratch:
            self._set_scratch()


    def _set_scratch(self):
        scratch_dir = None

        # get scratch directory in order:
        # (1) first try to get from program specific environment variable
        # different programs may need to use different scratch directory
        program_specific_enviornment_vars = os.path.expanduser(
            f"~/.chemsmart/{self.PROGRAM}/{self.PROGRAM}.envars"
        )
        if os.path.exists(program_specific_enviornment_vars):
            # extract any scratch export statement
            try:
                with open(program_specific_enviornment_vars) as f:
                    for line in f.readlines():
                        if "SCRATCH" in line:
                            scratch_dir = line.split("=")[-1]
            except OSError as e:
                raise RuntimeError(
                    f"Failed to read environment variable file: {program_specific_enviornment_vars}"
                ) from e

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
        return f"{self.__class__.__qualname__}<server={self.server}, job={self.job}>"

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
    def executables(self):
        """Subclasses to implement. Return None if no executables."""
        raise NotImplementedError

    def _prerun(self):
        # Subclasses can implement
        pass

    @abstractmethod
    def _run(self, **kwargs):
        raise NotImplementedError

    def _postrun(self):
        # Subclasses can implement
        pass

    @abstractmethod
    def get_command(self, **kwargs):
        raise NotImplementedError

    @abstractmethod
    def _create_process(self, command, env):
        raise NotImplementedError

    @abstractmethod
    def _get_command(self):
        raise NotImplementedError

    def run(self, **kwargs):
        self._prerun()
        self._get_command()
        self._create_process()
        self._run(**kwargs)
        self._postrun()

    def copy(self):
        return copy.copy(self)
