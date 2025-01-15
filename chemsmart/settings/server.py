import logging
import os
import subprocess
from functools import lru_cache
from chemsmart.utils.mixins import cached_property
from chemsmart.utils.mixins import RegistryMixin
from chemsmart.io.yaml import YAMLFile

from chemsmart.settings.user import ChemsmartUserSettings

user_settings = ChemsmartUserSettings()

logger = logging.getLogger(__name__)


class Server(RegistryMixin):
    def __init__(self, name, **kwargs):
        self.name = name
        self.kwargs = kwargs

    def __str__(self):
        return f"Server: {self.name}"

    def __eq__(self, other):
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)

    def __repr__(self):
        return f"Server(name={self.name})"

    @classmethod
    def from_dict(cls, d):
        return cls(**d)

    @classmethod
    def from_yaml(cls, name):
        if not name:
            raise ValueError("No yaml file provided.")
        server_yaml = YAMLFile(filename=name)
        return cls(name, **server_yaml.yaml_contents_dict["SERVER"])

    @cached_property
    def scheduler(self):
        return self.kwargs.get("SCHEDULER", None)

    @cached_property
    def queue_name(self):
        return self.kwargs.get("QUEUE_NAME", None)

    @cached_property
    def num_hours(self):
        return self.kwargs.get("NUM_HOURS", None)

    @cached_property
    def mem_gb(self):
        return self.kwargs.get("MEM_GB", 64)

    @cached_property
    def num_cores(self):
        return self.kwargs.get("NUM_CORES", 16)

    @cached_property
    def num_gpus(self):
        return self.kwargs.get("NUM_GPUS", 0)

    @cached_property
    def num_threads(self):
        return self.kwargs.get("NUM_THREADS", 16)

    @cached_property
    def submit_command(self):
        return self.kwargs.get("SUBMIT_COMMAND", None)

    @cached_property
    def scratch_dir(self):
        return self.kwargs.get("SCRATCH_DIR", None)

    @cached_property
    def scratch(self):
        return self.scratch_dir is not None

    @cached_property
    def use_hosts(self):
        return self.kwargs.get("USE_HOSTS", None)

    @cached_property
    def extra_commands(self):
        return self.kwargs.get("EXTRA_COMMANDS", None)

    def register(self):
        # if server already in registry, pass
        if self in Server._REGISTRY:
            return self
        Server._REGISTRY.append(self)
        return self

    @classmethod
    def current(cls):
        return cls.from_scheduler_type()

    @classmethod
    @lru_cache(maxsize=12)
    def from_scheduler_type(cls):
        """
        Create a Server instance based on the detected scheduler type.

        Returns:
            Server: An instance of the Server class matching the detected scheduler.
        """
        scheduler_type = cls.detect_server_scheduler()

        if scheduler_type == "Unknown Scheduler":
            logger.info("No scheduler detected. Using local server.")
            return cls.from_servename(servername="local")

        # Match scheduler type with available Server subclasses
        for server_cls in cls.subclasses():
            if getattr(server_cls, "SCHEDULER_TYPE", None) == scheduler_type:
                return server_cls.from_servename(scheduler_type)

        raise ValueError(
            f"No server class defined for scheduler type: {scheduler_type}. "
            f"Available servers: {cls.subclasses()}"
        )

    @staticmethod
    @lru_cache(maxsize=12)
    def detect_server_scheduler():
        """
        Detect the server's job scheduler system.

        Returns:
            str: The detected scheduler type (e.g., SLURM, PBS, LSF, SGE, HTCondor).
        """
        schedulers = [
            {
                "name": "SLURM",
                "env_vars": ["SLURM_JOB_ID", "SLURM_CLUSTER_NAME"],
                "commands": [["squeue"]],
            },
            {
                "name": "PBS",
                "env_vars": ["PBS_JOBID", "PBS_QUEUE"],
                "commands": [["qstat"]],
            },
            {
                "name": "LSF",
                "env_vars": ["LSB_JOBID", "LSB_MCPU_HOSTS"],
                "commands": [["bjobs"]],
            },
            {
                "name": "SGE",
                "env_vars": [],
                "commands": [["qstat"], ["qstat", "-help"]],
                "check_output": lambda output: "Grid Engine" in output
                or "Sun" in output,
            },
            {
                "name": "HTCondor",
                "env_vars": [],
                "commands": [["condor_q"]],
            },
        ]

        for scheduler in schedulers:
            # Check environment variables
            if any(env in os.environ for env in scheduler.get("env_vars", [])):
                logger.info(f"Detected scheduler: {scheduler['name']}")
                return scheduler["name"]

            # Check commands
            for command in scheduler.get("commands", []):
                try:
                    result = subprocess.run(
                        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                    )
                    if "check_output" in scheduler:
                        output = result.stdout.decode()
                        if scheduler["check_output"](output):
                            logger.info(
                                f"Detected scheduler: {scheduler['name']}"
                            )
                            return scheduler["name"]
                    else:
                        logger.info(f"Detected scheduler: {scheduler['name']}")
                        return scheduler["name"]
                except FileNotFoundError:
                    pass  # Command not found, move to the next scheduler

        # Default case: unknown scheduler
        logger.info("No scheduler detected.")
        return "Unknown Scheduler"

    @classmethod
    def from_servename(cls, servername):
        """Obtain server details from server name."""
        if servername is None:
            # by default return current server
            return cls.current()
        return cls._from_server_name(servername)

    @classmethod
    def _from_server_name(cls, server_name):
        """Get server settings from user directory .yaml file based on server name."""
        if server_name.endswith(".yaml"):
            server_name = server_name
        else:
            server_name = f"{server_name}.yaml"
        server_name_yaml_path = os.path.join(
            user_settings.user_server_dir, server_name
        )
        user_settings_manager = ServerSettingsManager(
            filename=server_name_yaml_path
        )
        server = cls._from_servers_manager(user_settings_manager)

        if server is not None:
            return server

        # could not find server settings
        templates_path = os.path.join(os.path.dirname(__file__), "templates")
        raise ValueError(
            f"No server implemented for {server_name}.\n\n"
            f"Place new server .yaml file in {user_settings.user_server_dir}.\n\n"
            f"Templates for server settings .yaml files are available at {templates_path}\n\n "
            f"Currently available servers: {user_settings.all_available_servers}"
        )

    @classmethod
    def _from_servers_manager(cls, manager):
        try:
            return manager.create()
        except FileNotFoundError:
            return None


class YamlServerSettings(Server):
    NAME = "yaml"

    def __init__(self, name, **kwargs):
        super().__init__(name, **kwargs)

    @classmethod
    def from_yaml(cls, filename):
        yaml_file = YAMLFile(filename)
        return cls(name=filename, **yaml_file.yaml_contents_dict["SERVER"])

    def __repr__(self):
        return f"YamlServerSettings(name={self.name})"

    def __str__(self):
        return f"YamlServerSettings: {self.name}"

    def __eq__(self, other):
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)

    def __call__(self):
        return self

    def register(self):
        return self


class ServerSettingsManager:
    """Manages server settings specified in the form of yaml files in a folder.

    Args:
        filename: yaml filename in the default servers folder.
    """

    def __init__(self, filename):
        if filename is None:
            raise ValueError("filename is not specified")
        self.filename = os.path.abspath(filename)

    def create(self):
        return YamlServerSettings.from_yaml(self.filename)


class SLURMServer(YamlServerSettings):
    NAME = "SLURM"
    SCHEDULER_TYPE = "SLURM"

    def __init__(self, **kwargs):
        super().__init__(filename=f"{self.NAME}.yaml", **kwargs)


class PBSServer(YamlServerSettings):
    NAME = "PBS"
    SCHEDULER_TYPE = "PBS"

    def __init__(self, **kwargs):
        super().__init__(self.NAME, **kwargs)


class LSFServer(YamlServerSettings):
    NAME = "LSF"
    SCHEDULER_TYPE = "LSF"

    def __init__(self, **kwargs):
        super().__init__(self.NAME, **kwargs)


class SGE_Server(YamlServerSettings):
    NAME = "SGE"
    SCHEDULER_TYPE = "SGE"

    def __init__(self, **kwargs):
        super().__init__(self.NAME, **kwargs)
