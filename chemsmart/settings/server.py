import logging
import os
import subprocess
from chemsmart.utils.mixins import cached_property
from chemsmart.utils.mixins import RegistryMixin
from chemsmart.io.yaml import YAMLFile

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

    def __call__(self):
        server_cls = [s for s in Server.subclasses() if self.name == s.NAME]
        if len(server_cls) == 0:
            raise ValueError(
                f"No server of defined name: {self.name}.\nAvailable servers: {Server.subclasses()}"
            )

        assert len(server_cls) == 1
        server_cls = server_cls[0]
        return server_cls(**self.kwargs)

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
        return self.kwargs.get("SCHEDULER")

    @cached_property
    def queue_name(self):
        return self.kwargs.get("QUEUE_NAME")

    @cached_property
    def num_hours(self):
        return self.kwargs.get("NUM_HOURS")

    @cached_property
    def mem_gb(self):
        return self.kwargs.get("MEM_GB")

    @cached_property
    def num_cores(self):
        return self.kwargs.get("NUM_CORES")

    @cached_property
    def num_gpus(self):
        return self.kwargs.get("NUM_GPUS")

    @cached_property
    def num_threads(self):
        return self.kwargs.get("NUM_THREADS")

    @cached_property
    def execution_type(self):
        return self.kwargs.get("EXECUTION_TYPE")

    @cached_property
    def scratch(self):
        return self.kwargs.get("SCRATCH")

    @cached_property
    def use_hosts(self):
        return self.kwargs.get("USE_HOSTS")

    @cached_property
    def extra_commands(self):
        return self.kwargs.get("EXTRA_COMMANDS")

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
    def from_name(cls, name):
        return cls(name)()

    @classmethod
    def from_scheduler_type(cls):
        """
        Create a Server instance based on the detected scheduler type.

        Returns:
            Server: An instance of the Server class matching the detected scheduler.
        """
        scheduler_type = cls.detect_server_scheduler()

        if scheduler_type == "Unknown Scheduler":
            raise ValueError("Could not detect a known scheduler type.")

        # Match scheduler type with available Server subclasses
        for server_cls in cls.subclasses():
            if (
                hasattr(server_cls, "SCHEDULER_TYPE")
                and server_cls.SCHEDULER_TYPE == scheduler_type
            ):
                return server_cls(name=scheduler_type)()

        raise ValueError(
            f"No server class defined for scheduler type: {scheduler_type}. "
            f"Available servers: {cls.subclasses()}"
        )

    @staticmethod
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
                            return scheduler["name"]
                    else:
                        return scheduler["name"]
                except FileNotFoundError:
                    pass  # Command not found, move to the next scheduler

        # Default case: unknown scheduler
        return "Unknown Scheduler"


class SLURMServer(Server):
    NAME = "SLURM"
    SCHEDULER_TYPE = "SLURM"

    def __init__(self, **kwargs):
        super().__init__(self.NAME, **kwargs)


class PBSServer(Server):
    NAME = "PBS"
    SCHEDULER_TYPE = "PBS"

    def __init__(self, **kwargs):
        super().__init__(self.NAME, **kwargs)


class LSFServer(Server):
    NAME = "LSF"
    SCHEDULER_TYPE = "LSF"

    def __init__(self, **kwargs):
        super().__init__(self.NAME, **kwargs)


class SGE_Server(Server):
    NAME = "SGE"
    SCHEDULER_TYPE = "SGE"

    def __init__(self, **kwargs):
        super().__init__(self.NAME, **kwargs)


