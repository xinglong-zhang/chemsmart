import logging

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
    def scratch_dir(self):
        return self.kwargs.get("SCRATCH_DIR")

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
