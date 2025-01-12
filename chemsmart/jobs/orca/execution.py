import logging
import os

import pyatoms.io.yaml
from pyatoms.io.yaml import YAMLCreatableMixin
from pyatoms.settings.server.executables import Executables
from pyatoms.settings.user import ConfigurationError

logger = logging.getLogger(__name__)

settings_registry = []


class OrcaExecutablesFactory(YAMLCreatableMixin):
    """Callable factory for creating OrcaExecutables.

    Usage:
        factory = OrcaExecutablesFactory.from_yaml('myserver.yaml')
        executable = factory().
    """

    YAML_KEY = "ORCA"

    def __init__(self, servername, **kwargs):
        self.servername = servername
        self.kwargs = kwargs

    def __call__(self):
        return CustomOrcaExecutables(servername=self.servername, **self.kwargs)

    @classmethod
    def _other_yaml_kwargs(cls, filename):
        name = os.path.splitext(os.path.basename(filename))[0]
        return {"servername": name}


class OrcaExecutables(Executables):
    EXECUTABLES = NotImplemented
    SUITABLE_SERVERS = NotImplemented
    ENVIRONMENT_VARIABLES = NotImplemented
    FACTORY: type = OrcaExecutablesFactory

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def get_executable(self):
        exe = os.path.normpath(os.path.expanduser(self.EXECUTABLES))
        if not os.path.exists(exe):
            raise FileNotFoundError(f"Executable: {exe} does not exist")
        return exe

    @classmethod
    def from_yaml(cls, filename):
        factory = OrcaExecutablesFactory.from_yaml(filename)
        return factory()

    @classmethod
    def from_servername(cls, servername, **kwargs):
        try:
            return super().from_servername(servername, **kwargs)
        except pyatoms.io.yaml.MalformedYamlError as e:
            raise ConfigurationError(
                f"\n\nOrca Executable for {servername} are not set. "
            ) from e


class PyatomsOrcaExecutables(OrcaExecutables):
    """For testing orca jobs within pyatoms package."""

    EXECUTABLES = "../../executable/fake/orca"
    SUITABLE_SERVERS = ["inpackage"]
    ENVIRONMENT_VARIABLES = {
        "set": {"OMP_NUM_THREADS": 1, "ORCAFOLDER": "../../executable"}
    }


class LocalOrcaExecutables(OrcaExecutables):
    EXECUTABLES = "~/programs/orca/orca"
    SUITABLE_SERVERS = ["local"]
    ENVIRONMENT_VARIABLES = {
        "set": {"OMP_NUM_THREADS": 1, "ORCAFOLDER": "~/programs/orca"}
    }


class CustomOrcaExecutables(OrcaExecutables):
    """ORCA execution settings for a custom server.

    If created from yaml file, servername will be the name of the yaml file without the suffix.
    """

    SUITABLE_SERVERS = []

    def __init__(self, servername, orcafolder, local_run=False):
        super().__init__(local_run=local_run)
        self.servername = servername
        self.orcafolder = orcafolder

    @property
    def SERVERNAME(self):
        return self.servername

    @property
    def EXECUTABLES(self):
        return os.path.join(self.orcafolder, "orca")

    @property
    def ENVIRONMENT_VARIABLES(self):
        return {"set": {"orca": self.orcafolder, "OMP_NUM_THREADS": 1}}


def execution_settings(cls):
    if cls.SERVERNAME in [p.SERVERNAME for p in settings_registry]:
        raise AssertionError(f"{cls} has repeated server name")
    settings_registry.append(cls)
    return cls


class OrcaExecutionSettings:
    EXECUTABLES = NotImplemented
    SERVERNAME = NotImplemented
    ENVIRONMENT_VARIABLES = NotImplemented
    # SCRATCH_FOLDER = NotImplemented

    def __init__(self, executable=None):
        self._executable = executable

    def get_executable(self, job, shell=False):
        exe = (
            self.EXECUTABLES["shell"]
            if shell
            else self.EXECUTABLES["non-shell"]
        )
        return os.path.normpath(os.path.expanduser(exe))

    @classmethod
    def from_servername(cls, servername, **kwargs):
        for setting in settings_registry:
            if servername == setting.SERVERNAME:
                return setting(**kwargs)
        raise ValueError(f"Settings for {servername} could not be found")
