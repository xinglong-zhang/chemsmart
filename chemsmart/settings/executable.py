import logging
import os.path
from typing import Optional

from chemsmart.io.yaml import YAMLFile
from chemsmart.settings.user import ChemsmartUserSettings
from chemsmart.utils.mixins import RegistryMixin
from chemsmart.utils.utils import strip_out_comments

user_settings = ChemsmartUserSettings()

logger = logging.getLogger(__name__)


class Executable(RegistryMixin):
    """Abstract base class for obtaining program executable.
    Given program type, the executable will be specified in server.yaml file.
    """

    PROGRAM: Optional[str] = None

    def __init__(
        self,
        executable_folder=None,
        local_run=False,
        conda_env=None,
        modules=None,
        scripts=None,
        envars=None,
    ):
        self.executable_folder = executable_folder
        self.local_run = local_run
        self.conda_env = conda_env
        self.modules = modules
        self.scripts = scripts
        self.envars = envars

    @classmethod
    def from_servername(cls, servername):
        if servername.endswith(".yaml"):
            server_yaml = servername
        else:
            server_yaml = f"{servername}.yaml"
        server_yaml_file = os.path.join(
            user_settings.user_server_dir, server_yaml
        )
        server_yaml = YAMLFile(filename=server_yaml_file)
        executable_folder = os.path.expanduser(
            server_yaml.yaml_contents_dict[cls.PROGRAM]["EXEFOLDER"]
        )
        local_run = server_yaml.yaml_contents_dict[cls.PROGRAM].get(
            "LOCAL_RUN", False
        )
        conda_env = server_yaml.yaml_contents_dict[cls.PROGRAM].get(
            "CONDA_ENV", None
        )
        modules = server_yaml.yaml_contents_dict[cls.PROGRAM].get(
            "MODULES", None
        )
        scripts = server_yaml.yaml_contents_dict[cls.PROGRAM].get(
            "SCRIPTS", None
        )
        envars = server_yaml.yaml_contents_dict[cls.PROGRAM].get(
            "ENVARS", None
        )
        if conda_env is not None:
            conda_env = strip_out_comments(conda_env)
        if modules is not None:
            modules = strip_out_comments(modules)
        if scripts is not None:
            scripts = strip_out_comments(scripts)
        if envars is not None:
            envars = strip_out_comments(envars)
        return cls(
            executable_folder=executable_folder,
            local_run=local_run,
            conda_env=conda_env,
            modules=modules,
            scripts=scripts,
            envars=envars,
        )

    @property
    def available_servers(self):
        return user_settings.all_available_servers

    @property
    def scratch_dir(self):
        if self.envars is not None:
            for line in self.envars.split("\n"):
                line = line.split("#")[0].strip()
                if "SCRATCH" in line:
                    return line.split("=")[1]
        return None

    @property
    def env(self):
        if self.envars is not None:
            env = {}
            for line in self.envars.split("\n"):
                if line.startswith("export"):
                    line = line.split("#")[0].strip()
                    line = line[7:]
                    key, value = line.split("=")
                    env[key] = value
            return env
        return None


class GaussianExecutable(Executable):
    PROGRAM = "GAUSSIAN"

    def __init__(self, executable_folder=None, **kwargs):
        super().__init__(executable_folder=executable_folder, **kwargs)

    def get_executable(self):
        if self.executable_folder is not None:
            executable_path = os.path.join(self.executable_folder, "g16")
            return executable_path


class ORCAExecutable(Executable):
    PROGRAM = "ORCA"

    def __init__(self, executable_folder=None, **kwargs):
        super().__init__(executable_folder=executable_folder, **kwargs)

    def get_executable(self):
        if self.executable_folder is not None:
            executable_path = os.path.join(self.executable_folder, "orca")
            return executable_path


class NCIPLOTExecutable(Executable):
    PROGRAM = "NCIPLOT"

    def __init__(self, executable_folder=None, **kwargs):
        super().__init__(executable_folder=executable_folder, **kwargs)

    def get_executable(self):
        if self.executable_folder is not None:
            executable_path = os.path.join(self.executable_folder, "nciplot")
            return executable_path
