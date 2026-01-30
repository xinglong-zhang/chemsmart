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
    """
    Abstract base class for obtaining program executable paths and configurations.

    This class provides a framework for managing executable configurations for
    different computational chemistry programs. It reads configuration from
    server YAML files and handles environment setup including conda environments,
    modules, scripts, and environment variables.
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
        """
        Initialize the Executable instance.

        Args:
            executable_folder (str, optional): Path to executable directory.
            local_run (bool): Whether to run locally. Defaults to False.
            conda_env (str, optional): Conda environment configuration.
            modules (str, optional): Module loading commands.
            scripts (str, optional): Additional script commands.
            envars (str, optional): Environment variable export commands.
        """
        self.executable_folder = executable_folder
        self.local_run = local_run
        self.conda_env = conda_env
        self.modules = modules
        self.scripts = scripts
        self.envars = envars

    @classmethod
    def from_servername(cls, servername):
        """
        Create an Executable instance from server configuration file.

        Reads configuration from a YAML file in the user's server directory
        and creates an instance with the appropriate settings for the specified
        computational chemistry program.

        Args:
            servername (str): Name of the server configuration file (with or
                            without .yaml extension).

        Returns:
            Executable: An instance configured with server-specific settings.
        """
        # Ensure .yaml extension is present
        if servername.endswith(".yaml"):
            server_yaml = servername
        else:
            server_yaml = f"{servername}.yaml"
        # Load server configuration from YAML file
        server_yaml_file = os.path.join(
            user_settings.user_server_dir, server_yaml
        )
        server_yaml = YAMLFile(filename=server_yaml_file)

        # Extract configuration for the specific program
        exe_folder_raw = server_yaml.yaml_contents_dict[cls.PROGRAM][
            "EXEFOLDER"
        ]
        executable_folder = (
            os.path.expanduser(exe_folder_raw) if exe_folder_raw else None
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

        # Strip comments from configuration strings
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
        """
        Get list of available server configurations.

        Returns:
            list: List of available server configuration names.
        """
        return user_settings.all_available_servers

    @property
    def scratch_dir(self):
        """
        Extract scratch directory path from environment variables.

        Parses the envars configuration to find SCRATCH directory definition.

        Returns:
            str or None: Path to scratch directory if defined, None otherwise.
        """
        if self.envars is not None:
            for line in self.envars.split("\n"):
                line = line.split("#")[0].strip()  # Remove comments
                if "SCRATCH" in line:
                    return line.split("=")[1]
        return None

    @property
    def env(self):
        """
        Parse environment variables from envars configuration.

        Extracts export statements from the envars string and returns them
        as a dictionary of environment variables.

        Returns:
            dict or None: Dictionary of environment variables if envars is set,
                         None otherwise.
        """
        if self.envars is not None:
            env = {}
            for line in self.envars.split("\n"):
                if line.startswith("export"):
                    line = line.split("#")[0].strip()  # Remove comments
                    line = line[7:]  # Remove 'export ' prefix
                    key, value = line.split("=")
                    env[key] = value
            return env
        return None


class GaussianExecutable(Executable):
    """
    Executable handler for Gaussian quantum chemistry software.

    This class provides specific implementation for managing Gaussian 16
    executable paths and configurations.
    """

    PROGRAM = "GAUSSIAN"

    def __init__(self, executable_folder=None, **kwargs):
        """
        Initialize GaussianExecutable instance.

        Args:
            executable_folder (str, optional): Path to Gaussian executable directory.
            **kwargs: Additional arguments passed to parent Executable class.
        """
        super().__init__(executable_folder=executable_folder, **kwargs)

    def get_executable(self):
        """
        Get the full path to the Gaussian executable.

        Returns:
            str or None: Full path to g16 executable if executable_folder is set,
                        None otherwise.
        """
        if self.executable_folder is not None:
            executable_path = os.path.join(self.executable_folder, "g16")
            return executable_path


class ORCAExecutable(Executable):
    """
    Executable handler for ORCA quantum chemistry software.

    This class provides specific implementation for managing ORCA
    executable paths and configurations.
    """

    PROGRAM = "ORCA"

    def __init__(self, executable_folder=None, **kwargs):
        """
        Initialize ORCAExecutable instance.

        Args:
            executable_folder (str, optional): Path to ORCA executable directory.
            **kwargs: Additional arguments passed to parent Executable class.
        """
        super().__init__(executable_folder=executable_folder, **kwargs)

    def get_executable(self):
        """
        Get the full path to the ORCA executable.

        Returns:
            str or None: Full path to orca executable if executable_folder is set,
                        None otherwise.
        """
        if self.executable_folder is not None:
            executable_path = os.path.join(self.executable_folder, "orca")
            return executable_path


class NCIPLOTExecutable(Executable):
    """
    Executable handler for NCIPLOT non-covalent interaction analysis software.

    This class provides specific implementation for managing NCIPLOT
    executable paths and configurations.
    """

    PROGRAM = "NCIPLOT"

    def __init__(self, executable_folder=None, **kwargs):
        """
        Initialize NCIPLOTExecutable instance.

        Args:
            executable_folder (str, optional): Path to NCIPLOT executable directory.
            **kwargs: Additional arguments passed to parent Executable class.
        """
        super().__init__(executable_folder=executable_folder, **kwargs)

    def get_executable(self):
        """
        Get the full path to the NCIPLOT executable.

        Returns:
            str or None: Full path to nciplot executable if executable_folder is set,
                        None otherwise.
        """
        if self.executable_folder is not None:
            executable_path = os.path.join(self.executable_folder, "nciplot")
            return executable_path


class XTBExecutable(Executable):
    """
    Executable handler for XTB semiempirical quantum chemistry software.

    This class provides specific implementation for managing XTB
    executable paths and configurations.
    """

    PROGRAM = "XTB"  # all CAPS as required by 
    # server_yaml.yaml_contents_dict[cls.PROGRAM]["EXEFOLDER"]
    

    def __init__(self, executable_folder=None, **kwargs):
        """
        Initialize XTBExecutable instance.

        Args:
            executable_folder (str, optional): Path to XTB executable directory.
            **kwargs: Additional arguments passed to parent Executable class.
        """
        super().__init__(executable_folder=executable_folder, **kwargs)

    def get_executable(self):
        """
        Get the full path to the XTB executable.

        Returns:
            str: Full path to xtb executable if executable_folder is set,
                 otherwise just "xtb" to use the one from conda environment.
        """
        if self.executable_folder is not None:
            executable_path = os.path.join(self.executable_folder, "xtb")
            return executable_path
        else:
            # Use xtb from conda environment (in PATH)
            return "xtb"
