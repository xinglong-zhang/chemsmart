import logging
import os.path
from abc import abstractmethod
from collections import ChainMap
from itertools import chain

from chemsmart.io.yaml import YAMLFile
from chemsmart.utils.mixins import RegistryMixin

logger = logging.getLogger(__name__)


class Executables(RegistryMixin):
    """Abstract base class for obtaining program executables.
    Given program type, the executable will be specified in server.yaml file."""

    PROGRAM = NotImplemented
    def __init__(self, executable_path=None):
        self.executable_path = executable_path

    def get_executable(self):
        if self.executable_path is not None:
            return os.path.abspath(self.executable_path)

    @classmethod
    def from_yaml(cls, filename):
        executable_yaml = YAMLFile(filename=filename)
        return cls(executable_yaml.yaml_contents_dict[cls.PROGRAM]["EXEFOLDER"])


class GaussianExecutable(Executables):
    PROGRAM = "GAUSSIAN"
    def __init__(self, executable_path=None):
        super().__init__(executable_path=executable_path)


class ORCAExecutable(Executables):
    PROGRAM = "ORCA"
    def __init__(self, executable_path=None):
        super().__init__(executable_path=executable_path)

