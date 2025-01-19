import os
import logging
from functools import cached_property
import glob
from chemsmart.io.yaml import YAMLFile

logger = logging.getLogger(__name__)


class ChemsmartUserSettings:
    USER_YAML_FILE = "usersettings.yaml"
    USER_CONFIG_DIR = os.path.expanduser("~/.chemsmart")

    def __init__(self):
        self.yaml = os.path.join(self.USER_CONFIG_DIR, self.USER_YAML_FILE)
        self.config_dir = self.USER_CONFIG_DIR
        try:
            self.data = YAMLFile(filename=self.yaml).yaml_contents_dict
        except FileNotFoundError:
            self.data = {}

    @property
    def user_server_dir(self):
        return os.path.join(self.config_dir, "server")

    @property
    def user_gaussian_settings_dir(self):
        return os.path.join(self.config_dir, "gaussian")

    @property
    def user_gaussian_envars(self):
        return os.path.join(self.user_gaussian_settings_dir, "gaussian.envars")

    @property
    def user_gaussian_modules(self):
        return os.path.join(
            self.user_gaussian_settings_dir, "gaussian.modules"
        )

    @property
    def user_gaussian_script(self):
        return os.path.join(self.user_gaussian_settings_dir, "gaussian.sh")

    @property
    def user_orca_envars(self):
        return os.path.join(self.user_orca_settings_dir, "orca.envars")

    @property
    def user_orca_modules(self):
        return os.path.join(self.user_orca_settings_dir, "orca.modules")

    @property
    def user_orca_script(self):
        return os.path.join(self.user_orca_settings_dir, "orca.sh")

    @property
    def user_orca_settings_dir(self):
        return os.path.join(self.config_dir, "orca")

    @cached_property
    def server_yaml_files(self):
        return glob.glob(os.path.join(self.user_server_dir, "*.yaml"))

    @cached_property
    def gaussian_project_yaml_files(self):
        return glob.glob(
            os.path.join(self.user_gaussian_settings_dir, "*.yaml")
        )

    @cached_property
    def orca_project_yaml_files(self):
        return glob.glob(os.path.join(self.user_orca_settings_dir, "*.yaml"))

    @cached_property
    def scratch(self):
        return self.data.get("SCRATCH", None)

    @cached_property
    def email(self):
        return self.data.get("EMAIL", None)

    @cached_property
    def project(self):
        return self.data.get("PROJECT", None)

    @cached_property
    def all_available_servers(self):
        return [
            os.path.basename(s).removesuffix(".yaml")
            for s in self.server_yaml_files
        ]

    @cached_property
    def all_available_gaussian_projects(self):
        gaussian_project_settings = [
            os.path.basename(g).removesuffix(".yaml")  # python 3.9+
            # os.path.basename(g).strip(".yaml") does not work since
            # m062x.yaml gets stripped to 062x
            for g in self.gaussian_project_yaml_files
        ]
        return gaussian_project_settings

    @cached_property
    def all_available_orca_projects(self):
        return [
            os.path.basename(o).removesuffix(".yaml")
            for o in self.orca_project_yaml_files
        ]
