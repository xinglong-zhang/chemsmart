import os
from functools import cached_property
import glob

from chemsmart.io.yaml import YAMLFile


class ChemsmartUserSettings:
    USER_YAML_FILE = "usersettings.yaml"
    USER_CONFIG_DIR = os.path.expanduser("~/.chemsmart")

    def __init__(self):
        self.yaml = os.path.join(self.USER_CONFIG_DIR, self.USER_YAML_FILE)
        self.config_dir = self.USER_CONFIG_DIR
        self.data = YAMLFile(filename=self.yaml).yaml_contents_dict

    @property
    def user_server_dir(self):
        return os.path.join(self.config_dir, "server")

    @property
    def user_gaussian_settings_dir(self):
        return os.path.join(self.config_dir, "gaussian")

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
            os.path.basename(s).strip(".yaml") for s in self.server_yaml_files
        ]

    @cached_property
    def all_available_gaussian_projects(self):
        return [
            os.path.basename(g).strip(".yaml")
            for g in self.gaussian_project_yaml_files
        ]

    @cached_property
    def all_available_orca_projects(self):
        return [
            os.path.basename(o).strip(".yaml")
            for o in self.orca_project_yaml_files
        ]

    def _yaml_file(self, project):
        file = os.path.join(self.path, f"{project}.yaml")
        if not os.path.exists(file):
            raise FileNotFoundError(f"Could not find file: {file}")
        return file
