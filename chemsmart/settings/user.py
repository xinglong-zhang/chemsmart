import os
from functools import cached_property
from chemsmart.io.yaml import YAMLFile


class ChemsmartUserSettings:
    USER_YAML_FILE = "usersettings.yaml"
    USER_CONFIG_DIR = os.path.expanduser("~/.chemsmart")

    def __init__(self):
        self.yaml = self.USER_YAML_FILE
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
    def scratch(self):
        return self.data.get("SCRATCH", None)

    @cached_property
    def email(self):
        return self.data.get("EMAIL", None)

    @cached_property
    def project(self):
        return self.data.get("PROJECT", None)
