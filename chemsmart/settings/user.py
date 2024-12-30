import os
from chemsmart.io.yaml import YAMLFile


class ChemsmartUserSettings:
    USER_YAML_FILE = "usersettings.yaml"
    USER_CONFIG_DIR = os.path.expanduser("~/.chemsmart")

    def __init__(self, yaml, config_dir, **kwargs):
        self.yaml = yaml
        self.config_dir = config_dir
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
