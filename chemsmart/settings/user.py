"""
User configuration management for ChemSmart computational chemistry software.

This module provides comprehensive user settings management including configuration
file handling, directory structure management, and access to user-specific settings
for computational chemistry software (Gaussian, ORCA) and server configurations.
Manages the ~/.chemsmart user configuration directory and associated files.

Classes:
    ChemsmartUserSettings: Main class for user configuration management

Dependencies:
    - chemsmart.io.yaml: YAML file handling utilities

Configuration Structure:
    ~/.chemsmart/
    ├── usersettings.yaml
    ├── server/
    ├── gaussian/
    │   ├── gaussian.envars
    │   ├── gaussian.modules
    │   └── gaussian.sh
    └── orca/
        ├── orca.envars
        └── orca.modules
"""

import glob
import logging
import os
from functools import cached_property

from chemsmart.io.yaml import YAMLFile

logger = logging.getLogger(__name__)


class ChemsmartUserSettings:
    """
    User configuration settings manager for ChemSmart.
    
    Manages user-specific configuration files, directories, and settings for
    computational chemistry software. Provides access to configuration paths,
    environment variables, and user-defined project settings.
    
    Attributes:
        USER_YAML_FILE (str): Name of the main user settings YAML file.
        USER_CONFIG_DIR (str): Path to the user configuration directory.
        yaml (str): Full path to the user settings YAML file.
        config_dir (str): User configuration directory path.
        data (dict): Loaded YAML configuration data.
    """
    
    USER_YAML_FILE = "usersettings.yaml"
    USER_CONFIG_DIR = os.path.expanduser("~/.chemsmart")

    def __init__(self):
        """
        Initialize user settings manager.
        
        Loads user configuration from YAML file if it exists, otherwise
        initializes with empty configuration.
        """
        self.yaml = os.path.join(self.USER_CONFIG_DIR, self.USER_YAML_FILE)
        self.config_dir = self.USER_CONFIG_DIR
        try:
            self.data = YAMLFile(filename=self.yaml).yaml_contents_dict
        except FileNotFoundError:
            self.data = {}

    @property
    def user_server_dir(self):
        """
        Get the user server configurations directory.
        
        Returns:
            str: Path to the directory containing user server configurations.
        """
        return os.path.join(self.config_dir, "server")

    @property
    def user_gaussian_settings_dir(self):
        """
        Get the user Gaussian settings directory.
        
        Returns:
            str: Path to the directory containing Gaussian-specific settings.
        """
        return os.path.join(self.config_dir, "gaussian")

    @property
    def user_gaussian_envars(self):
        """
        Get the path to Gaussian environment variables file.
        
        Returns:
            str: Path to the file containing Gaussian environment variables.
        """
        return os.path.join(self.user_gaussian_settings_dir, "gaussian.envars")

    @property
    def user_gaussian_modules(self):
        """
        Get the path to Gaussian modules file.
        
        Returns:
            str: Path to the file containing Gaussian module loading commands.
        """
        return os.path.join(
            self.user_gaussian_settings_dir, "gaussian.modules"
        )

    @property
    def user_gaussian_script(self):
        """
        Get the path to Gaussian execution script.
        
        Returns:
            str: Path to the Gaussian execution script file.
        """
        return os.path.join(self.user_gaussian_settings_dir, "gaussian.sh")

    @property
    def user_orca_envars(self):
        """
        Get the path to ORCA environment variables file.
        
        Returns:
            str: Path to the file containing ORCA environment variables.
        """
        return os.path.join(self.user_orca_settings_dir, "orca.envars")

    @property
    def user_orca_modules(self):
        """
        Get the path to ORCA modules file.
        
        Returns:
            str: Path to the file containing ORCA module loading commands.
        """
        return os.path.join(self.user_orca_settings_dir, "orca.modules")

    @property
    def user_orca_script(self):
        """
        Get the path to ORCA execution script.
        
        Returns:
            str: Path to the ORCA execution script file.
        """
        return os.path.join(self.user_orca_settings_dir, "orca.sh")

    @property
    def user_orca_settings_dir(self):
        """
        Get the user ORCA settings directory.
        
        Returns:
            str: Path to the directory containing ORCA-specific settings.
        """
        return os.path.join(self.config_dir, "orca")

    @cached_property
    def server_yaml_files(self):
        """
        Get list of server YAML configuration files.
        
        Returns:
            list: List of paths to server configuration YAML files.
        """
        return glob.glob(os.path.join(self.user_server_dir, "*.yaml"))

    @cached_property
    def gaussian_project_yaml_files(self):
        """
        Get list of Gaussian project YAML configuration files.
        
        Returns:
            list: List of paths to Gaussian project configuration YAML files.
        """
        return glob.glob(
            os.path.join(self.user_gaussian_settings_dir, "*.yaml")
        )

    @cached_property
    def orca_project_yaml_files(self):
        """
        Get list of ORCA project YAML configuration files.
        
        Returns:
            list: List of paths to ORCA project configuration YAML files.
        """
        return glob.glob(os.path.join(self.user_orca_settings_dir, "*.yaml"))

    @cached_property
    def scratch(self):
        """
        Get scratch directory configuration.
        
        Returns:
            str or None: Path to scratch directory or None if not configured.
        """
        return self.data.get("SCRATCH", None)

    @cached_property
    def email(self):
        """
        Get user email configuration.
        
        Returns:
            str or None: User email address or None if not configured.
        """
        return self.data.get("EMAIL", None)

    @cached_property
    def project(self):
        """
        Get default project configuration.
        
        Returns:
            str or None: Default project name or None if not configured.
        """
        return self.data.get("PROJECT", None)

    @cached_property
    def all_available_servers(self):
        """
        Get list of all available server configurations.
        
        Returns:
            list: List of server names (without .yaml extension) available
                  in the user server directory.
        """
        return [
            os.path.basename(s).removesuffix(".yaml")
            for s in self.server_yaml_files
        ]

    @cached_property
    def all_available_gaussian_projects(self):
        """
        Get list of all available Gaussian project configurations.
        
        Returns:
            list: List of Gaussian project names (without .yaml extension)
                  available in the user Gaussian settings directory.
        """
        gaussian_project_settings = [
            os.path.basename(g).removesuffix(".yaml")  # python 3.9+
            # os.path.basename(g).strip(".yaml") does not work since
            # m062x.yaml gets stripped to 062x
            for g in self.gaussian_project_yaml_files
        ]
        return gaussian_project_settings

    @cached_property
    def all_available_orca_projects(self):
        """
        Get list of all available ORCA project configurations.
        
        Returns:
            list: List of ORCA project names (without .yaml extension)
                  available in the user ORCA settings directory.
        """
        return [
            os.path.basename(o).removesuffix(".yaml")
            for o in self.orca_project_yaml_files
        ]
