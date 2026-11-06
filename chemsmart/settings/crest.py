import logging
import os

import yaml

from chemsmart.jobs.crest.settings import CRESTJobSettings
from chemsmart.settings.user import CHEMSMARTUserSettings
from chemsmart.utils.mixins import RegistryMixin

user_settings = CHEMSMARTUserSettings()
logger = logging.getLogger(__name__)


class CRESTProjectSettings(RegistryMixin):
    """Base CREST project settings."""

    PROJECT_NAME = "general"
    gfn_version = "gfn2"

    def main_settings(self):
        settings = CRESTJobSettings.default()
        settings.gfn_version = self.gfn_version
        return settings

    def conformer_settings(self):
        settings = self.main_settings().copy()
        settings.jobtype = "conformers"
        return settings

    @classmethod
    def from_project(cls, project):
        user_project_settings = cls._from_user_project_name(project)
        if user_project_settings is not None:
            return user_project_settings

        packaged_project_settings = cls._from_packaged_project_name(project)
        if packaged_project_settings is not None:
            return packaged_project_settings

        templates_path = os.path.join(os.path.dirname(__file__), "templates")
        raise FileNotFoundError(
            f"No CREST project settings implemented for {project}.\n\n"
            f"Place new CREST project settings .yaml file in "
            f"{user_settings.user_crest_settings_dir}.\n\n"
            f"Templates for such settings.yaml files are available at "
            f"{templates_path}\n\n "
            f"Currently available projects: "
            f"{user_settings.all_available_crest_projects}"
        )

    @classmethod
    def _from_projects_manager(cls, manager):
        try:
            return manager.create()
        except FileNotFoundError:
            return None

    @classmethod
    def _from_user_project_name(cls, project_name):
        if project_name is None:
            return None
        project_name_yaml_path = os.path.join(
            CHEMSMARTUserSettings().user_crest_settings_dir,
            f"{project_name}.yaml",
        )
        manager = CRESTProjectSettingsManager(filename=project_name_yaml_path)
        return cls._from_projects_manager(manager)

    @classmethod
    def _from_packaged_project_name(cls, project_name):
        if project_name is None:
            return None
        current_file_dir = os.path.dirname(os.path.abspath(__file__))
        project_name_yaml_path = os.path.join(
            current_file_dir,
            "templates",
            ".chemsmart",
            "crest",
            f"{project_name}.yaml",
        )
        manager = CRESTProjectSettingsManager(filename=project_name_yaml_path)
        return cls._from_projects_manager(manager)


class YamlCRESTProjectSettings(CRESTProjectSettings):
    """YAML-backed CREST project settings."""

    PROJECT_NAME = "yaml"

    def __init__(self, conformer_settings):
        self._conformer_settings = conformer_settings

    def conformer_settings(self):
        return self._conformer_settings.copy()

    @classmethod
    def from_yaml(cls, filename):
        return YamlCRESTProjectSettingsBuilder(filename=filename).build()


class YamlCRESTProjectSettingsBuilder:
    def __init__(self, filename):
        self.filename = os.path.abspath(filename)

    def build(self):
        project_settings = YamlCRESTProjectSettings(
            conformer_settings=self._settings_for_job("conformers"),
        )
        project_settings.PROJECT_NAME = self._parse_project_name()
        return project_settings

    def _read_config(self):
        with open(self.filename) as handle:
            return yaml.safe_load(handle) or {}

    def _settings_for_job(self, jobtype):
        config = self._read_config().get(jobtype, {})
        if not config:
            logger.warning(
                f"No CREST configuration found for {jobtype} in "
                f"{self.filename}. Using defaults."
            )
        config = dict(config)
        config.setdefault("jobtype", jobtype)
        return CRESTJobSettings.from_dict(config)

    def _parse_project_name(self):
        return os.path.basename(self.filename).removesuffix(".yaml")


class CRESTProjectSettingsManager:
    """Manages CREST project settings from YAML files."""

    def __init__(self, filename):
        if filename is None:
            raise ValueError("filename is not specified")
        self.filename = os.path.abspath(filename)

    def create(self):
        return YamlCRESTProjectSettings.from_yaml(self.filename)
