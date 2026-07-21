import logging
import os

import yaml

from chemsmart.jobs.xtb.settings import XTBJobSettings
from chemsmart.settings.user import ChemsmartUserSettings
from chemsmart.settings.workspace_project import workspace_project_path
from chemsmart.utils.mixins import RegistryMixin

user_settings = ChemsmartUserSettings()
logger = logging.getLogger(__name__)


class XTBProjectSettings(RegistryMixin):
    """Base xTB project settings."""

    PROJECT_NAME = "general"
    gfn_version = "gfn2"

    def main_settings(self):
        settings = XTBJobSettings.default()
        settings.gfn_version = self.gfn_version
        return settings

    def sp_settings(self):
        settings = self.main_settings().copy()
        settings.jobtype = "sp"
        return settings

    def opt_settings(self):
        settings = self.main_settings().copy()
        settings.jobtype = "opt"
        return settings

    def hess_settings(self):
        settings = self.main_settings().copy()
        settings.jobtype = "hess"
        return settings

    @classmethod
    def from_project(cls, project):
        # xTB needs no project YAML. Its method space is just the GFN
        # Hamiltonian (gfn0/gfn1/gfn2/gfnff) with a strong GFN2 default, unlike
        # the functional/basis/solvent choice a Gaussian or ORCA project must
        # encode. A bare ``chemsmart run xtb -f mol.xyz opt`` (no ``-p``)
        # therefore runs from XTBJobSettings.default() plus any CLI overrides;
        # this base class already exposes those defaults for every leaf.
        if project is None:
            return cls()

        # A named project still layers workspace -> user -> packaged, so a
        # project can pin non-default settings when a workflow wants them.
        workspace_project_settings = cls._from_workspace_project_name(project)
        if workspace_project_settings is not None:
            return workspace_project_settings

        user_project_settings = cls._from_user_project_name(project)
        if user_project_settings is not None:
            return user_project_settings

        packaged_project_settings = cls._from_packaged_project_name(project)
        if packaged_project_settings is not None:
            return packaged_project_settings

        # A project was explicitly named but not found: that is a user error,
        # not a request for defaults, so fail with a discoverable path.
        templates_path = os.path.join(os.path.dirname(__file__), "templates")
        raise FileNotFoundError(
            f"No xTB project settings implemented for {project}.\n\n"
            f"xTB does not require a project; omit -p to use GFN2 defaults, "
            f"or place a project .yaml file in "
            f"{user_settings.user_xtb_settings_dir}.\n\n"
            f"Templates for such settings.yaml files are available at "
            f"{templates_path}\n\n "
            f"Currently available projects: "
            f"{user_settings.all_available_xtb_projects}"
        )

    @classmethod
    def _from_projects_manager(cls, manager):
        try:
            return manager.create()
        except FileNotFoundError:
            return None

    @classmethod
    def _from_workspace_project_name(cls, project_name):
        if project_name is None:
            return None
        project_name_yaml_path = workspace_project_path(project_name, "xtb")
        manager = XTBProjectSettingsManager(
            filename=str(project_name_yaml_path)
        )
        return cls._from_projects_manager(manager)

    @classmethod
    def _from_user_project_name(cls, project_name):
        if project_name is None:
            return None
        project_name_yaml_path = os.path.join(
            ChemsmartUserSettings().user_xtb_settings_dir,
            f"{project_name}.yaml",
        )
        manager = XTBProjectSettingsManager(filename=project_name_yaml_path)
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
            "xtb",
            f"{project_name}.yaml",
        )
        manager = XTBProjectSettingsManager(filename=project_name_yaml_path)
        return cls._from_projects_manager(manager)


class YamlXTBProjectSettings(XTBProjectSettings):
    """YAML-backed xTB project settings."""

    PROJECT_NAME = "yaml"

    def __init__(self, sp_settings, opt_settings, hess_settings):
        self._sp_settings = sp_settings
        self._opt_settings = opt_settings
        self._hess_settings = hess_settings

    def sp_settings(self):
        return self._sp_settings.copy()

    def opt_settings(self):
        return self._opt_settings.copy()

    def hess_settings(self):
        return self._hess_settings.copy()

    @classmethod
    def from_yaml(cls, filename):
        return YamlXTBProjectSettingsBuilder(filename=filename).build()


class YamlXTBProjectSettingsBuilder:
    def __init__(self, filename):
        self.filename = os.path.abspath(filename)

    def build(self):
        project_settings = YamlXTBProjectSettings(
            sp_settings=self._settings_for_job("sp"),
            opt_settings=self._settings_for_job("opt"),
            hess_settings=self._settings_for_job("hess"),
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
                f"No xTB configuration found for {jobtype} in "
                f"{self.filename}. Using defaults."
            )
        config = dict(config)
        config.setdefault("jobtype", jobtype)
        return XTBJobSettings.from_dict(config)

    def _parse_project_name(self):
        return os.path.basename(self.filename).removesuffix(".yaml")


class XTBProjectSettingsManager:
    """Manages xTB project settings from YAML files."""

    def __init__(self, filename):
        if filename is None:
            raise ValueError("filename is not specified")
        self.filename = os.path.abspath(filename)

    def create(self):
        return YamlXTBProjectSettings.from_yaml(self.filename)
