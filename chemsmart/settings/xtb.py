import logging
import os

from chemsmart.jobs.xtb.settings import XTBJobSettings
from chemsmart.settings.user import ChemsmartUserSettings
from chemsmart.utils.mixins import RegistryMixin

user_settings = ChemsmartUserSettings()


logger = logging.getLogger(__name__)
project_settings_registry: list[str] = []
# Type annotation required for Python's type checker (e.g., mypy)
# which has stricter requirements when working in contexts where type inference
# is insufficient or ambiguous.
# This is not due to Python 3.9 itself, but type checking rules enforced by tools
# like mypy or stricter typing practices.


class XTBProjectSettings(RegistryMixin):
    """Most general xTB settings class with key defaults."""

    PROJECT_NAME = "general"
    gfn_version = "gfn2"

    def main_settings(self):
        """xTB main settings with key default values."""
        default_xtb_job_settings = XTBJobSettings.default()
        default_xtb_job_settings.gfn_version = self.gfn_version
        return default_xtb_job_settings

    def sp_settings(self):
        """xTB default settings for sp job."""
        settings = self.main_settings().copy()
        settings.job_type = "sp"
        return settings

    def opt_settings(self):
        """xTB default settings for opt job."""
        settings = self.main_settings().copy()
        settings.job_type = "opt"
        return settings

    def hess_settings(self):
        """xTB default settings for hess job."""
        settings = self.main_settings().copy()
        settings.job_type = "hess"
        return settings

    def md_settings(self):
        """xTB default settings for md job."""
        settings = self.main_settings().copy()
        settings.job_type = "md"
        return settings

    @classmethod
    def from_project(cls, project):
        """Get project settings based on project name."""
        user_project_settings = cls._from_user_project_name(project)
        if user_project_settings is not None:
            return user_project_settings
        else:
            chemsmart_test_project_settings = (
                cls._from_chemsmart_test_projects(project)
            )
            if chemsmart_test_project_settings is not None:
                return chemsmart_test_project_settings

        templates_path = os.path.join(os.path.dirname(__file__), "templates")
        raise FileNotFoundError(
            f"No project settings implemented for {project}.\n\n"
            f"Place new xTB project settings .yaml file in {user_settings.user_xtb_settings_dir}.\n\n"
            f"Templates for such settings.yaml files are available at {templates_path}\n\n "
            f"Currently available projects: {user_settings.all_available_xtb_projects}"
        )

    @classmethod
    def _from_projects_manager(cls, manager):
        try:
            return manager.create()
        except FileNotFoundError:
            return None

    @classmethod
    def _from_user_project_name(cls, project_name):
        """Get .yaml project settings file from user directory based on project name."""
        project_name_yaml_path = os.path.join(
            ChemsmartUserSettings().user_xtb_settings_dir,
            f"{project_name}.yaml",
        )
        user_settings_manager = XTBProjectSettingsManager(
            filename=project_name_yaml_path
        )
        settings = cls._from_projects_manager(user_settings_manager)

        if settings is not None:
            return settings


class YamlXTBProjectSettings(XTBProjectSettings):
    PROJECT_NAME = "yaml"

    def __init__(
        self,
        sp_settings,
        opt_settings,
        hess_settings,
        md_settings,
    ):
        self._sp_settings = sp_settings
        self._opt_settings = opt_settings
        self._hess_settings = hess_settings
        self._md_settings = md_settings

    def sp_settings(self):
        return self._sp_settings

    def opt_settings(self):
        return self._opt_settings

    def hess_settings(self):
        return self._hess_settings

    def md_settings(self):
        return self._md_settings

    @classmethod
    def from_yaml(cls, filename):
        builder = YamlXTBProjectSettingsBuilder(filename=filename)
        return builder.build()


class YamlXTBProjectSettingsBuilder:
    def __init__(self, filename):
        self.filename = filename

    def build(self):
        sp_settings = self._project_settings_for_job(job_type="sp")
        opt_settings = self._project_settings_for_job(job_type="opt")
        hess_settings = self._project_settings_for_job(job_type="hess")
        md_settings = self._project_settings_for_job(job_type="md")

        # Internal job_type promotion based on opt flag
        if hess_settings.opt:
            hess_settings.job_type = "ohess"

        if md_settings.opt:
            md_settings.job_type = "omd"

        project_settings = YamlXTBProjectSettings(
            sp_settings=sp_settings,
            opt_settings=opt_settings,
            hess_settings=hess_settings,
            md_settings=md_settings,
        )

        name = self._parse_project_name()
        project_settings.PROJECT_NAME = name
        return project_settings

    def _read_config(self):
        from chemsmart.jobs.settings import read_molecular_job_yaml

        return read_molecular_job_yaml(self.filename, program="xtb")

    def _project_settings_for_job(self, job_type):
        try:
            job_type_config = self._read_config().get(job_type)
            if job_type_config is not None:
                if "job_type" not in job_type_config:
                    job_type_config["job_type"] = job_type
                return XTBJobSettings.from_dict(job_type_config)
            # Fallback to default
            logger.warning(
                f"No specific configuration found for '{job_type}' in project settings. "
                f"Using default settings."
            )
            settings = XTBJobSettings.default()
            settings.job_type = job_type
            return settings

        except KeyError:
            # Should practically not happen with .get() but if read_config fails structure
            raise RuntimeError(f"Error reading settings for {job_type}")

    def _parse_project_name(self):
        return os.path.basename(self.filename).split(".")[0]


class XTBProjectSettingsManager:
    """Manages XTB project settings specified in the form of yaml files in a folder.

    Args:
        filename: yaml filename in the default XTB projects folder.
    """

    def __init__(self, filename):
        if filename is None:
            raise ValueError("filename is not specified")
        self.filename = os.path.abspath(filename)

    def create(self):
        return YamlXTBProjectSettings.from_yaml(self.filename)
