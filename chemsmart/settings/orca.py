import logging
import os

from chemsmart.jobs.orca.settings import (
    ORCAIRCJobSettings,
    ORCAJobSettings,
)
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


class ORCAProjectSettings(RegistryMixin):
    """Most general Gaussian settings class with key defaults."""

    PROJECT_NAME = "general"
    functional = None
    small_basis = None
    large_basis = None

    def main_settings(self):
        """Gaussian main settings with key default values."""
        default_orca_job_settings = ORCAJobSettings.default()
        default_orca_job_settings.functional = self.functional
        default_orca_job_settings.basis = self.small_basis
        return default_orca_job_settings

    def opt_settings(self):
        """Gaussian default settings for opt job."""
        settings = self.main_settings().copy()
        settings.job_type = "opt"
        return settings

    def modred_settings(self):
        """Gaussian default settings for modred job."""
        settings = self.main_settings().copy()
        settings.job_type = "modred"
        return settings

    def ts_settings(self):
        """Gaussian default settings for ts job."""
        settings = self.main_settings().copy()
        settings.job_type = "ts"
        return settings

    def irc_settings(self):
        """Gaussian default settings for irc job."""
        settings = self.main_settings().copy()
        settings = ORCAIRCJobSettings(
            **settings.__dict__
        )  # convert settings to GaussianIRCJobSettings
        settings.job_type = "irc"
        settings.freq = False
        return settings

    def scan_settings(self):
        """Gaussian default settings for scan job."""
        settings = self.main_settings().copy()
        settings.job_type = "scan"
        settings.freq = False
        return settings

    def nci_settings(self):
        """Gaussian default settings for nci job."""
        settings = self.main_settings().copy()
        settings.job_type = "nci"
        settings.freq = False
        return settings

    def wbi_settings(self):
        """Gaussian default settings for WBI job."""
        settings = self.main_settings().copy()
        settings.job_type = "wbi"
        settings.freq = False
        return settings

    def sp_settings(self):
        """Gaussian default settings for sp job."""
        settings = self.main_settings().copy()
        settings.job_type = "sp"
        settings.freq = False  # turn off freq calculation for sp job
        settings.basis = self.large_basis
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
            f"Place new ORCA project settings .yaml file in {user_settings.user_orca_settings_dir}.\n\n"
            f"Templates for such settings.yaml files are available at {templates_path}\n\n "
            f"Currently available projects: {user_settings.all_available_orca_projects}"
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
            ChemsmartUserSettings().user_orca_settings_dir,
            f"{project_name}.yaml",
        )
        user_settings_manager = ORCAProjectSettingsManager(
            filename=project_name_yaml_path
        )
        settings = cls._from_projects_manager(user_settings_manager)

        if settings is not None:
            return settings

    @classmethod
    def _from_chemsmart_test_projects(cls, project_name):
        """Get .yaml project settings file from chemsmart test projects."""
        current_file_dir = os.path.dirname(os.path.abspath(__file__))
        test_projects_dir = os.path.join(
            current_file_dir, "../../tests/data/ORCATests/project_yaml"
        )

        project_name_yaml_path = os.path.join(
            test_projects_dir, f"{project_name}.yaml"
        )
        project_settings_manager = ORCAProjectSettingsManager(
            filename=project_name_yaml_path
        )
        settings = cls._from_projects_manager(project_settings_manager)

        if settings is not None:
            return settings


class YamlORCAProjectSettings(ORCAProjectSettings):
    PROJECT_NAME = "yaml"

    def __init__(
        self,
        opt_settings,
        modred_settings,
        ts_settings,
        irc_settings,
        scan_settings,
        nci_settings,
        sp_settings,
        td_settings,
        wbi_settings,
    ):
        self._opt_settings = opt_settings
        self._modred_settings = modred_settings
        self._ts_settings = ts_settings
        self._irc_settings = irc_settings
        self._scan_settings = scan_settings
        self._nci_settings = nci_settings
        self._sp_settings = sp_settings
        self._td_settings = td_settings
        self._wbi_settings = wbi_settings

    def opt_settings(self):
        return self._opt_settings

    def modred_settings(self):
        return self._modred_settings

    def ts_settings(self):
        return self._ts_settings

    def irc_settings(self):
        return self._irc_settings

    def scan_settings(self):
        return self._scan_settings

    def nci_settings(self):
        return self._nci_settings

    def sp_settings(self):
        return self._sp_settings

    def td_settings(self):
        return self._td_settings

    def wbi_settings(self):
        return self._wbi_settings

    @classmethod
    def from_yaml(cls, filename):
        builder = YamlORCAProjectSettingsBuilder(filename=filename)
        return builder.build()


class YamlORCAProjectSettingsBuilder:
    def __init__(self, filename):
        self.filename = filename

    def build(self):
        opt_settings = self._project_settings_for_job(job_type="opt")
        modred_settings = self._project_settings_for_job(job_type="modred")
        ts_settings = self._project_settings_for_job(job_type="ts")
        irc_settings = self._project_settings_for_job(job_type="irc")
        scan_settings = self._project_settings_for_job(job_type="scan")
        nci_settings = self._project_settings_for_job(job_type="nci")
        sp_settings = self._project_settings_for_job(job_type="sp")
        td_settings = self._project_settings_for_job(job_type="td")
        wbi_settings = self._project_settings_for_job(job_type="wbi")

        project_settings = YamlORCAProjectSettings(
            opt_settings=opt_settings,
            modred_settings=modred_settings,
            ts_settings=ts_settings,
            irc_settings=irc_settings,
            scan_settings=scan_settings,
            nci_settings=nci_settings,
            sp_settings=sp_settings,
            td_settings=td_settings,
            wbi_settings=wbi_settings,
        )

        name = self._parse_project_name()
        project_settings.PROJECT_NAME = name
        return project_settings

    def _read_config(self):
        from chemsmart.jobs.settings import read_molecular_job_yaml

        return read_molecular_job_yaml(self.filename, program="orca")

    def _project_settings_for_job(self, job_type):
        # Define a dictionary to map job_type to corresponding settings class
        settings_mapping = {
            "irc": ORCAIRCJobSettings,
        }

        try:
            job_type_config = self._read_config().get(job_type)
            if job_type_config is not None:
                return settings_mapping.get(
                    job_type, ORCAJobSettings
                ).from_dict(job_type_config)
        except KeyError as e:
            raise RuntimeError(
                f"Gaussian settings for job {job_type} cannot be found!\n"
                f"Available Gaussian jobs with settings are: {self._read_config().keys()}"
            ) from e

    def _parse_project_name(self):
        return os.path.basename(self.filename).split(".")[0]


class ORCAProjectSettingsManager:
    """Manages Gaussian project settings specified in the form of yaml files in a folder.

    Args:
        filename: yaml filename in the default Gaussian projects folder.
    """

    def __init__(self, filename):
        if filename is None:
            raise ValueError("filename is not specified")
        self.filename = os.path.abspath(filename)

    def create(self):
        return YamlORCAProjectSettings.from_yaml(self.filename)
