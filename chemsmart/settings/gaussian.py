import logging
import os

from chemsmart.jobs.gaussian.settings import (
    GaussianIRCJobSettings,
    GaussianJobSettings,
    GaussianTDDFTJobSettings,
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


class GaussianProjectSettings(RegistryMixin):
    """
    Base class for Gaussian project settings with default configurations.

    This class provides a framework for managing Gaussian calculation settings
    across different job types. It defines default settings for common quantum
    chemistry calculations including geometry optimization, transition state
    searches, and various analysis methods.

    Attributes:
        PROJECT_NAME (str): Identifier for the project configuration.
        functional (str, optional): DFT functional to use for calculations.
        small_basis (str, optional): Basis set for initial/fast calculations.
        large_basis (str, optional): Basis set for high-accuracy calculations.
    """

    PROJECT_NAME = "general"
    functional = None
    small_basis = None
    large_basis = None

    def main_settings(self):
        """
        Get main Gaussian settings with key default values.

        Returns:
            GaussianJobSettings: Default job settings with functional and basis set.
        """
        default_gaussian_job_settings = GaussianJobSettings.default()
        default_gaussian_job_settings.functional = self.functional
        default_gaussian_job_settings.basis = self.small_basis
        return default_gaussian_job_settings

    def opt_settings(self):
        """
        Get default settings for geometry optimization jobs.

        Returns:
            GaussianJobSettings: Settings configured for structure optimization.
        """
        settings = self.main_settings().copy()
        settings.jobtype = "opt"
        return settings

    def modred_settings(self):
        """
        Get default settings for modified redundant coordinate optimization.

        Returns:
            GaussianJobSettings: Settings for constrained optimization jobs.
        """
        settings = self.main_settings().copy()
        settings.jobtype = "modred"
        return settings

    def ts_settings(self):
        """
        Get default settings for transition state optimization jobs.

        Returns:
            GaussianJobSettings: Settings for transition state searches.
        """
        settings = self.main_settings().copy()
        settings.jobtype = "ts"
        return settings

    def irc_settings(self):
        """
        Get default settings for Intrinsic Reaction Coordinate calculations.

        Returns:
            GaussianIRCJobSettings: IRC-specific settings with frequency disabled.
        """
        settings = self.main_settings().copy()
        # Convert to IRC-specific settings class
        settings = GaussianIRCJobSettings(**settings.__dict__)
        settings.jobtype = "irc"
        settings.freq = False
        return settings

    def scan_settings(self):
        """
        Get default settings for potential energy surface scan calculations.

        Returns:
            GaussianJobSettings: Settings for coordinate scanning with frequency disabled.
        """
        settings = self.main_settings().copy()
        settings.jobtype = "scan"
        settings.freq = False
        return settings

    def nci_settings(self):
        """
        Get default settings for Non-Covalent Interaction analysis.

        Returns:
            GaussianJobSettings: Settings for NCI analysis with frequency disabled.
        """
        settings = self.main_settings().copy()
        settings.jobtype = "nci"
        settings.freq = False
        return settings

    def wbi_settings(self):
        """
        Get default settings for Wiberg Bond Index calculations.

        Returns:
            GaussianJobSettings: Settings for WBI analysis with frequency disabled.
        """
        settings = self.main_settings().copy()
        settings.jobtype = "wbi"
        settings.freq = False
        return settings

    def sp_settings(self):
        """
        Get default settings for single point energy calculations.

        Uses the large basis set for higher accuracy and disables frequency
        calculations for computational efficiency.

        Returns:
            GaussianJobSettings: Settings for single point calculations.
        """
        settings = self.main_settings().copy()
        settings.jobtype = "sp"
        settings.freq = False  # Disable frequency calculation for efficiency
        settings.basis = self.large_basis  # Use high-accuracy basis set
        return settings

    @classmethod
    def from_project(cls, project):
        """
        Create project settings instance based on project name.

        Searches for project configuration in the following order:
        1. User-defined project settings directory
        2. Chemsmart test project configurations

        Args:
            project (str): Name of the project configuration to load.

        Returns:
            GaussianProjectSettings: Configured settings instance.

        Raises:
            FileNotFoundError: If no configuration is found for the specified project.
        """
        # First try user-defined project settings
        user_project_settings = cls._from_user_project_name(project)
        if user_project_settings is not None:
            return user_project_settings
        else:
            # Fall back to chemsmart test project settings
            chemsmart_test_project_settings = (
                cls._from_chemsmart_test_projects(project)
            )
            if chemsmart_test_project_settings is not None:
                return chemsmart_test_project_settings

        # Generate helpful error message with available options
        templates_path = os.path.join(os.path.dirname(__file__), "templates")
        raise FileNotFoundError(
            f"No project settings implemented for {project}.\n\n"
            f"Place new gaussian project settings .yaml file in {user_settings.user_gaussian_settings_dir}.\n\n"
            f"Templates for such settings.yaml files are available at {templates_path}\n\n "
            f"Currently available projects: {user_settings.all_available_gaussian_projects}"
        )

    @classmethod
    def _from_projects_manager(cls, manager):
        """
        Create settings from a project manager instance.

        Args:
            manager: Project settings manager instance.

        Returns:
            YamlGaussianProjectSettings or None: YAML-based settings instance if successful, None if failed.
        """
        try:
            return manager.create()
        except FileNotFoundError:
            return None

    @classmethod
    def _from_user_project_name(cls, project_name):
        """
        Load project settings from user directory based on project name.

        Args:
            project_name (str): Name of the project configuration file.

        Returns:
            YamlGaussianProjectSettings or None: YAML-based settings instance if found, None otherwise.
        """
        project_name_yaml_path = os.path.join(
            ChemsmartUserSettings().user_gaussian_settings_dir,
            f"{project_name}.yaml",
        )
        user_settings_manager = GaussianProjectSettingsManager(
            filename=project_name_yaml_path
        )
        settings = cls._from_projects_manager(user_settings_manager)

        if settings is not None:
            return settings

    @classmethod
    def _from_chemsmart_test_projects(cls, project_name):
        """
        Load project settings from chemsmart test projects directory.

        Args:
            project_name (str): Name of the test project configuration.

        Returns:
            YamlGaussianProjectSettings or None: YAML-based settings instance if found, None otherwise.
        """
        current_file_dir = os.path.dirname(os.path.abspath(__file__))
        test_projects_dir = os.path.join(
            current_file_dir, "../../tests/data/GaussianTests/project_yaml"
        )

        project_name_yaml_path = os.path.join(
            test_projects_dir, f"{project_name}.yaml"
        )
        project_settings_manager = GaussianProjectSettingsManager(
            filename=project_name_yaml_path
        )
        settings = cls._from_projects_manager(project_settings_manager)

        if settings is not None:
            return settings


class YamlGaussianProjectSettings(GaussianProjectSettings):
    """
    YAML-based implementation of Gaussian project settings.

    This class loads and manages Gaussian calculation settings from YAML
    configuration files. It provides specific settings for each job type
    as defined in the YAML configuration.

    Attributes:
        PROJECT_NAME (str): Set to "yaml" to identify YAML-based settings.
    """

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
        """
        Initialize YAML-based project settings.

        Args:
            opt_settings: Settings for optimization jobs.
            modred_settings: Settings for modified redundant coordinate jobs.
            ts_settings: Settings for transition state jobs.
            irc_settings: Settings for IRC calculations.
            scan_settings: Settings for coordinate scanning jobs.
            nci_settings: Settings for NCI analysis jobs.
            sp_settings: Settings for single point calculations.
            td_settings: Settings for TD-DFT calculations.
            wbi_settings: Settings for Wiberg bond index calculations.
        """
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
        builder = YamlGaussianProjectSettingsBuilder(filename=filename)
        return builder.build()


class YamlGaussianProjectSettingsBuilder:
    """
    Builder class for constructing YAML-based Gaussian project settings.

    This class reads YAML configuration files and builds appropriate
    GaussianJobSettings instances for each job type. It handles the
    mapping between YAML configuration and specific settings classes.

    Attributes:
        filename (str): Path to the YAML configuration file.
    """

    def __init__(self, filename):
        """
        Initialize the builder with a YAML configuration file.

        Args:
            filename (str): Path to the YAML configuration file.
        """
        self.filename = filename

    def build(self):
        """
        Build a complete YamlGaussianProjectSettings instance.

        Reads the YAML configuration and creates settings for all supported
        job types including optimization, transition states, IRC, scans, etc.

        Returns:
            YamlGaussianProjectSettings: Fully configured project settings.
        """
        # Build settings for each supported job type
        opt_settings = self._project_settings_for_job(jobtype="opt")
        modred_settings = self._project_settings_for_job(jobtype="modred")
        ts_settings = self._project_settings_for_job(jobtype="ts")
        irc_settings = self._project_settings_for_job(jobtype="irc")
        scan_settings = self._project_settings_for_job(jobtype="scan")
        nci_settings = self._project_settings_for_job(jobtype="nci")
        sp_settings = self._project_settings_for_job(jobtype="sp")
        td_settings = self._project_settings_for_job(jobtype="td")
        wbi_settings = self._project_settings_for_job(jobtype="wbi")

        # Create the project settings instance
        project_settings = YamlGaussianProjectSettings(
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

        # Set project name from filename and return
        name = self._parse_project_name()
        project_settings.PROJECT_NAME = name
        return project_settings

    def _read_config(self):
        """
        Read and parse the YAML configuration file.

        Returns:
            dict: Parsed YAML configuration data.
        """
        from chemsmart.jobs.settings import read_molecular_job_yaml

        return read_molecular_job_yaml(self.filename, program="gaussian")

    def _project_settings_for_job(self, jobtype):
        """
        Create job-specific settings from YAML configuration.

        Maps job types to appropriate settings classes and creates
        configured instances based on YAML data.

        Args:
            jobtype (str): Type of job (opt, ts, irc, td, etc.).

        Returns:
            GaussianJobSettings: Configured settings for the specified job type.

        Raises:
            RuntimeError: If configuration for the job type is not found.
        """
        # Map job types to their specific settings classes
        settings_mapping = {
            "irc": GaussianIRCJobSettings,
            "td": GaussianTDDFTJobSettings,
        }

        try:
            jobtype_config = self._read_config().get(jobtype)
            if jobtype_config is not None:
                # Use specific settings class if available, otherwise default
                return settings_mapping.get(
                    jobtype, GaussianJobSettings
                ).from_dict(jobtype_config)
        except KeyError as e:
            available_jobs = list(self._read_config().keys())
            raise RuntimeError(
                f"Gaussian settings for job {jobtype} cannot be found!\n"
                f"Available Gaussian jobs with settings are: {available_jobs}"
            ) from e

    def _parse_project_name(self):
        """
        Extract project name from the YAML filename.

        Returns:
            str: Project name derived from filename (without extension).
        """
        return os.path.basename(self.filename).split(".")[0]


class GaussianProjectSettingsManager:
    """
    Manager for Gaussian project settings from YAML configuration files.

    Provides management interface for loading Gaussian computational chemistry
    project settings from YAML files in a specified folder structure. Handles
    file validation and project settings creation.

    Attributes:
        filename (str): Absolute path to the YAML configuration file.
    """

    def __init__(self, filename):
        """
        Initialize the Gaussian project settings manager.

        Args:
            filename (str): Path to YAML configuration file containing
                Gaussian project settings.

        Raises:
            ValueError: If filename is None or not specified.
        """
        if filename is None:
            raise ValueError("filename is not specified")
        self.filename = os.path.abspath(filename)

    def create(self):
        """
        Create project settings from the specified YAML file.

        Loads and parses the YAML configuration file to create a complete
        Gaussian project settings instance with all job configurations.

        Returns:
            YamlGaussianProjectSettings: Configured project settings loaded
                from the YAML file.

        Raises:
            FileNotFoundError: If the specified YAML file does not exist.
            ValueError: If the YAML file is malformed or invalid.
        """
        return YamlGaussianProjectSettings.from_yaml(self.filename)
