import logging
import os

from chemsmart.jobs.orca.settings import (
    ORCAIRCJobSettings,
    ORCAJobSettings,
    ORCATSJobSettings,
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
    """
    Base class for ORCA project settings with default configurations.
    
    Provides fundamental settings and configurations for ORCA quantum chemistry
    calculations. Includes default values for computational parameters and
    methods for creating job-specific settings for various calculation types.
    
    Attributes:
        PROJECT_NAME (str): Name identifier for the project.
        functional (str): DFT functional to use for calculations.
        small_basis (str): Small basis set for initial calculations.
        large_basis (str): Large basis set for high-accuracy calculations.
    """

    PROJECT_NAME = "general"
    functional = None
    small_basis = None
    large_basis = None

    def main_settings(self):
        """
        Create main ORCA settings with default values.
        
        Provides the base settings configuration that serves as the foundation
        for all other job types. Includes functional and basis set assignments.
        
        Returns:
            ORCAJobSettings: Default ORCA job settings with functional and basis set.
        """
        default_orca_job_settings = ORCAJobSettings.default()
        default_orca_job_settings.functional = self.functional
        default_orca_job_settings.basis = self.small_basis
        return default_orca_job_settings

    def opt_settings(self):
        """
        ORCA settings for geometry optimization calculations.
        
        Creates settings for geometry optimization jobs that find
        minimum energy structures on the potential energy surface.
        
        Returns:
            ORCAJobSettings: Settings configured for geometry optimization.
        """
        settings = self.main_settings().copy()
        settings.job_type = "opt"
        return settings

    def modred_settings(self):
        """
        ORCA settings for modredundant coordinate optimization.
        
        Creates settings for constrained optimization using modredundant
        coordinates, allowing partial optimization of specific coordinates.
        
        Returns:
            ORCAJobSettings: Settings configured for modredundant optimization.
        """
        settings = self.main_settings().copy()
        settings.job_type = "modred"
        return settings

    def ts_settings(self):
        """
        ORCA settings for transition state optimization calculations.
        
        Creates settings for finding and optimizing transition states using
        specialized TS optimization algorithms in ORCA.
        
        Returns:
            ORCATSJobSettings: Settings configured for transition state optimization.
        """
        settings = self.main_settings().copy()
        settings = ORCATSJobSettings(
            **settings.__dict__
        )  # convert settings to ORCATSJobSettings
        settings.job_type = "ts"
        return settings

    def irc_settings(self):
        """
        ORCA settings for intrinsic reaction coordinate calculations.
        
        Creates settings for IRC calculations that trace reaction pathways
        from transition states to reactants and products. Frequency calculations
        are disabled by default for IRC jobs.
        
        Returns:
            ORCAIRCJobSettings: Settings configured for IRC calculations.
        """
        settings = self.main_settings().copy()
        settings = ORCAIRCJobSettings(
            **settings.__dict__
        )  # convert settings to ORCAIRCJobSettings
        settings.job_type = "irc"
        settings.freq = False
        return settings

    def scan_settings(self):
        """
        ORCA settings for potential energy surface scan calculations.
        
        Creates settings for relaxed coordinate scans to explore potential
        energy surfaces along specific reaction coordinates. Frequency
        calculations are disabled by default.
        
        Returns:
            ORCAJobSettings: Settings configured for coordinate scan calculations.
        """
        settings = self.main_settings().copy()
        settings.job_type = "scan"
        settings.freq = False
        return settings

    def nci_settings(self):
        """
        ORCA settings for non-covalent interaction analysis.
        
        Creates settings for NCI (Non-Covalent Interaction) analysis calculations
        that identify and characterize weak intermolecular interactions.
        Frequency calculations are disabled by default.
        
        Returns:
            ORCAJobSettings: Settings configured for NCI analysis.
        """
        settings = self.main_settings().copy()
        settings.job_type = "nci"
        settings.freq = False
        return settings

    def wbi_settings(self):
        """
        ORCA settings for Wiberg Bond Index calculations.
        
        Creates settings for WBI calculations that provide quantitative
        measures of bond orders and electron sharing between atoms.
        Frequency calculations are disabled by default.
        
        Returns:
            ORCAJobSettings: Settings configured for WBI calculations.
        """
        settings = self.main_settings().copy()
        settings.job_type = "wbi"
        settings.freq = False
        return settings

    def sp_settings(self):
        """
        ORCA settings for single point energy calculations.
        
        Creates settings for single point calculations at higher level of theory
        using large basis set. Frequency calculations are disabled and the
        large basis set is used for improved accuracy.
        
        Returns:
            ORCAJobSettings: Settings configured for single point calculations.
        """
        settings = self.main_settings().copy()
        settings.job_type = "sp"
        settings.freq = False  # turn off freq calculation for sp job
        settings.basis = self.large_basis
        return settings

    @classmethod
    def from_project(cls, project):
        """
        Get project settings based on project name.
        
        Loads project settings from various sources including user-defined
        settings and built-in test projects. Provides detailed error messages
        if the requested project is not found.
        
        Args:
            project (str): Name of the project to load settings for.
            
        Returns:
            YamlORCAProjectSettings: Configured YAML-based project settings instance.
            
        Raises:
            FileNotFoundError: If no project settings are found for the
                specified project name, with guidance on creating new settings.
        """
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
        """
        Load project settings using a project manager.
        
        Internal method for loading project settings through a manager instance.
        Handles exceptions and provides error logging for failed operations.
        
        Args:
            manager: Project settings manager instance.
            
        Returns:
            YamlORCAProjectSettings or None: YAML-based loaded settings or None if loading fails.
        """
        try:
            return manager.create()
        except FileNotFoundError:
            return None

    @classmethod
    def _from_user_project_name(cls, project_name):
        """
        Load YAML project settings from user directory.
        
        Searches for project settings files in the user's ORCA settings
        directory based on the project name.
        
        Args:
            project_name (str): Name of the project to load.
            
        Returns:
            YamlORCAProjectSettings or None: YAML-based loaded settings or None if not found.
        """
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
        """
        Load YAML project settings from test projects directory.
        
        Searches for project settings in the chemsmart test projects directory,
        typically used for built-in example configurations.
        
        Args:
            project_name (str): Name of the test project to load.
            
        Returns:
            YamlORCAProjectSettings or None: YAML-based loaded settings or None if not found.
        """
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
    """
    YAML-based ORCA project settings with job-specific configurations.
    
    This class loads and manages ORCA calculation settings from YAML
    configuration files, providing pre-configured settings for different
    job types. Each job type maintains its own settings instance.
    
    Attributes:
        PROJECT_NAME (str): Project name identifier for YAML settings.
        _opt_settings: Geometry optimization job settings.
        _modred_settings: Modredundant coordinate job settings.
        _ts_settings: Transition state optimization settings.
        _irc_settings: Intrinsic reaction coordinate settings.
        _scan_settings: Coordinate scan job settings.
        _nci_settings: Non-covalent interaction analysis settings.
        _sp_settings: Single point calculation settings.
        _td_settings: Time-dependent DFT settings.
        _wbi_settings: Wiberg bond index calculation settings.
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
        Initialize YAML-based ORCA project settings.
        
        Args:
            opt_settings: Settings for geometry optimization jobs.
            modred_settings: Settings for modredundant coordinate jobs.
            ts_settings: Settings for transition state optimization.
            irc_settings: Settings for IRC calculations.
            scan_settings: Settings for coordinate scan jobs.
            nci_settings: Settings for NCI analysis.
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
        """
        Get geometry optimization settings.
        
        Returns:
            ORCAJobSettings: Pre-configured optimization settings.
        """
        return self._opt_settings

    def modred_settings(self):
        """
        Get modredundant coordinate optimization settings.
        
        Returns:
            ORCAJobSettings: Pre-configured modredundant settings.
        """
        return self._modred_settings

    def ts_settings(self):
        """
        Get transition state optimization settings.
        
        Returns:
            ORCATSJobSettings: Pre-configured TS optimization settings.
        """
        return self._ts_settings

    def irc_settings(self):
        """
        Get intrinsic reaction coordinate calculation settings.
        
        Returns:
            ORCAIRCJobSettings: Pre-configured IRC calculation settings.
        """
        return self._irc_settings

    def scan_settings(self):
        """
        Get coordinate scan calculation settings.
        
        Returns:
            ORCAJobSettings: Pre-configured scan calculation settings.
        """
        return self._scan_settings

    def nci_settings(self):
        """
        Get non-covalent interaction analysis settings.
        
        Returns:
            ORCAJobSettings: Pre-configured NCI analysis settings.
        """
        return self._nci_settings

    def sp_settings(self):
        """
        Get single point calculation settings.
        
        Returns:
            ORCAJobSettings: Pre-configured single point settings.
        """
        return self._sp_settings

    def td_settings(self):
        """
        Get time-dependent DFT calculation settings.
        
        Returns:
            ORCAJobSettings: Pre-configured TD-DFT settings.
        """
        return self._td_settings

    def wbi_settings(self):
        """
        Get Wiberg bond index calculation settings.
        
        Returns:
            ORCAJobSettings: Pre-configured WBI calculation settings.
        """
        return self._wbi_settings

    @classmethod
    def from_yaml(cls, filename):
        """
        Create project settings from YAML configuration file.
        
        Args:
            filename (str): Path to YAML configuration file.
            
        Returns:
            YamlORCAProjectSettings: Configured project settings instance.
        """
        builder = YamlORCAProjectSettingsBuilder(filename=filename)
        return builder.build()


class YamlORCAProjectSettingsBuilder:
    """
    Builder class for constructing YAML-based ORCA project settings.
    
    This class reads YAML configuration files and builds appropriate
    ORCA project settings with job-specific configurations. Handles
    the mapping of YAML data to ORCA job settings objects.
    
    Attributes:
        filename (str): Path to the YAML configuration file.
    """
    
    def __init__(self, filename):
        """
        Initialize the YAML ORCA project settings builder.
        
        Args:
            filename (str): Path to YAML configuration file containing
                ORCA job definitions and settings.
        """
        self.filename = filename

    def build(self):
        """
        Build complete ORCA project settings from YAML configuration.
        
        Reads the YAML file and creates job-specific settings for all
        supported ORCA calculation types, constructing a complete
        project settings instance.
        
        Returns:
            YamlORCAProjectSettings: Complete project settings with all
                job types configured from YAML data.
                
        Raises:
            FileNotFoundError: If the YAML configuration file is not found.
            ValueError: If the YAML file contains invalid configuration.
        """
        opt_settings = self._project_settings_for_job(job_type="opt")
        modred_settings = self._project_settings_for_job(job_type="modred")
        ts_settings = self._project_settings_for_job(job_type="ts")
        irc_settings = self._project_settings_for_job(job_type="irc")
        scan_settings = self._project_settings_for_job(job_type="scan")
        nci_settings = self._project_settings_for_job(job_type="nci")
        sp_settings = self._project_settings_for_job(job_type="sp")
        td_settings = self._project_settings_for_job(job_type="td")
        wbi_settings = self._project_settings_for_job(job_type="wbi")

        # Create complete project settings with all job configurations
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

        return read_molecular_job_yaml(self.filename, program="orca")

    def _project_settings_for_job(self, job_type):
        """
        Create job-specific settings from YAML configuration.
        
        Maps job types to appropriate settings classes and creates
        configured instances based on YAML data.
        
        Args:
            job_type (str): Type of job (opt, ts, irc, scan, etc.).
            
        Returns:
            ORCAJobSettings or ORCAIRCJobSettings or ORCATSJobSettings: Configured 
                settings for the specified job type. Returns ORCAIRCJobSettings for 
                IRC jobs, ORCATSJobSettings for TS jobs, and ORCAJobSettings for 
                all other job types.
            
        Raises:
            RuntimeError: If configuration for the job type is not found.
        """
        # Map job types to their specific settings classes
        settings_mapping = {"irc": ORCAIRCJobSettings, "ts": ORCATSJobSettings}

        try:
            job_type_config = self._read_config().get(job_type)
            if job_type_config is not None:
                # Use specific settings class if available, otherwise default
                return settings_mapping.get(
                    job_type, ORCAJobSettings
                ).from_dict(job_type_config)
        except KeyError as e:
            available_jobs = list(self._read_config().keys())
            raise RuntimeError(
                f"ORCA settings for job {job_type} cannot be found!\n"
                f"Available ORCA jobs with settings are: {available_jobs}"
            ) from e

    def _parse_project_name(self):
        """
        Extract project name from the YAML filename.
        
        Returns:
            str: Project name derived from filename (without extension).
        """
        return os.path.basename(self.filename).split(".")[0]


class ORCAProjectSettingsManager:
    """
    Manager for ORCA project settings from YAML configuration files.
    
    Provides management interface for loading ORCA computational chemistry
    project settings from YAML files in a specified folder structure. Handles
    file validation and project settings creation.
    
    Attributes:
        filename (str): Absolute path to the YAML configuration file.
    """

    def __init__(self, filename):
        """
        Initialize the ORCA project settings manager.
        
        Args:
            filename (str): Path to YAML configuration file containing
                ORCA project settings.
                
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
        ORCA project settings instance with all job configurations.
        
        Returns:
            YamlORCAProjectSettings: Configured project settings loaded
                from the YAML file.
                
        Raises:
            FileNotFoundError: If the specified YAML file does not exist.
            ValueError: If the YAML file is malformed or invalid.
        """
        return YamlORCAProjectSettings.from_yaml(self.filename)
