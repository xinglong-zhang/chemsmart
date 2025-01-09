import logging
import os
from abc import abstractmethod
from functools import cache

from chemsmart.jobs.gaussian.settings import (
    GaussianJobSettings,
    GaussianIRCJobSettings,
    GaussianTDDFTJobSettings,
)
from chemsmart.settings.user import ChemsmartUserSettings

logger = logging.getLogger(__name__)
project_settings_registry = []


class GaussianProjectSettings:
    """Most general Gaussian settings class with key defaults."""

    PROJECT_NAME = "general"
    functional = None
    small_basis = None
    large_basis = None

    def main_settings(self):
        """Gaussian main settings with key default values."""
        default_gaussian_job_settings = GaussianJobSettings.default()
        default_gaussian_job_settings.functional = self.functional
        default_gaussian_job_settings.basis = self.small_basis
        return default_gaussian_job_settings

    def opt_settings(self):
        """Gaussian default settings for opt job."""
        settings = self.main_settings().copy()
        settings.job_type = "opt"
        return settings

    def modred_settings(self):
        """Gaussian default settings for modredundant job."""
        settings = self.main_settings().copy()
        settings.job_type = "modredundant"
        return settings

    def ts_settings(self):
        """Gaussian default settings for ts job."""
        settings = self.main_settings().copy()
        settings.job_type = "ts"
        return settings

    def irc_settings(self):
        """Gaussian default settings for irc job."""
        settings = self.main_settings().copy()
        settings = GaussianIRCJobSettings(
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
        settings.freq = False
        settings.basis = self.large_basis
        return settings


class YamlGaussianProjectSettings(GaussianProjectSettings):
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
        builder = YamlGaussianProjectSettingsBuilder(filename=filename)
        return builder.build()


class YamlGaussianProjectSettingsBuilder:
    def __init__(self, filename):
        self.filename = filename

    def build(self):
        opt_settings = self._project_settings_for_job(job_type="opt")
        modred_settings = self._project_settings_for_job(
            job_type="modredundant"
        )
        ts_settings = self._project_settings_for_job(job_type="ts")
        irc_settings = self._project_settings_for_job(job_type="irc")
        scan_settings = self._project_settings_for_job(job_type="scan")
        nci_settings = self._project_settings_for_job(job_type="nci")
        sp_settings = self._project_settings_for_job(job_type="sp")
        td_settings = self._project_settings_for_job(job_type="td")
        wbi_settings = self._project_settings_for_job(job_type="wbi")

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

        name = self._parse_project_name()
        project_settings.PROJECT_NAME = name
        return project_settings

    def _read_config(self):
        from chemsmart.jobs.settings import read_molecular_job_yaml

        return read_molecular_job_yaml(self.filename)

    def _project_settings_for_job(self, job_type):
        # Define a dictionary to map job_type to corresponding settings class
        settings_mapping = {
            "irc": GaussianIRCJobSettings,
            "td": GaussianTDDFTJobSettings,
        }

        try:
            job_type_config = self._read_config().get(job_type)
            if job_type_config is not None:
                return settings_mapping.get(
                    job_type, GaussianJobSettings
                ).from_dict(job_type_config)
        except KeyError as e:
            raise RuntimeError(
                f"Gaussian settings for job {job_type} cannot be found!\n"
                f"Available Gaussian jobs with settings are: {self._read_config().keys()}"
            ) from e

    def _parse_project_name(self):
        return os.path.basename(self.filename).split(".")[0]


class GaussianProjectSettingsManager:
    """Manages Gaussian project settings specified in the form of yaml files in a folder.

    Args:
        filename: yaml filename in the default Gaussian projects folder.
    """

    def __init__(self, filename):
        if filename is None:
            raise ValueError("filename is not specified")
        self.filename = os.path.abspath(filename)

    def create(self):
        return YamlGaussianProjectSettings.from_yaml(self.filename)
