"""PySCF project settings.

A PySCF project YAML uses the same ``gas:`` / ``solv:`` shape as Gaussian and
ORCA, because PySCF's functional, basis and phase genuinely vary per job --
which is the reason that split exists. ``gas:`` feeds ``opt`` and ``hess``;
``solv:`` feeds ``sp``. That is the standard workflow: optimise and take
frequencies at one level, refine the single point in solvent.

Unlike Gaussian and ORCA, PySCF passes its own jobtype lists to
``read_molecular_job_yaml`` so that exactly three configs are produced rather
than the fifteen the shared default emits.
"""

import hashlib
import logging
import os

from chemsmart.jobs.pyscf.settings import PySCFJobSettings
from chemsmart.settings.user import CHEMSMARTUserSettings
from chemsmart.utils.mixins import RegistryMixin

user_settings = CHEMSMARTUserSettings()

logger = logging.getLogger(__name__)

#: Jobtypes fed by the ``gas:`` section of a PySCF project YAML.
PYSCF_GAS_PHASE_JOBS = ["opt", "hess"]

#: Jobtypes fed by the ``solv:`` section.
PYSCF_SP_JOBS = ["sp"]

#: Every jobtype a PySCF project YAML produces.
PYSCF_JOBTYPES = PYSCF_GAS_PHASE_JOBS + PYSCF_SP_JOBS


class PySCFProjectSettings(RegistryMixin):
    """Base PySCF project settings.

    PySCF has no defensible default functional or basis -- unlike xTB, whose
    GFN2 default is a complete method -- so a project YAML is required. The
    base class exists to define the interface and to fail with a discoverable
    message when a named project is missing.
    """

    PROJECT_NAME = "general"
    functional = None
    basis = None

    def main_settings(self):
        """Return the base settings all jobtypes derive from."""
        settings = PySCFJobSettings.default()
        settings.functional = self.functional
        settings.basis = self.basis
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
        """Resolve a named project to settings.

        Order: user config dir, then the packaged test projects, then raise
        with the list of projects that were actually found.
        """
        user_project_settings = cls._from_user_project_name(project)
        if user_project_settings is not None:
            return user_project_settings

        test_project_settings = cls._from_chemsmart_test_projects(project)
        if test_project_settings is not None:
            return test_project_settings

        templates_path = os.path.join(os.path.dirname(__file__), "templates")
        raise FileNotFoundError(
            f"No PySCF project settings implemented for {project}.\n\n"
            f"Place new PySCF project settings .yaml file in "
            f"{user_settings.user_pyscf_settings_dir}.\n\n"
            f"Templates for such settings.yaml files are available at "
            f"{templates_path}\n\n "
            f"Currently available projects: "
            f"{user_settings.all_available_pyscf_projects}"
        )

    @classmethod
    def _from_projects_manager(cls, manager):
        try:
            return manager.create()
        except FileNotFoundError:
            return None

    @classmethod
    def _from_user_project_name(cls, project_name):
        project_name_yaml_path = os.path.join(
            CHEMSMARTUserSettings().user_pyscf_settings_dir,
            f"{project_name}.yaml",
        )
        manager = PySCFProjectSettingsManager(filename=project_name_yaml_path)
        return cls._from_projects_manager(manager)

    @classmethod
    def _from_chemsmart_test_projects(cls, project_name):
        current_file_dir = os.path.dirname(os.path.abspath(__file__))
        test_projects_dir = os.path.join(
            current_file_dir, "../../tests/data/PySCFTests/project_yaml"
        )
        project_name_yaml_path = os.path.join(
            test_projects_dir, f"{project_name}.yaml"
        )
        manager = PySCFProjectSettingsManager(filename=project_name_yaml_path)
        return cls._from_projects_manager(manager)


class YamlPySCFProjectSettings(PySCFProjectSettings):
    """PySCF project settings backed by a project YAML."""

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
        return YamlPySCFProjectSettingsBuilder(filename=filename).build()


class YamlPySCFProjectSettingsBuilder:
    """Build :class:`YamlPySCFProjectSettings` from a project YAML."""

    def __init__(self, filename):
        self.filename = filename

    def build(self):
        project_settings = YamlPySCFProjectSettings(
            sp_settings=self._project_settings_for_job("sp"),
            opt_settings=self._project_settings_for_job("opt"),
            hess_settings=self._project_settings_for_job("hess"),
        )
        project_settings.PROJECT_NAME = self._parse_project_name()
        return project_settings

    def _read_config(self):
        from chemsmart.jobs.settings import read_molecular_job_yaml

        return read_molecular_job_yaml(
            self.filename,
            program="pyscf",
            gas_phase_jobs=PYSCF_GAS_PHASE_JOBS,
            sp_jobs=PYSCF_SP_JOBS,
        )

    def _project_settings_for_job(self, jobtype):
        config = self._read_config()
        jobtype_config = config.get(jobtype)
        if jobtype_config is None:
            raise RuntimeError(
                f"PySCF settings for job {jobtype} cannot be found!\n"
                f"Available PySCF jobs with settings are: "
                f"{list(config.keys())}"
            )
        settings = PySCFJobSettings.from_dict(jobtype_config)
        settings._project_yaml_digest = self._project_yaml_digest()
        return settings

    def _project_yaml_digest(self):
        """Return the SHA-256 identity of the resolved project artifact."""
        digest = hashlib.sha256()
        with open(self.filename, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()

    def _parse_project_name(self):
        return os.path.basename(self.filename).split(".")[0]


class PySCFProjectSettingsManager:
    """Load PySCF project settings from a YAML file."""

    def __init__(self, filename):
        if filename is None:
            raise ValueError("filename is not specified")
        self.filename = os.path.abspath(filename)

    def create(self):
        return YamlPySCFProjectSettings.from_yaml(self.filename)
