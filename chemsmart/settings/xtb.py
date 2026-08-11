"""Project-YAML settings for ChemSmart's strict xTB execution surface."""

import logging
import os
from collections.abc import Mapping

import yaml

from chemsmart.jobs.xtb.settings import XTBJobSettings
from chemsmart.settings.user import CHEMSMARTUserSettings
from chemsmart.utils.mixins import RegistryMixin

logger = logging.getLogger(__name__)


class XTBProjectSettings(RegistryMixin):
    """Default or YAML-backed settings for xTB ``sp``, ``opt`` and ``hess``."""

    PROJECT_NAME = "general"
    gfn_version = "gfn2"
    source_file = None

    def explicit_fields(self, jobtype):
        """Return fields explicitly declared for *jobtype* by this project.

        Built-in defaults are deliberately not explicit scientific state.
        This distinction lets an input artifact retain its charge and
        multiplicity unless an approved project section declares both.
        """

        if jobtype not in ("sp", "opt", "hess"):
            raise ValueError(f"Unknown xTB jobtype: {jobtype!r}.")
        return frozenset()

    def main_settings(self):
        settings = XTBJobSettings.default()
        settings.gfn_version = self.gfn_version
        return settings

    def sp_settings(self):
        settings = self.main_settings().copy()
        settings.jobtype = "sp"
        return settings.validate(expected_jobtype="sp")

    def opt_settings(self):
        settings = self.main_settings().copy()
        settings.jobtype = "opt"
        return settings.validate(expected_jobtype="opt")

    def hess_settings(self):
        settings = self.main_settings().copy()
        settings.jobtype = "hess"
        return settings.validate(expected_jobtype="hess")

    @classmethod
    def from_project(cls, project):
        """Resolve an explicit file, user project, or packaged project.

        Omitting ``-p`` is intentional for xTB and returns ChemSmart's GFN2
        defaults.  A named but unresolved project is never silently replaced by
        those defaults.
        """
        if project is None:
            return cls()

        project = os.fspath(project)
        if os.path.isfile(project):
            return XTBProjectSettingsManager(filename=project).create()

        for resolver in (
            cls._from_user_project_name,
            cls._from_packaged_project_name,
        ):
            project_settings = resolver(project)
            if project_settings is not None:
                return project_settings

        user_dir = cls._user_xtb_settings_dir()
        templates_path = os.path.join(os.path.dirname(__file__), "templates")
        raise FileNotFoundError(
            f"No xTB project settings implemented for {project!r}.\n\n"
            "Omit -p to use ChemSmart's GFN2 defaults, pass an explicit YAML "
            f"path, or place <name>.yaml in {user_dir}.\n\n"
            f"Templates are available at {templates_path}.\n\n"
            f"Available projects: {cls.available_projects()}"
        )

    @classmethod
    def _user_xtb_settings_dir(cls):
        # Keep this slice independent of the coordinator-owned user-settings
        # registry addition while honoring CHEMSMART_CONFIG_DIR.
        return os.path.join(CHEMSMARTUserSettings.resolve_config_dir(), "xtb")

    @classmethod
    def _packaged_xtb_settings_dir(cls):
        return os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "templates",
            ".chemsmart",
            "xtb",
        )

    @classmethod
    def available_projects(cls):
        names = set()
        for directory in (
            cls._user_xtb_settings_dir(),
            cls._packaged_xtb_settings_dir(),
        ):
            if not os.path.isdir(directory):
                continue
            for filename in os.listdir(directory):
                if filename.endswith((".yaml", ".yml")):
                    names.add(os.path.splitext(filename)[0])
        return sorted(names)

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
        for suffix in (".yaml", ".yml"):
            path = os.path.join(
                cls._user_xtb_settings_dir(), f"{project_name}{suffix}"
            )
            settings = cls._from_projects_manager(
                XTBProjectSettingsManager(path)
            )
            if settings is not None:
                return settings
        return None

    @classmethod
    def _from_packaged_project_name(cls, project_name):
        if project_name is None:
            return None
        for suffix in (".yaml", ".yml"):
            path = os.path.join(
                cls._packaged_xtb_settings_dir(), f"{project_name}{suffix}"
            )
            settings = cls._from_projects_manager(
                XTBProjectSettingsManager(path)
            )
            if settings is not None:
                return settings
        return None


class YamlXTBProjectSettings(XTBProjectSettings):
    """Immutable-by-copy view of a validated xTB project YAML."""

    PROJECT_NAME = "yaml"

    def __init__(
        self,
        sp_settings,
        opt_settings,
        hess_settings,
        explicit_fields=None,
    ):
        self._sp_settings = sp_settings
        self._opt_settings = opt_settings
        self._hess_settings = hess_settings
        self._explicit_fields = {
            jobtype: frozenset(fields)
            for jobtype, fields in (explicit_fields or {}).items()
        }

    def explicit_fields(self, jobtype):
        if jobtype not in ("sp", "opt", "hess"):
            raise ValueError(f"Unknown xTB jobtype: {jobtype!r}.")
        return self._explicit_fields.get(jobtype, frozenset())

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
    """Fail-closed loader for the exact three-section xTB YAML dialect."""

    SECTIONS = ("sp", "opt", "hess")

    def __init__(self, filename):
        self.filename = os.path.abspath(filename)
        self._config = None

    def build(self):
        config = self._read_config()
        unknown_sections = sorted(set(config) - set(self.SECTIONS))
        if unknown_sections:
            raise ValueError(
                "Unknown xTB project section(s): "
                + ", ".join(unknown_sections)
            )
        resolved = {}
        explicit_fields = {}
        for jobtype in self.SECTIONS:
            resolved[jobtype], explicit_fields[jobtype] = (
                self._settings_for_job(jobtype)
            )
        project_settings = YamlXTBProjectSettings(
            sp_settings=resolved["sp"],
            opt_settings=resolved["opt"],
            hess_settings=resolved["hess"],
            explicit_fields=explicit_fields,
        )
        project_settings.PROJECT_NAME = self._parse_project_name()
        project_settings.source_file = self.filename
        return project_settings

    def _read_config(self):
        if self._config is None:
            with open(self.filename) as handle:
                config = yaml.safe_load(handle)
            if config is None:
                config = {}
            if not isinstance(config, Mapping):
                raise TypeError("An xTB project YAML must contain a mapping.")
            self._config = dict(config)
        return self._config

    def _settings_for_job(self, jobtype):
        raw = self._read_config().get(jobtype)
        if raw is None:
            logger.warning(
                f"No xTB configuration found for {jobtype} in "
                f"{self.filename}; using ChemSmart defaults."
            )
            config = {}
        elif not isinstance(raw, Mapping):
            raise TypeError(
                f"xTB project section {jobtype!r} must be a mapping or null."
            )
        else:
            config = dict(raw)

        explicit_fields = frozenset(config)
        explicit_state = explicit_fields & {"charge", "multiplicity"}
        if explicit_state and explicit_state != {"charge", "multiplicity"}:
            raise ValueError(
                f"xTB project section {jobtype!r} must declare charge and "
                "multiplicity together, or omit both to preserve source "
                "electronic state."
            )
        if jobtype != "opt" and "optimization_level" in explicit_fields:
            raise ValueError(
                "xTB optimization_level is valid only in the opt project "
                f"section, not {jobtype!r}."
            )

        declared_jobtype = config.get("jobtype")
        if declared_jobtype is not None:
            if not isinstance(declared_jobtype, str):
                raise TypeError(
                    f"xTB {jobtype}.jobtype must be a string when provided."
                )
            if declared_jobtype.strip().lower() != jobtype:
                raise ValueError(
                    f"xTB project section {jobtype!r} declares incompatible "
                    f"jobtype {declared_jobtype!r}."
                )
        config["jobtype"] = jobtype
        settings = XTBJobSettings.from_dict(config).validate(
            expected_jobtype=jobtype
        )
        return settings, explicit_fields

    def _parse_project_name(self):
        return os.path.basename(self.filename).removesuffix(".yaml").removesuffix(
            ".yml"
        )


class XTBProjectSettingsManager:
    """Create a validated xTB project-settings object from a YAML file."""

    def __init__(self, filename):
        if filename is None:
            raise ValueError("filename is not specified")
        self.filename = os.path.abspath(filename)

    def create(self):
        return YamlXTBProjectSettings.from_yaml(self.filename)
