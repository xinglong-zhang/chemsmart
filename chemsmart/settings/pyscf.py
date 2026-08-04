"""PySCF project settings with a canonical stage-specific representation.

New project files use explicit ``sp``, ``opt``, ``hess`` and experimental
``td`` sections.  The
historical ChemSmart ``gas``/``solv`` dialect remains a supported migration
input: ``gas`` feeds ``opt`` and ``hess`` while ``solv`` feeds ``sp``.  A
project may temporarily contain both representations, but duplicate
definitions must describe the same effective calculation.

Only the stage requested by the CLI is materialized and validated.  This is
important for practical projects: a valid single-point section must not be
rejected because a future optimization section is incomplete.  Rendering is
always canonical and therefore never emits ``gas`` or ``solv``.
"""

import hashlib
import logging
import os
from collections.abc import Mapping

import yaml

from chemsmart.jobs.pyscf.settings import PySCFJobSettings
from chemsmart.settings.user import CHEMSMARTUserSettings
from chemsmart.utils.mixins import RegistryMixin

logger = logging.getLogger(__name__)

PYSCF_JOBTYPES = ("sp", "opt", "hess", "td")
PYSCF_MIGRATION_SECTIONS = ("gas", "solv")
PYSCF_ALLOWED_SECTIONS = PYSCF_JOBTYPES + PYSCF_MIGRATION_SECTIONS
PYSCF_STAGE_SOURCES = {
    "sp": ("sp", "solv"),
    "opt": ("opt", "gas"),
    "hess": ("hess", "gas"),
    "td": ("td",),
}

# Loader-owned defaults that materially affect the generated PySCF script.
# Canonical project output records these values even when the source proposal
# omitted them, so the promoted YAML is a complete view of the settings the
# host will apply.  Method and basis remain mandatory user/model proposals;
# charge and multiplicity remain owned by the input artifact unless both were
# explicitly supplied.
PYSCF_CANONICAL_EFFECTIVE_FIELDS = (
    "density_fit",
    "engine",
    "freq",
    "scf_maxiter",
    "scf_tol",
)
PYSCF_CANONICAL_OPT_EFFECTIVE_FIELDS = (
    "opt_maxsteps",
    "opt_solver",
)


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
        settings.freq = True
        return settings

    def td_settings(self):
        """Return preview-only vertical-excitation settings."""

        settings = self.main_settings().copy()
        settings.jobtype = "td"
        return settings

    def explicit_fields(self, jobtype):
        """Return project fields explicitly supplied for *jobtype*.

        Programmatic project classes have no YAML provenance, so their base
        implementation reports no explicit fields.  YAML-backed projects
        override this method and use it to preserve source electronic state.
        """

        if jobtype not in PYSCF_JOBTYPES:
            raise ValueError(f"Unknown PySCF jobtype: {jobtype!r}.")
        return frozenset()

    @classmethod
    def from_yaml(cls, filename):
        """Load an explicit PySCF project artifact through the strict loader."""

        return YamlPySCFProjectSettingsBuilder(filename=filename).build()

    @classmethod
    def from_project(cls, project):
        """Resolve an exact file, user project, or packaged project.

        PySCF has no scientifically complete default. Omission and unresolved
        names therefore fail closed rather than silently selecting a fixture.
        """
        if project is not None:
            project = os.fspath(project)
            if os.path.isfile(project):
                return PySCFProjectSettingsManager(filename=project).create()

        for resolver in (
            cls._from_user_project_name,
            cls._from_packaged_project_name,
        ):
            project_settings = resolver(project)
            if project_settings is not None:
                return project_settings

        templates_path = os.path.join(os.path.dirname(__file__), "templates")
        raise FileNotFoundError(
            f"No PySCF project settings implemented for {project!r}.\n\n"
            "Pass an explicit YAML path or place <name>.yaml in "
            f"{cls._user_pyscf_settings_dir()}.\n\n"
            f"Templates are available at {templates_path}.\n\n"
            f"Available projects: {cls.available_projects()}"
        )

    @classmethod
    def _user_pyscf_settings_dir(cls):
        return os.path.join(
            CHEMSMARTUserSettings.resolve_config_dir(), "pyscf"
        )

    @classmethod
    def _packaged_pyscf_settings_dir(cls):
        return os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "templates",
            ".chemsmart",
            "pyscf",
        )

    @classmethod
    def available_projects(cls):
        names = set()
        for directory in (
            cls._user_pyscf_settings_dir(),
            cls._packaged_pyscf_settings_dir(),
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
                cls._user_pyscf_settings_dir(), f"{project_name}{suffix}"
            )
            settings = cls._from_projects_manager(
                PySCFProjectSettingsManager(path)
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
                cls._packaged_pyscf_settings_dir(),
                f"{project_name}{suffix}",
            )
            settings = cls._from_projects_manager(
                PySCFProjectSettingsManager(path)
            )
            if settings is not None:
                return settings
        return None


class YamlPySCFProjectSettings(PySCFProjectSettings):
    """Lazy, immutable-by-copy view of a PySCF project YAML."""

    PROJECT_NAME = "yaml"

    def __init__(self, builder):
        self._builder = builder
        self._settings_cache = {}
        self._explicit_fields_cache = {}

    def _settings_for_job(self, jobtype):
        if jobtype not in self._settings_cache:
            settings, explicit_fields = self._builder.settings_for_job(jobtype)
            self._settings_cache[jobtype] = settings
            self._explicit_fields_cache[jobtype] = explicit_fields
        return self._settings_cache[jobtype].copy()

    def explicit_fields(self, jobtype):
        self._settings_for_job(jobtype)
        return self._explicit_fields_cache[jobtype]

    def sp_settings(self):
        return self._settings_for_job("sp")

    def opt_settings(self):
        return self._settings_for_job("opt")

    def hess_settings(self):
        return self._settings_for_job("hess")

    def td_settings(self):
        return self._settings_for_job("td")

    def canonical_sections(self, jobtypes=None):
        """Return stage-keyed settings suitable for canonical YAML output."""

        return self._builder.canonical_sections(jobtypes=jobtypes)

    def render_canonical_yaml(self, jobtypes=None):
        """Render deterministic stage sections, never migration aliases."""

        return yaml.safe_dump(
            self.canonical_sections(jobtypes=jobtypes),
            sort_keys=True,
            default_flow_style=False,
            allow_unicode=True,
        )

    @classmethod
    def from_yaml(cls, filename):
        return YamlPySCFProjectSettingsBuilder(filename=filename).build()


class YamlPySCFProjectSettingsBuilder:
    """Build a lazy dual-dialect :class:`YamlPySCFProjectSettings`."""

    def __init__(self, filename=None, *, sections=None):
        if (filename is None) == (sections is None):
            raise ValueError(
                "Provide exactly one PySCF project filename or section map."
            )
        self.filename = (
            os.path.abspath(filename) if filename is not None else None
        )
        self._input_sections = (
            {str(name): dict(values) for name, values in sections.items()}
            if sections is not None
            else None
        )
        self._config = None
        self._source_bytes = None
        self._resolved = {}

    def build(self):
        # Validate only the document envelope here.  Scientific settings are
        # resolved lazily by the invoked stage.
        self._read_config()
        project_settings = YamlPySCFProjectSettings(builder=self)
        project_settings.PROJECT_NAME = self._parse_project_name()
        project_settings.source_file = self.filename
        return project_settings

    def _read_config(self):
        if self._config is not None:
            return self._config
        if self.filename is not None:
            with open(self.filename, "rb") as handle:
                self._source_bytes = handle.read()
            config = yaml.safe_load(self._source_bytes)
        else:
            config = self._input_sections
            self._source_bytes = yaml.safe_dump(
                config,
                sort_keys=True,
                default_flow_style=False,
                allow_unicode=True,
            ).encode("utf-8")
        if not isinstance(config, Mapping):
            raise ValueError("PySCF project YAML must contain a mapping.")
        unknown = sorted(set(config).difference(PYSCF_ALLOWED_SECTIONS))
        if unknown:
            raise ValueError(
                "Unsupported PySCF top-level section(s): "
                f"{', '.join(map(str, unknown))}; use canonical sp, opt, "
                "hess, and td sections or migration-only gas and solv "
                "sections."
            )
        if not any(name in config for name in PYSCF_ALLOWED_SECTIONS):
            raise ValueError(
                "PySCF project YAML contains no supported calculation section."
            )
        for name in PYSCF_ALLOWED_SECTIONS:
            if name in config and not isinstance(config[name], Mapping):
                raise ValueError(
                    f"PySCF section {name!r} must contain a settings mapping."
                )
        self._config = {
            name: dict(config[name])
            for name in PYSCF_ALLOWED_SECTIONS
            if name in config
        }
        return self._config

    def available_jobtypes(self):
        """Return stages that have either a canonical or migration source."""

        config = self._read_config()
        return tuple(
            jobtype
            for jobtype in PYSCF_JOBTYPES
            if any(source in config for source in PYSCF_STAGE_SOURCES[jobtype])
        )

    @staticmethod
    def _semantic_settings(settings):
        return {
            key: value
            for key, value in vars(settings).items()
            if not key.startswith("_")
        }

    def _candidate_for_section(self, section, jobtype):
        raw = dict(self._read_config()[section])
        explicit_fields = frozenset(raw)
        explicit_state = explicit_fields & {"charge", "multiplicity"}
        if explicit_state and explicit_state != {"charge", "multiplicity"}:
            raise ValueError(
                f"PySCF section {section!r} must declare charge and "
                "multiplicity together, or omit both to preserve the source "
                "electronic state."
            )

        declared_jobtype = raw.get("jobtype")
        if declared_jobtype is not None:
            if not isinstance(declared_jobtype, str):
                raise TypeError(
                    f"PySCF {section}.jobtype must be a string when provided."
                )
            if declared_jobtype.strip().lower() != jobtype:
                raise ValueError(
                    f"PySCF section {section!r} cannot provide {jobtype!r}; "
                    f"it declares jobtype {declared_jobtype!r}."
                )

        raw["jobtype"] = jobtype
        if jobtype == "opt":
            if section == "opt" and raw.get("freq") is True:
                raise ValueError(
                    "PySCF canonical opt section requires freq=False; add a "
                    "separate hess stage so the optimized-geometry handoff "
                    "remains explicit."
                )
            # Historical gas.freq requested an in-process Hessian.  Migration
            # keeps the settings readable but materializes the modern explicit
            # opt node; the hess node is available independently from gas.
            raw["freq"] = False
        elif jobtype == "hess":
            if section == "hess" and raw.get("freq") is False:
                raise ValueError(
                    "PySCF canonical hess section cannot declare freq=False."
                )
            # A standalone hess command is itself the explicit frequency
            # request.  In migration ``gas.freq`` only governed whether opt
            # appended a Hessian, so it must not disable the hess command.
            raw["freq"] = True

        settings = PySCFJobSettings.from_dict(raw)
        settings._project_yaml_digest = self._project_yaml_digest()
        settings.validate()
        canonical_raw = dict(raw)
        canonical_raw.pop("jobtype", None)
        return settings, explicit_fields, canonical_raw

    def _resolve_job(self, jobtype):
        if jobtype not in PYSCF_JOBTYPES:
            raise ValueError(f"Unknown PySCF jobtype: {jobtype!r}.")
        if jobtype in self._resolved:
            settings, explicit_fields, canonical_raw = self._resolved[jobtype]
            return settings.copy(), explicit_fields, dict(canonical_raw)

        config = self._read_config()
        source_names = [
            source
            for source in PYSCF_STAGE_SOURCES[jobtype]
            if source in config
        ]
        if not source_names:
            raise RuntimeError(
                f"PySCF settings for job {jobtype} cannot be found!\n"
                f"Available PySCF jobs with settings are: "
                f"{list(self.available_jobtypes())}"
            )

        candidates = {
            source: self._candidate_for_section(source, jobtype)
            for source in source_names
        }
        if len(candidates) == 2:
            canonical_name, migration_name = PYSCF_STAGE_SOURCES[jobtype]
            canonical = candidates[canonical_name][0]
            migration = candidates[migration_name][0]
            canonical_semantics = self._semantic_settings(canonical)
            migration_semantics = self._semantic_settings(migration)
            if canonical_semantics != migration_semantics:
                differing = sorted(
                    key
                    for key in set(canonical_semantics)
                    | set(migration_semantics)
                    if canonical_semantics.get(key)
                    != migration_semantics.get(key)
                )
                raise ValueError(
                    f"Conflicting PySCF definitions for {jobtype!r}: "
                    f"canonical section {canonical_name!r} and migration "
                    f"section {migration_name!r} differ in "
                    f"{', '.join(differing)}."
                )

        selected_name = (
            jobtype if jobtype in candidates else source_names[0]
        )
        selected_settings, _, canonical_raw = candidates[selected_name]
        explicit_fields = frozenset().union(
            *(item[1] for item in candidates.values())
        )
        self._resolved[jobtype] = (
            selected_settings.copy(),
            explicit_fields,
            dict(canonical_raw),
        )
        return selected_settings.copy(), explicit_fields, dict(canonical_raw)

    def settings_for_job(self, jobtype):
        settings, explicit_fields, _ = self._resolve_job(jobtype)
        return settings, explicit_fields

    def canonical_sections(self, jobtypes=None):
        """Materialize canonical, loader-effective stage sections.

        When *jobtypes* is omitted only stages represented by the source file
        are emitted.  Supplying a subset preserves lazy stage validation.
        Explicit values are retained, while scientifically relevant loader
        defaults are added by the host so a promoted project does not hide
        settings that will affect the generated script.
        """

        if jobtypes is None:
            jobtypes = self.available_jobtypes()
        elif isinstance(jobtypes, str):
            jobtypes = (jobtypes,)
        else:
            jobtypes = tuple(jobtypes)

        unknown = sorted(set(jobtypes).difference(PYSCF_JOBTYPES))
        if unknown:
            raise ValueError(
                "Unknown PySCF canonical stage(s): " + ", ".join(unknown)
            )
        rendered = {}
        for jobtype in PYSCF_JOBTYPES:
            if jobtype not in jobtypes:
                continue
            settings, _, canonical_raw = self._resolve_job(jobtype)
            fields = PYSCF_CANONICAL_EFFECTIVE_FIELDS
            if jobtype == "opt":
                fields += PYSCF_CANONICAL_OPT_EFFECTIVE_FIELDS
            for field in fields:
                canonical_raw.setdefault(field, getattr(settings, field))
            rendered[jobtype] = canonical_raw
        return rendered

    def _project_yaml_digest(self):
        """Return the SHA-256 identity of the resolved project artifact."""
        self._read_config()
        return hashlib.sha256(self._source_bytes).hexdigest()

    def _parse_project_name(self):
        if self.filename is None:
            return "materialized"
        return os.path.splitext(os.path.basename(self.filename))[0]


def materialize_canonical_project_sections(sections):
    """Return host-owned canonical PySCF sections from a typed candidate.

    This is the project-renderer hook used by the agent.  It shares the exact
    dual-dialect loader and validation path with file-backed projects, but is
    compute-free and does not create an intermediate file.
    """

    if not isinstance(sections, Mapping):
        raise ValueError("PySCF project sections must contain a mapping.")
    return YamlPySCFProjectSettingsBuilder(
        sections=sections
    ).canonical_sections()


class PySCFProjectSettingsManager:
    """Load PySCF project settings from a YAML file."""

    def __init__(self, filename):
        if filename is None:
            raise ValueError("filename is not specified")
        self.filename = os.path.abspath(filename)

    def create(self):
        return YamlPySCFProjectSettings.from_yaml(self.filename)
