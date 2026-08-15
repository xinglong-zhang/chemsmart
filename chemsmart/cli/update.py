import copy
import io
import logging
import os
import re
import shutil
import subprocess
import tempfile
from abc import ABC, abstractmethod
from collections.abc import Mapping, MutableMapping, Sequence
from dataclasses import dataclass
from importlib import resources
from pathlib import Path

import click
import tomlkit
import yaml

from chemsmart.utils.logger import create_logger

# Initialize logger
logger = logging.getLogger(__name__)
create_logger(debug=True, stream=True)

RESERVED_TOP_LEVEL_FIELDS = {"SERVER"}


class _ExplicitYamlNull:
    """Marker serialized as an explicit, unquoted YAML ``Null`` value."""


def _represent_explicit_yaml_null(representer, _value):
    """Represent the interactive null marker as a YAML null scalar."""
    return representer.represent_scalar("tag:yaml.org,2002:null", "Null")


class ConfigUpdateError(Exception):
    """Raised when a config update cannot be planned or written safely."""


class ProjectUpdateError(Exception):
    """Raised when project templates cannot be added safely."""


@dataclass(frozen=True)
class ConfigUpdateReport:
    """Result of processing one server YAML file.

    Attributes:
        path: Server YAML path that was checked.
        added_programs: Top-level program sections copied from the template.
        changed: Whether the in-memory YAML content changed.
    """

    path: Path
    added_programs: tuple[str, ...] = ()
    changed: bool = False

    @property
    def is_up_to_date(self) -> bool:
        """Return True when the file needed no update."""
        return not self.changed and not self.added_programs


@dataclass
class _PreparedConfig:
    """User config paired with the bundled template before writing.

    Attributes:
        path: Server YAML path being updated.
        original_text: Raw file content before any in-memory changes.
        data: Round-trip YAML mapping loaded from the user file.
        template_data: Round-trip YAML mapping loaded from ``server.yaml``.
        missing_programs: Template program sections absent from the user file.
    """

    path: Path
    original_text: str
    data: MutableMapping
    template_data: Mapping
    missing_programs: tuple[str, ...] = ()


class Updater(ABC):
    """Base class for CHEMSMART update operations."""

    @abstractmethod
    def update(self):
        """Perform the configured update. Subclasses can implement."""
        pass


class DependencyUpdater(Updater):
    """
    Handles dependency updates in pyproject.toml.
    """

    def __init__(self):
        self._package_path = Path(__file__).resolve().parent.parent.parent
        self._pyproject_path = self._package_path / "pyproject.toml"

    def update(self):
        return self.update_pyproject_toml()

    @property
    def package_path(self) -> Path:
        logger.debug(f"Package path: {self._package_path}")
        return self._package_path

    @property
    def pyproject_path(self) -> Path:
        return self._pyproject_path

    def update_pyproject_toml(self):
        """
        Update dependencies in pyproject.toml.
        """
        requirements_path = self._generate_requirements()
        if requirements_path:
            self._update_toml(requirements_path)

    def _generate_requirements(self, ignore_dirs=None) -> Path:
        """
        Run pipreqs and extract dependencies.
        """
        ignore_dirs = ignore_dirs or [".git", ".github"]
        ignore_arg = ",".join(ignore_dirs)

        with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
            tmp_path = Path(tmp_file.name)
            logger.debug(f"Temporary requirements file: {tmp_path}")

        cmd = [
            "pipreqs",
            str(self.package_path),
            "--force",
            "--ignore",
            ignore_arg,
            "--savepath",
            str(tmp_path),
        ]
        logger.info(f"Running pipreqs: {' '.join(cmd)}")

        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            logger.error(f"Error running pipreqs: {stderr.decode()}")
            return None
        logger.info(stdout.decode())
        return tmp_path

    def _get_existing_dependencies(self):
        """
        Extract dependencies from pyproject.toml and environment.yml.
        """
        dependencies = set()
        dependencies.update(self._get_existing_dependencies_from_pyproject())
        dependencies.update(
            self._get_existing_dependencies_from_environment_yml()
        )
        return dependencies

    def _get_existing_dependencies_from_pyproject(self):
        """
        Extract dependencies from pyproject.toml.
        """
        if not self.pyproject_path.exists():
            logger.error("pyproject.toml not found.")
            return set()

        with self.pyproject_path.open("r") as f:
            pyproject_toml = tomlkit.parse(f.read())

        # Safely handle missing 'project' or 'dependencies' keys
        return set(pyproject_toml.get("project", {}).get("dependencies", []))

    def _get_existing_dependencies_from_environment_yml(self):
        """
        Extract dependencies from environment.yml.
        """
        env_yml_path = Path("environment.yml")
        if not env_yml_path.exists():
            logger.info(
                "environment.yml not found, skipping Conda dependencies."
            )
            return set()

        with env_yml_path.open("r") as f:
            env_data = yaml.safe_load(f)

        dependencies = set()
        conda_deps = env_data.get("dependencies", [])
        for dep in conda_deps:
            if isinstance(dep, str):
                dependencies.add(dep)
            elif isinstance(dep, dict) and "pip" in dep:
                dependencies.update(dep["pip"])

        return dependencies

    def _get_missing_dependencies(self, requirements_path):
        """
        Identify dependencies missing from pyproject.toml, preserving
        versions.
        """
        if not requirements_path.exists():
            logger.error("Generated requirements file not found.")
            return set()

        # Read detected dependencies from pipreqs output (with versions)
        with requirements_path.open("r") as f:
            detected_deps = {line.strip() for line in f if line.strip()}

        # Extract existing dependencies
        existing_deps = self._get_existing_dependencies()

        # Package name normalization mapping (for comparison)
        package_mapping = {
            "pymol-open-source": "pymol",  # Treat pymol-open-source as
            # pymol for matching
        }
        reverse_mapping = {
            v: k for k, v in package_mapping.items()
        }  # For display

        # Function to extract package name (for comparison)
        # Normalizes according to PEP 503: lowercase and replace [-_.] with -
        def extract_pkg_name(dep):
            pkg_name = re.split(r"[=<>~!]", dep, maxsplit=1)[0].strip()
            # PEP 503 normalization: lowercase first,
            # then replace runs of [-_.] with single dash
            pkg_name = pkg_name.lower()
            pkg_name = re.sub(r"[-_.]+", "-", pkg_name)
            return package_mapping.get(pkg_name, pkg_name)

        # Normalize package names for comparison
        detected_pkgs = {extract_pkg_name(dep) for dep in detected_deps}
        existing_pkgs = {extract_pkg_name(dep) for dep in existing_deps}

        # Create display version of detected_pkgs with reverse mapping
        detected_pkgs_display = {
            reverse_mapping.get(pkg, pkg) for pkg in detected_pkgs
        }
        logger.debug(f"Detected packages: {detected_pkgs_display}")
        logger.debug(f"Existing packages: {existing_pkgs}")

        # Identify missing dependencies with their full strings
        missing_dependencies = {
            dep
            for dep in detected_deps
            if extract_pkg_name(dep) not in existing_pkgs
        }
        return missing_dependencies

    def _update_toml(self, requirements_path):
        """
        Update pyproject.toml with missing dependencies, including versions,
        in lowercase.
        """
        missing_dependencies = self._get_missing_dependencies(
            requirements_path
        )
        if not missing_dependencies:
            logger.info("No new dependencies to add.")
            return

        with self.pyproject_path.open("r") as f:
            pyproject_toml = tomlkit.parse(f.read())

        # Function to normalize package name to lowercase while preserving
        # version
        def normalize_dep(dep):
            pkg_name, *version = re.split(r"([=<>~!])", dep, maxsplit=1)
            if version:  # If there's a version specifier
                return pkg_name.lower() + "".join(version)
            return pkg_name.lower()  # No version, just lowercase the name

        # Adjust and normalize dependencies
        adjusted_deps = {
            normalize_dep(dep).replace("==", ">=")
            for dep in missing_dependencies
        }
        for dep in adjusted_deps:
            pyproject_toml["project"]["dependencies"].append(dep)

        with self.pyproject_path.open("w") as f:
            f.write(tomlkit.dumps(pyproject_toml))

        logger.info(f"Added missing dependencies: {', '.join(adjusted_deps)}")


class VersionUpdater(Updater):
    """Handles package version updates across project metadata files."""

    def __init__(self, version_number: str):
        package_path = Path(__file__).resolve().parent.parent.parent
        self.version_number = version_number
        self._pyproject_path = package_path / "pyproject.toml"
        self._version_file_path = package_path / "chemsmart" / "VERSION"
        self._docs_conf_file_path = (
            package_path / "docs" / "source" / "conf.py"
        )

    def update(self):
        return self.update_version_number(self.version_number)

    @property
    def pyproject_path(self) -> Path:
        return self._pyproject_path

    @property
    def version_file_path(self) -> Path:
        return self._version_file_path

    @property
    def docs_conf_file_path(self) -> Path:
        return self._docs_conf_file_path

    def update_version_number(self, version_number):
        """Update version number in chemsmart/VERSION, pyproject.toml,
        and docs/source/conf.py"""
        from chemsmart.utils.repattern import release_pattern, version_pattern

        logger.info(f"Updating to version number: {version_number}")

        # 1) Update chemsmart/VERSION
        if not self.version_file_path.exists():
            logger.error(
                f"Version file not found at {self.version_file_path}."
            )
        else:
            logger.debug(
                f"Writing version {version_number} to {self.version_file_path}"
            )
            self.version_file_path.write_text(
                f"{version_number}\n", encoding="utf-8"
            )

        # 2) Update pyproject.toml
        # version = "0.1.9"
        if not self.pyproject_path.exists():
            logger.error(f"pyproject.toml not found at: {self.pyproject_path}")
        else:
            logger.debug(f"Updating pyproject.toml at {self.pyproject_path}")
            pyproject_text = self.pyproject_path.read_text(encoding="utf-8")

            # Replace: version = "0.1.9" (preserve spacing and comments)
            version_replacement = r"\g<1>" + version_number + r"\g<2>"

            new_pyproject_text, n_subs = re.subn(
                version_pattern,
                version_replacement,
                pyproject_text,
                count=1,
            )

            if n_subs == 0:
                logger.warning(
                    "No 'version = \"...\"' assignment found in pyproject.toml; "
                    "project version not updated."
                )
            else:
                self.pyproject_path.write_text(
                    new_pyproject_text, encoding="utf-8"
                )
                logger.info(
                    f"Updated project version to {version_number} "
                    f"in {self.pyproject_path}"
                )

        # 3) Update docs/source/conf.py: release = "x.y.z"
        if not self.docs_conf_file_path.exists():
            logger.warning(
                f"Sphinx conf.py not found at: {self.docs_conf_file_path} "
                "- skipping docs version update."
            )
        else:
            logger.debug(
                f"Updating Sphinx release in {self.docs_conf_file_path}"
            )
            text = self.docs_conf_file_path.read_text(encoding="utf-8")

            # Replace: release = "0.1.9"
            replacement = r"\g<1>" + version_number + r"\g<2>"

            new_text, n_subs = re.subn(
                release_pattern, replacement, text, count=1, flags=re.MULTILINE
            )

            if n_subs == 0:
                logger.warning(
                    "No 'release = \"...\"' assignment found in conf.py; "
                    "Sphinx version not updated."
                )
            else:
                self.docs_conf_file_path.write_text(new_text, encoding="utf-8")
                logger.info(
                    f"Updated Sphinx release to {version_number} "
                    f"in {self.docs_conf_file_path}"
                )

        logger.info(f"Version update to {version_number} completed.")


class ProjectsUpdater(Updater):
    """Add missing project template directories and files."""

    def update(self) -> tuple[str, ...]:
        """Copy missing project template content into the user config.

        Existing files are never overwritten. Top-level template directories
        are discovered automatically, except ``server``, which is managed by
        :class:`ConfigurationUpdater`.

        Returns:
            Names of project directories that received missing content.
        """
        config_dir = self._config_dir()
        updated_projects: list[str] = []

        try:
            self._ensure_directory(config_dir)
            with resources.as_file(
                self._template_config_dir()
            ) as template_dir:
                project_templates = sorted(
                    (
                        path
                        for path in template_dir.iterdir()
                        if path.name != "server" and path.is_dir()
                    ),
                    key=lambda path: path.name.casefold(),
                )
                for project_template in project_templates:
                    destination = config_dir / project_template.name
                    if self._copy_missing(project_template, destination):
                        updated_projects.append(project_template.name)
        except OSError as exc:
            message = exc.strerror or str(exc)
            raise ProjectUpdateError(
                f"Could not update project templates: {message}"
            ) from exc

        return tuple(updated_projects)

    @staticmethod
    def _config_dir() -> Path:
        """Return the resolved user configuration directory."""
        from chemsmart.settings.user import CHEMSMARTUserSettings

        configured_dir = CHEMSMARTUserSettings.resolve_config_dir()
        return Path(os.path.abspath(configured_dir))

    @staticmethod
    def _template_config_dir():
        """Return the bundled ``.chemsmart`` template directory."""
        return (
            resources.files("chemsmart.settings") / "templates" / ".chemsmart"
        )

    @staticmethod
    def _ensure_directory(path: Path) -> None:
        """Create a missing directory and reject unsafe path conflicts."""
        if path.is_symlink():
            raise ProjectUpdateError(f"Refusing to update symlink: {path}")
        if path.exists() and not path.is_dir():
            raise ProjectUpdateError(f"Expected a directory: {path}")
        path.mkdir(parents=True, exist_ok=True)

    def _copy_missing(self, source: Path, destination: Path) -> bool:
        """Recursively copy missing content without replacing destinations."""
        if source.is_symlink() or destination.is_symlink():
            raise ProjectUpdateError(
                f"Refusing to update symlink: {destination}"
            )

        if source.is_dir():
            changed = not destination.exists()
            self._ensure_directory(destination)
            for child in sorted(
                source.iterdir(), key=lambda path: path.name.casefold()
            ):
                changed = (
                    self._copy_missing(child, destination / child.name)
                    or changed
                )
            return changed

        if not source.is_file():
            raise ProjectUpdateError(
                f"Unsupported project template entry: {source}"
            )
        if destination.exists():
            if not destination.is_file():
                raise ProjectUpdateError(f"Expected a file: {destination}")
            return False

        created = False
        try:
            with (
                source.open("rb") as source_file,
                destination.open("xb") as destination_file,
            ):
                created = True
                shutil.copyfileobj(source_file, destination_file)
            shutil.copystat(source, destination)
        except BaseException:
            if created:
                try:
                    destination.unlink()
                except FileNotFoundError:
                    pass
            raise
        return True


class ConfigurationUpdater(Updater):
    """Handles server YAML configuration updates."""

    def __init__(self, server: str | Sequence[str] | None = None):
        if server is None:
            self.servers = ()
        elif isinstance(server, str):
            self.servers = (server,)
        else:
            self.servers = tuple(server)

    def update(self) -> list[ConfigUpdateReport]:
        """Update server YAML files with missing top-level program sections.

        When ``servers`` is empty, all existing ``*.yaml`` files in the
        resolved server config directory are scanned. Otherwise, only the
        selected existing files are checked; each name may be passed with or
        without ``.yaml``. Missing server YAML files are not created.

        In an interactive terminal, one optional ``EXEFOLDER`` override is
        requested for each unique missing program whose template defines that
        field. Blank input keeps the bundled template value. Overrides are
        applied only to newly copied program sections.
        """
        self._require_yaml_dependency()
        server_dir = self._server_config_dir()
        yaml_files = self._select_server_yaml_files(server_dir, self.servers)
        template_data = self._load_server_template()

        prepared_configs: list[_PreparedConfig] = []
        for yaml_file in yaml_files:
            original_text, data = self._load_user_yaml(yaml_file)
            prepared_configs.append(
                _PreparedConfig(
                    path=yaml_file,
                    original_text=original_text,
                    data=data,
                    template_data=template_data,
                )
            )

        self._plan_missing_programs(prepared_configs)
        program_exefolders = self._prompt_program_exefolders(prepared_configs)
        changed_reports = self._apply_update_plan(
            prepared_configs, program_exefolders
        )
        self._write_changed_files(prepared_configs, changed_reports)
        return changed_reports

    @staticmethod
    def _require_yaml_dependency():
        from ruamel.yaml import YAML
        from ruamel.yaml.error import YAMLError

        return YAML, YAMLError

    @staticmethod
    def _server_config_dir() -> Path:
        """Return the resolved ``server`` config directory."""
        from chemsmart.settings.user import CHEMSMARTUserSettings

        configured_dir = CHEMSMARTUserSettings.resolve_config_dir()
        return Path(os.path.abspath(configured_dir)) / "server"

    def _select_server_yaml_files(
        self, server_dir: Path, servers: Sequence[str]
    ) -> list[Path]:
        """Select existing server YAML files to process."""
        if servers:
            selected_files: list[Path] = []
            seen: set[str] = set()
            for server in servers:
                server_file = self._normalize_server_filename(server)
                if server_file in seen:
                    continue
                server_path = server_dir / server_file
                if not server_path.is_file():
                    raise ConfigUpdateError(
                        f"Server YAML not found in {server_dir}: {server_file}"
                    )
                seen.add(server_file)
                selected_files.append(server_path)
            return selected_files

        if not server_dir.is_dir():
            raise ConfigUpdateError(
                f"Server config directory not found: {server_dir}"
            )

        try:
            return sorted(server_dir.glob("*.yaml"))
        except OSError as exc:
            raise ConfigUpdateError(
                f"Could not read server config directory {server_dir}: "
                f"{self._format_error(exc)}"
            ) from exc

    @staticmethod
    def _normalize_server_filename(server: str) -> str:
        """Normalize and validate a server file name."""
        if not server:
            raise ConfigUpdateError("--server must not be empty.")
        server_path = Path(server)
        if server_path.name != server:
            raise ConfigUpdateError(
                "--server only accepts a YAML file name from the config "
                "server directory, not a path."
            )
        if server_path.suffix and server_path.suffix != ".yaml":
            raise ConfigUpdateError("--server must name a .yaml file.")
        return server if server.endswith(".yaml") else f"{server}.yaml"

    @staticmethod
    def _yaml():
        """Create the round-trip YAML parser/emitter used for config files."""
        from ruamel.yaml import YAML

        yaml_parser = YAML()
        yaml_parser.representer.add_representer(
            _ExplicitYamlNull, _represent_explicit_yaml_null
        )
        yaml_parser.preserve_quotes = True
        yaml_parser.indent(mapping=4, sequence=4, offset=2)
        yaml_parser.width = 4096
        return yaml_parser

    def _load_server_template(self) -> Mapping:
        """Load the bundled ``server.yaml`` update template."""
        from ruamel.yaml.error import YAMLError

        template = (
            resources.files("chemsmart.settings") / "templates" / "server.yaml"
        )
        yaml_parser = self._yaml()
        try:
            data = yaml_parser.load(template.read_text(encoding="utf-8"))
        except YAMLError as exc:
            raise ConfigUpdateError(
                f"Bundled template {template.name} is not valid YAML: {exc}"
            ) from exc
        except (OSError, UnicodeError) as exc:
            raise ConfigUpdateError(
                f"Could not read bundled template {template.name}: "
                f"{self._format_error(exc)}"
            ) from exc
        if not isinstance(data, Mapping):
            raise ConfigUpdateError(
                f"Bundled template {template.name} must contain a mapping."
            )
        return data

    def _load_user_yaml(self, path: Path) -> tuple[str, MutableMapping]:
        """Load a user server YAML file as round-trip text and mapping."""
        from ruamel.yaml.error import YAMLError

        try:
            original_text = path.read_text(encoding="utf-8")
            data = self._yaml().load(original_text)
        except YAMLError as exc:
            raise ConfigUpdateError(
                f"{path.name} is not valid YAML: {exc}"
            ) from exc
        except (OSError, UnicodeError) as exc:
            raise ConfigUpdateError(
                f"Could not read {path.name}: {self._format_error(exc)}"
            ) from exc

        if not isinstance(data, MutableMapping):
            raise ConfigUpdateError(
                f"{path.name} must contain a YAML mapping at the top level."
            )
        return original_text, data

    @staticmethod
    def _template_programs(template_data: Mapping) -> tuple[str, ...]:
        """Return top-level template sections treated as program configs."""
        return tuple(
            str(field)
            for field in template_data.keys()
            if str(field).upper() not in RESERVED_TOP_LEVEL_FIELDS
        )

    def _plan_missing_programs(
        self, prepared_configs: list[_PreparedConfig]
    ) -> None:
        """Record each config's missing template program sections."""
        for prepared in prepared_configs:
            prepared.missing_programs = tuple(
                program
                for program in self._template_programs(prepared.template_data)
                if program not in prepared.data
            )

    @staticmethod
    def _programs_with_exefolder(
        prepared_configs: list[_PreparedConfig],
    ) -> tuple[str, ...]:
        """Return unique missing programs whose template has EXEFOLDER."""
        programs: list[str] = []
        seen: set[str] = set()
        for prepared in prepared_configs:
            for program in prepared.missing_programs:
                section = prepared.template_data[program]
                if (
                    program not in seen
                    and isinstance(section, Mapping)
                    and "EXEFOLDER" in section
                ):
                    seen.add(program)
                    programs.append(program)
        return tuple(programs)

    @staticmethod
    def _is_interactive() -> bool:
        """Return whether stdin is an interactive terminal."""
        return click.get_text_stream("stdin").isatty()

    def _prompt_program_exefolders(
        self, prepared_configs: list[_PreparedConfig]
    ) -> dict[str, str | _ExplicitYamlNull]:
        """Collect one normalized EXEFOLDER override per missing program.

        Case-insensitive ``null`` input becomes an explicit YAML null value;
        all other non-empty input has its user-home marker expanded.
        """
        if not self._is_interactive():
            return {}

        overrides: dict[str, str | _ExplicitYamlNull] = {}
        for program in self._programs_with_exefolder(prepared_configs):
            folder = click.prompt(
                f"{program} EXEFOLDER "
                "(press Enter to use the bundled template value)",
                default="",
                show_default=False,
            ).strip()
            if folder:
                overrides[program] = (
                    _ExplicitYamlNull()
                    if folder.casefold() == "null"
                    else os.path.expanduser(folder)
                )
        return overrides

    def _apply_update_plan(
        self,
        prepared_configs: list[_PreparedConfig],
        program_exefolders: Mapping[str, str | _ExplicitYamlNull],
    ) -> list[ConfigUpdateReport]:
        """Apply planned additions and optional paths to prepared configs."""
        reports: list[ConfigUpdateReport] = []
        for prepared in prepared_configs:
            added_programs: list[str] = []
            for program in prepared.missing_programs:
                section = copy.deepcopy(prepared.template_data[program])
                if (
                    isinstance(section, MutableMapping)
                    and "EXEFOLDER" in section
                ):
                    if program in program_exefolders:
                        section["EXEFOLDER"] = program_exefolders[program]
                    elif section["EXEFOLDER"] is None:
                        section["EXEFOLDER"] = _ExplicitYamlNull()
                prepared.data[program] = section
                added_programs.append(program)
            reports.append(
                ConfigUpdateReport(
                    path=prepared.path,
                    added_programs=tuple(added_programs),
                    changed=bool(added_programs),
                )
            )
        return reports

    def _write_changed_files(
        self,
        prepared_configs: list[_PreparedConfig],
        reports: list[ConfigUpdateReport],
    ) -> None:
        """Write changed configs using per-file atomic replacement."""
        from ruamel.yaml.error import YAMLError

        reports_by_path = {report.path: report for report in reports}
        pending_writes: list[tuple[Path, str]] = []
        for prepared in prepared_configs:
            report = reports_by_path[prepared.path]
            if not report.changed:
                continue
            if prepared.path.is_symlink():
                raise ConfigUpdateError(
                    "Refusing to update symlinked server YAML: "
                    f"{prepared.path}"
                )
            try:
                updated_text = self._dump_yaml(prepared.data)
            except (OSError, UnicodeError, YAMLError) as exc:
                raise ConfigUpdateError(
                    f"Could not update {prepared.path.name}: "
                    f"{self._format_error(exc)}"
                ) from exc
            if updated_text != prepared.original_text:
                pending_writes.append((prepared.path, updated_text))

        for path, updated_text in pending_writes:
            self._atomic_write_text(path, updated_text)

    def _dump_yaml(self, data: Mapping) -> str:
        """Serialize a round-trip YAML mapping to text."""
        stream = io.StringIO()
        self._yaml().dump(data, stream)
        return stream.getvalue()

    def _atomic_write_text(self, path: Path, text: str) -> None:
        """Replace one file with same-directory temporary output."""
        if path.is_symlink():
            raise ConfigUpdateError(
                f"Refusing to update symlinked server YAML: {path}"
            )

        tmp_name: str | None = None
        try:
            stat_result = path.stat()
            fd, tmp_name = tempfile.mkstemp(
                prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
            )
            with os.fdopen(fd, "w", encoding="utf-8") as tmp_file:
                tmp_file.write(text)
            os.chmod(tmp_name, stat_result.st_mode)
            os.replace(tmp_name, path)
        except OSError as exc:
            if tmp_name is not None:
                try:
                    os.unlink(tmp_name)
                except FileNotFoundError:
                    pass
            raise ConfigUpdateError(
                f"Could not update {path.name}: {self._format_error(exc)}"
            ) from exc

    @staticmethod
    def _format_error(exc: BaseException) -> str:
        """Return a concise message for expected filesystem-style errors."""
        if isinstance(exc, OSError) and exc.strerror:
            return exc.strerror
        return str(exc)


@click.group(name="update")
def update():
    """
    Manage updates in the chemsmart package.
    """


@update.command()
def deps():
    """
    Automatically update dependencies in pyproject.toml.
    """
    logger.info("Updating dependencies...")
    DependencyUpdater().update()
    logger.info("Update complete.")


@update.command()
@click.option(
    "-v",
    "--version-number",
    type=str,
    required=True,
    help="Version number to be updated to.",
)
def version(version_number: str):
    """
    Automatically update chemsmart version in chemsmart/VERSION,
    pyproject.toml and docs/source/conf.py.
    """
    logger.info("Updating version number...")
    VersionUpdater(version_number).update()
    logger.info("Update complete.")


@update.command(name="projects")
def projects():
    """Add missing project directories and files from bundled templates."""

    try:
        updated_projects = ProjectsUpdater().update()
    except ProjectUpdateError as exc:
        raise click.ClickException(str(exc)) from exc

    if updated_projects:
        click.echo("Updated project templates:")
        for project in updated_projects:
            click.echo(f"  added missing content: {project}")
    else:
        click.echo("Project templates are already up to date.")


@update.command(name="configs")
@click.option(
    "-s",
    "--server",
    type=str,
    multiple=True,
    help=(
        "Existing server YAML to update, with or without .yaml. "
        "Repeat to select multiple files."
    ),
)
def configs(server: tuple[str, ...]):
    """Update existing server YAML program sections from server.yaml."""

    try:
        reports = ConfigurationUpdater(server).update()
    except ModuleNotFoundError as exc:
        if exc.name in {"ruamel", "ruamel.yaml"}:
            raise click.ClickException(
                "ruamel.yaml is required for 'chemsmart update configs'.\n"
                "Install it with `pip install 'ruamel.yaml>=0.18.16'`\n"
                "or rerun `make env`."
            ) from exc
        raise
    except ConfigUpdateError as exc:
        raise click.ClickException(str(exc)) from exc

    for report in reports:
        click.echo(f"{report.path.name}:")
        if report.added_programs:
            for program in report.added_programs:
                click.echo(f"  added program: {program}")
        if report.is_up_to_date:
            click.echo("  already up to date")


if __name__ == "__main__":
    update()
