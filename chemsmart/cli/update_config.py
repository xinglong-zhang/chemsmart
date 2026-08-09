"""Update existing server YAML files from bundled server templates.

This module adds missing top-level program sections from matched bundled
templates while leaving ``SERVER`` and existing program sections untouched. It
does not recurse into existing program sections, and it uses round-trip YAML
plus same-directory temporary replacement to preserve formatting where possible
and write each changed file safely.
"""

import copy
import io
import os
import tempfile
from collections.abc import Mapping, MutableMapping
from dataclasses import dataclass
from importlib import resources
from pathlib import Path

from ruamel.yaml import YAML
from ruamel.yaml.error import YAMLError

from chemsmart.settings.user import CHEMSMARTUserSettings

RESERVED_TOP_LEVEL_FIELDS = {"SERVER"}
LOCAL_SCHEDULER_VALUES = {"local", "none", "null"}


class ConfigUpdateError(Exception):
    """Raised when a config update cannot be planned or written safely."""


@dataclass(frozen=True)
class ConfigUpdateReport:
    """Result of processing one server YAML file.

    Attributes:
        path: Server YAML path that was checked.
        added_programs: Top-level program sections copied from the template.
        skipped_reason: Reason the file was skipped, if no template matched.
        changed: Whether the in-memory YAML content changed.
    """

    path: Path
    added_programs: tuple[str, ...] = ()
    skipped_reason: str | None = None
    changed: bool = False

    @property
    def is_up_to_date(self) -> bool:
        """Return True when the file needed no update and was not skipped."""
        return (
            not self.changed
            and not self.added_programs
            and self.skipped_reason is None
        )


@dataclass
class _PreparedConfig:
    """User config paired with the matched template before writing.

    Attributes:
        path: Server YAML path being updated.
        original_text: Raw file content before any in-memory changes.
        data: Round-trip YAML mapping loaded from the user file.
        template_data: Round-trip YAML mapping loaded from the matched template.
    """

    path: Path
    original_text: str
    data: MutableMapping
    template_data: Mapping


def update_server_configs(
    server: str | None = None,
) -> list[ConfigUpdateReport]:
    """Update server YAML files with missing top-level program sections.

    When ``server`` is ``None``, all existing ``*.yaml`` files in the resolved
    server config directory are scanned. When ``server`` is provided, only that
    existing file is checked; the name may be passed with or without ``.yaml``.
    Missing server YAML files are not created.

    The update only copies program sections that exist in the matched bundled
    template but are absent from the user YAML. It does not modify ``SERVER``,
    overwrite existing program sections, or inspect/fill nested fields inside
    existing program sections.

    Args:
        server: Optional server YAML file name from the config directory.

    Returns:
        A report for each checked or skipped YAML file.

    Raises:
        ConfigUpdateError: If files, templates, YAML, or writes cannot be
            handled safely.
    """
    server_dir = _server_config_dir()
    yaml_files = _select_server_yaml_files(server_dir, server)
    templates = _load_server_templates()

    prepared_configs: list[_PreparedConfig] = []
    reports: list[ConfigUpdateReport] = []
    for yaml_file in yaml_files:
        original_text, data = _load_user_yaml(yaml_file)
        template_name = _match_template_name(
            yaml_file.name, data, templates, strict=server is not None
        )
        if template_name is None:
            reports.append(
                ConfigUpdateReport(
                    path=yaml_file,
                    skipped_reason="could not match a bundled server template",
                )
            )
            continue
        prepared_configs.append(
            _PreparedConfig(
                path=yaml_file,
                original_text=original_text,
                data=data,
                template_data=templates[template_name],
            )
        )

    changed_reports = _build_update_reports(prepared_configs)
    _write_changed_files(prepared_configs, changed_reports)
    return reports + changed_reports


def _server_config_dir() -> Path:
    """Return the resolved ``server`` config directory."""
    configured_dir = os.environ.get(
        "CHEMSMART_CONFIG_DIR", CHEMSMARTUserSettings.USER_CONFIG_DIR
    )
    return Path(os.path.abspath(os.path.expanduser(configured_dir))) / "server"


def _select_server_yaml_files(
    server_dir: Path, server: str | None
) -> list[Path]:
    """Select existing server YAML files to process.

    Args:
        server_dir: Resolved directory containing user server YAML files.
        server: Optional server file name, with or without ``.yaml``.

    Returns:
        Existing YAML paths to check, sorted for full-directory scans.

    Raises:
        ConfigUpdateError: If the directory cannot be read or the requested
            server file is invalid or missing.
    """
    if server is not None:
        server_file = _normalize_server_filename(server)
        server_path = server_dir / server_file
        if not server_path.is_file():
            raise ConfigUpdateError(
                f"Server YAML not found in {server_dir}: {server_file}"
            )
        return [server_path]

    if not server_dir.is_dir():
        raise ConfigUpdateError(
            f"Server config directory not found: {server_dir}"
        )

    try:
        return sorted(server_dir.glob("*.yaml"))
    except OSError as exc:
        raise ConfigUpdateError(
            f"Could not read server config directory {server_dir}: "
            f"{_format_error(exc)}"
        ) from exc


def _normalize_server_filename(server: str) -> str:
    """Normalize and validate a server file name.

    Args:
        server: User-provided server name or ``.yaml`` file name.

    Returns:
        A file name ending in ``.yaml``.

    Raises:
        ConfigUpdateError: If the value is empty, names a path, or uses a
            non-YAML suffix.
    """
    if not server:
        raise ConfigUpdateError("--server must not be empty.")
    server_path = Path(server)
    if server_path.name != server:
        raise ConfigUpdateError(
            "--server only accepts a YAML file name from the config server "
            "directory, not a path."
        )
    if server_path.suffix and server_path.suffix != ".yaml":
        raise ConfigUpdateError("--server must name a .yaml file.")
    return server if server.endswith(".yaml") else f"{server}.yaml"


def _yaml() -> YAML:
    """Create the round-trip YAML parser/emitter used for config files."""
    yaml = YAML()
    yaml.preserve_quotes = True
    yaml.indent(mapping=4, sequence=4, offset=2)
    yaml.width = 4096
    return yaml


def _load_server_templates() -> dict[str, Mapping]:
    """Load bundled server templates keyed by template file name.

    Returns:
        Mapping from template file name to top-level YAML mapping.

    Raises:
        ConfigUpdateError: If bundled templates cannot be read, parsed, or do
            not contain top-level mappings.
    """
    template_dir = (
        resources.files("chemsmart.settings")
        / "templates"
        / ".chemsmart"
        / "server"
    )
    templates: dict[str, Mapping] = {}
    yaml = _yaml()
    try:
        template_files = list(template_dir.iterdir())
    except OSError as exc:
        raise ConfigUpdateError(
            f"Could not read bundled server templates: {_format_error(exc)}"
        ) from exc

    for template in template_files:
        if not template.is_file() or not template.name.endswith(".yaml"):
            continue
        try:
            data = yaml.load(template.read_text(encoding="utf-8"))
        except YAMLError as exc:
            raise ConfigUpdateError(
                f"Bundled template {template.name} is not valid YAML: {exc}"
            ) from exc
        except (OSError, UnicodeError) as exc:
            raise ConfigUpdateError(
                f"Could not read bundled template {template.name}: "
                f"{_format_error(exc)}"
            ) from exc
        if not isinstance(data, Mapping):
            raise ConfigUpdateError(
                f"Bundled template {template.name} must contain a mapping."
            )
        templates[template.name] = data
    return templates


def _load_user_yaml(path: Path) -> tuple[str, MutableMapping]:
    """Load a user server YAML file as round-trip text and mapping.

    Args:
        path: Existing user server YAML path.

    Returns:
        The original file text and parsed top-level mapping.

    Raises:
        ConfigUpdateError: If the file cannot be read, is invalid YAML, or does
            not contain a top-level mapping.
    """
    try:
        original_text = path.read_text(encoding="utf-8")
        data = _yaml().load(original_text)
    except YAMLError as exc:
        raise ConfigUpdateError(
            f"{path.name} is not valid YAML: {exc}"
        ) from exc
    except (OSError, UnicodeError) as exc:
        raise ConfigUpdateError(
            f"Could not read {path.name}: {_format_error(exc)}"
        ) from exc

    if not isinstance(data, MutableMapping):
        raise ConfigUpdateError(
            f"{path.name} must contain a YAML mapping at the top level."
        )
    return original_text, data


def _match_template_name(
    filename: str,
    data: Mapping,
    templates: Mapping[str, Mapping],
    strict: bool,
) -> str | None:
    """Choose the bundled template for a user YAML file.

    The file name is matched first. Custom names fall back to
    ``SERVER.SCHEDULER`` values: SLURM, PBS, or local/null/none. In non-strict
    scans an unmatched file is skipped by returning ``None``; strict mode raises
    an error for the explicitly selected file.

    Args:
        filename: User YAML file name.
        data: Parsed user YAML mapping.
        templates: Available bundled templates keyed by file name.
        strict: Whether unmatched custom files should raise an error.

    Returns:
        Matched template file name, or ``None`` when skipped in non-strict mode.

    Raises:
        ConfigUpdateError: If no template can be matched in strict mode.
    """
    if filename in templates:
        return filename

    scheduler = _scheduler_value(data)
    if scheduler == "slurm" and "SLURM.yaml" in templates:
        return "SLURM.yaml"
    if scheduler == "pbs" and "PBS.yaml" in templates:
        return "PBS.yaml"
    if scheduler in LOCAL_SCHEDULER_VALUES and "local.yaml" in templates:
        return "local.yaml"

    if strict:
        raise ConfigUpdateError(
            f"Could not match a bundled template for {filename}. "
            "Set SERVER.SCHEDULER to SLURM, PBS, or Null/local, or name the "
            "file after an existing template."
        )
    return None


def _scheduler_value(data: Mapping) -> str | None:
    """Return normalized ``SERVER.SCHEDULER`` or ``None`` if unavailable."""
    server_data = data.get("SERVER")
    if not isinstance(server_data, Mapping):
        return None
    if "SCHEDULER" not in server_data:
        return None
    scheduler = server_data["SCHEDULER"]
    if scheduler is None:
        return "null"
    return str(scheduler).strip().lower()


def _template_programs(template_data: Mapping) -> tuple[str, ...]:
    """Return top-level template sections treated as program configs."""
    return tuple(
        str(field)
        for field in template_data.keys()
        if str(field).upper() not in RESERVED_TOP_LEVEL_FIELDS
    )


def _build_update_reports(
    prepared_configs: list[_PreparedConfig],
) -> list[ConfigUpdateReport]:
    """Copy missing program sections into prepared user configs.

    Args:
        prepared_configs: User YAML mappings paired with matched templates.

    Returns:
        Reports describing the in-memory changes for each prepared file.
    """
    reports: list[ConfigUpdateReport] = []
    for prepared in prepared_configs:
        added_programs: list[str] = []
        for program in _template_programs(prepared.template_data):
            if program not in prepared.data:
                prepared.data[program] = copy.deepcopy(
                    prepared.template_data[program]
                )
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
    prepared_configs: list[_PreparedConfig],
    reports: list[ConfigUpdateReport],
) -> None:
    """Write changed configs using per-file atomic replacement.

    The function first renders changed mappings and collects pending writes,
    then writes each file in sequence. This avoids writing files that are known
    symlinks, but it does not provide a transaction across multiple files.

    Args:
        prepared_configs: User YAML mappings after in-memory updates.
        reports: Reports produced for the prepared configs.

    Raises:
        ConfigUpdateError: If a changed file is a symlink, cannot be rendered,
            or cannot be written safely.
    """
    reports_by_path = {report.path: report for report in reports}
    pending_writes: list[tuple[Path, str]] = []
    for prepared in prepared_configs:
        report = reports_by_path[prepared.path]
        if not report.changed:
            continue
        if prepared.path.is_symlink():
            raise ConfigUpdateError(
                f"Refusing to update symlinked server YAML: {prepared.path}"
            )
        try:
            updated_text = _dump_yaml(prepared.data)
        except (OSError, UnicodeError, YAMLError) as exc:
            raise ConfigUpdateError(
                f"Could not update {prepared.path.name}: {_format_error(exc)}"
            ) from exc
        if updated_text != prepared.original_text:
            pending_writes.append((prepared.path, updated_text))

    for path, updated_text in pending_writes:
        _atomic_write_text(path, updated_text)


def _dump_yaml(data: Mapping) -> str:
    """Serialize a round-trip YAML mapping to text."""
    stream = io.StringIO()
    _yaml().dump(data, stream)
    return stream.getvalue()


def _atomic_write_text(path: Path, text: str) -> None:
    """Replace one file with same-directory temporary output.

    The replacement is atomic for the single target file via ``os.replace``.
    Symlinks are refused only when this function is asked to write them.

    Args:
        path: Target file to replace.
        text: Complete replacement content.

    Raises:
        ConfigUpdateError: If the path is a symlink, a temporary file cannot be
            created, permissions cannot be copied, or replacement fails.
    """
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
            f"Could not update {path.name}: {_format_error(exc)}"
        ) from exc


def _format_error(exc: BaseException) -> str:
    """Return a concise message for expected filesystem-style errors."""
    if isinstance(exc, OSError) and exc.strerror:
        return exc.strerror
    return str(exc)
