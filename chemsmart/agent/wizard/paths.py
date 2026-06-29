"""Path and file-permission helpers for wizard persistence."""

from __future__ import annotations

import os
import re
from pathlib import Path

_SERVER_NAME_PATTERN = re.compile(r"^[a-zA-Z0-9_][a-zA-Z0-9_.-]{0,63}$")


def server_config_dir() -> Path:
    """Return the per-user wizard/server configuration directory."""

    return Path.home() / ".chemsmart" / "server"


def validate_server_name(name: str) -> str:
    """Validate a wizard server basename."""

    if not isinstance(name, str) or not _SERVER_NAME_PATTERN.fullmatch(name):
        raise ValueError(f"Invalid server name: {name!r}")
    return name


def server_yaml_path(name: str) -> Path:
    """Return the contained server YAML path for a validated name."""

    validate_server_name(name)
    target = server_config_dir() / f"{name}.yaml"
    _assert_server_path_contained(target)
    return target


def server_cache_path(name: str) -> Path:
    """Return the contained server cache path for a validated name."""

    validate_server_name(name)
    target = server_config_dir() / f"{name}.cache.json"
    _assert_server_path_contained(target)
    return target


def ensure_private_directory(path: Path) -> None:
    """Create a user-private directory if needed."""

    path.mkdir(mode=0o700, parents=True, exist_ok=True)
    os.chmod(path, 0o700)


def write_private_text(
    target: Path,
    text: str,
    *,
    overwrite: bool,
) -> str:
    """Write UTF-8 text using explicit private permissions."""

    _assert_server_path_contained(target)
    ensure_private_directory(target.parent)

    flags = os.O_WRONLY | os.O_CREAT
    if overwrite:
        flags |= os.O_TRUNC
    else:
        flags |= os.O_EXCL

    fd = os.open(target, flags, 0o600)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(text)
    finally:
        os.chmod(target, 0o600)
    return str(target)


def _assert_server_path_contained(target: Path) -> None:
    base_dir = server_config_dir().resolve(strict=False)
    resolved_target = target.resolve(strict=False)
    if not resolved_target.is_relative_to(base_dir):
        raise ValueError(
            f"Path escapes ~/.chemsmart/server containment: {target}"
        )
