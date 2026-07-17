"""Canonical session loading and explicit legacy-session migration."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import uuid
from pathlib import Path

from chemsmart.agent.models import SessionState, utc_now_iso

CURRENT_STATE_NAME = "session.json"
LEGACY_STATE_NAME = "state.json"
MIGRATION_MANIFEST_NAME = "migration_manifest.json"
SESSION_SCHEMA_VERSION = 2


class LegacySessionFormatError(RuntimeError):
    """Raised when a legacy session requires an explicit migration."""


class SessionMigrationError(RuntimeError):
    """Raised when a session cannot be migrated without data loss."""


def current_session_dirs(
    session_root: str | os.PathLike[str],
) -> list[Path]:
    """Return canonical session directories, newest name first."""
    root = Path(session_root)
    if not root.is_dir():
        return []
    sessions = [
        path
        for path in root.iterdir()
        if path.is_dir()
        and not path.name.startswith(".")
        and (path / CURRENT_STATE_NAME).is_file()
    ]
    return sorted(sessions, reverse=True)


def load_current_session_state(
    session_dir: str | os.PathLike[str],
    *,
    required: bool = False,
) -> SessionState | None:
    """Load canonical state without silently interpreting legacy artifacts."""
    directory = Path(session_dir)
    current_path = directory / CURRENT_STATE_NAME
    if current_path.is_file():
        return SessionState.load(current_path)
    if (directory / LEGACY_STATE_NAME).is_file():
        raise LegacySessionFormatError(_legacy_migration_message(directory))
    if required:
        raise FileNotFoundError(current_path)
    return None


def migrate_legacy_session(
    source: str | os.PathLike[str],
    destination: str | os.PathLike[str] | None = None,
) -> dict[str, object]:
    """Copy a legacy session into a new canonical session directory."""
    source_path = Path(source).expanduser().resolve()
    if not source_path.is_dir():
        raise SessionMigrationError(
            f"Legacy session directory does not exist: {source_path}"
        )
    if (source_path / CURRENT_STATE_NAME).exists():
        raise SessionMigrationError(
            f"Session already uses {CURRENT_STATE_NAME}: {source_path}"
        )
    legacy_path = source_path / LEGACY_STATE_NAME
    if not legacy_path.is_file():
        raise SessionMigrationError(
            f"Legacy session has no {LEGACY_STATE_NAME}: {source_path}"
        )

    target_path = _migration_destination(source_path, destination)
    _validate_destination(source_path, target_path)
    source_hashes = _tree_hashes(source_path)
    legacy_state = _load_legacy_state(legacy_path)
    legacy_session_id = legacy_state.session_id
    legacy_state.session_id = target_path.name

    target_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = target_path.with_name(
        f".{target_path.name}.migrating-{uuid.uuid4().hex[:8]}"
    )
    try:
        shutil.copytree(source_path, temporary, symlinks=True)
        legacy_state.save(temporary / CURRENT_STATE_NAME)
        legacy_state.save(temporary / LEGACY_STATE_NAME)
        manifest = {
            "schema_version": SESSION_SCHEMA_VERSION,
            "migrated_at": utc_now_iso(),
            "source": str(source_path),
            "destination": str(target_path),
            "source_session_id": legacy_session_id,
            "session_id": legacy_state.session_id,
            "source_artifact_sha256": source_hashes,
        }
        (temporary / MIGRATION_MANIFEST_NAME).write_text(
            json.dumps(manifest, indent=2, sort_keys=True),
            encoding="utf-8",
        )
        temporary.replace(target_path)
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise

    return manifest


def resolve_session_source(
    value: str,
    session_root: str | os.PathLike[str],
) -> Path:
    """Resolve either an explicit path or a session id under the root."""
    explicit = Path(value).expanduser()
    if explicit.exists():
        return explicit.resolve()
    return (Path(session_root) / value).resolve()


def _load_legacy_state(path: Path) -> SessionState:
    try:
        return SessionState.load(path)
    except Exception as exc:
        raise SessionMigrationError(
            f"Legacy session state is malformed: {path}"
        ) from exc


def _migration_destination(
    source: Path,
    destination: str | os.PathLike[str] | None,
) -> Path:
    if destination is None:
        return source.with_name(f"{source.name}-migrated-v2")
    return Path(destination).expanduser().resolve()


def _validate_destination(source: Path, destination: Path) -> None:
    if destination.exists():
        raise SessionMigrationError(
            f"Migration destination already exists: {destination}"
        )
    if destination == source or source in destination.parents:
        raise SessionMigrationError(
            "Migration destination must not be the source or its descendant."
        )


def _tree_hashes(root: Path) -> dict[str, str]:
    hashes: dict[str, str] = {}
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            payload = f"symlink:{os.readlink(path)}".encode()
        elif path.is_file():
            payload = path.read_bytes()
        else:
            continue
        hashes[str(path.relative_to(root))] = hashlib.sha256(
            payload
        ).hexdigest()
    return hashes


def _legacy_migration_message(session_dir: Path) -> str:
    return (
        "Legacy session format detected. Migrate it explicitly with: "
        f"chemsmart agent migrate-session {session_dir}"
    )


__all__ = [
    "CURRENT_STATE_NAME",
    "LEGACY_STATE_NAME",
    "LegacySessionFormatError",
    "MIGRATION_MANIFEST_NAME",
    "SESSION_SCHEMA_VERSION",
    "SessionMigrationError",
    "current_session_dirs",
    "load_current_session_state",
    "migrate_legacy_session",
    "resolve_session_source",
]
