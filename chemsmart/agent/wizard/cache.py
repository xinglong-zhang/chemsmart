"""Sidecar cache persistence for wizard node-refresh data."""

from __future__ import annotations

import json
import logging
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timedelta, timezone

UTC = timezone.utc
from pathlib import Path
from typing import Literal

from chemsmart.agent.wizard.paths import server_cache_path, write_private_text

logger = logging.getLogger(__name__)

CacheStatus = Literal["fresh", "stale", "error"]
_VALID_STATUSES = {"fresh", "stale", "error"}


@dataclass(frozen=True)
class CacheEntry:
    server_name: str
    host: str | None
    mode: str
    scheduler: str | None
    probed_at: str
    source_commands: dict[str, str]
    partitions: list[dict]
    node_summary: dict
    program_candidates: dict
    status: CacheStatus
    last_error: str | None


def cache_path(server_name: str) -> Path:
    """Return the sidecar cache path for a wizard server name."""

    return server_cache_path(server_name)


def load_cache(name: str) -> CacheEntry | None:
    """Load a cache entry, returning ``None`` on missing or corrupt data."""

    path = cache_path(name)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return None
    except (OSError, json.JSONDecodeError) as exc:
        logger.warning("Failed to load wizard cache %s: %s", path, exc)
        return None

    try:
        return _entry_from_payload(payload)
    except (TypeError, ValueError) as exc:
        logger.warning("Ignoring corrupt wizard cache %s: %s", path, exc)
        return None


def write_cache(entry: CacheEntry) -> str:
    """Persist a cache entry as JSON and return its absolute path."""

    path = cache_path(entry.server_name)
    return write_private_text(
        path,
        json.dumps(asdict(entry), indent=2, sort_keys=True) + "\n",
        overwrite=True,
    )


def is_stale(entry: CacheEntry, ttl_hours: int = 24) -> bool:
    """Return ``True`` when the cache entry is older than the TTL."""

    try:
        probed_at = _parse_iso8601_utc(entry.probed_at)
    except ValueError:
        logger.warning(
            "Cache entry for %s has invalid probed_at %r; treating as stale",
            entry.server_name,
            entry.probed_at,
        )
        return True
    age = datetime.now(UTC) - probed_at
    return age > timedelta(hours=ttl_hours)


def mark_status(
    entry: CacheEntry,
    status: CacheStatus,
    last_error: str | None = None,
) -> CacheEntry:
    """Return a new cache entry with an updated status/error field."""

    _validate_status(status)
    return replace(entry, status=status, last_error=last_error)


def _entry_from_payload(payload: object) -> CacheEntry:
    if not isinstance(payload, dict):
        raise TypeError("cache payload must be a JSON object")

    status = payload.get("status")
    _validate_status(status)

    entry = CacheEntry(
        server_name=_require_str(payload, "server_name"),
        host=_optional_str(payload.get("host")),
        mode=_require_str(payload, "mode"),
        scheduler=_optional_str(payload.get("scheduler")),
        probed_at=_require_str(payload, "probed_at"),
        source_commands=_require_dict(payload, "source_commands"),
        partitions=_require_list(payload, "partitions"),
        node_summary=_require_dict(payload, "node_summary"),
        program_candidates=_require_dict(payload, "program_candidates"),
        status=status,
        last_error=_optional_str(payload.get("last_error")),
    )
    _parse_iso8601_utc(entry.probed_at)
    return entry


def _validate_status(status: object) -> None:
    if status not in _VALID_STATUSES:
        raise ValueError(f"invalid cache status: {status!r}")


def _require_str(payload: dict, field_name: str) -> str:
    value = payload.get(field_name)
    if not isinstance(value, str):
        raise TypeError(f"{field_name} must be a string")
    return value


def _optional_str(value: object) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise TypeError("optional string field must be a string or null")
    return value


def _require_dict(payload: dict, field_name: str) -> dict:
    value = payload.get(field_name)
    if not isinstance(value, dict):
        raise TypeError(f"{field_name} must be an object")
    return value


def _require_list(payload: dict, field_name: str) -> list:
    value = payload.get(field_name)
    if not isinstance(value, list):
        raise TypeError(f"{field_name} must be a list")
    return value


def _parse_iso8601_utc(value: str) -> datetime:
    normalized = value.replace("Z", "+00:00")
    parsed = datetime.fromisoformat(normalized)
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=UTC)
    return parsed.astimezone(UTC)
