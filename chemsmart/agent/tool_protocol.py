from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from chemsmart.agent.permissions import RuntimePermissionMode


@dataclass(frozen=True)
class RuntimeToolMetadata:
    read_only: bool = False
    edit_safe: bool = False
    requires_mode: frozenset[RuntimePermissionMode] = frozenset()
    ui_summary_template: str | None = None
    side_effect: str | None = None


def is_read_only(spec: Any) -> bool:
    return _metadata_for(spec).read_only


def is_allowed_in_mode(
    spec: Any,
    mode: RuntimePermissionMode,
) -> bool:
    metadata = _metadata_for(spec)
    if metadata.requires_mode and mode not in metadata.requires_mode:
        return False
    if mode == RuntimePermissionMode.PLAN:
        return False
    if mode == RuntimePermissionMode.BYPASS:
        return True
    if metadata == RuntimeToolMetadata():
        return True
    if mode == RuntimePermissionMode.READ_ONLY:
        return metadata.read_only
    if mode == RuntimePermissionMode.ACCEPT_EDITS:
        return metadata.read_only or metadata.edit_safe
    return False


def _metadata_for(spec: Any) -> RuntimeToolMetadata:
    metadata = getattr(spec, "metadata", None)
    if isinstance(metadata, RuntimeToolMetadata):
        return metadata
    return RuntimeToolMetadata()
