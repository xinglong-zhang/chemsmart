"""Generic JSON-safe projection shared by agent persistence services."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable


def generic_json_safe(
    value: Any,
    *,
    recurse: Callable[[Any], Any] | None = None,
) -> Any:
    project = recurse or (lambda item: generic_json_safe(item))
    if isinstance(value, dict):
        return {str(key): project(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [project(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, bytes):
        return {"type": "bytes", "length": len(value)}
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return {"type": value.__class__.__name__, "repr": repr(value)}


__all__ = ["generic_json_safe"]
