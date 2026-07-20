"""Shared normalization for nested chemistry intent values."""

from __future__ import annotations

import ast
import re
from typing import Any


def structured_sequence(
    value: Any,
    *,
    allow_numeric_text: bool = False,
) -> tuple[Any, ...] | None:
    if isinstance(value, str):
        try:
            value = ast.literal_eval(value)
        except (SyntaxError, ValueError):
            if not allow_numeric_text:
                return None
            tokens = re.findall(r"-?\d+(?:\.\d+)?", value)
            if not tokens:
                return None
            value = [float(token) for token in tokens]
    if not isinstance(value, (list, tuple)):
        return None
    normalized: list[Any] = []
    for item in value:
        if isinstance(item, (list, tuple)):
            nested = structured_sequence(
                item,
                allow_numeric_text=allow_numeric_text,
            )
            if nested is None:
                return None
            normalized.append(nested)
        elif isinstance(item, bool):
            normalized.append(item)
        elif isinstance(item, (int, float)):
            normalized.append(float(item))
        else:
            normalized.append(str(item).strip().lower())
    return tuple(normalized)


def unwrap_singleton_group(
    value: tuple[Any, ...] | None,
) -> tuple[Any, ...] | None:
    if value is not None and len(value) == 1 and isinstance(value[0], tuple):
        return value[0]
    return value


__all__ = ["structured_sequence", "unwrap_singleton_group"]
