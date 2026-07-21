"""Shared numeric equivalence for generated scan directives."""

from __future__ import annotations

from typing import Any


def scan_value_matches(key: str, expected: Any, observed: Any) -> bool:
    if key == "num_steps":
        return int(expected) == int(observed)
    return abs(float(expected) - float(observed)) <= 1e-9


__all__ = ["scan_value_matches"]
