"""Coordinate, atom-range, and input-structure command contracts."""

from __future__ import annotations

import ast
from pathlib import Path
from typing import Any

from chemsmart.agent.harness.command_rules.models import (
    CommandContractIssue,
    reject,
)
from chemsmart.agent.harness.command_rules.tokens import option_value
from chemsmart.utils.utils import get_list_from_string_range


def parse_coordinate_literal(
    value: str | object,
) -> tuple[list[list[int]] | None, str | None]:
    if not isinstance(value, str):
        return None, "coordinates must have a value"
    try:
        parsed = ast.literal_eval(value)
    except (SyntaxError, ValueError):
        return None, "coordinates must be a Python/JSON list literal"
    if not isinstance(parsed, list) or not parsed:
        return None, "coordinates must be a non-empty list"

    if all(
        isinstance(item, int) and not isinstance(item, bool) for item in parsed
    ):
        groups: list[list[Any]] = [parsed]
    elif all(isinstance(item, list) for item in parsed):
        groups = parsed
    else:
        return None, "coordinates cannot mix atom indices and nested lists"

    normalized: list[list[int]] = []
    for group in groups:
        if not 2 <= len(group) <= 4:
            return None, "each coordinate must contain 2, 3, or 4 atom indices"
        if not all(
            isinstance(item, int) and not isinstance(item, bool) and item > 0
            for item in group
        ):
            return None, "coordinate atom indices must be positive integers"
        if len(group) != len(set(group)):
            return None, "a coordinate cannot repeat the same atom index"
        normalized.append(list(group))
    canonical_groups = {
        min(tuple(group), tuple(reversed(group))) for group in normalized
    }
    if len(normalized) != len(canonical_groups):
        return (
            None,
            "coordinate groups must not be duplicated or reversed duplicates",
        )
    return normalized, None


def parse_numeric_sequence(
    value: str | object,
    *,
    integers: bool,
) -> tuple[list[int | float] | None, str | None]:
    if not isinstance(value, str):
        return None, "must have a value"
    try:
        parsed = ast.literal_eval(value)
    except (SyntaxError, ValueError):
        return None, "must be a numeric scalar or list literal"
    values = parsed if isinstance(parsed, list) else [parsed]
    if not values:
        return None, "cannot be an empty list"
    if any(isinstance(item, bool) for item in values):
        return None, "must contain numeric values"
    if integers:
        if not all(isinstance(item, int) and item > 0 for item in values):
            return None, "must contain positive integers"
        return list(values), None
    if not all(isinstance(item, (int, float)) for item in values):
        return None, "must contain numeric values"
    return list(values), None


def range_literal_error(value: str | object) -> str | None:
    if not isinstance(value, str) or not value.strip():
        return "atom indices must have a value"
    if "[" in value or "]" in value or any(char.isspace() for char in value):
        return "atom indices must use a flat comma/range form without brackets or spaces"
    try:
        indices = get_list_from_string_range(value)
    except (TypeError, ValueError):
        return "atom indices must use a valid comma/range form such as 1-6,9"
    if not indices or any(index <= 0 for index in indices):
        return "atom indices must be positive and 1-indexed"
    if len(indices) != len(set(indices)):
        return "atom indices must not contain duplicates or overlapping ranges"
    return None


def atom_bound_issues(
    groups: list[list[int]],
    program_tokens: list[str],
    *,
    cwd: str | Path | None,
    rule_id: str,
    label: str,
) -> list[CommandContractIssue]:
    atom_count = input_atom_count(program_tokens, cwd=cwd)
    if atom_count is None:
        return []
    invalid = sorted(
        {
            index
            for group in groups
            for index in group
            if index < 1 or index > atom_count
        }
    )
    if not invalid:
        return []
    return [
        reject(
            rule_id,
            f"{label} contains atom indices outside the input structure",
            {
                "atom_count": atom_count,
                "invalid_indices": invalid,
                "input": option_value(program_tokens, ("-f", "--filename")),
            },
            (f"atom indices in the inclusive range 1-{atom_count}",),
        )
    ]


def input_atom_count(
    program_tokens: list[str],
    *,
    cwd: str | Path | None,
) -> int | None:
    source = option_value(program_tokens, ("-f", "--filename"))
    if not isinstance(source, str):
        return None
    path = Path(source).expanduser()
    if not path.is_absolute():
        path = Path(cwd or Path.cwd()) / path
    if path.suffix.lower() != ".xyz" or not path.is_file():
        return None
    try:
        first_line = path.read_text(encoding="utf-8").splitlines()[0].strip()
        count = int(first_line)
    except (OSError, UnicodeDecodeError, IndexError, ValueError):
        return None
    return count if count > 0 else None


__all__ = [
    "atom_bound_issues",
    "input_atom_count",
    "parse_coordinate_literal",
    "parse_numeric_sequence",
    "range_literal_error",
]
