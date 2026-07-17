"""Token and scalar helpers shared by command-contract rules."""

from __future__ import annotations

import ast


MISSING = object()


def token_index(tokens: list[str], target: str) -> int | None:
    try:
        return tokens.index(target)
    except ValueError:
        return None


def has_option(tokens: list[str], aliases: tuple[str, ...]) -> bool:
    return option_value(tokens, aliases) is not MISSING


def option_value(tokens: list[str], aliases: tuple[str, ...]) -> str | object:
    for index, token in enumerate(tokens):
        for alias in aliases:
            if token == alias:
                if index + 1 < len(tokens):
                    return tokens[index + 1]
                return None
            if alias.startswith("--") and token.startswith(f"{alias}="):
                return token.split("=", 1)[1]
    return MISSING


def is_positive_integer_literal(value: str | object) -> bool:
    try:
        parsed = ast.literal_eval(value) if isinstance(value, str) else value
    except (SyntaxError, ValueError):
        return False
    return isinstance(parsed, int) and not isinstance(parsed, bool) and parsed > 0


def is_unit_interval_literal(value: str | object) -> bool:
    try:
        parsed = ast.literal_eval(value) if isinstance(value, str) else value
    except (SyntaxError, ValueError):
        return False
    return (
        isinstance(parsed, (int, float))
        and not isinstance(parsed, bool)
        and 0 < float(parsed) <= 1
    )


__all__ = [
    "MISSING",
    "has_option",
    "is_positive_integer_literal",
    "is_unit_interval_literal",
    "option_value",
    "token_index",
]
