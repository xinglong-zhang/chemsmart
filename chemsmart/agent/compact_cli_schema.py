"""Compact, lossless-enough signatures for request-scoped CLI synthesis.

The full Click schema is the runtime source of truth, but its prose, defaults,
and repeated parameter names are poor prompt material.  This module projects a
pruned schema into positional command paths and the minimum option contract the
model needs to emit a legal command.  Validation still uses the full schema.
"""

from __future__ import annotations

from typing import Any

JsonDict = dict[str, Any]


def compact_cli_signature(schema: JsonDict) -> JsonDict:
    """Return deterministic command signatures without help/default bloat."""

    commands: list[JsonDict] = []
    root_name = str(schema.get("name") or "chemsmart")
    _walk(schema, prefix=(), commands=commands, include_node=False)
    return {"root": root_name, "commands": commands}


def compact_signature_paths(signature: JsonDict) -> set[str]:
    """Return command paths represented by a compact signature."""

    return {
        str(command.get("path") or "")
        for command in signature.get("commands") or []
        if isinstance(command, dict) and command.get("path")
    }


def _walk(
    node: JsonDict,
    *,
    prefix: tuple[str, ...],
    commands: list[JsonDict],
    include_node: bool,
) -> None:
    if include_node:
        entry: JsonDict = {"path": " ".join(prefix)}
        options = [
            compact
            for option in node.get("options") or []
            if isinstance(option, dict)
            and (compact := _compact_parameter(option)) is not None
        ]
        if options:
            entry["args"] = options
        commands.append(entry)

    children = node.get("subcommands") or {}
    if not isinstance(children, dict):
        return
    for name in sorted(children):
        child = children[name]
        if isinstance(child, dict):
            _walk(
                child,
                prefix=(*prefix, str(name)),
                commands=commands,
                include_node=True,
            )


def _compact_parameter(parameter: JsonDict) -> str | None:
    """Encode one Click parameter without repeating JSON field names.

    Syntax: ``flags[:type][:{choices}][!][[]]`` where ``!`` means required
    and ``[]`` means repeatable.  Plain string values omit ``:str``.
    """

    opts = [str(item) for item in parameter.get("opts") or [] if item]
    name = str(parameter.get("name") or "").strip()
    if not opts and not name:
        return None

    compact = "|".join(opts) if opts else f"<{name}>"
    value_type = _compact_type(parameter.get("type"))
    choices = parameter.get("choices") or (
        value_type.get("choices") if isinstance(value_type, dict) else None
    )
    if isinstance(choices, list) and choices:
        compact += ":{" + "|".join(str(item) for item in choices) + "}"
    elif not isinstance(value_type, str) or value_type not in {"str", "value"}:
        compact += ":" + _type_token(value_type)
    if parameter.get("required"):
        compact += "!"
    if parameter.get("multiple"):
        compact += "[]"
    return compact


def _type_token(value: Any) -> str:
    if not isinstance(value, dict):
        return str(value)
    kind = str(value.get("type") or "value")
    bounds = [
        f"{key}={value[key]}"
        for key in ("min", "max")
        if value.get(key) is not None
    ]
    return kind + ("(" + ",".join(bounds) + ")" if bounds else "")


def _compact_type(value: Any) -> Any:
    if not isinstance(value, dict):
        return value or "str"
    kind = value.get("type") or "value"
    choices = value.get("choices")
    if isinstance(choices, list) and choices:
        return {"type": kind, "choices": choices}
    bounds = {
        key: value[key]
        for key in ("min", "max")
        if value.get(key) is not None
    }
    return {"type": kind, **bounds} if bounds else kind
