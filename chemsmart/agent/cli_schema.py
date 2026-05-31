"""Build a deterministic JSON schema for the ChemSmart click CLI tree."""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path
from typing import Any

import click

from chemsmart import __version__

JsonDict = dict[str, Any]


def build_chemsmart_cli_schema() -> JsonDict:
    """Return a stable, JSON-serializable schema for the ChemSmart CLI.

    The schema starts at :mod:`chemsmart.cli.main`'s ``entry_point`` command
    and recursively resolves click groups, including ChemSmart's lazy
    ``DeferredGroup`` proxies. Each command node contains its own click
    parameters and recursively nested subcommands.
    """

    numexpr_logger = logging.getLogger("numexpr.utils")
    previous_numexpr_level = numexpr_logger.level
    numexpr_logger.setLevel(logging.WARNING)
    from chemsmart.cli.main import entry_point

    numexpr_logger.setLevel(previous_numexpr_level)
    with click.Context(entry_point, info_name="chemsmart") as ctx:
        return _command_schema(entry_point, ctx)


def dump_schema_to_json(path: str | Path) -> JsonDict:
    """Write the ChemSmart CLI schema plus metadata to ``path``.

    Args:
        path: Destination JSON file path.

    Returns:
        The metadata-wrapped schema that was written to disk.
    """

    document = schema_with_metadata(build_chemsmart_cli_schema())
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return document


def schema_with_metadata(schema: JsonDict) -> JsonDict:
    """Return ``schema`` with ChemSmart version and stable content hash."""

    body = json.dumps(schema, sort_keys=True, separators=(",", ":"))
    schema_hash = hashlib.sha256(body.encode("utf-8")).hexdigest()
    document = dict(schema)
    document["_meta"] = {
        "chemsmart_version": __version__,
        "schema_hash": schema_hash,
    }
    return document


def _command_schema(command: click.Command, ctx: click.Context) -> JsonDict:
    resolved_command = _resolve_deferred_group(command)
    if resolved_command is not command:
        command = resolved_command
        ctx = click.Context(
            command, info_name=ctx.info_name, parent=ctx.parent
        )

    schema: JsonDict = {
        "name": ctx.info_name or command.name,
        "description": _command_description(command),
        "options": [_parameter_schema(param) for param in command.params],
        "subcommands": {},
    }

    if isinstance(command, click.Group):
        subcommands: dict[str, JsonDict] = {}
        for command_name in command.list_commands(ctx):
            child = command.get_command(ctx, command_name)
            if child is None or getattr(child, "hidden", False):
                continue
            child = _resolve_deferred_group(child)
            child_ctx = click.Context(
                child, info_name=command_name, parent=ctx
            )
            subcommands[command_name] = _command_schema(child, child_ctx)
        schema["subcommands"] = subcommands

    return schema


def _resolve_deferred_group(command: click.Command) -> click.Command:
    """Resolve ChemSmart lazy click groups without mutating their cache."""

    # DeferredGroup._load() stores the loaded command on the process-global
    # entry_point tree. Schema generation must be side-effect free so other
    # tests and callers still observe the lightweight help placeholders.
    loader = getattr(command, "_loader", None)
    if callable(loader):
        loaded = loader()
        if isinstance(loaded, click.Command):
            return loaded
    return command


def _command_description(command: click.Command) -> str | None:
    help_text = command.help or command.short_help
    if help_text is None:
        return None
    return inspect_cleandoc(help_text)


def inspect_cleandoc(text: str) -> str:
    """Normalize click help text without importing ``inspect`` at call sites."""

    import inspect

    return inspect.cleandoc(text)


def _parameter_schema(param: click.Parameter) -> JsonDict:
    default_info = _serialize_default(param.default)
    schema: JsonDict = {
        "name": param.name,
        "opts": _parameter_opts(param),
        "type": _serialize_type(param.type),
        "choices": _choices(param.type),
        "default": default_info["default"],
        "help": _parameter_help(param),
        "is_flag": bool(getattr(param, "is_flag", False)),
        "multiple": bool(param.multiple),
        "required": bool(param.required),
        "nargs": param.nargs,
    }
    if default_info["default_is_callable"]:
        schema["default_is_callable"] = True
    return schema


def _parameter_opts(param: click.Parameter) -> list[str]:
    if not isinstance(param, click.Option):
        return []
    opts = [*param.opts, *param.secondary_opts]
    filtered = [opt for opt in opts if opt not in {"--help", "-h"}]
    return sorted(dict.fromkeys(filtered), key=_option_sort_key)


def _option_sort_key(option: str) -> tuple[int, str]:
    return (1 if option.startswith("--") else 0, option)


def _parameter_help(param: click.Parameter) -> str | None:
    help_text = getattr(param, "help", None)
    if help_text is None:
        return None
    return inspect_cleandoc(help_text)


def _serialize_type(param_type: click.ParamType) -> str | JsonDict:
    if isinstance(param_type, click.Choice):
        return {"type": "choice", "choices": list(param_type.choices)}
    if isinstance(param_type, click.Path):
        return {
            "type": "path",
            "exists": bool(param_type.exists),
            "file_okay": bool(param_type.file_okay),
            "dir_okay": bool(param_type.dir_okay),
        }
    if isinstance(param_type, click.IntRange):
        return {
            "type": "int",
            "min": param_type.min,
            "max": param_type.max,
        }
    if isinstance(param_type, click.FloatRange):
        return {
            "type": "float",
            "min": param_type.min,
            "max": param_type.max,
        }
    if param_type is click.STRING or isinstance(
        param_type, click.types.StringParamType
    ):
        return "str"
    if param_type is click.INT or isinstance(
        param_type, click.types.IntParamType
    ):
        return "int"
    if param_type is click.FLOAT or isinstance(
        param_type, click.types.FloatParamType
    ):
        return "float"
    if param_type is click.BOOL or isinstance(
        param_type, click.types.BoolParamType
    ):
        return "bool"
    return getattr(param_type, "name", None) or param_type.__class__.__name__


def _choices(param_type: click.ParamType) -> list[str] | None:
    if isinstance(param_type, click.Choice):
        return list(param_type.choices)
    return None


def _serialize_default(default: Any) -> JsonDict:
    if callable(default):
        return {"default": None, "default_is_callable": True}
    normalized = _jsonable_value(default)
    try:
        json.dumps(normalized)
    except (TypeError, ValueError):
        return {"default": None, "default_is_callable": True}
    return {"default": normalized, "default_is_callable": False}


def _jsonable_value(value: Any) -> Any:
    if isinstance(value, tuple):
        return [_jsonable_value(item) for item in value]
    if isinstance(value, list):
        return [_jsonable_value(item) for item in value]
    if isinstance(value, dict):
        return {
            str(key): _jsonable_value(item)
            for key, item in sorted(
                value.items(), key=lambda pair: str(pair[0])
            )
        }
    return value
