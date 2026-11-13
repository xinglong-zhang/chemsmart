"""Discovery and validation of CHEMSMART CLI plugins."""

from __future__ import annotations

import importlib.metadata
import logging
from collections import Counter
from functools import lru_cache, wraps
from typing import Iterable

import click

from chemsmart.utils.cli import _add_subcommand_info_to_ctx
from chemsmart.utils.mixins import RegistryMixin

logger = logging.getLogger(__name__)

JOB_COMMANDS = "chemsmart.job_commands"
COMMANDS = "chemsmart.commands"


class PluginLoadError(RuntimeError):
    """Raised when a plugin cannot be loaded or validated."""


def _record_job_command(command):
    """Record ordinary Click invocations without replacing plugin classes."""
    if getattr(command, "_chemsmart_records_invocation", False):
        return command
    invoke = command.invoke
    get_params = command.get_params

    @wraps(get_params)
    def recording_get_params(ctx):
        params = get_params(ctx)
        for param in params:
            _record_parameter_value(param)
        return params

    command.get_params = recording_get_params

    @wraps(invoke)
    def record_and_invoke(ctx):
        record = _add_subcommand_info_to_ctx(ctx)
        record["parameters"] = tuple(ctx.command.get_params(ctx))
        record["raw_values"] = getattr(ctx, "_chemsmart_raw_values", {})
        return invoke(ctx)

    command.invoke = record_and_invoke
    if isinstance(command, click.Group):
        get_command = command.get_command

        @wraps(get_command)
        def recording_get_command(ctx, name):
            child = get_command(ctx, name)
            return _record_job_command(child) if child is not None else None

        # Intercept lookup so lazy and dynamically added children also work.
        command.get_command = recording_get_command
    command._chemsmart_records_invocation = True
    return command


def _record_parameter_value(param):
    """Capture resolved defaults/env values before conversion or callbacks."""
    if getattr(param, "_chemsmart_records_value", False):
        return
    process_value = param.process_value

    @wraps(process_value)
    def record_and_process(ctx, value):
        if not hasattr(ctx, "_chemsmart_raw_values"):
            ctx._chemsmart_raw_values = {}
        # Copy parser-owned sequences: callbacks may mutate their input.
        raw = tuple(value) if isinstance(value, (list, tuple)) else value
        ctx._chemsmart_raw_values[param.name] = raw
        return process_value(ctx, value)

    param.process_value = record_and_process
    param._chemsmart_records_value = True


def _provider(entry_point: importlib.metadata.EntryPoint) -> str:
    distribution = getattr(entry_point, "dist", None)
    name = getattr(distribution, "name", None) if distribution else None
    return f"{entry_point.name} = {entry_point.value}" + (
        f" ({name})" if name else ""
    )


def _declared_classes(command):
    """Validate the whole declaration before activating any plugin class."""
    classes = getattr(command, "chemsmart_classes", ())
    if not isinstance(classes, tuple):
        raise PluginLoadError("command.chemsmart_classes must be a tuple.")
    for cls in classes:
        if (
            not isinstance(cls, type)
            or cls is RegistryMixin
            or not issubclass(cls, RegistryMixin)
            or not cls.REGISTERABLE
        ):
            raise PluginLoadError(
                "command.chemsmart_classes must contain registerable "
                "RegistryMixin subclasses."
            )
    return classes


@lru_cache(maxsize=None)
def discover(
    group: str, strict: bool = False, reserved: frozenset[str] = frozenset()
) -> tuple[click.Command, ...]:
    """Load validated commands for an entry-point group once per process."""
    entries = sorted(
        importlib.metadata.entry_points().select(group=group),
        key=lambda entry: (entry.name, entry.value),
    )
    commands: list[click.Command] = []
    counts = Counter(entry.name for entry in entries)
    for entry in entries:
        try:
            if entry.name in reserved:
                raise PluginLoadError(
                    f"{group} entry point {_provider(entry)} conflicts with "
                    "an existing command."
                )
            if counts[entry.name] > 1:
                raise PluginLoadError(
                    f"{group} has multiple entry points named {entry.name!r}."
                )
            with RegistryMixin.defer_plugin_registration():
                command = entry.load()
            if not isinstance(command, click.Command):
                raise PluginLoadError(
                    f"{group} entry point {_provider(entry)} is not a "
                    "Click Command."
                )
            if command.name != entry.name:
                raise PluginLoadError(
                    f"{group} entry point {_provider(entry)} has name "
                    f"{entry.name!r}, but command.name is {command.name!r}."
                )
            classes = _declared_classes(command)
            if group == JOB_COMMANDS:
                _record_job_command(command)
            for cls in classes:
                if cls not in RegistryMixin._REGISTRY:
                    RegistryMixin._REGISTRY.append(cls)
            commands.append(command)
        except Exception as error:
            if strict:
                if isinstance(error, PluginLoadError):
                    raise
                raise PluginLoadError(
                    f"Failed to load {group} entry point {_provider(entry)}: "
                    f"{error}"
                ) from error
            logger.warning(
                "Skipping broken %s entry point %s: %s",
                group,
                _provider(entry),
                error,
            )
    return tuple(commands)


def add_commands(
    group: str,
    command_group: click.Group,
    builtins: Iterable[click.Command],
    *,
    strict: bool = False,
) -> None:
    """Register built-ins followed by validated plugin commands."""
    names: dict[str, str] = {
        name: "existing command" for name in command_group.commands
    }
    for command in builtins:
        if command.name is not None:
            if command.name in names:
                raise PluginLoadError(
                    f"{group} command name {command.name!r} conflicts with "
                    f"{names[command.name]}."
                )
            names[command.name] = "built-in"
        command_group.add_command(command)
    for command in discover(group, strict=strict, reserved=frozenset(names)):
        if command.name in names:
            error = PluginLoadError(
                f"{group} command name {command.name!r} conflicts between "
                f"{names[command.name]} and a plugin entry point."
            )
            if strict:
                raise error
            logger.warning("Skipping plugin: %s", error)
            continue
        names[command.name] = "plugin"
        command_group.add_command(command)
