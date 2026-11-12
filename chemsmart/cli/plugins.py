"""Discovery and validation of CHEMSMART CLI plugins."""

from __future__ import annotations

import importlib.metadata
import logging
from collections import Counter
from functools import lru_cache
from typing import Iterable

import click

from chemsmart.utils.cli import MyCommand, MyGroup
from chemsmart.utils.mixins import RegistryMixin

logger = logging.getLogger(__name__)

JOB_COMMANDS = "chemsmart.job_commands"
COMMANDS = "chemsmart.commands"


class PluginLoadError(RuntimeError):
    """Raised when a plugin cannot be loaded or validated."""


def _validate_job_command(command):
    """Ensure every job command records arguments for submission scripts."""
    if isinstance(command, click.Group):
        if not isinstance(command, MyGroup):
            raise PluginLoadError("Job command groups must use MyGroup.")
        for child in command.commands.values():
            _validate_job_command(child)
    elif not isinstance(command, MyCommand):
        raise PluginLoadError("Job commands must use MyCommand.")


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
            if group == JOB_COMMANDS:
                _validate_job_command(command)
            classes = _declared_classes(command)
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
