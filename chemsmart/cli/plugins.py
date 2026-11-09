"""Discovery and validation of CHEMSMART CLI plugins."""

from __future__ import annotations

import importlib.metadata
import logging
import sys
from collections import Counter
from functools import lru_cache
from typing import Iterable

import click

from chemsmart.utils.mixins import RegistryMixin

logger = logging.getLogger(__name__)

JOB_COMMANDS = "chemsmart.job_commands"
COMMANDS = "chemsmart.commands"


class PluginLoadError(RuntimeError):
    """Raised when a plugin cannot be loaded or validated."""


def _provider(entry_point: importlib.metadata.EntryPoint) -> str:
    distribution = getattr(entry_point, "dist", None)
    name = getattr(distribution, "name", None) if distribution else None
    return f"{entry_point.name} = {entry_point.value}" + (
        f" ({name})" if name else ""
    )


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
        registry_before = list(RegistryMixin._REGISTRY)
        modules_before = set(sys.modules)
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
            commands.append(command)
        except Exception as error:
            # Imports may register classes before failing validation. Restore
            # the registry and discard newly imported modules that defined
            # those classes so a later import can register them again.
            added = [
                cls for cls in RegistryMixin._REGISTRY
                if cls not in registry_before
            ]
            RegistryMixin._REGISTRY[:] = registry_before
            module_names = {cls.__module__ for cls in added}
            module_names.add(entry.value.split(":", 1)[0])
            for name in module_names - modules_before:
                sys.modules.pop(name, None)
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
