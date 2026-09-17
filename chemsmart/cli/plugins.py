"""Discovery and validation of CHEMSMART CLI plugins."""

from __future__ import annotations

import importlib.metadata
import logging
from functools import lru_cache
from typing import Iterable

import click

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
def discover(group: str, strict: bool = False) -> tuple[click.Command, ...]:
    """Load validated commands for an entry-point group once per process."""
    entries = sorted(
        importlib.metadata.entry_points().select(group=group),
        key=lambda entry: (entry.name, entry.value),
    )
    commands: list[click.Command] = []
    seen: dict[str, str] = {}
    for entry in entries:
        try:
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
            if command.name in seen:
                raise PluginLoadError(
                    f"{group} command name {command.name!r} conflicts between "
                    f"{seen[command.name]} and {_provider(entry)}."
                )
            seen[command.name] = _provider(entry)
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
    names: dict[str, str] = {}
    for command in builtins:
        if command.name is not None:
            names[command.name] = "built-in"
        command_group.add_command(command)
    for command in discover(group, strict=strict):
        if command.name in names:
            raise PluginLoadError(
                f"{group} command name {command.name!r} conflicts between "
                f"{names[command.name]} and a plugin entry point."
            )
        names[command.name] = "plugin"
        command_group.add_command(command)
