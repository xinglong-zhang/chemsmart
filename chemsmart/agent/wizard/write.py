"""Persist wizard-rendered server YAML files."""

from __future__ import annotations

from chemsmart.agent.wizard.paths import server_yaml_path, write_private_text


def write_server_yaml(
    name: str,
    yaml_text: str,
    overwrite: bool = False,
) -> str:
    """Write a server YAML under ``~/.chemsmart/server`` and return its path."""

    target = server_yaml_path(name)
    return write_private_text(target, yaml_text, overwrite=overwrite)
