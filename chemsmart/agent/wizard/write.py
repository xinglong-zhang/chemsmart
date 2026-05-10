"""Persist wizard-rendered server YAML files."""

from __future__ import annotations

from pathlib import Path


def write_server_yaml(
    name: str,
    yaml_text: str,
    overwrite: bool = False,
) -> str:
    """Write a server YAML under ``~/.chemsmart/server`` and return its path."""

    target = Path.home() / ".chemsmart" / "server" / f"{name}.yaml"
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and not overwrite:
        raise FileExistsError(f"Server YAML already exists: {target}")
    target.write_text(yaml_text, encoding="utf-8")
    return str(target)
