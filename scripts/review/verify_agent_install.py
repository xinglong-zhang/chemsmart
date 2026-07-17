#!/usr/bin/env python3
"""Verify ChemSmart's installed core and optional agent boundaries."""

from __future__ import annotations

import argparse
import json
import sys
from importlib import metadata, resources
from pathlib import Path
from typing import Any


def _distribution_root() -> Path:
    import chemsmart

    return Path(chemsmart.__file__).resolve().parent


def _verify_core() -> dict[str, Any]:
    from chemsmart.agent.harness.basis_sets.catalog import load_basis_catalog
    from chemsmart.agent.harness.command_semantics import (
        evaluate_command_semantics,
    )
    from chemsmart.agent.model_command_parser import parse_model_command
    from chemsmart.agent.project_yaml import validate_project_yaml
    from chemsmart.agent.v8_adapter import adapt

    package = resources.files("chemsmart.agent.harness.basis_sets")
    catalog_path = package.joinpath("bse_basis_catalog.json")
    catalog = load_basis_catalog()
    assert catalog_path.is_file()
    assert catalog.get("basis_sets")
    for symbol in (
        adapt,
        evaluate_command_semantics,
        parse_model_command,
        validate_project_yaml,
    ):
        assert callable(symbol)
    return {
        "basis_catalog_entries": len(catalog["basis_sets"]),
        "deterministic_imports": 4,
    }


def _verify_agent() -> dict[str, Any]:
    from chemsmart.agent import AgentSession
    from chemsmart.agent.cli import agent
    from chemsmart.agent.provider_config import load_active_provider_config
    from chemsmart.agent.providers import get_provider

    for symbol in (AgentSession, agent, get_provider, load_active_provider_config):
        assert callable(symbol)
    return {"agent_imports": 4}


def _verify_tui() -> dict[str, Any]:
    from chemsmart.agent.tui import launch_tui
    from chemsmart.agent.tui.screens.chat import ChatScreen

    assert callable(launch_tui)
    assert isinstance(ChatScreen.__name__, str)
    return {"tui_imports": 2}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--profile",
        choices=("core", "agent", "tui"),
        required=True,
    )
    parser.add_argument(
        "--require-installed",
        action="store_true",
        help="Reject an import resolved directly from the source checkout.",
    )
    args = parser.parse_args()

    root = _distribution_root()
    if args.require_installed and "site-packages" not in root.as_posix():
        raise RuntimeError(f"chemsmart resolved outside site-packages: {root}")

    receipt: dict[str, Any] = {
        "distribution": metadata.version("chemsmart"),
        "profile": args.profile,
        "python": ".".join(map(str, sys.version_info[:3])),
        "root": str(root),
    }
    receipt.update(_verify_core())
    if args.profile in {"agent", "tui"}:
        receipt.update(_verify_agent())
    if args.profile == "tui":
        receipt.update(_verify_tui())
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
