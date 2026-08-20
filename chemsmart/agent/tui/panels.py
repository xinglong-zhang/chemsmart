"""Chrome renderables: wordmark, phase chip, and (later) live panels."""

from __future__ import annotations

from importlib import metadata

from rich.text import Text


def product_version() -> str:
    try:
        return metadata.version("chemsmart")
    except metadata.PackageNotFoundError:  # pragma: no cover - dev checkout
        return "dev"


def wordmark() -> Text:
    return Text.assemble(
        ("ChemSmart", "bold cyan"),
        (" Agent", "bold"),
        (f"  v{product_version()}", "dim"),
    )


def phase_chip(phase: str, hint: str) -> Text:
    return Text.assemble(
        ("● ", "bold cyan"),
        (phase, "bold"),
        ("  ·  ", "dim"),
        (hint, "dim"),
    )


__all__ = ["phase_chip", "product_version", "wordmark"]
