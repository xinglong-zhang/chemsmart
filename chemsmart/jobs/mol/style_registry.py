"""PyMOL style names that are safe to import without PyMOL installed.

The actual template imports :mod:`pymol` when PyMOL executes it.  CLI command
registration only needs these public style identifiers, so it must not import
that execution-only template during ``chemsmart --help``.
"""

from __future__ import annotations

SCIENTIFIC_STYLE_COMMAND_NAMES = (
    "glossy",
    "comic",
    "soft_cartoon",
    "editorial_minimal",
    "soft_ceramic",
    "neon_coordination_core",
    "matte_clay",
    "xray_wire",
    "steric_surface",
    "quasi_chemdraw_bold",
)

PYMOL_SCIENTIFIC_STYLE_COMMANDS = {
    name: name for name in SCIENTIFIC_STYLE_COMMAND_NAMES
}
PYMOL_VISUALIZE_STYLE_CLI_CHOICES = [
    style.replace("_", "-") for style in SCIENTIFIC_STYLE_COMMAND_NAMES
]

_SCIENTIFIC_STYLE_TEMPLATE = "zhang_group_scientific_styles.py"
PYMOL_STYLE_TEMPLATES = {
    "pymol": "zhang_group_pymol_style.py",
    "cylview": "zhang_group_pymol_style.py",
    "cylview_flat": "zhang_group_pymol_style.py",
    **dict.fromkeys(
        SCIENTIFIC_STYLE_COMMAND_NAMES, _SCIENTIFIC_STYLE_TEMPLATE
    ),
}


def normalize_pymol_style(style):
    """Return a supported PyMOL style keyword."""

    normalized = (style or "pymol").lower().replace("-", "_")
    if normalized not in PYMOL_STYLE_TEMPLATES:
        raise ValueError(f"The style {style} is not available!")
    return normalized


def is_pymol_derived_style(style):
    """Return whether a style uses the scientific PyMOL template."""

    if style is None:
        return False
    return normalize_pymol_style(style) in PYMOL_SCIENTIFIC_STYLE_COMMANDS


__all__ = [
    "PYMOL_SCIENTIFIC_STYLE_COMMANDS",
    "PYMOL_STYLE_TEMPLATES",
    "PYMOL_VISUALIZE_STYLE_CLI_CHOICES",
    "SCIENTIFIC_STYLE_COMMAND_NAMES",
    "is_pymol_derived_style",
    "normalize_pymol_style",
]
