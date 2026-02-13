"""Utility functions for iterate module."""

import logging
import os
import re

import click

from chemsmart.utils.repattern import safe_label_pattern
from chemsmart.utils.utils import get_list_from_string_range

logger = logging.getLogger(__name__)

# Define allowed keys for configuration (TOML format) validation
ALLOWED_TOP_LEVEL_KEYS = {"skeletons", "substituents"}
ALLOWED_SKELETON_KEYS = {
    "file_path",
    "label",
    "link_index",
    "skeleton_indices",
}
ALLOWED_SUBSTITUENT_KEYS = {
    "file_path",
    "label",
    "link_index",
}

ITERATE_TEMPLATE = """# Chemsmart Iterate Configuration Template
# =========================================
# This template defines skeletons and substituents for iterative structure generation.
#
# Structure Sources:
#   - file_path: Path to molecular structure file (.xyz, .com, .log, .sdf, etc.)
#
# Required Fields:
#   - label: Unique identifier for this structure
#   - link_index: Atom index (1-based) where substitution occurs
#
# Optional Fields:
#   - skeleton_indices: Atom indices to keep, format: "1, 3-10, 15" (1-based)
#                       If not specified, all atoms except link_index are kept
#                       IMPORTANT: If specified, link_index MUST be included in skeleton_indices, otherwise an error will occur.
#
# Note: Index values may be provided as integers, lists, or strings (ranges/comma-separated).

# ==============================================================================
# SKELETONS
# ==============================================================================
# Define base molecular scaffolds that will be modified with substituents.

# Example 1: Skeleton from file
[[skeletons]]
file_path = "/path/to/skeleton1.xyz"
label = "benzene_core"
link_index = "6,7,8"
skeleton_indices = "1-5"

# Example 2: Another skeleton
# [[skeletons]]
# file_path = "..."
# label = "..."
# ...


# ==============================================================================
# SUBSTITUENTS
# ==============================================================================
# Define functional groups or fragments to attach to skeletons.

# Example 1: Substituent from file
[[substituents]]
file_path = "/path/to/methyl.xyz"
label = "methyl"
link_index = "1"

# Add more substituents as needed...
# [[substituents]]
# ...


# ==============================================================================
# USAGE NOTES
# ==============================================================================
# 1. Atom indices are 1-based (first atom is index 1)
# 2. link_index specifies the atom that will form the new bond
# 3. For skeletons, the atom at link_index is typically replaced/removed
# 4. For substituents, the atom at link_index bonds to the skeleton
# 5. skeleton_indices format: "1, 3-10, 15, 20-25" (ranges and individual indices)
"""


def generate_template(
    output_path: str = "iterate_template.toml", overwrite: bool = False
) -> str:
    """
    Generate a template configuration file for iterate jobs.

    Parameters
    ----------
    output_path : str
        Path to write the template file. Default is 'iterate_template.toml'.
    overwrite : bool
        If True, overwrite existing file. Default is False.

    Returns
    -------
    str
        Path to the generated template file.

    Raises
    ------
    FileExistsError
        If file exists and overwrite is False.
    """
    # Add .toml extension if not present
    if not output_path.endswith(".toml"):
        output_path = f"{output_path}.toml"

    # Check if file exists
    if os.path.exists(output_path) and not overwrite:
        raise FileExistsError(
            f"File '{output_path}' already exists. "
            f"Use --overwrite or specify a different filename."
        )

    # Write template
    with open(output_path, "w") as f:
        f.write(ITERATE_TEMPLATE)

    logger.info(f"Generated template: {output_path}")
    return output_path


def _parse_index_string(
    value, entry_type: str, idx: int, field_name: str
) -> list[int] | None:
    """
    Parse index string into a list of integers using utility function.
    Wrapper to handle empty values and single integers.

    Parameters
    ----------
    value : int, str, or None
        The value to parse
    entry_type : str
        "skeleton" or "substituent" (for error messages context, though currently unused for errors)
    idx : int
        Entry index
    field_name : str
        Field name

    Returns
    -------
    list[int] | None
        List of 1-based indices, or None if value is empty/null
    """
    # Handle empty/null values
    if value is None or value == "" or value == "null":
        return None

    # Handle single integer
    if isinstance(value, int):
        if value <= 0:
            raise click.BadParameter(
                f"{entry_type.capitalize()} entry {idx + 1}: Found invalid index {value} in '{field_name}'. "
                f"Index must be positive (1-based).",
                param_hint="'-f' / '--filename'",
            )
        return [value]

    # Use utility function for string parsing
    try:
        parsed_indices = None
        if isinstance(value, str):
            parsed_indices = get_list_from_string_range(value)
        elif isinstance(value, list) and all(
            isinstance(x, int) for x in value
        ):
            # Already a list of ints (though TOML parser usually handles this, sometimes flexible)
            parsed_indices = value

        if parsed_indices is not None:
            # S2 Check: Validate not empty
            if len(parsed_indices) == 0:
                raise click.BadParameter(
                    f"{entry_type.capitalize()} entry {idx + 1}: Found empty list in '{field_name}'. "
                    f"At least one index must be provided.",
                    param_hint="'-f' / '--filename'",
                )

            # S2 Check: Validate positive non-zero indices
            if any(i <= 0 for i in parsed_indices):
                raise click.BadParameter(
                    f"{entry_type.capitalize()} entry {idx + 1}: Found invalid index <= 0 in '{field_name}'. "
                    f"All indices must be positive (1-based). Found: {parsed_indices}",
                    param_hint="'-f' / '--filename'",
                )
            return parsed_indices

    except click.BadParameter:
        raise
    except Exception:
        raise click.BadParameter(
            f"{entry_type.capitalize()} entry {idx + 1}: Invalid format '{value}' in '{field_name}'. "
            f"Expected integer, comma-separated list, or range (e.g. '1-5').",
            param_hint="'-f' / '--filename'",
        )

    # If it's not None, not Int, and not String (or string parsing failed silently elsewhere)
    raise click.BadParameter(
        f"{entry_type.capitalize()} entry {idx + 1}: '{field_name}' has invalid type {type(value).__name__}.",
        param_hint="'-f' / '--filename'",
    )


def validate_config(config: dict, filename: str) -> dict:
    """
    Validate configuration (from .toml file) and normalize values.

    Parameters
    ----------
    config : dict
        Raw configuration dictionary from parsed file
    filename : str
        Path to configuration file (for error messages)

    Returns
    -------
    dict
        Validated and normalized configuration

    Raises
    ------
    click.BadParameter
        If configuration contains invalid keys or structure
    """
    # Check top-level keys
    unknown_top_keys = set(config.keys()) - ALLOWED_TOP_LEVEL_KEYS
    if unknown_top_keys:
        raise click.BadParameter(
            f"Unknown top-level key(s) in configuration: {unknown_top_keys}. "
            f"Allowed keys: {ALLOWED_TOP_LEVEL_KEYS}",
            param_hint="'-f' / '--filename'",
        )

    validated_config = {}

    # Validate and normalize skeletons
    if "skeletons" in config:
        skeletons = config["skeletons"]
        if skeletons is None:
            validated_config["skeletons"] = []
        elif not isinstance(skeletons, list):
            raise click.BadParameter(
                f"'skeletons' must be a list, got {type(skeletons).__name__}",
                param_hint="'-f' / '--filename'",
            )
        else:
            validated_config["skeletons"] = [
                _validate_skeleton_entry(entry, idx, filename)
                for idx, entry in enumerate(skeletons)
            ]
    else:
        validated_config["skeletons"] = []

    # Validate and normalize substituents
    if "substituents" in config:
        substituents = config["substituents"]
        if substituents is None:
            validated_config["substituents"] = []
        elif not isinstance(substituents, list):
            raise click.BadParameter(
                f"'substituents' must be a list, got {type(substituents).__name__}",
                param_hint="'-f' / '--filename'",
            )
        else:
            validated_config["substituents"] = [
                _validate_substituent_entry(entry, idx, filename)
                for idx, entry in enumerate(substituents)
            ]
    else:
        validated_config["substituents"] = []

    # Business Logic Validation
    _validate_business_logic(
        validated_config["skeletons"],
        validated_config["substituents"],
        filename,
    )

    return validated_config


def _validate_business_logic(
    skeletons: list, substituents: list, filename: str
):
    """
    Orchestrate business logic checks for all entities.
    """
    _validate_skeletons_logic(skeletons, filename)

    _validate_substituents_logic(substituents, filename)


def _validate_skeletons_logic(skeletons: list, filename: str):
    """
    Validate business rules for all skeleton entries.
    """
    for idx, entry in enumerate(skeletons):
        _validate_single_skeleton_rules(entry, idx, filename)


def _validate_single_skeleton_rules(entry: dict, idx: int, filename: str):
    """
    Collection of all business rules for a single skeleton.
    """
    # Rule 1: Check if link_index is included in skeleton_indices
    _rule_skeleton_indices_must_contain_link(entry, idx, filename)

    # Rule 2: file_path is required
    _rule_file_path_required(entry, idx, filename)

    # Rule 3: label must be safe for filenames
    _rule_label_syntax(entry, idx, "Skeleton", filename)


def _rule_label_syntax(entry: dict, idx: int, entry_type: str, filename: str):
    """
    Rule: Label must only contain safe characters [a-zA-Z0-9_-.].
    """
    label = entry.get("label")
    if label:
        if not re.match(safe_label_pattern, label):
            raise click.BadParameter(
                f"{entry_type} entry {idx + 1} label '{label}': "
                f"Contains invalid characters. Allowed characters: a-z, A-Z, 0-9, _, -, .",
                param_hint=filename,
            )


def _rule_file_path_required(entry: dict, idx: int, filename: str):
    """
    Rule: file_path must be provided and not empty.
    """
    if not entry.get("file_path"):
        raise click.BadParameter(
            f"Skeleton entry {idx + 1} (label='{entry.get('label', 'unnamed')}'): "
            f"Missing required field 'file_path'.",
            param_hint=filename,
        )


def _rule_skeleton_indices_must_contain_link(
    entry: dict, idx: int, filename: str
):
    """
    Rule: If skeleton_indices is specified, link_index must be included in it.
    """
    skeleton_indices = entry.get("skeleton_indices")
    link_indices = entry.get("link_index")

    if skeleton_indices is not None and link_indices:
        skel_set = set(skeleton_indices)
        missing = [i for i in link_indices if i not in skel_set]

        if missing:
            raise click.BadParameter(
                f"Skeleton entry {idx + 1} (label='{entry.get('label')}'): "
                f"The link_index {missing} is not included in 'skeleton_indices'. "
                f"When 'skeleton_indices' is specified, it must contain the link atom.",
                param_hint=filename,
            )


def _validate_substituents_logic(substituents: list, filename: str):
    """
    Validate business rules for all substituent entries.
    """
    for idx, entry in enumerate(substituents):
        _validate_single_substituent_rules(entry, idx, filename)


def _validate_single_substituent_rules(entry: dict, idx: int, filename: str):
    """
    Collection of all business rules for a single substituent.
    """
    # Rule 1: file_path is required
    if not entry.get("file_path"):
        raise click.BadParameter(
            f"Substituent entry {idx + 1} (label='{entry.get('label', 'unnamed')}'): "
            f"Missing required field 'file_path'.",
            param_hint=filename,
        )

    # Rule 2: link_index must be single value
    link_index = entry.get("link_index")
    if link_index and len(link_index) > 1:
        raise click.BadParameter(
            f"Substituent entry {idx + 1} (label='{entry.get('label', 'unnamed')}'): "
            f"Multiple values found in 'link_index': {link_index}. "
            f"Substituents must have exactly one link atom.",
            param_hint=filename,
        )

    # Rule 3: label must be safe for filenames
    _rule_label_syntax(entry, idx, "Substituent", filename)


def _validate_skeleton_entry(entry: dict, idx: int, filename: str) -> dict:
    """
    Validate a single skeleton entry.

    Parameters
    ----------
    entry : dict
        Skeleton entry from TOML
    idx : int
        Index of the entry (for error messages)
    filename : str
        Path to TOML file (for error messages)

    Returns
    -------
    dict
        Validated and normalized skeleton entry
    """
    if entry is None:
        raise click.BadParameter(
            f"Skeleton entry {idx + 1} is empty/null",
            param_hint="'-f' / '--filename'",
        )

    if not isinstance(entry, dict):
        raise click.BadParameter(
            f"Skeleton entry {idx + 1} must be a dictionary, got {type(entry).__name__}",
            param_hint="'-f' / '--filename'",
        )

    # Check for unknown keys
    unknown_keys = set(entry.keys()) - ALLOWED_SKELETON_KEYS
    if unknown_keys:
        raise click.BadParameter(
            f"Unknown key(s) in skeleton entry {idx + 1}: {unknown_keys}. "
            f"Allowed keys: {ALLOWED_SKELETON_KEYS}",
            param_hint="'-f' / '--filename'",
        )

    # Normalize entry: convert empty/null values to None
    normalized = {}
    for key in ALLOWED_SKELETON_KEYS:
        value = entry.get(key)

        # Special handling for link_index and skeleton_indices: parse to list[int]
        if key == "link_index":
            normalized[key] = _parse_index_string(
                value, "skeleton", idx, "link_index"
            )
        elif key == "skeleton_indices":
            normalized[key] = _parse_index_string(
                value, "skeleton", idx, "skeleton_indices"
            )
        else:
            # Convert empty string, null, or missing to None
            if value is None or value == "" or value == "null":
                normalized[key] = None
            else:
                normalized[key] = value

    return normalized


def _validate_substituent_entry(entry: dict, idx: int, filename: str) -> dict:
    """
    Validate a single substituent entry.

    Parameters
    ----------
    entry : dict
        Substituent entry from TOML
    idx : int
        Index of the entry (for error messages)
    filename : str
        Path to TOML file (for error messages)

    Returns
    -------
    dict
        Validated and normalized substituent entry
    """
    if entry is None:
        raise click.BadParameter(
            f"Substituent entry {idx + 1} is empty/null",
            param_hint="'-f' / '--filename'",
        )

    if not isinstance(entry, dict):
        raise click.BadParameter(
            f"Substituent entry {idx + 1} must be a dictionary, got {type(entry).__name__}",
            param_hint="'-f' / '--filename'",
        )

    # Check for unknown keys
    unknown_keys = set(entry.keys()) - ALLOWED_SUBSTITUENT_KEYS
    if unknown_keys:
        raise click.BadParameter(
            f"Unknown key(s) in substituent entry {idx + 1}: {unknown_keys}. "
            f"Allowed keys: {ALLOWED_SUBSTITUENT_KEYS}",
            param_hint="'-f' / '--filename'",
        )

    # Normalize entry: convert empty/null values to None
    normalized = {}
    for key in ALLOWED_SUBSTITUENT_KEYS:
        value = entry.get(key)

        # Special handling for link_index: parse to list[int]
        if key == "link_index":
            normalized[key] = _parse_index_string(
                value, "substituent", idx, "link_index"
            )
        else:
            # Convert empty string, null, or missing to None
            if value is None or value == "" or value == "null":
                normalized[key] = None
            else:
                normalized[key] = value

    return normalized
