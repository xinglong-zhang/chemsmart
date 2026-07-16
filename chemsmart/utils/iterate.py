"""Utility functions for iterate module."""

import logging
import os
import re

import click

from chemsmart.utils.repattern import safe_label_pattern
from chemsmart.utils.utils import get_list_from_string_range

logger = logging.getLogger(__name__)

# Define allowed keys for configuration validation
ALLOWED_TOP_LEVEL_KEYS = {"skeletons", "substituents", "algorithm"}
ALLOWED_ALGORITHM_KEYS = {"name", "options"}
ALLOWED_SKELETON_KEYS = {
    "file_path",
    "label",
    "link_index",
    "skeleton_indices",
    "slots",
}
ALLOWED_SUBSTITUENT_KEYS = {
    "file_path",
    "label",
    "link_index",
    "groups",
}
ALLOWED_EMBEDDED_SLOT_KEYS = {
    "group",
    "link_indices",
}

ITERATE_YAML_TEMPLATE = """\
# Chemsmart Iterate Configuration Template (YAML)
# ================================================
#
# The configuration is split into two independent layers:
#   1. Input format layer  -> chosen via the CLI subcommand ('yaml', 'cdxml').
#      This file *is* the YAML input; it lists skeletons and substituents.
#   2. Algorithm layer      -> how substituent positions are optimized.
#      Configured via the optional top-level 'algorithm' block below, and/or
#      via a CLI algorithm subcommand ('yaml lagrange', 'yaml etkdg').
#
# ALGORITHM PRIORITY (low -> high):
#   built-in default (lagrange_multipliers)
#     < 'algorithm' block in this file
#     < CLI algorithm subcommand name
#     < CLI algorithm options explicitly passed on the command line.
# If the CLI selects a different algorithm than this file, the options in
# this file (which belong to the other algorithm) are ignored.
#
# Each skeleton contributes to a global contiguous 1-based group sequence:
#   - Skeleton with "link_index": occupies one implicit group per skeleton.
#   - Skeleton with "slots":      occupies one group per slot.
#
# Substituents declare which groups they can fill via the "groups" field.
#
# RULES:
#   - A skeleton must have EITHER link_index OR slots, not both.
#   - Group numbers must be globally contiguous starting from 1.
#   - All substituents must specify "groups".
#   - Atom indices are 1-based (first atom is index 1).
#   - If skeleton_indices is specified, it must include all link atoms.

# === Algorithm layer (optional; default is lagrange_multipliers) ===
algorithm:
  name: lagrange_multipliers      # or 'etkdg'
  options:
    sphere_direction_samples_num: 96
    axial_rotations_sample_num: 6
  # --- ETKDG (RDKit ETKDGv3; local by default) ---
  # name: etkdg
  # options:
  #   use_global_optimization: false  # true re-embeds every atom
  #   num_conformers: 10              # best (lowest-energy) one is kept
  #   random_seed: 42
  #   max_iterations: 2000
  #   force_field: none               # or uff / mmff94 / mmff94s

skeletons:
  # === link_index shorthand (auto-assigned group 1) ===
  - file_path: "/path/to/skeleton1.xyz"
    label: "benzene_core"
    link_index: "6,7,8"
    skeleton_indices: "1-5"    # optional: atoms to keep

  # === Explicit slots (groups 2, 3) ===
  # - file_path: "/path/to/complex_skeleton.xyz"
  #   label: "complex_skel"
  #   skeleton_indices: "1-30"  # optional
  #   slots:
  #     - group: 2
  #       link_indices: "10,11"
  #     - group: 3
  #       link_indices: "25"

substituents:
  - file_path: "/path/to/methyl.xyz"
    label: "methyl"
    link_index: 1
    groups: [1]          # available for group 1 (skeleton1)

  # - file_path: "/path/to/ethyl.xyz"
  #   label: "ethyl"
  #   link_index: 1
  #   groups: [1, 2, 3]  # available for groups 1, 2, and 3
"""


def generate_yaml_template(
    output_path: str = "iterate_template.yaml", overwrite: bool = False
) -> str:
    """
    Generate a YAML template configuration file for iterate jobs.

    Parameters
    ----------
    output_path : str
        Path to write the template file. Default is 'iterate_template.yaml'.
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
    if not output_path.endswith((".yaml", ".yml")):
        output_path = f"{output_path}.yaml"

    if os.path.exists(output_path) and not overwrite:
        raise FileExistsError(
            f"File '{output_path}' already exists. "
            f"Use --overwrite or specify a different filename."
        )

    with open(output_path, "w") as f:
        f.write(ITERATE_YAML_TEMPLATE)

    logger.info(f"Generated YAML template: {output_path}")
    return output_path


def _parse_index_string(
    value, context: str, field_name: str
) -> list[int] | None:
    """
    Parse index string into a list of integers using utility function.
    Wrapper to handle empty values and single integers.

    Parameters
    ----------
    value : int, str, or None
        The value to parse
    context : str
        Human-readable caller context for error messages,
        e.g. "Skeleton entry 1" or "Skeleton 2, slot 3"
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
                f"{context}: "
                f"Found invalid index {value} in '{field_name}'. "
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
            # Already a list of ints
            parsed_indices = value

        if parsed_indices is not None:
            # S2 Check: Validate not empty
            if len(parsed_indices) == 0:
                raise click.BadParameter(
                    f"{context}: "
                    f"Found empty list in '{field_name}'. "
                    f"At least one index must be provided.",
                    param_hint="'-f' / '--filename'",
                )

            # S2 Check: Validate positive non-zero indices
            if any(i <= 0 for i in parsed_indices):
                raise click.BadParameter(
                    f"{context}: "
                    f"Found invalid index <= 0 in '{field_name}'. "
                    f"All indices must be positive (1-based). Found: {parsed_indices}",
                    param_hint="'-f' / '--filename'",
                )
            return parsed_indices

    except click.BadParameter:
        raise
    except Exception:
        raise click.BadParameter(
            f"{context}: "
            f"Invalid format '{value}' in '{field_name}'. "
            f"Expected integer, comma-separated list, or range (e.g. '1-5').",
            param_hint="'-f' / '--filename'",
        )

    # If it's not None, not Int, and not String
    # (or string parsing failed silently elsewhere)
    raise click.BadParameter(
        f"{context}: "
        f"'{field_name}' has invalid type {type(value).__name__}.",
        param_hint="'-f' / '--filename'",
    )


def _resolve_file_path(file_path, config_dir: str):
    """Resolve a molecule file path against the config file's directory.

    Relative ``file_path`` values are resolved relative to the directory of
    the YAML configuration file (not the current working directory), so
    molecule files can live next to their config and do not depend on where
    the command is invoked. Absolute paths and empty values are returned
    unchanged.
    """
    if not file_path:
        return file_path
    if os.path.isabs(file_path):
        return file_path
    return os.path.normpath(os.path.join(config_dir, file_path))


def validate_yaml_config(config: dict, filename: str) -> dict:
    """
    Validate YAML configuration and normalize values.
    Supports link_index (shorthand) and slots (explicit) formats.
    All skeletons participate in global contiguous group numbering.

    Parameters
    ----------
    config : dict
        Raw configuration dictionary from parsed YAML file
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
            f"Unknown top-level key(s) in configuration: "
            f"{unknown_top_keys}. "
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
                f"'skeletons' must be a list, "
                f"got {type(skeletons).__name__}",
                param_hint="'-f' / '--filename'",
            )
        else:
            validated_config["skeletons"] = [
                _validate_skeleton_entry(entry, skel_idx)
                for skel_idx, entry in enumerate(skeletons)
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
                f"'substituents' must be a list, "
                f"got {type(substituents).__name__}",
                param_hint="'-f' / '--filename'",
            )
        else:
            validated_config["substituents"] = [
                _validate_substituent_entry(entry, sub_idx)
                for sub_idx, entry in enumerate(substituents)
            ]
    else:
        validated_config["substituents"] = []

    # Resolve relative molecule paths against the config file's directory so
    # that file_path values never depend on the current working directory.
    # The original (as-written) value is preserved for diagnostics.
    config_dir = os.path.dirname(os.path.abspath(filename))
    for entry in (
        validated_config["skeletons"] + validated_config["substituents"]
    ):
        raw_path = entry.get("file_path")
        entry["file_path_raw"] = raw_path
        entry["file_path"] = _resolve_file_path(raw_path, config_dir)

    # Validate and normalize the optional algorithm block
    validated_config["algorithm"] = _validate_algorithm_entry(
        config.get("algorithm"), filename
    )

    # Enforce link_index XOR slots for each skeleton
    for idx, skel in enumerate(validated_config["skeletons"]):
        has_link = skel.get("link_index") is not None
        has_slots = skel.get("slots") is not None
        if has_link and has_slots:
            raise click.BadParameter(
                f"Skeleton entry {idx + 1} ('{skel.get('label', 'unnamed')}'): "
                f"Cannot have both 'link_index' and 'slots'.",
                param_hint="'-f' / '--filename'",
            )
        if not has_link and not has_slots:
            raise click.BadParameter(
                f"Skeleton entry {idx + 1} ('{skel.get('label', 'unnamed')}'): "
                f"Must have either 'link_index' or 'slots'.",
                param_hint="'-f' / '--filename'",
            )

    # Common Business Logic Validation
    _validate_business_logic(
        validated_config["skeletons"],
        validated_config["substituents"],
        filename,
    )

    # YAML-specific cross-validation (global group numbering)
    _validate_yaml_multi_site_logic(validated_config, filename)

    return validated_config


def _validate_algorithm_entry(entry, filename: str):
    """
    Validate and normalize the optional top-level ``algorithm`` block.

    Parameters
    ----------
    entry : dict or None
        Raw ``algorithm`` block from the parsed YAML file.
    filename : str
        Path to configuration file (for error messages).

    Returns
    -------
    dict or None
        Normalized ``{"name": <canonical>, "options": {...}}`` mapping, or
        ``None`` when no algorithm block is present.

    Raises
    ------
    click.BadParameter
        If the block is malformed, references an unknown algorithm, or
        contains options not supported by that algorithm.
    """
    # Imported lazily to keep importing this module lightweight (the
    # algorithm registry lives alongside the job settings).
    from chemsmart.jobs.iterate.settings import (
        available_algorithm_names,
        normalize_algorithm_name,
        validate_algorithm_options,
    )

    if entry is None:
        return None

    if not isinstance(entry, dict):
        raise click.BadParameter(
            f"'algorithm' must be a mapping, got {type(entry).__name__}.",
            param_hint="'-f' / '--filename'",
        )

    unknown_keys = set(entry.keys()) - ALLOWED_ALGORITHM_KEYS
    if unknown_keys:
        raise click.BadParameter(
            f"Unknown key(s) in 'algorithm' block: {unknown_keys}. "
            f"Allowed keys: {ALLOWED_ALGORITHM_KEYS}",
            param_hint="'-f' / '--filename'",
        )

    name = entry.get("name")
    if name is None:
        raise click.BadParameter(
            "'algorithm' block is missing required field 'name'. "
            f"Available algorithms: {available_algorithm_names()}.",
            param_hint="'-f' / '--filename'",
        )

    try:
        canonical_name = normalize_algorithm_name(name)
    except ValueError as exc:
        raise click.BadParameter(str(exc), param_hint="'-f' / '--filename'")

    options = entry.get("options")
    if options is None:
        options = {}
    elif not isinstance(options, dict):
        raise click.BadParameter(
            f"'algorithm.options' must be a mapping, "
            f"got {type(options).__name__}.",
            param_hint="'-f' / '--filename'",
        )

    try:
        validated_options = validate_algorithm_options(canonical_name, options)
    except ValueError as exc:
        raise click.BadParameter(str(exc), param_hint="'-f' / '--filename'")

    return {"name": canonical_name, "options": validated_options}


def _validate_business_logic(
    skeletons: list, substituents: list, filename: str
):
    """
    Orchestrate business logic checks for all entities.
    """
    _validate_skeletons_logic(skeletons, filename)

    _validate_substituents_logic(substituents, filename)

    _validate_unique_labels(skeletons, substituents, filename)


def _validate_unique_labels(
    skeletons: list, substituents: list, filename: str
):
    """Reject duplicate skeleton or substituent labels.

    Labels are used to build combination labels and output filenames, so a
    duplicate would silently overwrite another structure's result/file.
    """

    def _check(entries: list, kind: str):
        seen: dict = {}
        for idx, entry in enumerate(entries):
            label = entry.get("label")
            if label is None:
                continue
            if label in seen:
                raise click.BadParameter(
                    f"{kind} entry {idx + 1}: duplicate label '{label}' "
                    f"(already used by {kind.lower()} entry "
                    f"{seen[label] + 1}). Labels must be unique to avoid "
                    f"output collisions.",
                    param_hint=filename,
                )
            seen[label] = idx

    _check(skeletons, "Skeleton")
    _check(substituents, "Substituent")


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

    # Rule 4: link_index values must not repeat within one skeleton
    _rule_no_duplicate_link_index(entry, idx, filename)

    # Rule 5: a connection site must not be shared across slots
    _rule_no_overlapping_slot_links(entry, idx, filename)


def _rule_no_duplicate_link_index(entry: dict, idx: int, filename: str):
    """Rule: link_index values must be unique within a skeleton."""
    link_indices = entry.get("link_index")
    if not link_indices:
        return
    seen = set()
    duplicates = sorted({i for i in link_indices if i in seen or seen.add(i)})
    if duplicates:
        raise click.BadParameter(
            f"Skeleton entry {idx + 1} "
            f"(label='{entry.get('label', 'unnamed')}'): "
            f"Duplicate link_index value(s) {duplicates}. "
            f"Each connection site must be listed once.",
            param_hint=filename,
        )


def _rule_no_overlapping_slot_links(entry: dict, idx: int, filename: str):
    """Rule: a link index must belong to at most one slot."""
    slots = entry.get("slots")
    if not slots:
        return
    owner: dict = {}
    for slot in slots:
        for link in slot["link_indices"]:
            if link in owner:
                raise click.BadParameter(
                    f"Skeleton entry {idx + 1} "
                    f"(label='{entry.get('label', 'unnamed')}'): "
                    f"link index {link} is used by multiple slots "
                    f"(groups {owner[link]} and {slot['group']}). "
                    f"Each connection site may belong to only one slot.",
                    param_hint=filename,
                )
            owner[link] = slot["group"]


def _rule_label_syntax(entry: dict, idx: int, entry_type: str, filename: str):
    """
    Rule: Label must only contain safe characters [a-zA-Z0-9_-.].
    """
    label = entry.get("label")
    if label:
        if not re.match(safe_label_pattern, label):
            raise click.BadParameter(
                f"{entry_type} entry {idx + 1} label '{label}': "
                f"Contains invalid characters. "
                f"Allowed characters: a-z, A-Z, 0-9, _, -, .",
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


def _validate_skeleton_entry(entry: dict, skel_idx: int) -> dict:
    """
    Validate a single skeleton entry.

    Parameters
    ----------
    entry : dict
        Skeleton entry from config
    skel_idx : int
        Index of the entry (for error messages)

    Returns
    -------
    dict
        Validated and normalized skeleton entry
    """
    if entry is None:
        raise click.BadParameter(
            f"Skeleton entry {skel_idx + 1} is empty/null",
            param_hint="'-f' / '--filename'",
        )

    if not isinstance(entry, dict):
        raise click.BadParameter(
            f"Skeleton entry {skel_idx + 1} must be a dictionary, "
            f"got {type(entry).__name__}",
            param_hint="'-f' / '--filename'",
        )

    # Check for unknown keys
    unknown_keys = set(entry.keys()) - ALLOWED_SKELETON_KEYS
    if unknown_keys:
        raise click.BadParameter(
            f"Unknown key(s) in skeleton entry {skel_idx + 1}: {unknown_keys}. "
            f"Allowed keys: {ALLOWED_SKELETON_KEYS}",
            param_hint="'-f' / '--filename'",
        )

    # Normalize common fields
    normalized = {}
    for key in ("file_path", "label"):
        value = entry.get(key)
        if value is None or value == "" or value == "null":
            normalized[key] = None
        else:
            normalized[key] = value

    # Parse link_index
    skel_context = f"Skeleton entry {skel_idx + 1}"
    normalized["link_index"] = _parse_index_string(
        entry.get("link_index"), skel_context, "link_index"
    )

    # Parse skeleton_indices
    normalized["skeleton_indices"] = _parse_index_string(
        entry.get("skeleton_indices"), skel_context, "skeleton_indices"
    )

    # Handle slots field
    slots_raw = entry.get("slots")
    if slots_raw is not None:
        if not isinstance(slots_raw, list):
            raise click.BadParameter(
                f"Skeleton entry {skel_idx + 1}: 'slots' must be a list.",
                param_hint="'-f' / '--filename'",
            )
        if len(slots_raw) == 0:
            raise click.BadParameter(
                f"Skeleton entry {skel_idx + 1}: 'slots' cannot be empty.",
                param_hint="'-f' / '--filename'",
            )
        normalized["slots"] = [
            _validate_embedded_slot_entry(slot, skel_idx, slot_idx)
            for slot_idx, slot in enumerate(slots_raw)
        ]
    else:
        normalized["slots"] = None

    return normalized


def _validate_substituent_entry(entry: dict, sub_idx: int) -> dict:
    """
    Validate a single substituent entry.

    Parameters
    ----------
    entry : dict
        Substituent entry from config
    sub_idx : int
        Index of the entry (for error messages)

    Returns
    -------
    dict
        Validated and normalized substituent entry
    """
    if entry is None:
        raise click.BadParameter(
            f"Substituent entry {sub_idx + 1} is empty/null",
            param_hint="'-f' / '--filename'",
        )

    if not isinstance(entry, dict):
        raise click.BadParameter(
            f"Substituent entry {sub_idx + 1} must be a dictionary, "
            f"got {type(entry).__name__}",
            param_hint="'-f' / '--filename'",
        )

    # Check for unknown keys
    unknown_keys = set(entry.keys()) - ALLOWED_SUBSTITUENT_KEYS
    if unknown_keys:
        raise click.BadParameter(
            f"Unknown key(s) in substituent entry {sub_idx + 1}: {unknown_keys}. "
            f"Allowed keys: {ALLOWED_SUBSTITUENT_KEYS}",
            param_hint="'-f' / '--filename'",
        )

    # Normalize common fields
    normalized = {}
    for key in ("file_path", "label"):
        value = entry.get(key)
        if value is None or value == "" or value == "null":
            normalized[key] = None
        else:
            normalized[key] = value

    # Parse link_index
    normalized["link_index"] = _parse_index_string(
        entry.get("link_index"),
        f"Substituent entry {sub_idx + 1}",
        "link_index",
    )

    # Handle groups field
    groups_raw = entry.get("groups")
    if groups_raw is not None:
        if not isinstance(groups_raw, list):
            raise click.BadParameter(
                f"Substituent entry {sub_idx + 1}: 'groups' must be a list.",
                param_hint="'-f' / '--filename'",
            )
        if len(groups_raw) == 0:
            raise click.BadParameter(
                f"Substituent entry {sub_idx + 1}: 'groups' cannot be empty.",
                param_hint="'-f' / '--filename'",
            )
        if not all(isinstance(g, int) and g > 0 for g in groups_raw):
            raise click.BadParameter(
                f"Substituent entry {sub_idx + 1}: all 'groups' entries "
                f"must be positive integers.",
                param_hint="'-f' / '--filename'",
            )
        normalized["groups"] = groups_raw
    else:
        normalized["groups"] = None

    return normalized


def _validate_embedded_slot_entry(entry, skel_idx: int, slot_idx: int) -> dict:
    """
    Validate a single embedded slot entry within a skeleton.

    Parameters
    ----------
    entry : dict
        Slot entry (should have 'group' and 'link_indices')
    skel_idx : int
        Parent skeleton index (for error messages)
    slot_idx : int
        Slot index within the skeleton (for error messages)

    Returns
    -------
    dict
        Validated slot with {"group": int, "link_indices": list[int]}
    """
    prefix = f"Skeleton {skel_idx + 1}, slot {slot_idx + 1}"

    if entry is None:
        raise click.BadParameter(
            f"{prefix}: is empty/null.",
            param_hint="'-f' / '--filename'",
        )

    if not isinstance(entry, dict):
        raise click.BadParameter(
            f"{prefix}: must be a dictionary, " f"got {type(entry).__name__}.",
            param_hint="'-f' / '--filename'",
        )

    unknown = set(entry.keys()) - ALLOWED_EMBEDDED_SLOT_KEYS
    if unknown:
        raise click.BadParameter(
            f"{prefix}: unknown key(s) {unknown}. "
            f"Allowed: {ALLOWED_EMBEDDED_SLOT_KEYS}",
            param_hint="'-f' / '--filename'",
        )

    # group: required, positive integer
    group = entry.get("group")
    if group is None:
        raise click.BadParameter(
            f"{prefix}: missing required field 'group'.",
            param_hint="'-f' / '--filename'",
        )
    if not isinstance(group, int) or group <= 0:
        raise click.BadParameter(
            f"{prefix}: 'group' must be a positive integer, "
            f"got {group!r}.",
            param_hint="'-f' / '--filename'",
        )

    # link_indices: required, parse to list[int]
    link_indices = _parse_index_string(
        entry.get("link_indices"), prefix, "link_indices"
    )
    if not link_indices:
        raise click.BadParameter(
            f"{prefix}: 'link_indices' is required "
            f"and must contain at least one index.",
            param_hint="'-f' / '--filename'",
        )

    seen: set = set()
    duplicates = sorted({i for i in link_indices if i in seen or seen.add(i)})
    if duplicates:
        raise click.BadParameter(
            f"{prefix}: 'link_indices' contains duplicate value(s) "
            f"{duplicates}.",
            param_hint="'-f' / '--filename'",
        )

    return {"group": group, "link_indices": link_indices}


def _validate_yaml_multi_site_logic(validated: dict, filename: str):
    """
    Cross-validate YAML configuration with global group numbering.

    ALL skeletons participate in a global contiguous 1-based group
    numbering sequence:
    - Skeletons with ``link_index`` (no slots) each occupy exactly
      one implicit group number.
    - Skeletons with ``slots`` occupy as many group numbers as they
      have slots, and those numbers must match the expected range.

    E.g. skeleton A (link_index) → group 1;
         skeleton B (3 slots)    → groups {2, 3, 4};
         skeleton C (link_index) → group 5.

    Rules:
    1. Group numbers form a contiguous 1-based sequence across all skeletons.
    2. All substituents must have a ``groups`` field.
    3. Substituent ``groups`` values must reference existing group numbers.
    4. For skeletons with slots: if skeleton_indices is specified,
       it must contain all slot link_indices.
    """
    skeletons = validated["skeletons"]
    substituents = validated["substituents"]

    # --- Compute global group sequence ---
    all_groups = {}  # group_number -> skeleton_label
    next_expected_start = 1

    for idx, skel in enumerate(skeletons):
        skel_label = skel.get("label") or f"skeleton entry {idx + 1}"
        slots = skel.get("slots")

        if slots is None:
            # No-slots skeleton: occupies one implicit group
            implicit_group = next_expected_start
            all_groups[implicit_group] = skel_label
            next_expected_start += 1
        else:
            # Slots skeleton: check global uniqueness + contiguous range
            for slot in slots:
                group = slot["group"]
                if group in all_groups:
                    raise click.BadParameter(
                        f"Skeleton '{skel_label}': "
                        f"Group number {group} is already used by "
                        f"skeleton '{all_groups[group]}'. "
                        f"Group numbers must be globally unique.",
                        param_hint=filename,
                    )

            skel_groups = sorted(slot["group"] for slot in slots)
            n_slots = len(skel_groups)
            expected = list(
                range(next_expected_start, next_expected_start + n_slots)
            )

            if skel_groups != expected:
                raise click.BadParameter(
                    f"Skeleton '{skel_label}': slot group numbers must be "
                    f"contiguous 1-based integers. "
                    f"Expected {expected}, got {skel_groups}. "
                    f"Groups are assigned sequentially across all skeletons "
                    f"(skeletons without slots each occupy one group number). "
                    f"Previous skeleton ended at group "
                    f"{next_expected_start - 1}.",
                    param_hint=filename,
                )

            for slot in slots:
                all_groups[slot["group"]] = skel_label
            next_expected_start += n_slots

    # --- Validate all substituents have 'groups' ---
    for idx, sub in enumerate(substituents):
        groups = sub.get("groups")
        if groups is None:
            raise click.BadParameter(
                f"Substituent entry {idx + 1} "
                f"('{sub.get('label', 'unnamed')}'): "
                f"'groups' field is required in YAML format. "
                f"Available groups: {sorted(all_groups.keys())}",
                param_hint=filename,
            )
        for g in groups:
            if g not in all_groups:
                raise click.BadParameter(
                    f"Substituent entry {idx + 1} "
                    f"('{sub.get('label', 'unnamed')}'): "
                    f"Group {g} does not reference any skeleton group. "
                    f"Available groups: {sorted(all_groups.keys())}",
                    param_hint=filename,
                )

    # --- Validate skeleton_indices for slots-skeletons ---
    for idx, skel in enumerate(skeletons):
        slots = skel.get("slots")
        if slots is None:
            continue
        skeleton_indices = skel.get("skeleton_indices")
        if skeleton_indices is not None:
            skel_set = set(skeleton_indices)
            for slot in slots:
                missing = [
                    i for i in slot["link_indices"] if i not in skel_set
                ]
                if missing:
                    raise click.BadParameter(
                        f"Skeleton entry {idx + 1} "
                        f"('{skel.get('label', 'unnamed')}'): "
                        f"Slot group {slot['group']} link_indices {missing} "
                        f"are not in 'skeleton_indices'.",
                        param_hint=filename,
                    )
