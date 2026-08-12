"""Direct command-line input adapter for Iterate jobs."""

import ast
import os

import click

from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.iterate import validate_iterate_config

from .algorithms import register_iterate_algorithm_commands
from .builder import build_iterate_job
from .options import (
    click_iterate_combination_options,
    click_iterate_execution_options,
    validate_iterate_output_options,
)

DIRECT_PARAM_HINT = "'iterate direct'"


def _literal_list(value: str, option_name: str) -> list:
    try:
        parsed = ast.literal_eval(value)
    except (SyntaxError, ValueError) as exc:
        raise click.BadParameter(
            f"{option_name} must be a Python list literal: {exc}.",
            param_hint=option_name,
        )

    if not isinstance(parsed, list):
        raise click.BadParameter(
            f"{option_name} must be a Python list literal.",
            param_hint=option_name,
        )
    if len(parsed) == 0:
        raise click.BadParameter(
            f"{option_name} cannot be an empty list.",
            param_hint=option_name,
        )
    return parsed


def _parse_flat_list(value: str, option_name: str) -> list:
    parsed = _literal_list(value, option_name)
    if any(isinstance(item, list) for item in parsed):
        raise click.BadParameter(
            f"{option_name} must be a flat list.",
            param_hint=option_name,
        )
    return parsed


def _parse_skeleton_groups(
    value: str,
) -> tuple[list | None, list[list] | None]:
    parsed = _literal_list(value, "-skg / --skeleton-groups")
    has_nested = any(isinstance(item, list) for item in parsed)
    has_flat = any(not isinstance(item, list) for item in parsed)

    if has_nested and has_flat:
        raise click.BadParameter(
            "-skg / --skeleton-groups cannot mix integers and nested lists.",
            param_hint="-skg / --skeleton-groups",
        )

    if not has_nested:
        return parsed, None

    slots = []
    for slot_idx, slot in enumerate(parsed, start=1):
        if not isinstance(slot, list):
            raise click.BadParameter(
                "-skg / --skeleton-groups cannot mix integers and nested lists.",
                param_hint="-skg / --skeleton-groups",
            )
        if len(slot) == 0:
            raise click.BadParameter(
                f"-skg / --skeleton-groups slot {slot_idx} cannot be empty.",
                param_hint="-skg / --skeleton-groups",
            )
        if any(isinstance(item, list) for item in slot):
            raise click.BadParameter(
                "-skg / --skeleton-groups cannot be nested more than two levels.",
                param_hint="-skg / --skeleton-groups",
            )
        slots.append(slot)
    return None, slots


def _none_to_default(value: str | None) -> str | None:
    if value is None or value.lower() == "none":
        return None
    return value


def _require_count(
    name: str, values: tuple, expected: int | None = None
) -> None:
    count = len(values)
    if expected is None:
        if count == 0:
            raise click.UsageError(f"At least one {name} entry is required.")
        return
    if count != expected:
        raise click.UsageError(
            f"{name} count mismatch: expected {expected}, got {count}."
        )


def _optional_tuple(values: tuple, expected: int, name: str) -> tuple:
    if len(values) == 0:
        return (None,) * expected
    if len(values) != expected:
        raise click.UsageError(
            f"{name} count mismatch: expected {expected}, got {len(values)}."
        )
    return values


def build_iterate_direct_config(
    *,
    skeleton_files: tuple[str, ...],
    skeleton_groups: tuple[str, ...],
    skeleton_indices: tuple[str, ...],
    skeleton_labels: tuple[str, ...],
    substituent_files: tuple[str, ...],
    substituent_indices: tuple[int, ...],
    substituent_groups: tuple[str, ...],
    substituent_labels: tuple[str, ...],
    path_base_dir: str,
) -> dict:
    """Convert ordered direct-input tuples into normalized standard config."""
    _require_count("skeleton file", skeleton_files)
    _require_count("skeleton group", skeleton_groups, len(skeleton_files))
    _require_count("substituent file", substituent_files)
    _require_count(
        "substituent index", substituent_indices, len(substituent_files)
    )
    _require_count(
        "substituent group", substituent_groups, len(substituent_files)
    )

    skeleton_indices = _optional_tuple(
        skeleton_indices, len(skeleton_files), "skeleton indices"
    )
    skeleton_labels = _optional_tuple(
        skeleton_labels, len(skeleton_files), "skeleton label"
    )
    substituent_labels = _optional_tuple(
        substituent_labels, len(substituent_files), "substituent label"
    )

    skeletons = []
    next_group = 1
    for file_path, groups_text, indices_text, label_text in zip(
        skeleton_files, skeleton_groups, skeleton_indices, skeleton_labels
    ):
        link_index, slot_indices = _parse_skeleton_groups(groups_text)
        entry = {
            "file_path": file_path,
            "label": _none_to_default(label_text),
            "skeleton_indices": (
                None
                if indices_text is None or indices_text.lower() == "none"
                else _parse_flat_list(
                    indices_text, "-ski / --skeleton-indices"
                )
            ),
        }
        if slot_indices is None:
            entry["link_index"] = link_index
            entry["slots"] = None
            next_group += 1
        else:
            entry["link_index"] = None
            entry["slots"] = []
            for indices in slot_indices:
                entry["slots"].append(
                    {"group": next_group, "link_indices": indices}
                )
                next_group += 1
        skeletons.append(entry)

    substituents = []
    for file_path, link_index, groups_text, label_text in zip(
        substituent_files,
        substituent_indices,
        substituent_groups,
        substituent_labels,
    ):
        sub_groups = _parse_flat_list(
            groups_text, "-subg / --substituent-groups"
        )
        substituents.append(
            {
                "file_path": file_path,
                "label": _none_to_default(label_text),
                "link_index": [link_index],
                "groups": sub_groups,
            }
        )

    return validate_iterate_config(
        {"skeletons": skeletons, "substituents": substituents},
        param_hint=DIRECT_PARAM_HINT,
        path_base_dir=path_base_dir,
    )


def _build_iterate_direct_job(ctx, cli_algorithm_name=None, cli_options=None):
    """Build an Iterate job after converting stored direct-input tuples."""
    data = ctx.obj["iterate"]
    if data.get("config") is None:
        data["config"] = build_iterate_direct_config(
            skeleton_files=data["skeleton_files"],
            skeleton_groups=data["skeleton_groups"],
            skeleton_indices=data["skeleton_indices"],
            skeleton_labels=data["skeleton_labels"],
            substituent_files=data["substituent_files"],
            substituent_indices=data["substituent_indices"],
            substituent_groups=data["substituent_groups"],
            substituent_labels=data["substituent_labels"],
            path_base_dir=data["path_base_dir"],
        )
    return build_iterate_job(
        ctx,
        cli_algorithm_name=cli_algorithm_name,
        cli_options=cli_options,
    )


@click.group(name="direct", cls=MyGroup, invoke_without_command=True)
@click.option(
    "-skf",
    "--skeleton-file",
    "skeleton_file",
    multiple=True,
    type=str,
    help="Skeleton molecule file. Repeat with matching -skg entries.",
)
@click.option(
    "-skg",
    "--skeleton-groups",
    "skeleton_groups",
    multiple=True,
    type=str,
    help="Skeleton attachment groups as [1,2,3] shorthand or [[1,2],[3,4]] slots.",
)
@click.option(
    "-ski",
    "--skeleton-indices",
    "skeleton_indices",
    multiple=True,
    type=str,
    help="Optional skeleton core indices as a flat list; use none for a placeholder.",
)
@click.option(
    "-skl",
    "--skeleton-label",
    "skeleton_label",
    multiple=True,
    type=str,
    help="Optional skeleton label; use none for the default skeletonN label.",
)
@click.option(
    "-subf",
    "--substituent-file",
    "substituent_file",
    multiple=True,
    type=str,
    help="Substituent molecule file. Repeat with matching -subi and -subg entries.",
)
@click.option(
    "-subi",
    "--substituent-index",
    "substituent_index",
    multiple=True,
    type=int,
    help="Single 1-based substituent link atom.",
)
@click.option(
    "-subg",
    "--substituent-groups",
    "substituent_groups",
    multiple=True,
    type=str,
    help="Global skeleton groups this substituent can fill, e.g. [1,3].",
)
@click.option(
    "-subl",
    "--substituent-label",
    "substituent_label",
    multiple=True,
    type=str,
    help="Optional substituent label; use none for the default substituentN label.",
)
@click_iterate_execution_options
@click_iterate_combination_options
@click_job_options
@click.pass_context
def direct_cmd(
    ctx,
    skeleton_file,
    skeleton_groups,
    skeleton_indices,
    skeleton_label,
    substituent_file,
    substituent_index,
    substituent_groups,
    substituent_label,
    timeout,
    outputfile,
    directory,
    separate_outputs,
    combination_mode,
    max_substituted_sites,
    skip_completed,
):
    """
    Run iterate jobs from ordered command-line skeleton/substituent entries.

    ``-skg [1,2,3]`` is shorthand for one implicit group. ``-skg [[1,2]]``
    keeps one explicit slot. Shells such as zsh may require quoting these
    list arguments or prefixing the command with ``noglob``.
    """
    ctx.ensure_object(dict)
    directory = validate_iterate_output_options(
        ctx,
        outputfile=outputfile,
        directory=directory,
        separate_outputs=separate_outputs,
    )
    ctx.obj["iterate"] = {
        "config": None,
        "config_file": None,
        "skeleton_files": skeleton_file,
        "skeleton_groups": skeleton_groups,
        "skeleton_indices": skeleton_indices,
        "skeleton_labels": skeleton_label,
        "substituent_files": substituent_file,
        "substituent_indices": substituent_index,
        "substituent_groups": substituent_groups,
        "substituent_labels": substituent_label,
        "path_base_dir": os.getcwd(),
        "timeout": timeout,
        "outputfile": outputfile,
        "directory": directory,
        "separate_outputs": separate_outputs,
        "combination_mode": combination_mode,
        "max_substituted_sites": max_substituted_sites,
        "skip_completed": skip_completed,
    }

    if ctx.invoked_subcommand is None:
        return _build_iterate_direct_job(ctx)


register_iterate_algorithm_commands(direct_cmd, _build_iterate_direct_job)
