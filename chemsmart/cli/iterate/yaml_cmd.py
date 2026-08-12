"""
CLI subcommands for running iterate jobs from YAML configuration files.

The ``yaml`` command represents the *input format* layer (it reads a YAML
configuration file). The position-optimization *algorithm* is a separate
layer, selected via an optional algorithm subcommand:

\b
    yaml etkdg   RDKit ETKDGv3 algorithm (local by default; --global).
    yaml jlgo    Joint Lagrange Geometry Optimization (multi-substituent).

When no algorithm subcommand is given, the algorithm declared in the YAML
``algorithm`` block (or the built-in default ``etkdg``) is
used. Algorithm parameters live only on their respective subcommands so that
different algorithms cannot share the same option namespace.
"""

import functools
import logging
import os

import click
import yaml

from chemsmart.cli.job import click_filename_options, click_job_options
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.iterate import (
    generate_yaml_template,
    validate_yaml_config,
)

from .algorithms import register_iterate_algorithm_commands
from .builder import build_iterate_job
from .options import (
    click_iterate_combination_options,
    click_iterate_execution_options,
    validate_iterate_output_options,
)

logger = logging.getLogger(__name__)


def click_yaml_common_options(f):
    """
    Common (algorithm-agnostic) Click options for the ``yaml`` command group.

    Algorithm-specific parameters are intentionally *not* included here; they
    are defined on the individual algorithm subcommands instead.
    """

    @click.option(
        "-g",
        "--generate-template",
        is_flag=False,
        flag_value="iterate_template.yaml",
        default=None,
        type=str,
        help="Generate a template configuration file and exit. "
        "Optionally specify output path (default: iterate_template.yaml).",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def _load_yaml_config(filename: str) -> dict:
    """Load, validate, and normalize an Iterate YAML configuration."""
    # Validate filename (deferred here so that subcommand '--help' works
    # without requiring '-f').
    if not filename:
        raise click.BadParameter(
            "A configuration file is required.",
            param_hint="'-f' / '--filename'",
        )

    if not os.path.exists(filename):
        raise click.BadParameter(
            f"File '{filename}' does not exist.",
            param_hint="'-f' / '--filename'",
        )

    if not filename.endswith((".yaml", ".yml")):
        raise click.BadParameter(
            f"File '{filename}' must be a YAML file "
            f"(ending with .yaml or .yml).",
            param_hint="'-f' / '--filename'",
        )

    # Load and validate the YAML configuration file. A parse error or a
    # non-mapping top level is a usage/configuration error -> Click exit 2.
    try:
        with open(filename, "r") as f:
            raw_config = yaml.safe_load(f)
    except yaml.YAMLError as exc:
        raise click.BadParameter(
            f"File '{filename}' is not valid YAML: {exc}",
            param_hint="'-f' / '--filename'",
        )

    if raw_config is None:
        raw_config = {}

    if not isinstance(raw_config, dict):
        raise click.BadParameter(
            f"File '{filename}' must contain a YAML mapping at the top "
            f"level, got {type(raw_config).__name__}.",
            param_hint="'-f' / '--filename'",
        )

    config = validate_yaml_config(raw_config, filename)
    logger.debug(f"Loaded YAML configuration from '{filename}'")
    return config


def _is_sub_invocation(ctx: click.Context) -> bool:
    """Return whether the current command is nested below ``sub``."""
    current = ctx
    while current is not None:
        if current.info_name == "sub":
            return True
        current = current.parent
    return False


def _build_yaml_iterate_job(ctx, cli_algorithm_name=None, cli_options=None):
    """Load YAML config if needed, then build the Iterate job."""
    data = ctx.obj["iterate"]
    if data.get("config") is None:
        filename = data.get("filename")
        data["config"] = _load_yaml_config(filename)
        data["config_file"] = filename
    return build_iterate_job(
        ctx,
        cli_algorithm_name=cli_algorithm_name,
        cli_options=cli_options,
    )


@click.group(name="yaml", cls=MyGroup, invoke_without_command=True)
@click_yaml_common_options
@click_iterate_execution_options
@click_iterate_combination_options
@click_job_options
@click_filename_options
@click.pass_context
def yaml_cmd(
    ctx,
    filename,
    timeout,
    outputfile,
    generate_template,
    directory,
    separate_outputs,
    combination_mode,
    max_substituted_sites,
    skip_completed,
    **kwargs,
):
    """
    Run iterate jobs from a YAML configuration file.

    The YAML file defines skeletons and substituents. The optimization
    algorithm can be declared in the YAML ``algorithm`` block or selected
    via an algorithm subcommand (``jlgo``, ``etkdg``). When both are
    given, the CLI takes precedence.

    All skeletons participate in global contiguous group numbering;
    skeletons with link_index occupy one implicit group each, while
    skeletons with slots occupy one group per slot.

    Examples:

    \b
    chemsmart run iterate yaml -f config.yaml
    chemsmart run iterate yaml -f config.yaml jlgo
    chemsmart run iterate yaml -f config.yaml jlgo \\
        --no-adaptive-sampling \\
        --link-sphere-samples 48 \\
        --axial-samples 4
    chemsmart run iterate yaml -f config.yaml etkdg \\
        --num-conformers 50 --random-seed 1
    chemsmart run iterate yaml -g
    chemsmart run iterate yaml -g my_config.yaml
    """
    ctx.ensure_object(dict)

    # Handle -g option: generate template and exit
    if generate_template is not None:
        if _is_sub_invocation(ctx):
            raise click.UsageError(
                "--generate-template is a local utility and cannot be "
                "submitted to a scheduler."
            )
        template_path = generate_yaml_template(
            generate_template, overwrite=False
        )
        click.echo(f"Generated template: {template_path}")
        ctx.exit(0)

    directory = validate_iterate_output_options(
        ctx,
        outputfile=outputfile,
        directory=directory,
        separate_outputs=separate_outputs,
    )

    # Store shared parameters for the algorithm subcommands / shared executor.
    # Filename validation and YAML loading are deferred to the executor so
    # that 'yaml <algorithm> --help' works without requiring '-f'.
    ctx.obj["iterate"] = {
        "filename": filename,
        "config": None,
        "config_file": filename,
        "timeout": timeout,
        "outputfile": outputfile,
        "directory": directory,
        "separate_outputs": separate_outputs,
        "combination_mode": combination_mode,
        "max_substituted_sites": max_substituted_sites,
        "skip_completed": skip_completed,
    }

    # No algorithm subcommand: run with the YAML/default algorithm.
    if ctx.invoked_subcommand is None:
        return _build_yaml_iterate_job(ctx)


register_iterate_algorithm_commands(yaml_cmd, _build_yaml_iterate_job)
