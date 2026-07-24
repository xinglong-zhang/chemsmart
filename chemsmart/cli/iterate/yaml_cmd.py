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
import sys

import click
import yaml

from chemsmart.cli.job import click_filename_options, click_job_options
from chemsmart.jobs.iterate.job import IterateJob
from chemsmart.jobs.iterate.settings import (
    IterateJobSettings,
    resolve_algorithm_config,
)
from chemsmart.utils.cli import MyCommand, MyGroup
from chemsmart.utils.iterate import (
    generate_yaml_template,
    validate_yaml_config,
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
    @click.option(
        "--separate-outputs/--no-separate-outputs",
        default=False,
        show_default=True,
        help="Save each structure as a separate XYZ file.",
    )
    @click.option(
        "-t",
        "--timeout",
        default=120,
        type=click.IntRange(min=1),
        show_default=True,
        help="Timeout in seconds for each worker process.",
    )
    @click.option(
        "-cm",
        "--combination-mode",
        default="independent",
        type=click.Choice(
            ["independent", "global"],
            case_sensitive=False,
        ),
        show_default=True,
        help="Combination strategy for skeleton slots. "
        "Each slot specifies a group number; only substituents belonging to "
        "that group are candidates for the slot. "
        "Each position includes a 'None' (keep original) option. "
        "'independent' (default): each slot is expanded separately. "
        "E.g. slot R1 (group 1), slot R2 (group 2); "
        "sub A, B in group 1, sub C in group 2: "
        "R1 → mol(R1=A), mol(R1=B); R2 → mol(R2=C) → 3 structures. "
        "'global': all slots combined via single Cartesian product. "
        "Same setup → mol(R1=A), mol(R1=B), mol(R2=C), "
        "mol(R1=A,R2=C), mol(R1=B,R2=C) → 5 structures.",
    )
    @click.option(
        "-d",
        "--directory",
        default=None,
        type=click.Path(file_okay=False, dir_okay=True),
        help="Directory to save output files. Use only with --separate-outputs.",
    )
    @click.option(
        "-o",
        "--outputfile",
        default=None,
        type=str,
        help="Output filename (without .xyz extension) for generated structures. "
        "Defaults to <configuration_stem>_iterate.xyz. "
        "Use only with --no-separate-outputs.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def _collect_cli_options(options: dict) -> dict:
    """Drop unset CLI options before merging them over YAML settings."""
    return {name: value for name, value in options.items() if value is not None}


def _is_sub_invocation(ctx: click.Context) -> bool:
    """Return whether the current command is nested below ``sub``."""
    current = ctx
    while current is not None:
        if current.info_name == "sub":
            return True
        current = current.parent
    return False


def _build_iterate_job(
    ctx, cli_algorithm_name=None, cli_options=None
) -> IterateJob:
    """
    Shared Job builder for the ``yaml`` group and its algorithm subcommands.

    Resolve the effective algorithm configuration and return an Iterate job.

    Parameters
    ----------
    ctx : click.Context
        Context whose ``obj["iterate"]`` holds the common parameters
        (populated by the ``yaml`` group). The YAML configuration file is
        loaded and validated here.
    cli_algorithm_name : str, optional
        Algorithm name selected via a subcommand (``None`` for the group's
        no-subcommand path).
    cli_options : dict, optional
        Algorithm options explicitly passed on the command line.
    """
    data = ctx.obj["iterate"]
    filename = data["filename"]
    nprocs = ctx.obj["jobrunner"].num_cores
    if nprocs < 1:
        raise click.BadParameter(
            "must be at least 1",
            param_hint="'-n' / '--num-cores'",
        )
    timeout = data["timeout"]
    outputfile = data["outputfile"]
    directory = data["directory"]
    separate_outputs = data["separate_outputs"]
    combination_mode = data["combination_mode"]

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
    logger.debug(f"  Skeletons: {len(config['skeletons'])}")
    logger.debug(f"  Substituents: {len(config['substituents'])}")

    try:
        algorithm_config = resolve_algorithm_config(
            yaml_algorithm=config.get("algorithm"),
            cli_algorithm_name=cli_algorithm_name,
            cli_options=cli_options,
        )
    except ValueError as exc:
        # An unsupported algorithm/option or an out-of-range option value is a
        # usage/configuration error -> Click exit code 2.
        raise click.BadParameter(str(exc))

    logger.debug(
        f"Using algorithm '{algorithm_config.name}' "
        f"with options {algorithm_config.options}"
    )
    logger.debug(f"  Combination mode: {combination_mode}")

    # Create job settings
    job_settings = IterateJobSettings(
        config_file=filename,
        algorithm_config=algorithm_config,
        combination_mode=combination_mode,
    )
    job_settings.skeleton_list = config["skeletons"]
    job_settings.substituent_list = config["substituents"]

    job = IterateJob(
        settings=job_settings,
        jobrunner=ctx.obj.get("jobrunner"),
        nprocs=nprocs,
        timeout=timeout,
        outputfile=outputfile,
        separate_outputs=separate_outputs,
        output_directory=directory,
        command_line=" ".join(sys.argv),
        show_worker_logs=logger.isEnabledFor(logging.DEBUG),
        skip_completed=data["skip_completed"],
    )

    logger.debug(f"Created IterateJob with {nprocs} process(es)")
    return job


@click.group(name="yaml", cls=MyGroup, invoke_without_command=True)
@click_yaml_common_options
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

    # Validate arguments based on separate_outputs flag
    source_output = ctx.get_parameter_source("outputfile")
    source_directory = ctx.get_parameter_source("directory")

    if separate_outputs:
        if source_output == click.core.ParameterSource.COMMANDLINE:
            raise click.UsageError(
                "Option '-o' / '--outputfile' is not allowed when '--separate-outputs' "
                "is enabled. Please use '-d' / '--directory' to specify the output location."
            )
        if directory is None:
            directory = os.getcwd()
    else:
        if source_directory == click.core.ParameterSource.COMMANDLINE:
            raise click.UsageError(
                "Option '-d' / '--directory' is not allowed when '--no-separate-outputs' "
                "(default) is active. Please use '-o' / '--outputfile' to specify the output file."
            )

    # Store shared parameters for the algorithm subcommands / shared executor.
    # Filename validation and YAML loading are deferred to the executor so
    # that 'yaml <algorithm> --help' works without requiring '-f'.
    ctx.obj["iterate"] = {
        "filename": filename,
        "timeout": timeout,
        "outputfile": outputfile,
        "directory": directory,
        "separate_outputs": separate_outputs,
        "combination_mode": combination_mode,
        "skip_completed": skip_completed,
    }

    # No algorithm subcommand: run with the YAML/default algorithm.
    if ctx.invoked_subcommand is None:
        return _build_iterate_job(ctx)


@yaml_cmd.command(name="jlgo", cls=MyCommand)
@click.option(
    "--adaptive-sampling/--no-adaptive-sampling",
    default=None,
    help="Run a fixed coarse sampling stage first; the six full-stage "
    "sampling/pruning options (--link-sphere-samples, "
    "--orientation-sphere-samples, --axial-samples, --candidate-pool-size, "
    "--preselect, --beam-width) only take effect when the coarse stage does "
    "not produce an acceptable optimized structure. Enabled by default. "
    "--max-starts and "
    "--slsqp-maxiter always apply. Use --no-adaptive-sampling to always "
    "apply the full sampling parameters.",
)
@click.option(
    "--link-sphere-samples",
    default=None,
    type=int,
    help="Full-stage number of linking-atom bond-sphere position samples "
    "(default: 48).",
)
@click.option(
    "--orientation-sphere-samples",
    default=None,
    type=int,
    help="Full-stage number of substituent principal-axis direction samples "
    "(default: 24).",
)
@click.option(
    "--axial-samples",
    default=None,
    type=int,
    help="Number of axial rotations per orientation direction (default: 4).",
)
@click.option(
    "--candidate-pool-size",
    default=None,
    type=int,
    help="Per-substituent candidate pool size kept after region exclusion "
    "(default: 20).",
)
@click.option(
    "--preselect",
    default=None,
    type=int,
    help="Top joint combinations fed into greedy start selection "
    "(default: 48).",
)
@click.option(
    "--beam-width",
    default=None,
    type=int,
    help="Beam width retained per layer during feasible-domain pruning "
    "(default: 4096).",
)
@click.option(
    "--max-starts",
    default=None,
    type=int,
    help="Maximum number of 6K-dimensional joint starts handed to SLSQP "
    "(default: 8).",
)
@click.option(
    "--slsqp-maxiter",
    default=None,
    type=int,
    help="Maximum SLSQP iterations per start (default: 200).",
)
@click.pass_context
def jlgo(
    ctx,
    adaptive_sampling,
    link_sphere_samples,
    orientation_sphere_samples,
    axial_samples,
    candidate_pool_size,
    preselect,
    beam_width,
    max_starts,
    slsqp_maxiter,
):
    """
    Run Joint Lagrange Geometry Optimization (JLGO).

    Attaches one or more substituents to the skeleton in a single joint
    (6K-dimensional) optimization. Options passed here override the matching
    values in the YAML ``algorithm`` block; options left unset keep their YAML
    value.

    With adaptive sampling enabled (the default), a fixed coarse stage runs
    first and the six full-stage sampling/pruning options
    (--link-sphere-samples, --orientation-sphere-samples, --axial-samples,
    --candidate-pool-size, --preselect, --beam-width) only take effect when
    the coarse stage does not produce an acceptable optimized structure;
    --max-starts and --slsqp-maxiter always apply. Use
    --no-adaptive-sampling to always apply the full sampling parameters.

    Examples:

    \b
    chemsmart run iterate yaml -f config.yaml jlgo
    chemsmart run iterate yaml -f config.yaml jlgo \\
        --no-adaptive-sampling \\
        --max-starts 16 \\
        --slsqp-maxiter 300
    """
    cli_options = _collect_cli_options(
        {
            "use_adaptive_sampling": adaptive_sampling,
            "n_link_sphere": link_sphere_samples,
            "n_orientation_sphere": orientation_sphere_samples,
            "n_axial": axial_samples,
            "candidate_pool_size": candidate_pool_size,
            "preselect": preselect,
            "beam_width": beam_width,
            "max_starts": max_starts,
            "slsqp_maxiter": slsqp_maxiter,
        },
    )
    return _build_iterate_job(
        ctx,
        cli_algorithm_name="jlgo",
        cli_options=cli_options,
    )


@yaml_cmd.command(name="etkdg", cls=MyCommand)
@click.option(
    "--global/--local",
    default=None,
    help="Embedding mode. 'local' (default) keeps the skeleton fixed and "
    "only re-embeds the substituent; 'global' re-embeds every atom.",
)
@click.option(
    "--num-conformers",
    default=None,
    type=click.IntRange(min=1),
    help="Number of ETKDG conformers to try per attachment; the "
    "lowest-energy one is kept (default: 10).",
)
@click.option(
    "--random-seed",
    default=None,
    type=int,
    help="Base RDKit random seed (-1 for a non-reproducible random seed; "
    "default: 42).",
)
@click.option(
    "--max-iterations",
    default=None,
    type=click.IntRange(min=0),
    help="Maximum ETKDG embedding iterations (0 uses the RDKit default; "
    "default: 2000).",
)
@click.option(
    "--random-coords/--no-random-coords",
    default=None,
    help="Start embedding from random coordinates (usually more robust; "
    "enabled by default).",
)
@click.option(
    "--enforce-chirality/--no-enforce-chirality",
    default=None,
    help="Enforce the input chirality during embedding (disabled by default).",
)
@click.option(
    "--force-field",
    default=None,
    type=click.Choice(
        ["none", "uff", "mmff94", "mmff94s"],
        case_sensitive=False,
    ),
    help="Optional force-field post-optimization after embedding "
    "(default: none).",
)
@click.pass_context
def etkdg(
    ctx,
    num_conformers,
    random_seed,
    max_iterations,
    random_coords,
    enforce_chirality,
    force_field,
    **cli_values,
):
    """
    Optimize substituent positions with the RDKit ETKDGv3 algorithm.

    By default this runs in ``local`` mode: the skeleton atoms are held
    fixed and only the substituent is re-embedded. Pass ``--global`` to
    re-embed the whole molecule. All options override the matching values
    in the YAML ``algorithm`` block; options left unset keep their YAML
    value.

    Examples:

    \b
    chemsmart run iterate yaml -f config.yaml etkdg
    chemsmart run iterate yaml -f config.yaml etkdg --global
    chemsmart run iterate yaml -f config.yaml etkdg \\
        --num-conformers 50 --random-seed 1
    """
    # Click derives the parameter name ``global`` directly from ``--global``.
    # Because ``global`` is a Python keyword, it must be received through
    # **cli_values instead of appearing as a formal function parameter.
    cli_options = _collect_cli_options(
        {
            "use_global_optimization": cli_values["global"],
            "num_conformers": num_conformers,
            "random_seed": random_seed,
            "max_iterations": max_iterations,
            "use_random_coordinates": random_coords,
            "enforce_chirality": enforce_chirality,
            "force_field": force_field,
        },
    )
    return _build_iterate_job(
        ctx,
        cli_algorithm_name="etkdg",
        cli_options=cli_options,
    )
