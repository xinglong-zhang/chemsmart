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

from chemsmart.cli.job import click_filename_options
from chemsmart.jobs.iterate.job import IterateJob
from chemsmart.jobs.iterate.report import ERROR_CODE_DESCRIPTIONS
from chemsmart.jobs.iterate.runner import IterateJobRunner
from chemsmart.jobs.iterate.settings import (
    IterateJobSettings,
    resolve_algorithm_config,
)
from chemsmart.utils.iterate import (
    generate_yaml_template,
    validate_yaml_config,
)

logger = logging.getLogger(__name__)


class _ProgressReporter:
    """CLI-owned, display-only progress reporter for iterate runs.

    Instances are callable with ``(completed, total)`` and are passed to
    :meth:`IterateJob.run` as ``progress_callback``. This keeps all terminal
    presentation in the CLI layer; the runner never depends on Click.

    On the first call a one-line header is printed. On an interactive
    terminal an in-place ASCII progress bar is drawn and refreshed; on a
    non-interactive stream (e.g. output redirected to a file) the dynamic bar
    is suppressed to avoid carriage-return spam while the header is still
    shown.
    """

    _BAR_WIDTH = 30

    def __init__(self, nprocs: int, enabled: bool):
        self._nprocs = nprocs
        self._enabled = enabled
        self._started = False
        self._bar_active = False

    def __call__(self, completed: int, total: int) -> None:
        if not self._started:
            self._started = True
            word = "process" if self._nprocs == 1 else "processes"
            click.echo(
                f"Running {total} combinations with {self._nprocs} {word}"
            )
            if self._enabled:
                click.echo("")

        if not self._enabled:
            return

        span = max(total, 1)
        filled = min(self._BAR_WIDTH, int(self._BAR_WIDTH * completed / span))
        bar = "█" * filled + " " * (self._BAR_WIDTH - filled)
        click.echo(f"\rIterate  [{bar}]  {completed}/{total}", nl=False)
        self._bar_active = True
        if completed >= total:
            click.echo("")
            self._bar_active = False

    def close(self) -> None:
        """Terminate an in-progress bar line (e.g. on early exit)."""
        if self._bar_active:
            click.echo("")
            self._bar_active = False


def click_yaml_common_options(f):
    """
    Common (algorithm-agnostic) Click options for the ``yaml`` command group.

    Algorithm-specific parameters are intentionally *not* included here; they
    are defined on the individual algorithm subcommands instead.
    """

    @click.option(
        "-g",
        "--generate-template",
        "generate_template_path",
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
        "-np",
        "--nprocs",
        default=1,
        type=click.IntRange(min=1),
        show_default=True,
        help="Number of processes for parallel execution.",
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
        "combination_mode",
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
        default="iterate_out",
        type=str,
        show_default=True,
        help="Output filename (without .xyz extension) for generated structures. "
        "Use only with --no-separate-outputs.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def _collect_explicit_options(ctx, options: dict) -> dict:
    """
    Return only the options the user explicitly passed on the command line.

    Options left at their Click default (i.e. not typed by the user) are
    excluded so that they never override values coming from the YAML
    ``algorithm`` block. This relies on ``ctx.get_parameter_source``.
    """
    explicit = {}
    for name, value in options.items():
        source = ctx.get_parameter_source(name)
        if source == click.core.ParameterSource.COMMANDLINE:
            explicit[name] = value
    return explicit


def _execute_iterate_job(
    ctx, cli_algorithm_name=None, cli_options=None
) -> None:
    """
    Shared execution path for the ``yaml`` group and its algorithm subcommands.

    Resolves the effective algorithm configuration (default < YAML < CLI),
    builds the job settings, and runs the iterate job.

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
    nprocs = data["nprocs"]
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

    # The existing ``chemsmart run --debug`` logging option also selects the
    # verbose Iterate display: no progress bar and no worker/RDKit suppression.
    debug_mode = logger.isEnabledFor(logging.DEBUG)

    # Create job runner
    jobrunner = IterateJobRunner(show_worker_logs=debug_mode)

    # Create job
    job = IterateJob(
        settings=job_settings,
        jobrunner=jobrunner,
        nprocs=nprocs,
        timeout=timeout,
        outputfile=outputfile,
        separate_outputs=separate_outputs,
        output_directory=directory,
        command_line=" ".join(sys.argv),
    )

    logger.debug(f"Created IterateJob with {nprocs} process(es)")

    # Run the job
    logger.debug("Running iterate job to generate molecular structures.")

    progress = _ProgressReporter(
        nprocs=nprocs,
        enabled=bool(
            getattr(sys.stdout, "isatty", None) and sys.stdout.isatty()
        ),
    )
    try:
        summary = job.run(progress_callback=None if debug_mode else progress)
    except KeyboardInterrupt:
        # User interrupt (SIGINT): outputs are written only after all
        # combinations finish, so in-memory results may not have reached disk.
        progress.close()
        click.echo(
            "Interrupted by user. Existing output files were not deleted.\n"
            "Partial in-memory results may not have been written.",
            err=True,
        )
        ctx.exit(130)
    except Exception as e:
        progress.close()
        logger.error(f"Error running iterate job: {e}")
        raise click.ClickException(str(e))

    progress.close()

    # Combination-level statistics (shown only when combinations were run).
    if summary.total > 0:
        successful = summary.structures_written
        failed = summary.failed + summary.timed_out + summary.write_failed
        num_width = max(4, len(str(summary.total)))
        click.echo("")
        click.echo(f"{'Total combinations:':<24}{summary.total:>{num_width}}")
        click.echo(
            f"{'Successful combinations:':<24}{successful:>{num_width}}"
        )
        click.echo(f"{'Failed combinations:':<24}{failed:>{num_width}}")
        click.echo("")

    # Surface the run report location, or a clear write error.
    if summary.summary_path:
        click.echo(f"Report: {summary.summary_path}")
    elif summary.summary_write_error:
        click.echo(
            f"Error: run report could not be written: "
            f"{summary.summary_write_error}",
            err=True,
        )

    # Output location: a single line for the merged file, or the directory and
    # file count for separate per-structure files.
    if summary.structures_written > 0:
        if separate_outputs:
            out_dir = directory or os.getcwd()
            click.echo(
                f"Output directory: {out_dir} "
                f"({summary.structures_written} XYZ files)"
            )
        elif summary.output_paths:
            click.echo(f"Output: {summary.output_paths[0]}")

    # The runner decides the exit code from the full contract. Exit 0 only for
    # a completely clean run whose report is on disk; any error means exit 1.
    if summary.exit_code == 0:
        return

    # Non-zero: surface the Gaussian-style error codes (or the report-write
    # failure) without depending on any top-level status string.
    if summary.summary_write_error:
        message = (
            f"The run report could not be written: "
            f"{summary.summary_write_error}"
        )
    elif summary.error_codes:
        details = "; ".join(
            f"{code} ({ERROR_CODE_DESCRIPTIONS.get(code, 'error')})"
            for code in summary.error_codes
        )
        message = (
            f"Iterate error termination [{details}]. "
            f"See the run report for details."
        )
    else:
        message = "Iterate error termination. See the run report for details."
    raise click.ClickException(message)


@click.group(name="yaml", invoke_without_command=True)
@click_yaml_common_options
@click_filename_options
@click.pass_context
def yaml_cmd(
    ctx,
    filename,
    nprocs,
    timeout,
    outputfile,
    generate_template_path,
    directory,
    separate_outputs,
    combination_mode,
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
    if generate_template_path is not None:
        template_path = generate_yaml_template(
            generate_template_path, overwrite=False
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
        "nprocs": nprocs,
        "timeout": timeout,
        "outputfile": outputfile,
        "directory": directory,
        "separate_outputs": separate_outputs,
        "combination_mode": combination_mode,
    }

    # No algorithm subcommand: run with the YAML/default algorithm.
    if ctx.invoked_subcommand is None:
        _execute_iterate_job(ctx)


@yaml_cmd.command(name="jlgo")
@click.option(
    "--adaptive-sampling/--no-adaptive-sampling",
    "use_adaptive_sampling",
    default=True,
    show_default=True,
    help="Run a fixed coarse sampling stage first; the six full-stage "
    "sampling/pruning options (--link-sphere-samples, "
    "--orientation-sphere-samples, --axial-samples, --candidate-pool-size, "
    "--preselect, --beam-width) only take effect when the coarse stage does "
    "not produce an acceptable optimized structure. --max-starts and "
    "--slsqp-maxiter always apply. Use --no-adaptive-sampling to always "
    "apply the full sampling parameters.",
)
@click.option(
    "--link-sphere-samples",
    "n_link_sphere",
    default=48,
    type=int,
    show_default=True,
    help="Full-stage number of linking-atom bond-sphere position samples.",
)
@click.option(
    "--orientation-sphere-samples",
    "n_orientation_sphere",
    default=24,
    type=int,
    show_default=True,
    help="Full-stage number of substituent principal-axis direction samples.",
)
@click.option(
    "--axial-samples",
    "n_axial",
    default=4,
    type=int,
    show_default=True,
    help="Number of axial rotations per orientation direction.",
)
@click.option(
    "--candidate-pool-size",
    "candidate_pool_size",
    default=20,
    type=int,
    show_default=True,
    help="Per-substituent candidate pool size kept after region exclusion.",
)
@click.option(
    "--preselect",
    "preselect",
    default=48,
    type=int,
    show_default=True,
    help="Top joint combinations fed into greedy start selection.",
)
@click.option(
    "--beam-width",
    "beam_width",
    default=4096,
    type=int,
    show_default=True,
    help="Beam width retained per layer during feasible-domain pruning.",
)
@click.option(
    "--max-starts",
    "max_starts",
    default=8,
    type=int,
    show_default=True,
    help="Maximum number of 6K-dimensional joint starts handed to SLSQP.",
)
@click.option(
    "--slsqp-maxiter",
    "slsqp_maxiter",
    default=200,
    type=int,
    show_default=True,
    help="Maximum SLSQP iterations per start.",
)
@click.pass_context
def jlgo(
    ctx,
    use_adaptive_sampling,
    n_link_sphere,
    n_orientation_sphere,
    n_axial,
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
    cli_options = _collect_explicit_options(
        ctx,
        {
            "use_adaptive_sampling": use_adaptive_sampling,
            "n_link_sphere": n_link_sphere,
            "n_orientation_sphere": n_orientation_sphere,
            "n_axial": n_axial,
            "candidate_pool_size": candidate_pool_size,
            "preselect": preselect,
            "beam_width": beam_width,
            "max_starts": max_starts,
            "slsqp_maxiter": slsqp_maxiter,
        },
    )
    _execute_iterate_job(
        ctx,
        cli_algorithm_name="jlgo",
        cli_options=cli_options,
    )


@yaml_cmd.command(name="etkdg")
@click.option(
    "--global/--local",
    "use_global_optimization",
    default=False,
    show_default=True,
    help="Embedding mode. 'local' (default) keeps the skeleton fixed and "
    "only re-embeds the substituent; 'global' re-embeds every atom.",
)
@click.option(
    "--num-conformers",
    "num_conformers",
    default=10,
    type=click.IntRange(min=1),
    show_default=True,
    help="Number of ETKDG conformers to try per attachment; the "
    "lowest-energy one is kept.",
)
@click.option(
    "--random-seed",
    "random_seed",
    default=42,
    type=int,
    show_default=True,
    help="Base RDKit random seed (-1 for a non-reproducible random seed).",
)
@click.option(
    "--max-iterations",
    "max_iterations",
    default=2000,
    type=click.IntRange(min=0),
    show_default=True,
    help="Maximum ETKDG embedding iterations (0 uses the RDKit default).",
)
@click.option(
    "--random-coords/--no-random-coords",
    "use_random_coordinates",
    default=True,
    show_default=True,
    help="Start embedding from random coordinates (usually more robust).",
)
@click.option(
    "--enforce-chirality/--no-enforce-chirality",
    "enforce_chirality",
    default=False,
    show_default=True,
    help="Enforce the input chirality during embedding.",
)
@click.option(
    "--force-field",
    "force_field",
    default="none",
    type=click.Choice(
        ["none", "uff", "mmff94", "mmff94s"],
        case_sensitive=False,
    ),
    show_default=True,
    help="Optional force-field post-optimization after embedding.",
)
@click.pass_context
def etkdg(
    ctx,
    use_global_optimization,
    num_conformers,
    random_seed,
    max_iterations,
    use_random_coordinates,
    enforce_chirality,
    force_field,
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
    cli_options = _collect_explicit_options(
        ctx,
        {
            "use_global_optimization": use_global_optimization,
            "num_conformers": num_conformers,
            "random_seed": random_seed,
            "max_iterations": max_iterations,
            "use_random_coordinates": use_random_coordinates,
            "enforce_chirality": enforce_chirality,
            "force_field": force_field,
        },
    )
    _execute_iterate_job(
        ctx,
        cli_algorithm_name="etkdg",
        cli_options=cli_options,
    )
