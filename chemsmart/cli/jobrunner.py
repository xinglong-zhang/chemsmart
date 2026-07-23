import functools
import logging

import click

logger = logging.getLogger(__name__)

_DEPRECATED_MAX_TASKS_WARNINGS = {
    "n": "-N is deprecated; use -M/--max-tasks.",
    "num_nodes": (
        "--num-nodes and --number-of-nodes are deprecated; use -M/--max-tasks. "
        "They do not set server.num_nodes."
    ),
}
_deprecated_max_tasks_warnings_emitted: set[str] = set()

PARALLEL_HELP_SUB = (
    "Max concurrent scheduler array tasks when submitting a batch (SLURM, "
    "PBS, or LSF), and expand nestable jobs (crest/QRC/dias/traj) into one "
    "array task per child. Default is serial: top-level batches use one "
    "array task at a time (%1), and nestable jobs submit as a single parent "
    "with nested serial children. Pass --run-in-parallel to enable "
    "concurrent array tasks / nestable expansion; cap concurrency with "
    "-M/--max-tasks."
)

MAX_TASKS_HELP_SUB = (
    "With chemsmart sub --run-in-parallel: max concurrent scheduler array "
    "tasks (-M for the %%M cap in SLURM ``--array=1-N%%M``; also PBS "
    "``-J 1-N%%M``, LSF ``-J name[1-N%%M]``). Ignored under the serial "
    "default. Each array task still uses one node (--nodes=1 on SLURM); "
    "-n sets cores per task."
)

MAX_TASKS_N_ALIAS_HELP = "Deprecated; use -M/--max-tasks."

NUM_NODES_MAX_TASKS_ALIAS_HELP = "Deprecated; use -M/--max-tasks (does not request multiple nodes per task)."

MAX_TASKS_HELP_RUN = (
    "Unused for chemsmart run array throttling (run is serial). "
    "Accepted for CLI parity with sub; does not set server node count."
)

PARALLEL_HELP_RUN = (
    "Batch children always run serially with full resources; "
    "use chemsmart sub for cluster concurrency. "
    "--no-run-in-parallel is the default behavior."
)

_MAX_TASKS_CLI_OPTIONS = (
    ("primary", "array_concurrency", None),
    ("n", "max_tasks_n_alias", "n"),
    ("num_nodes", "max_tasks_num_nodes_alias", "num_nodes"),
)


def warn_deprecated_max_tasks_alias(warning_key: str) -> None:
    """Emit a one-time warning for a deprecated max-tasks alias."""
    if warning_key in _deprecated_max_tasks_warnings_emitted:
        return
    message = _DEPRECATED_MAX_TASKS_WARNINGS.get(warning_key)
    if message is None:
        return
    _deprecated_max_tasks_warnings_emitted.add(warning_key)
    logger.warning(message)


def merge_max_tasks_options(
    ctx: click.Context,
    array_concurrency: int | None,
    max_tasks_n_alias: int | None,
    max_tasks_num_nodes_alias: int | None,
) -> int | None:
    """Resolve max-tasks CLI options, warning on deprecated aliases."""
    values_by_source = {
        "primary": array_concurrency,
        "n": max_tasks_n_alias,
        "num_nodes": max_tasks_num_nodes_alias,
    }
    cli_values: dict[str, int | None] = {}
    for source_key, param_name, warning_key in _MAX_TASKS_CLI_OPTIONS:
        if (
            ctx.get_parameter_source(param_name)
            != click.core.ParameterSource.COMMANDLINE
        ):
            continue
        if warning_key is not None:
            warn_deprecated_max_tasks_alias(warning_key)
        cli_values[source_key] = values_by_source[source_key]

    if not cli_values:
        return array_concurrency

    unique_values = set(cli_values.values())
    if len(unique_values) > 1:
        raise click.UsageError(
            "Conflicting max-task options were set. Use -M/--max-tasks "
            "(accepted alias: --max-concurrency)."
        )
    return next(iter(unique_values))


def click_jobrunner_options(f=None, *, entry_point=None):
    """Common job runner configuration options.

    Pass ``entry_point="run"`` or ``entry_point="sub"`` so
    ``--run-in-parallel/--no-run-in-parallel`` help matches that command.
    """
    if entry_point == "sub":
        parallel_help = PARALLEL_HELP_SUB
        max_tasks_help = MAX_TASKS_HELP_SUB
    elif entry_point == "run":
        parallel_help = PARALLEL_HELP_RUN
        max_tasks_help = MAX_TASKS_HELP_RUN
    else:
        parallel_help = f"{PARALLEL_HELP_SUB} {PARALLEL_HELP_RUN}"
        max_tasks_help = f"{MAX_TASKS_HELP_SUB} {MAX_TASKS_HELP_RUN}"

    def decorator(func):
        @click.option(
            "-s",
            "--server",
            type=str,
            default=None,
            help="Server name. If not specified, will try to automatically "
            "determine and use the current server.",
        )
        @click.option(
            "-n",
            "--num-cores",
            type=int,
            help="Number of cores for each job.",
        )
        @click.option(
            "-g",
            "--num-gpus",
            type=int,
            default=None,
            help="Number of GPUs per node. Defaults to number of GPUs on "
            "the specified server if None.",
        )
        @click.option(
            "-m", "--mem-gb", type=int, default=None, help="Memory in GB."
        )
        @click.option(
            "-M",
            "--max-tasks",
            "--max-concurrency",
            "array_concurrency",
            type=int,
            default=None,
            help=max_tasks_help,
        )
        @click.option(
            "-N",
            "max_tasks_n_alias",
            type=int,
            default=None,
            help=MAX_TASKS_N_ALIAS_HELP,
        )
        @click.option(
            "--num-nodes",
            "--number-of-nodes",
            "max_tasks_num_nodes_alias",
            type=int,
            default=None,
            help=NUM_NODES_MAX_TASKS_ALIAS_HELP,
        )
        @click.option(
            "--fake/--no-fake",
            default=False,
            type=bool,
            help="If true, fake job runners will be used.",
        )
        @click.option(
            "--scratch/--no-scratch",
            default=None,
            type=bool,
            help=(
                "Run in scratch mode or without a scratch folder. "
                "Omit both flags to use program SCRATCH from the server YAML "
                "when set, otherwise the job runner's SCRATCH class default; "
                "use --scratch or --no-scratch to override explicitly."
            ),
        )
        @click.option(
            "--delete-scratch/--no-delete-scratch",
            default=False,
            type=bool,
            help="If job was run in scratch, delete the scratch folder after the job "
            "is completed successfully.",
        )
        @click.option(
            "--run-in-parallel/--no-run-in-parallel",
            default=False,
            type=bool,
            help=parallel_help,
        )
        @functools.wraps(func)
        def wrapper_common_options(*args, **kwargs):
            ctx = click.get_current_context()
            kwargs["array_concurrency"] = merge_max_tasks_options(
                ctx,
                kwargs.pop("array_concurrency"),
                kwargs.pop("max_tasks_n_alias"),
                kwargs.pop("max_tasks_num_nodes_alias"),
            )
            return func(*args, **kwargs)

        return wrapper_common_options

    if f is not None:
        return decorator(f)
    return decorator
