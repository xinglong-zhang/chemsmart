"""Reusable Click options for Iterate input commands."""

import functools
import os

import click


def click_iterate_execution_options(f):
    """Common execution options for Iterate input commands."""

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
        "'global': all slots combined via single Cartesian product.",
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
        "Defaults to <configuration_stem>_iterate.xyz for YAML input and "
        "iterate_iterate.xyz for direct input. Use only with --no-separate-outputs.",
    )
    @functools.wraps(f)
    def wrapper_iterate_execution_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_iterate_execution_options


def click_iterate_combination_options(f):
    """Common combination-generation options for Iterate inputs."""

    @click.option(
        "-ms",
        "--max-substituted-sites",
        default=None,
        type=click.IntRange(min=0),
        help="Maximum number of substituted skeleton attachment sites per "
        "combination. Unset or 0 means unlimited.",
    )
    @functools.wraps(f)
    def wrapper_iterate_combination_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_iterate_combination_options


def validate_iterate_output_options(
    ctx,
    *,
    outputfile: str | None,
    directory: str | None,
    separate_outputs: bool,
) -> str | None:
    """Validate mutually exclusive Iterate output options."""
    source_output = ctx.get_parameter_source("outputfile")
    source_directory = ctx.get_parameter_source("directory")

    if separate_outputs:
        if source_output == click.core.ParameterSource.COMMANDLINE:
            raise click.UsageError(
                "Option '-o' / '--outputfile' is not allowed when '--separate-outputs' "
                "is enabled. Please use '-d' / '--directory' to specify the output location."
            )
        if directory is None:
            return os.getcwd()
        return directory

    if source_directory == click.core.ParameterSource.COMMANDLINE:
        raise click.UsageError(
            "Option '-d' / '--directory' is not allowed when '--no-separate-outputs' "
            "(default) is active. Please use '-o' / '--outputfile' to specify the output file."
        )
    return directory
