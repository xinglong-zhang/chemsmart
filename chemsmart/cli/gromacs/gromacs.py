"""
GROMACS Command Line Interface

This module provides the main CLI interface for GROMACS workflows.
"""

import functools
import logging

import click

from chemsmart.cli.job import click_filename_options
from chemsmart.utils.cli import MyGroup

logger = logging.getLogger(__name__)


def click_gromacs_options(f):
    """
    Common click options decorator for GROMACS jobs.
    """

    @click.option(
        "--project", "-p", type=str, default=None, help="Project settings."
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


@click.group(cls=MyGroup)
@click_gromacs_options
@click_filename_options
@click.pass_context
def gromacs(ctx, project, filename):
    """
    Main CLI command group for running GROMACS jobs.
    """

    logger.debug(f"GROMACS project: {project}")
    logger.debug(f"GROMACS filename: {filename}")

    ctx.ensure_object(dict)
    ctx.obj["project"] = project
    ctx.obj["filename"] = filename
    ctx.obj["jobrunner"] = None


@gromacs.result_callback()
@click.pass_context
def gromacs_process_pipeline(ctx, *args, **kwargs):
    """
    Result callback function for processing GROMACS command pipeline.
    """
    kwargs.update({"subcommand": ctx.invoked_subcommand})
    ctx.obj[ctx.info_name] = kwargs
    logger.debug(
        f"Pipeline completed for subcommand: {ctx.invoked_subcommand}"
    )
    return args[0] if args else None
from chemsmart.cli.gromacs.em import em