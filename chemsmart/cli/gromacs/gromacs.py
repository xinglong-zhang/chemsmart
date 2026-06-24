"""
GROMACS Command Line Interface.

This module provides the main CLI group for GROMACS workflows.

The group-level responsibility is intentionally small:
1. collect project / input options;
2. load GromacsProjectSettings when a project YAML is provided;
3. store shared values in ctx.obj for leaf commands such as em.
"""

import functools
import logging
from pathlib import Path

import click

from chemsmart.cli.job import click_filename_options
from chemsmart.settings.gromacs import GromacsProjectSettings
from chemsmart.utils.cli import MyGroup

logger = logging.getLogger(__name__)


def click_gromacs_options(f):
    """
    Common click options decorator for GROMACS jobs.
    """

    @click.option(
        "--project",
        "-p",
        type=str,
        default=None,
        help=(
            "GROMACS project YAML file or project name. "
            "For now, if this value points to an existing file, "
            "it will be treated as a YAML file path."
        ),
    )
    @click.option(
        "--project-yaml",
        "project_yaml",
        type=click.Path(
            exists=True,
            dir_okay=False,
            resolve_path=True,
        ),
        default=None,
        help=(
            "Explicit GROMACS project YAML file. "
            "This has priority over --project."
        ),
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def _load_project_settings(project=None, project_yaml=None):
    """
    Load GROMACS project settings from CLI project options.

    project_yaml has priority because click already validates that it exists.
    project is kept for compatibility with ChemSmart-style -p/--project usage.
    """
    if project_yaml is not None:
        return GromacsProjectSettings.from_yaml(project_yaml)

    if project is not None:
        project_path = Path(project)

        if project_path.exists() and project_path.is_file():
            return GromacsProjectSettings.from_yaml(project_path)

        raise FileNotFoundError(
            "GROMACS --project currently expects a project YAML file path. "
            f"Received: {project}"
        )

    return None


@click.group(cls=MyGroup)
@click_gromacs_options
@click_filename_options
@click.pass_context
def gromacs(ctx, project, project_yaml, filename):
    """
    Main CLI command group for running GROMACS jobs.
    """

    logger.debug("GROMACS project: %s", project)
    logger.debug("GROMACS project YAML: %s", project_yaml)
    logger.debug("GROMACS filename: %s", filename)

    project_settings = _load_project_settings(
        project=project,
        project_yaml=project_yaml,
    )

    ctx.ensure_object(dict)
    ctx.obj["project"] = project
    ctx.obj["project_yaml"] = project_yaml
    ctx.obj["project_settings"] = project_settings
    ctx.obj["filename"] = filename


@gromacs.result_callback()
@click.pass_context
def gromacs_process_pipeline(ctx, *args, **kwargs):
    """
    Result callback function for processing GROMACS command pipeline.
    """
    kwargs.update({"subcommand": ctx.invoked_subcommand})
    ctx.obj[ctx.info_name] = kwargs

    logger.debug(
        "Pipeline completed for subcommand: %s",
        ctx.invoked_subcommand,
    )

    return args[0] if args else None
