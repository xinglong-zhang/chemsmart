"""
CLI subcommand for running iterate jobs from CDXML files.

This subcommand will support reading skeleton/substituent definitions
directly from ChemDraw CDXML files in a future release.
"""

import logging

import click

from chemsmart.cli.job import click_filename_options

logger = logging.getLogger(__name__)


@click.command(name="cdxml")
@click_filename_options
@click.pass_context
def cdxml(ctx, **kwargs):
    """
    Run iterate jobs from a CDXML file (coming soon).

    This command will generate new molecular structures by reading
    skeleton and substituent definitions directly from ChemDraw CDXML files.

    This feature is under development.
    """
    raise click.ClickException(
        "The 'iterate cdxml' command is not yet implemented. "
        "Please define the skeletons and substituents in a YAML configuration "
        "file and use 'chemsmart run iterate yaml -f config.yaml'."
    )
