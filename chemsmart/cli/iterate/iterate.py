"""
CLI group for iterate commands.

Provides subcommands for generating molecular structures by attaching
substituents to skeleton molecules. Input can be provided via YAML
configuration files or CDXML files (coming soon).
"""

import logging

import click

from chemsmart.utils.cli import MyGroup

from .cdxml import cdxml
from .yaml_cmd import yaml_cmd

logger = logging.getLogger(__name__)


@click.group(cls=MyGroup)
@click.pass_context
def iterate(ctx, **kwargs):
    """
    Generate molecular structures by attaching substituents to skeletons.

    Use a subcommand to specify the input format:

    \b
    chemsmart run iterate yaml  -f config.yaml   YAML configuration file
    chemsmart run iterate cdxml -f scheme.cdxml   CDXML file (coming soon)
    """
    ctx.ensure_object(dict)


iterate.add_command(yaml_cmd, name="yaml")
iterate.add_command(cdxml)
