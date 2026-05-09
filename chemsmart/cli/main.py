# ruff: noqa: E402
"""
CLI interface for chemsmart project.

This module provides the main entry point for the chemsmart command-line
interface, organizing various subcommands and providing the ASCII art
banner display.
"""

import logging

import click

logging.getLogger("numexpr.utils").setLevel(logging.WARNING)

from chemsmart import __version__
from chemsmart.utils.cli import MyGroup

from .config import config
from .run import run
from .sub import sub
from .update import update

try:
    from chemsmart.agent.cli import agent
except ImportError as exc:
    _AGENT_IMPORT_ERROR = str(exc)

    @click.group(name="agent", invoke_without_command=True)
    def agent():
        """AI-scientist agent commands (install with `pip install -e .[agent-tui]`)."""
        click.echo("agent support is not installed. Run:", err=True)
        click.echo('  pip install -e ".[agent-tui]"', err=True)
        click.echo(f"(import error: {_AGENT_IMPORT_ERROR})", err=True)
        raise click.exceptions.Exit(1)


@click.group(cls=MyGroup)
@click.pass_context
@click.version_option(version=__version__, prog_name="CHEMSMART")
@click.option("--verbose", is_flag=True, default=True)
def entry_point(ctx, verbose):
    """
    Main entry point for the chemsmart CLI.
    """
    if verbose:
        debug = True
        stream = True
    else:
        debug = False
        stream = False
    if ctx.invoked_subcommand == "agent":
        return

    # Set up logging
    from chemsmart.utils.logger import create_logger

    logger = create_logger(debug=debug, stream=stream)

    # ASCII Arts for CHEMSMART
    logger.info("\n")
    logger.info(
        "   "
        + " " * 25
        + "  ____ _   _ _____ __  __ ____  __  __    _    ____ _____ "
    )
    logger.info(
        "   "
        + " " * 25
        + r" / ___| | | | ____|  \/  / ___||  \/  |  / \  |  _ \_   _|"
    )
    logger.info(
        "   "
        + " " * 25
        + r"| |   | |_| |  _| | |\/| \___ \| |\/| | / _ \ | |_) || |  "
    )
    logger.info(
        "   "
        + " " * 25
        + r"| |___|  _  | |___| |  | |___) | |  | |/ ___ \|  _ < | |  "
    )
    logger.info(
        "   "
        + " " * 25
        + r" \____|_| |_|_____|_|  |_|____/|_|  |_/_/   \_\_| \_\|_|  "
        + "\n"
    )


entry_point.add_command(run)
entry_point.add_command(sub)
entry_point.add_command(config)
entry_point.add_command(update)
entry_point.add_command(agent)


def main():  # pragma: no cover
    """
    The main function executes on commands:
    ``python -m chemsmart`` and ``$ chemsmart``.

    This is your program's entry point.

    You can change this function to do whatever you want.
    Examples:
        * Run a test suite
        * Run a server
        * Do some other stuff
        * Run a command line application (Click, Typer, ArgParse)
        * List all available tasks
        * Run an application (Flask, FastAPI, Django, etc.)
    """
    obj = {}
    entry_point(obj=obj)


if __name__ == "__main__":
    main()
