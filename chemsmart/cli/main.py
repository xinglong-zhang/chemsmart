"""CLI interface for chemsmart project."""

import logging

import click

from chemsmart.utils.cli import MyGroup

from .config import config
from .run import run
from .sub import sub
from .update import update

logger = logging.getLogger(__name__)


@click.group(cls=MyGroup)
@click.pass_context
@click.option(
    "--verbose/--no-verbose",
    is_flag=True,
    type=bool,
    default=False,
    help="Turn on verbose logging.",
)
def entry_point(ctx, verbose):
    if verbose:
        from chemsmart.utils.logger import create_logger

        create_logger(debug=True, stream=True)
        logging.info("Verbose mode activated.")


entry_point.add_command(run)
entry_point.add_command(sub)
entry_point.add_command(config)
entry_point.add_command(update)


def main():  # pragma: no cover
    """
    The main function executes on commands:
    `python -m chemsmart` and `$ chemsmart `.

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
