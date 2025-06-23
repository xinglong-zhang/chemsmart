"""CLI interface for chemsmart project."""

import click

from chemsmart.utils.cli import MyGroup

from .config import config
from .run import run
from .sub import sub
from .update import update


@click.group(cls=MyGroup)
@click.pass_context
@click.option("--verbose", is_flag=True, default=True)
def entry_point(ctx, verbose):
    if verbose:
        debug = True
        stream = True
    else:
        debug = False
        stream = False
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
        + " / ___| | | | ____|  \/  / ___||  \/  |  / \  |  _ \_   _|"
    )
    logger.info(
        "   "
        + " " * 25
        + "| |   | |_| |  _| | |\/| \___ \| |\/| | / _ \ | |_) || |  "
    )
    logger.info(
        "   "
        + " " * 25
        + "| |___|  _  | |___| |  | |___) | |  | |/ ___ \|  _ < | |  "
    )
    logger.info(
        "   "
        + " " * 25
        + " \____|_| |_|_____|_|  |_|____/|_|  |_/_/   \_\_| \_\|_|  \n"
    )


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
