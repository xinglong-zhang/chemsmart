"""CLI interface for chemsmart project.

Be creative! do whatever you want!

- Install click or typer and create a CLI app
- Use builtin argparse
- Start a web application
- Import things from your .base module
"""

import click
import logging
from .run import run
from .sub import sub
from chemsmart.utils.cli import MyGroup


@click.group(cls=MyGroup)
@click.pass_context
@click.option("--verbose", is_flag=True, default=False)
def entry_point(ctx, verbose):
    if verbose:
        from chemsmart.utils.logger import create_logger

        create_logger(debug=True)
        logging.basicConfig(level=logging.INFO)
        logging.info("Verbose mode activated.")


entry_point.add_command(run)


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
