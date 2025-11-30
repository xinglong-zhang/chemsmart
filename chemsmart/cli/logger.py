import functools
import logging

import click

logger = logging.getLogger(__name__)


def logger_options(f):
    """Logging configuration options."""

    @click.option(
        "-d",
        "--debug/--no-debug",
        default=False,
        help="Turn on debug logging.",
    )
    @click.option(
        "--stream/--no-stream",
        default=True,
        help="Turn on logging to stdout.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options
