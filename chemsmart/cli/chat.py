import contextlib
import logging
import platform
from multiprocessing import set_start_method

import click

from chemsmart.chatbot.chatbot import Chatbot

system_type = platform.system()

if system_type in ("Darwin", "Windows"):
    with contextlib.suppress(RuntimeError):
        set_start_method("fork")

logger = logging.getLogger(__name__)


@click.group(name="chat")
@click.pass_context
def chat(ctx):
    """Start the chatbot."""
    chatbot = Chatbot()
    ctx.ensure_object(
        dict
    )  # Initialize the Click context object if not already initialized
    ctx.obj["chatbot"] = chatbot

    if ctx.invoked_subcommand is None:
        chatbot.chats()


@chat.command()
@click.pass_context
def train(ctx):
    """Train the chatbot.

    Examples:
        chemsmart chat train
    """
    chatbot = ctx.obj["chatbot"]
    chatbot.train()


if __name__ == "__main__":
    try:
        chat()
    except KeyboardInterrupt as e:
        raise e
