"""CLI interface for ChemSmart chatbot."""

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


@click.group(name="chat", invoke_without_command=True)
@click.option(
    "--auto-execute/--no-auto-execute",
    default=False,
    help="Automatically execute generated commands",
)
@click.pass_context
def chat(ctx, auto_execute):
    """Start the ChemSmart chatbot for interactive command generation.

    The chatbot helps you create chemistry simulation job commands
    using natural language. Describe what you want to do, and it
    will generate the appropriate chemsmart command.

    Examples:
        chemsmart chat
        chemsmart chat --auto-execute
    """
    ctx.ensure_object(dict)
    chatbot = Chatbot(auto_execute=auto_execute)
    ctx.obj["chatbot"] = chatbot

    if ctx.invoked_subcommand is None:
        chatbot.chats()


@chat.command()
@click.argument("description", nargs=-1, required=True)
@click.option(
    "--execute/--no-execute",
    default=False,
    help="Execute the generated command",
)
@click.pass_context
def generate(ctx, description, execute):
    """Generate a command from a natural language description.

    Examples:
        chemsmart chat generate "submit optimization job for molecule.xyz"
        chemsmart chat generate --execute "run single point on test.com"
    """
    chatbot = ctx.obj.get("chatbot")
    if not chatbot:
        chatbot = Chatbot(auto_execute=execute)

    description_text = " ".join(description)
    command = chatbot.generate_command(description_text)

    if command:
        click.echo(f"Generated command: {command}")
        if execute:
            click.echo("Executing...")
            output = chatbot.executor.execute(command)
            click.echo(f"Output:\n{output}")
    else:
        click.echo("Could not generate a command from the description.")
        click.echo("Please try a more specific description.")


@chat.command()
@click.pass_context
def train(ctx):
    """Train the chatbot model (requires additional dependencies).

    This command fine-tunes the LLM on chemsmart commands.
    Requires: transformers, torch, datasets

    Examples:
        chemsmart chat train
    """
    chatbot = ctx.obj.get("chatbot")
    if not chatbot:
        chatbot = Chatbot()
    chatbot.train()


if __name__ == "__main__":
    try:
        chat()
    except KeyboardInterrupt:
        pass
