# ruff: noqa: E402
"""
CLI interface for chemsmart project.

This module provides the main entry point for the chemsmart command-line
interface, organizing various subcommands and providing the ASCII art
banner display.
"""

from __future__ import annotations

import importlib
import logging
from collections.abc import Callable

import click

logging.getLogger("numexpr.utils").setLevel(logging.WARNING)

from chemsmart import __version__
from chemsmart.utils.cli import MyGroup


class DeferredGroup(click.Group):
    """Lightweight command proxy that imports the real group on demand."""

    def __init__(
        self,
        *,
        loader: Callable[[], click.Group],
        name: str,
        help_text: str,
        short_help: str,
        invoke_without_command: bool = False,
        help_only_params: list[click.Parameter] | None = None,
        help_only_subcommands: dict[str, str] | None = None,
        use_placeholder_for_help: bool = False,
    ) -> None:
        super().__init__(
            name=name,
            help=help_text,
            short_help=short_help,
            invoke_without_command=invoke_without_command,
        )
        self._loader = loader
        self._loaded: click.Group | None = None
        self._use_placeholder_for_help = use_placeholder_for_help
        for param in help_only_params or []:
            self.params.append(param)
        for command_name, command_help in (
            help_only_subcommands or {}
        ).items():
            self.add_command(
                click.Command(
                    name=command_name,
                    help=command_help,
                    short_help=command_help,
                )
            )

    def _load(self) -> click.Group:
        if self._loaded is None:
            self._loaded = self._loader()
        return self._loaded

    def _use_placeholder_context(self, args: list[str]) -> bool:
        return self._use_placeholder_for_help and args in (
            ["--help"],
            ["-h"],
        )

    def make_context(self, info_name, args, parent=None, **extra):
        if self._use_placeholder_context(list(args)):
            return super().make_context(
                info_name,
                args,
                parent=parent,
                **extra,
            )
        return self._load().make_context(
            info_name, args, parent=parent, **extra
        )

    def invoke(self, ctx):
        return self._load().invoke(ctx)

    def list_commands(self, ctx):
        if self._loaded is None and self._use_placeholder_for_help:
            return sorted(self.commands)
        return self._load().list_commands(ctx)

    def get_command(self, ctx, cmd_name):
        if self._loaded is None and self._use_placeholder_for_help:
            return self.commands.get(cmd_name)
        return self._load().get_command(ctx, cmd_name)


def _load_group(module_path: str, attr_name: str) -> click.Group:
    module = importlib.import_module(module_path)
    return getattr(module, attr_name)


def _missing_agent_group(import_error: str) -> click.Group:
    @click.group(name="agent", invoke_without_command=True)
    def agent():
        """AI-scientist agent commands (install with `pip install -e .[agent-tui]`)."""
        click.echo("agent support is not installed. Run:", err=True)
        click.echo('  pip install -e ".[agent-tui]"', err=True)
        click.echo(f"(import error: {import_error})", err=True)
        raise click.exceptions.Exit(1)

    return agent


def _load_agent_group() -> click.Group:
    try:
        return _load_group("chemsmart.agent.cli", "agent")
    except ImportError as exc:
        return _missing_agent_group(str(exc))


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


entry_point.add_command(
    DeferredGroup(
        loader=lambda: _load_group("chemsmart.cli.run", "run"),
        name="run",
        help_text="Main command for running chemsmart jobs.",
        short_help="Run chemsmart jobs.",
        invoke_without_command=False,
    )
)
entry_point.add_command(
    DeferredGroup(
        loader=lambda: _load_group("chemsmart.cli.sub", "sub"),
        name="sub",
        help_text="Submit chemsmart jobs to remote schedulers.",
        short_help="Submit chemsmart jobs.",
        invoke_without_command=False,
    )
)
entry_point.add_command(
    DeferredGroup(
        loader=lambda: _load_group("chemsmart.cli.config", "config"),
        name="config",
        help_text="Set up configuration files and environment variables.",
        short_help="Configure chemsmart.",
        invoke_without_command=True,
    )
)
entry_point.add_command(
    DeferredGroup(
        loader=lambda: _load_group("chemsmart.cli.update", "update"),
        name="update",
        help_text="Maintenance utilities for chemsmart sources.",
        short_help="Update chemsmart metadata.",
        invoke_without_command=False,
    )
)
entry_point.add_command(
    DeferredGroup(
        loader=_load_agent_group,
        name="agent",
        help_text="AI-scientist agent commands for chemsmart.",
        short_help="Run the chemsmart agent.",
        invoke_without_command=True,
        help_only_params=[
            click.Option(
                ["--plain"],
                is_flag=True,
                default=False,
                help=(
                    "Run the agent TUI in plain inline mode for "
                    "conservative terminals."
                ),
            ),
            click.Option(
                ["--verbose", "-v"],
                is_flag=True,
                default=False,
                help=(
                    "Keep console logging at DEBUG instead of writing it "
                    "only to the agent log file."
                ),
            ),
        ],
        help_only_subcommands={
            "doctor": "Validate api.env, AI_PROVIDER, and provider connectivity.",
            "run": "Plan and execute an agent workflow.",
            "ask": "Run a one-shot dry-run request and stream Rich output.",
            "resume": "Resume an agent session.",
            "sessions": "List saved agent sessions.",
            "tools": "List registered agent tools.",
            "wizard": "Probe a target server and render or write a wizard YAML config.",
            "wizard-verify": "Verify wizard/server transport wiring for an existing YAML.",
            "wizard-refresh": "Refresh or reuse the wizard node sidecar cache for a server.",
        },
        use_placeholder_for_help=True,
    )
)


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
