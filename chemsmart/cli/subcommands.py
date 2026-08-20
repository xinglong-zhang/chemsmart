"""Program command families for ``chemsmart run`` and ``chemsmart sub``.

The eleven program CLIs are heavy (they pull program adapters and their
scientific dependencies), so the run/sub groups load them lazily: a
command name resolves to its module only when that program is actually
invoked or its help is rendered.
"""

import importlib

import click

from chemsmart.utils.cli import MyGroup

#: Program command name -> defining module. The attribute name equals the
#: command name in every module.
PROGRAM_COMMAND_MODULES = {
    "database": "chemsmart.cli.database",
    "gaussian": "chemsmart.cli.gaussian",
    "grouper": "chemsmart.cli.grouper",
    "iterate": "chemsmart.cli.iterate",
    "mol": "chemsmart.cli.mol",
    "nciplot": "chemsmart.cli.nciplot",
    "orca": "chemsmart.cli.orca",
    "pka": "chemsmart.cli.pka",
    "pyscf": "chemsmart.cli.pyscf",
    "thermochemistry": "chemsmart.cli.thermochemistry",
    "xtb": "chemsmart.cli.xtb",
}


def load_program_command(name: str) -> click.Command:
    """Import and return one program command family by name."""

    module = importlib.import_module(PROGRAM_COMMAND_MODULES[name])
    return getattr(module, name)


class LazyProgramGroup(click.Group):
    """A group whose program children import on first use.

    Tree walkers (the agent CLI schema reads ``.commands`` directly) must
    always see the complete surface, so reading the mapping materialises
    every program module; dispatching one command imports only that one.
    """

    def __init__(self, *args, **kwargs):
        self._program_commands: dict[str, click.Command] = {}
        self._materialized = False
        super().__init__(*args, **kwargs)

    @property
    def commands(self) -> dict[str, click.Command]:
        if not self._materialized:
            self._materialized = True
            for name in PROGRAM_COMMAND_MODULES:
                if name not in self._program_commands:
                    self._program_commands[name] = load_program_command(name)
        return self._program_commands

    @commands.setter
    def commands(self, value) -> None:
        self._program_commands.update(value)

    def list_commands(self, ctx):
        return sorted(
            set(self._program_commands) | set(PROGRAM_COMMAND_MODULES)
        )

    def get_command(self, ctx, cmd_name):
        command = self._program_commands.get(cmd_name)
        if command is not None:
            return command
        if cmd_name in PROGRAM_COMMAND_MODULES:
            command = load_program_command(cmd_name)
            self._program_commands[cmd_name] = command
            return command
        return None


class LazyProgramMyGroup(LazyProgramGroup, MyGroup):
    """Lazy program loading plus MyGroup's invocation recording."""
