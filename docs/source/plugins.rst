CLI plugins
===========

CHEMSMART discovers commands from standard Python entry points.  Plugins may
provide job commands in ``chemsmart.job_commands`` or local administration and
reporting commands in ``chemsmart.commands``:

.. code-block:: toml

   [project.entry-points."chemsmart.job_commands"]
   mlip = "chemsmart_mlip.cli.jobs:mlip"

   [project.entry-points."chemsmart.commands"]
   mlip = "chemsmart_mlip.cli.project:mlip"

Each entry point must expose a Click command or group, and its entry-point
name must match ``command.name``. Job commands are available under both
``chemsmart run`` and ``chemsmart sub``.
Job command groups must use ``chemsmart.utils.cli.MyGroup`` and leaf
commands must use ``chemsmart.utils.cli.MyCommand``, including nested
commands. These classes record arguments for submission scripts.

In normal mode, duplicate plugin names are all skipped, and plugins colliding
with built-in commands are skipped. Broken plugins are also skipped. Each
failure is logged; strict discovery or registration raises an error instead.
Job plugins that introduce a new program must load an ``Executable`` subclass
with a matching ``PROGRAM``, so the submitter can resolve its server
configuration. Multiple matches are an error.

A minimal top-level plugin command (``chemsmart.commands``) is:

.. code-block:: python

   import click

   @click.command()
   def status():
       click.echo("ready")

Job commands must also declare the runtime classes they use on the entry-point
command. For example, in a separate plugin package:

.. code-block:: python

   import click
   from chemsmart.utils.cli import MyGroup
   from .jobs import ExampleJob, ExampleJobRunner, ExampleExecutable

   @click.group(cls=MyGroup)
   def example():
       pass

   # Define the job subcommands here, using cls=MyCommand.
   example.chemsmart_classes = (
       ExampleJob, ExampleJobRunner, ExampleExecutable,
   )

``chemsmart_classes`` is a tuple of registerable ``RegistryMixin`` subclasses.
Omit it for commands that introduce no runtime classes. Attach the declaration
to the root entry-point command, including classes needed by nested commands.
This also works for shared classes already imported by another plugin.

During discovery, external classes do not automatically enter CHEMSMART's
registry. The host validates the command and the entire declaration before
activating declared classes, once per class object. A failed plugin activates
none of its declarations. Core ``chemsmart.*`` classes continue to register
normally, including lazy core imports.

The loader preserves Python's normal module cache; it does not unload shared
dependencies or restore a global registry snapshot. Plugin imports must avoid
other side effects such as writing files or directly editing the registry:
those actions cannot be undone by command validation. Classes imported directly
outside plugin discovery retain the normal automatic registration behavior.

The plugin distribution must be installed in every environment, including
each compute node, that executes a generated ``chemsmart run ...`` command.
