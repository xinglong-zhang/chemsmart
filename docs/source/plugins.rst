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
with a matching ``PROGRAM`` when their command is imported, so the submitter
can resolve its server configuration. Multiple matches are an error.

A minimal plugin command is:

.. code-block:: python

   import click

   @click.command()
   def status():
       click.echo("ready")

The plugin distribution must be installed in every environment, including
each compute node, that executes a generated ``chemsmart run ...`` command.
