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

Duplicate names, including collisions with built-in commands, are errors.
Broken optional plugins are warned about and skipped, while strict discovery
used by validation raises the loading error.

A minimal plugin command is:

.. code-block:: python

   import click

   @click.command()
   def status():
       click.echo("ready")

The plugin distribution must be installed in every environment, including
each compute node, that executes a generated ``chemsmart run ...`` command.
