#############
 CLI plugins
#############

CHEMSMART discovers commands from standard Python entry points. Plugins may provide job commands in
``chemsmart.job_commands`` or local administration and reporting commands in ``chemsmart.commands``:

.. code:: toml

   [project.entry-points."chemsmart.job_commands"]
   mlip = "chemsmart_mlip.cli.jobs:mlip"

   [project.entry-points."chemsmart.commands"]
   mlip = "chemsmart_mlip.cli.project:mlip"

Each entry point must expose a Click command or group, and its entry-point name must match ``command.name``. Job
commands are available under both ``chemsmart run`` and ``chemsmart sub``. Job plugins use ordinary Click commands and
groups. The host records their invocations and input parameter values for submission scripts, including nested or lazy
commands. Use ordinary Click command classes in plugins; the host owns argument recording. ``MyGroup`` and ``MyCommand``
are used by core commands and are not part of the plugin contract. Option spellings, positional arguments, repeated
options, and flags are reconstructed from the Click parameter definitions. Values are captured before type conversion
and callbacks, so replay does not apply a callback's transformation twice. Resolved defaults and environment values are
included; file inputs and paths with ``resolve_path=True`` are made absolute for execution on compute nodes.

In normal mode, duplicate plugin names are all skipped, and plugins colliding with built-in commands are skipped. Broken
plugins are also skipped. Each failure is logged; strict discovery or registration raises an error instead. Job plugins
that introduce a new program must load an ``Executable`` subclass with a matching ``PROGRAM``, so the submitter can
resolve its server configuration. Multiple matches are an error.

A minimal top-level plugin command (``chemsmart.commands``) is:

.. code:: python

   import click


   @click.command()
   def status():
       click.echo("ready")

Job commands must also declare the runtime classes they use on the entry-point command. For example, in a separate
plugin package:

.. code:: python

   import click
   from .jobs import ExampleJob, ExampleJobRunner, ExampleExecutable


   @click.group()
   def example():
       pass


   @example.command()
   def train():
       return ExampleJob()


   example.chemsmart_classes = (
       ExampleJob,
       ExampleJobRunner,
       ExampleExecutable,
   )

``chemsmart_classes`` is a tuple of registerable ``RegistryMixin`` subclasses. Omit it for commands that introduce no
runtime classes. Attach the declaration to the root entry-point command, including classes needed by nested commands.
This also works for shared classes already imported by another plugin. Job commands return a host ``Job`` (or a list of
jobs), and declare the job, runner and executable adapters they introduce. The host owns parameter recording,
scheduling, and script generation; the plugin owns its domain logic.

During discovery, external classes do not automatically enter CHEMSMART's registry. The host validates the command and
the entire declaration before activating declared classes, once per class object. A failed plugin activates none of its
declarations. Core ``chemsmart.*`` classes continue to register normally, including lazy core imports.

The loader preserves Python's normal module cache; it does not unload shared dependencies or restore a global registry
snapshot. Plugin imports must avoid other side effects such as writing files or directly editing the registry: those
actions cannot be undone by command validation. Classes imported directly outside plugin discovery retain the normal
automatic registration behavior.

The plugin distribution must be installed in every environment, including each compute node, that executes a generated
``chemsmart run ...`` command.
