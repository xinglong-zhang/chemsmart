#################################
 Command Line Interface Overview
#################################

CHEMSMART provides a comprehensive command-line interface for quantum chemistry calculations and molecular analysis.
This guide covers the fundamental CLI structure, execution modes, and common options.

*************************
 Basic Command Structure
*************************

CHEMSMART offers two main execution modes:

-  **Local execution**: Use ``chemsmart run`` to execute tasks on the current terminal.
-  **HPC submission**: Use ``chemsmart sub`` to submit jobs to high-performance computing clusters.

The provider-neutral ``chemsmart agent`` command plans through these same
project loaders and Click commands. ``agent plan`` and ``agent run`` are
preview-only; ``agent tui`` adds the visible human review and one explicit
``/approve`` action for release-qualified real execution. See
:doc:`agent-workflows`.

The basic command structure is:

.. code:: bash

   chemsmart run/sub [OPTIONS] <CMD> [CMD_OPTIONS] <SUBCMD> [SUBCMD_OPTIONS]

*****************************
 Common Options for All Jobs
*****************************

Server and Resource Options
===========================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-s, --server``
      -  string
      -  Server name from ``~/.chemsmart/server/*.yaml`` (auto-detected if not specified)

   -  -  ``-n, --num-cores``
      -  int
      -  Number of cores per job

   -  -  ``-g, --num-gpus``
      -  int
      -  Number of GPUs per node (defaults to server configuration)

   -  -  ``-m, --mem-gb``
      -  int
      -  Memory allocation in gigabytes

.. note::

   The ``-s`` option takes the server name without the ``.yaml`` extension. The ``-n``, ``-g``, and ``-m`` options
   override the server defaults.

Execution Control Options
=========================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-S/-R, --skip-completed/--no-skip-completed``
      -  bool
      -  Skip or rerun completed jobs (default: skip)

   -  -  ``--fake/--no-fake``
      -  bool
      -  Enable simulation mode with fake job runners (default: disabled)

   -  -  ``--scratch/--no-scratch``
      -  bool or None
      -  Force scratch on (``--scratch``) or off (``--no-scratch``). Omit both (``None``) to use program YAML
         ``SCRATCH`` when set, otherwise the documented program default (see :doc:`configuration-server-settings`)

.. note::

   Use ``-R`` at the end of the command to rerun a completed job.

.. note::

   **CLI (``chemsmart run`` / ``chemsmart sub``)**

   When both ``--scratch`` and ``--no-scratch`` are omitted, scratch mode is decided before execution:

   #. Explicit ``--scratch`` or ``--no-scratch`` wins.
   #. Else program ``SCRATCH`` in server YAML when the selected runner has a registered executable/library configuration
      (including Gaussian, ORCA, NCIPLOT, PySCF, and xTB).
   #. Else the program default (enabled for Gaussian, ORCA, and NCIPLOT; disabled for operations that do not need
      scratch storage).

   Scratch **path** (when mode is on) is resolved separately from program ``ENVARS``, then ``SERVER.SCRATCH_DIR``, then
   user settings. See :ref:`scratch-behavior` in :doc:`configuration-server-settings`.

.. note::

   ``--fake`` automatically selects the program-matched fake runner based on the command group:

   -  ``chemsmart run --fake gaussian ...`` / ``chemsmart sub --fake gaussian ...`` uses the Gaussian fake runner.
   -  ``chemsmart run --fake orca ...`` / ``chemsmart sub --fake orca ...`` uses the ORCA fake runner.
   -  ``chemsmart run --fake pyscf ...`` / ``chemsmart sub --fake pyscf ...`` generates a PySCF preview artifact without
      importing or running PySCF.
   -  ``chemsmart run --fake xtb ...`` / ``chemsmart sub --fake xtb ...`` generates an xTB preview without invoking xTB.

   In these fake modes, executable-path checks for the corresponding real program are not required and the corresponding
   fake runner will be used without needing to specify its path.

Debugging and Logging Options
=============================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-d, --debug/--no-debug``
      -  bool
      -  Enable debug logging (default: disabled)

   -  -  ``--stream/--no-stream``
      -  bool
      -  Enable logging to stdout (default: enabled)

*****************************
 Submission-Specific Options
*****************************

These options are only available with ``chemsmart sub``:

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-t, --time-hours``
      -  float
      -  Maximum job runtime in hours

   -  -  ``-q, --queue``
      -  string
      -  HPC queue name

   -  -  ``-v, --verbose/--no-verbose``
      -  bool
      -  Enable verbose output and debug logging (default: disabled)

   -  -  ``--test/--no-test``
      -  bool
      -  Generate scripts without submitting (default: disabled)

   -  -  ``--print-command/--no-print-command``
      -  bool
      -  Print the generated command (default: disabled)

********************
 Available Commands
********************

Use ``chemsmart run --help`` and ``chemsmart sub --help`` as the live command
inventory. Both families register Gaussian, ORCA, PySCF, xTB, NCIPLOT,
molecular/PyMOL operations, pKa and thermochemistry analysis, molecular
databases, structure grouping, and ITERATE workflows. Program-specific pages
describe the available leaves and project settings.

Top-level ``chemsmart agent`` exposes public ``plan``, preview-only ``run``,
and interactive ``tui`` commands. In the TUI, the human reviews the displayed
YAML/CLI DAG and enters ``/approve`` once for a real release-qualified run. No
hash or approval-file token is part of that human interface.

************
 Next Steps
************

For specific job types, see the detailed tutorials:

-  :doc:`gaussian-cli-options`
-  :doc:`orca-cli-options`
-  :doc:`pyscf-cli-options`
-  :doc:`xtb-cli-options`
-  :doc:`pymol-cli-options`
-  :doc:`thermochemistry-analysis`
-  :doc:`grouper-cli-options`

.. note::

   CHEMSMART checks job name uniqueness. If a job with the same name is already running, submission will be blocked. Use
   ``-a`` (append label) or ``-l`` (label) options to create unique job names.
