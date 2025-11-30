#################################
 Command Line Interface Overview
#################################

Chemsmart provides a comprehensive command-line interface for quantum chemistry calculations and molecular analysis.
This guide covers the fundamental CLI structure, execution modes, and common options.

*************************
 Basic Command Structure
*************************

Chemsmart offers two main execution modes:

-  **Local execution**: Use ``chemsmart run`` to execute tasks on the current terminal.
-  **HPC submission**: Use ``chemsmart sub`` to submit jobs to high-performance computing clusters.

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
      -  bool
      -  Run in scratch space or working directory (default: scratch)

.. note::

   Use ``-R`` at the end of the command to rerun a completed job.

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

-  ``gaussian``: Run or submit Gaussian jobs
-  ``orca``: Run or submit ORCA jobs
-  ``mol``: Run PyMOL visualization and analysis jobs
-  ``thermochemistry``: Run thermochemistry analysis jobs

************
 Next Steps
************

For specific job types, see the detailed tutorials:

-  :doc:`gaussian-cli-options`
-  :doc:`orca-cli-options`
-  :doc:`pymol-cli-options`
-  :doc:`thermochemistry-analysis`

.. note::

   Chemsmart checks job name uniqueness. If a job with the same name is already running, submission will be blocked. Use
   ``-a`` (append label) or ``-l`` (label) options to create unique job names.
