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
         ``SCRATCH`` when set, otherwise the job-runner class default (see :doc:`configuration-server-settings`)

   -  -  ``--run-in-parallel/--no-run-in-parallel``

      -  bool

      -  With ``chemsmart sub``: allow concurrent scheduler array tasks and expand nestable jobs (QRC/DIAS/CREST/traj)
         into array tasks. Default is serial: top-level batches use one array task at a time (``%1``); nestable jobs
         submit as a single parent with nested serial children. Local ``run`` keeps children serial.

   -  -  ``-M, --max-tasks``
      -  int
      -  Max concurrent array tasks (``%M`` in SLURM/PBS/LSF). Ignored under serial mode. Deprecated aliases: ``-N``,
         ``--num-nodes``

.. note::

   Use ``-R`` at the end of the command to rerun a completed job.

.. note::

   **CLI (``chemsmart run`` / ``chemsmart sub``)**

   When both ``--scratch`` and ``--no-scratch`` are omitted, scratch mode is decided in ``JobRunner.from_job`` before
   the typed runner is built:

   #. Explicit ``--scratch`` or ``--no-scratch`` wins.
   #. Else program ``SCRATCH`` in server YAML (Gaussian, ORCA, NCIPLOT only).
   #. Else the job-runner class default (``True`` for Gaussian/ORCA/NCIPLOT; ``False`` for PyMOL, thermochemistry,
      etc.).

   **Programmatic API (direct constructor)**

   If you construct a runner yourself with ``scratch=None``, server YAML is **not** read—you get the class ``SCRATCH``
   default only. That can differ from the CLI path when YAML would override the class default.

   Scratch **path** (when mode is on) is resolved separately from program ``ENVARS``, then ``SERVER.SCRATCH_DIR``, then
   user settings. See :ref:`scratch-behavior` in :doc:`configuration-server-settings`.

.. note::

   ``--fake`` automatically selects the program-matched fake runner based on the command group:

   -  ``chemsmart run --fake gaussian ...`` / ``chemsmart sub --fake gaussian ...`` uses the Gaussian fake runner.
   -  ``chemsmart run --fake orca ...`` / ``chemsmart sub --fake orca ...`` uses the ORCA fake runner.

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

.. _cli-batch-array-submission:

*********************************
 Batch jobs and scheduler arrays
*********************************

A top-level ``BatchJob`` (for example a multi-molecule pKa table, or other multi-target CLI batches) is submitted as
**one** scheduler array with **one task per child**. Nestable jobs (QRC, DIAS, CREST, traj) are different: see
``--run-in-parallel`` above.

**Default (serial array)**

Without ``--run-in-parallel``, batch arrays use one concurrent task (``%1`` on SLURM/PBS/LSF). Children still share a
single array job, but only one runs at a time.

**Concurrent tasks**

Pass ``--run-in-parallel`` and optionally ``-M`` / ``--max-tasks`` to cap concurrency (``%M`` in SLURM
``--array=1-N%M``, PBS ``#PBS -J 1-N%M``, or LSF ``#BSUB -J "name[1-N%M]"``). Each array task still uses one node;
``-M`` is **not** a node count. Without ``-M``, all tasks may run at once unless ``CHEMSMART_MAX_SUBMITTERS`` (or the
server / jobrunner max-submitters setting) limits concurrency.

.. code:: bash

   chemsmart sub --run-in-parallel -M 4 gaussian -p my_project -f batch_input.csv ...

**Dry-run scripts**

.. code:: bash

   chemsmart sub --test gaussian -p my_project -f batch_input.csv ...

This writes submit and run scripts without queueing. Add ``--print-command`` to print the reconstructed ``chemsmart
run`` arguments.

Batch array submission requires **SLURM**, **PBS/Torque**, or **LSF**. See :doc:`configuration-server-settings` for
``SCHEDULER`` and related server options.

********************
 Available Commands
********************

-  ``gaussian``: Run or submit Gaussian jobs
-  ``orca``: Run or submit ORCA jobs
-  ``mol``: Run PyMOL visualization and analysis jobs
-  ``thermochemistry``: Run thermochemistry analysis jobs
-  ``grouper``: Run structure grouping jobs

************
 Next Steps
************

For specific job types, see the detailed tutorials:

-  :doc:`gaussian-cli-options`
-  :doc:`orca-cli-options`
-  :doc:`pymol-cli-options`
-  :doc:`thermochemistry-analysis`
-  :doc:`grouper-cli-options`

.. note::

   CHEMSMART checks job name uniqueness. If a job with the same name is already running, submission will be blocked. Use
   ``-a`` (append label) or ``-l`` (label) options to create unique job names.
