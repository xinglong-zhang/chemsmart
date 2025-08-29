Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to contact our team if you have questions or feedback.

Chemsmart Command Line Interface Tutorial
==========================================

ChemSmart provides a comprehensive command-line interface for quantum chemistry calculations and molecular analysis. This tutorial covers the fundamental CLI structure, execution modes, and common options that apply to all ChemSmart jobs. Whether you're running calculations locally or submitting to HPC clusters, ChemSmart offers flexible tools for computational chemistry workflows.

Basic Command Structure
^^^^^^^^^^^^^^^^^^^^^^^^

Chemsmart offers two main execution modes:

| **Local execution**: Use ``chemsmart run`` to execute tasks on current terminal.
| **HPC submission**: Use ``chemsmart sub`` to submit jobs to high-performance computing clusters.


*   The basic command structure is:

    .. code-block:: console

        chemsmart run/sub [OPTIONS] <CMD> [CMD_OPTIONS] <SUBCMD> [SUBCMD_OPTIONS]


OPTIONS for All Chemsmart Jobs
^^^^^^^^^^^^^^^^

.. list-table:: Server and Resource Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-s, --server``
     - string
     - If not specified, will try to automatically determine and use the current server
   * - ``-n, --num-cores``
     - int
     - Number of cores for each job
   * - ``-g, --num-gpus``
     - int
     - Number of GPUs per node. Defaults to number of GPUs on specified server if None
   * - ``-m, --mem-gb``
     - int
     - Memory allocation in gigabytes

.. note::

    |``-s`` is followed by the one of the servers specified in ``~/.chemsmart/server/*.yaml`` files, without ``.yaml`` extension.
    |``-n``, ``-g``, and ``-m`` can be used to override the default server.yaml.

.. list-table:: Execution Control Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-S/-R, --skip-completed/--no-skip-completed``
     - bool
     - To run completed job again. Use ``-R`` to rerun completed job (default=True)
   * - ``--fake/--no-fake``
     - bool
     - Enable simulation mode. If true, fake jobrunners will be used (default=False)
   * - ``--scratch/--no-scratch``
     - bool
     - Run in temporary scratch space or working directory (default=True)

.. note::

    ``-R`` need to add to the end of the command to rerun completed job！


.. list-table:: Debugging and Logging Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-d, --debug/--no-debug``
     - bool
     - Turns on debug logging for troubleshooting job failures (default=False)
   * - ``--stream/--no-stream``
     - bool
     - Turns on logging to stdout (default=True)


''chemsmart sub''-Specific Options
^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-t, --time-hours``
     - float
     - Maximum job runtime in hours
   * - ``-q, --queue``
     - string
     - Specify the HPC queue name
   * - ``-v, --verbose/--no-verbose``
     - bool
     - Turns on logging to stream output and debug logging (default=False)
   * - ``--test/--no-test``
     - bool
     - If true, job will not be submitted; only run and submit scripts will be written (default=False)
   * - ``--print-command/--no-print-command``
     - bool
     - Print the command generated (default=False)

COMMAND for Jobs
^^^^^^^^^^^^^^^^

*   use ``gaussian`` to run/sub a gaussian job
*   use ``orca`` to run/sub an ORCA job
*   use ``mol`` to run/sub a Pymol analysis job
*   use ``thermochemistry`` to run/sub a thermochemistry analysis job

In subsequent tutorials and examples, we will use "run" or "sub" based on common scenarios.

Next Steps
^^^^^^^^^^^^^^^^
For specific job types, see the detailed tutorials.

.. note::

    | Chemsmart checks the uniqueness of job names. If a job with the same name is already running, the *chemsmart sub* will be blocked.
    | Users can use the ``-a`` and ``-l`` options in the following tutorials to modify the job name and label, ensuring each one is unique.