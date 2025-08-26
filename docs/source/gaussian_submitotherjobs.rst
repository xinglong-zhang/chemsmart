Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to contact our team if you have questions or feedback.

Submit Other Gaussian Jobs
===========================

ChemSmart provides additional Gaussian job capabilities including multi-step link jobs, custom user-defined calculations, and direct execution of Gaussian input files.

Link Jobs
---------

Run multi-step Gaussian calculations with linked job steps.

.. code-block:: console

    chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] link [SUBCMD_OPTIONS]

Link Job Specific OPTIONS
^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Link Job Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-j, --jobtype <string>``
     - string
     - Gaussian job type. Options: ["opt", "ts", "modred", "scan", "sp"]
   * - ``-st, --stable <string>``
     - string
     - Gaussian stability test options (default="opt")
   * - ``-g, --guess <string>``
     - string
     - Gaussian guess options (default="mix")
   * - ``-r, --route <string>``
     - string
     - Route for link section (default=None)

Basic Usage
^^^^^^^^^^^

* **Link job with specific job type**:

    .. code-block:: console

        chemsmart sub gaussian -p link_opt -f molecule.xyz link -j opt


Examples
^^^^^^^^

Kc8 project opt    and DNA project sp


Custom User Jobs
----------------

Generally, if a user wants to run job that is currently not present in our package, one can run custom job

.. code-block:: console

    chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] userjob [SUBCMD_OPTIONS]

Custom Job Specific OPTIONS
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Custom Job Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-r, --route <string>``
     - string
     - User-defined route for Gaussian calculation (required)
   * - ``-a, --append-info <string>``
     - string
     - Information to be appended at the end of the file (default=None)

Basic Usage
^^^^^^^^^^^

**Custom job with user-defined route**:

*   to create an input file named ``user_defined_job.com`` with user-specified route ``mnr functional/basis solvent`` etc and ``B 1 2 F\nA 1 2 3 F`` at the end of the input file after the specification of coordinates, run

    .. code-block:: console

        chemsmart sub -s shared gaussian -p test -f test.com -l user_defined_job userjob -r 'mnr functional/basis solvent etc' -a 'B 1 2 F\nA 1 2 3 F'


Direct Input File Execution
----------------------------

If a user wants to run a job with pre-prepared Gaussian input file directly, one can run the job directly without modifications.

.. code-block:: console

        chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] com

Basic Usage
^^^^^^^^^^^

**Direct execution of Gaussian input file**:

    .. code-block:: console

        chemsmart sub -s share gaussian -p test -f input_file.com com

    or for input file with .gjf extension
    .. code-block:: console

        chemsmart sub -s share gaussian -p test -f input_file.gjf com
