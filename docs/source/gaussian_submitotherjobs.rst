Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to contact our team if you have questions or feedback.

Submit Other Gaussian Jobs
===========================

ChemSmart provides additional Gaussian job capabilities including multi-step link jobs, custom user-defined calculations, and direct execution of Gaussian input files.

Link Jobs
---------

Run multi-step Gaussian calculations with linked job steps.

Link Job Specific Options
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

* **Basic link job**:

    .. code-block:: console

        chemsmart sub gaussian -p link_calc -f molecule.xyz link

* **Link job with specific job type**:

    .. code-block:: console

        chemsmart sub gaussian -p link_opt -f molecule.xyz link -j opt

* **Link job with stability test**:

    .. code-block:: console

        chemsmart sub gaussian -p link_stable -f molecule.xyz link -st opt

* **Link job with custom guess and route**:

    .. code-block:: console

        chemsmart sub gaussian -p link_custom -f molecule.xyz link -g core -r "scf=tight"

**Example applications**: Kc8 project and DNA project workflows

Custom User Jobs
----------------

Run user-defined Gaussian calculations with custom route specifications.

Custom Job Specific Options
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

* **Custom job with user-defined route**:

    .. code-block:: console

        chemsmart sub gaussian -p custom_calc -f molecule.xyz userjob -r "B3LYP/6-31G(d) opt freq"

* **Custom job with additional information**:

    .. code-block:: console

        chemsmart sub gaussian -p custom_detailed -f molecule.xyz userjob -r "MP2/cc-pVDZ" -a "Additional basis set information"

* **Advanced custom route example**:

    .. code-block:: console

        chemsmart sub gaussian -p advanced_custom -f complex.xyz userjob -r "wB97XD/6-311++G(d,p) scrf=(smd,solvent=water) opt"

Direct Input File Execution
----------------------------

Run Gaussian input files (.com or .gjf) directly without modifications.

Basic Usage
^^^^^^^^^^^

* **Run .com file directly**:

    .. code-block:: console

        chemsmart sub gaussian -f input_file.com com

* **Run .gjf file directly**:

    .. code-block:: console

        chemsmart sub gaussian -f calculation.gjf com

**Note**: This method uses the input file exactly as provided, without any modifications from ChemSmart.
