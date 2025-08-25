Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to contact our team if you have questions or feedback.

General CLI Options for All Mol Jobs
====================================

ChemSmart provides comprehensive command-line options for molecular visualization and analysis using PyMOL. Use ``chemsmart run mol --help`` for detailed information about all available options.

Basic Command Structure
^^^^^^^^^^^^^^^^^^^^^^^

The basic command structure for Mol jobs is:

.. code-block:: console

    chemsmart run [GENERAL_OPTIONS] mol [MOL_OPTIONS] <SUBCOMMAND> [SUBCOMMAND_OPTIONS]

MOL_OPTIONS
^^^^^^^^^^^
Works for all Mol jobs

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-f, --filename <string>``
     - string
     - Filename from which new Gaussian input is prepared (required)
   * - ``-l, --label <string>``
     - string
     - Write user input filename for the job (without extension). Will overwrite your original filename (default=None)
   * - ``-a, --append-label <string>``
     - string
     - Name to be appended to file for the job. Will append name to current filename (default=None)
   * - ``-i, --index <string>``
     - string
     - Index of molecules to use; 1-based indices. Default to the last molecule structure (default="-1")
   * - ``--pubchem <string>``
     - string
     - Get molecule structure from PubChem database using identifier (default=None)
   * - ``-o, --overwrite``
     - bool
     - Overwrite existing files (default=False)

SUBCOMMANDS for Different Mol Jobs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Basic Visualization
   :header-rows: 1
   :widths: 15 85

   * - Subcommand
     - Description
   * - ``visualize``
     - CLI for running automatic PyMOL visualization and saving as pse file
   * - ``movie``
     - CLI for generating automatic PyMOL movie for rotating molecules

.. list-table:: Reaction Analysis
   :header-rows: 1
   :widths: 15 85

   * - Subcommand
     - Description
   * - ``irc``
     - CLI for generating automatic PyMOL IRC movie and saving the trajectory

.. list-table:: Electronic Structure Analysis
   :header-rows: 1
   :widths: 15 85

   * - Subcommand
     - Description
   * - ``mo``
     - CLI for generating molecular orbitals (MOs) and saving as pse file
   * - ``spin``
     - CLI for generating spin density and saving as pse file

.. list-table:: Interaction Analysis
   :header-rows: 1
   :widths: 15 85

   * - Subcommand
     - Description
   * - ``nci``
     - CLI for generating automatic PyMOL NCI plot and saving as pse file


Next Steps
^^^^^^^^^^^^^^^^^^^^^^^

For specific visualization types, see the detailed tutorials:

* **Basic Visualization**: Standard molecular structure visualization and movie generation
* **Electronic Structure Analysis**: Molecular orbital and spin density visualization
* **Reaction Analysis**: IRC trajectory and transition state visualization
* **Interaction Analysis**: Non-covalent interaction (NCI) analysis and visualization
