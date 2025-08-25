Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to contact our team if you have questions or feedback.

General CLI Options for All Gaussian Jobs
==========================================

ChemSmart provides comprehensive command-line options for Gaussian quantum chemistry calculations. Use ``chemsmart sub gaussian --help`` for detailed information about all available options.

Basic Command Structure
^^^^^^^^^^^^^^^^

The basic command structure for Gaussian jobs is:

.. code-block:: console

    chemsmart run/sub [GENERAL_OPTIONS] gaussian [GAUSSIAN_OPTIONS] <SUBCOMMAND> [SUBCOMMAND_OPTIONS] [SOLVENT_OPTIONS]

GAUSSIAN_OPTIONS
^^^^^^^^^^^^^^^^
Works for all Gaussian jobs

.. list-table:: Project and File Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-p, --project``
     - string
     - Project settings (required)
   * - ``-f, --filename``
     - string
     - Filename from which new Gaussian input is prepared (required)
   * - ``-l, --label``
     - string
     - Write user input filename for the job (without extension). Will overwrite your original filename
   * - ``-a, --append-label``
     - string
     - Name to be appended to file for the job. Will append name to current filename
   * - ``-i, --index``
     - string
     - Index of molecules to use; 1-based indices. Default to the last molecule structure (-1)
   * - ``-P, --pubchem``
     - string
     - Queries structure from PubChem using name, smiles, cid and conformer information


.. list-table:: Molecular Properties Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-c, --charge``
     - int
     - Charge of the molecule
   * - ``-m, --multiplicity``
     - int
     - Multiplicity of the molecule
   * - ``-t, --title``
     - string
     - Gaussian job title

.. list-table:: Method and Basis Set Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-x, --functional``
     - string
     - New functional to run
   * - ``-b, --basis``
     - string
     - New basis set to run
   * - ``-s, --semiempirical``
     - string
     - Semiempirical method to run. Overwrite the -p yaml

.. list-table:: Route and Calculation Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-o, --additional-opt-options``
     - string
     - Additional opt options
   * - ``-r, --additional-route-parameters``
     - string
     - Additional route parameters
   * - ``-A, --append-additional-info``
     - string
     - Additional information to be appended at the end of the input file. E.g., scrf=read
   * - ``-C, --custom-solvent``
     - string
     - Additional information to be appended at the end of the input file. E.g., scrf=read
   * - ``-d, --dieze-tag``
     - string
     - Dieze tag for gaussian job; possible options include "n", "p", "t" to get "#n", "#p", "#t", respectively
   * - ``--forces/--no-forces``
     - bool
     - Whether to calculate forces (default=False)


SOLVENT-OPTIONS for Gaussian Jobs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Please set solvent parameters at the end when needed.

.. list-table:: Solvent Model Options
   :header-rows: 1
   :widths: 35 15 50

   * - Option
     - Type
     - Description
   * - ``--remove-solvent/--no-remove-solvent``
     - bool
     - Whether to use solvent model in the job. Defaults to project settings
   * - ``-sm, --solvent-model <string>``
     - string
     - Solvent model to be used for single point (default = none, use project settings)
   * - ``-si, --solvent-id <string>``
     - string
     - Solvent ID to be used for single point (default = none, use project settings)
   * - ``-so, --solvent-options <string>``
     - string
     - Additional solvent options in scrf=() route. default=None. E.g., iterative in scrf=(smd,water,iterative)

SUBCOMMAND for Different Gaussian Jobs
^^^^^^^^^^^^^^^^

.. list-table:: Conformational Sampling & Dynamics
   :header-rows: 1
   :widths: 15 85

   * - Subcommand
     - Description
   * - ``crest``
     - for running Gaussian CREST jobs
   * - ``traj``
     - for running Gaussian trajectory jobs

.. list-table:: Structure Optimization
   :header-rows: 1
   :widths: 15 85

   * - Subcommand
     - Description
   * - ``opt``
     - for optimization calculation for Gaussian
   * - ``crestopt``
     - for CREST-optimization calculation for Gaussian
   * - ``modred``
     - for running Gaussian modred jobs

.. list-table:: Transition State Search
   :header-rows: 1
   :widths: 15 85

   * - Subcommand
     - Description
   * - ``ts``
     - for transition state calculation for Gaussian
   * - ``irc``
     - for running Gaussian IRC jobs
   * - ``scan``
     - for running Gaussian scan jobs

.. list-table:: Electronic Structure Properties & Analyses
   :header-rows: 1
   :widths: 15 85

   * - Subcommand
     - Description
   * - ``sp``
     - for single point calculation for Gaussian
   * - ``nci``
     - for NCI for Gaussian
   * - ``dias``
     - for Distortion-Interaction/Activation-Strain analysis
   * - ``resp``
     - for RESP for Gaussian
   * - ``td``
     - for time-dependent DFT calculation for Gaussian
   * - ``wbi``
     - for WBI jobs

.. list-table:: Other Jobs
   :header-rows: 1
   :widths: 15 85

   * - Subcommand
     - Description
   * - ``com``
     - CLI for running Gaussian input file as is
   * - ``link``
     - CLI for Gaussian link jobs
   * - ``userjob``
     - CLI for running Gaussian custom jobs


Next Steps
^^^^^^^^^^^^^^^^

For specific calculation types, see the detailed tutorials:

*   Submit Conformational Sampling & Dynamics Jobs

*   Submit Structure Optimization Jobs

*   Submit Transition State Search Jobs

*   Submit Electronic Structure Properties & Analyses Jobs

*   Submit Other Jobs

