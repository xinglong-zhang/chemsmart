Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to
contact our team if you have questions or feedback.

###########################################
 General CLI Options for All Gaussian Jobs
###########################################

ChemSmart provides comprehensive command-line options for Gaussian quantum chemistry calculations. Use ``chemsmart sub
gaussian --help`` for detailed information about all available options.

*************************
 Basic Command Structure
*************************

The basic command structure for Gaussian jobs is:

.. code:: console

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] <SUBCMD> [SUBCMD_OPTIONS]

******************
 GAUSSIAN_OPTIONS
******************

Works for all Gaussian jobs

.. list-table:: Project and File Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-p, --project``
      -  string
      -  Choose the project settings

   -  -  ``-f, --filename``
      -  string
      -  Filename from which new Gaussian input is prepared

   -  -  ``-l, --label``
      -  string
      -  Write user input filename for the job (without extension). Labels should avoid containing special characters
         such as apostrophe (``'``), comma (`,`), asterisks (``*``), etc.

   -  -  ``-a, --append-label``
      -  string
      -  Name to be appended to file for the job. Will append name to current filename

   -  -  ``-i, --index``
      -  string
      -  Index of molecules to use; 1-based indices. Default to the last molecule structure (-1)

   -  -  ``-t, --title``
      -  string
      -  Gaussian job title

   -  -  ``-P, --pubchem``
      -  string
      -  Queries structure from PubChem using name, smiles, cid and conformer information

.. note::

   |  ``-p`` is followed by the one of the projects specified in ``~/.chemsmart/project/*.yaml`` files, without
      ``.yaml`` extension.
   |  ``-f`` is followed by an input file the user wishes to run job on. Note that this input file can be any format,
      such as ``.xyz``, Gaussian ``.com``, ``.gjf`` or ``.log`` file or ORCA ``.inp`` or ``.out`` file.

**Specify Name of the File**

-  Users can specify the name of the file to be created for the job, without file extension, they want to run by using
   the option ``-l``, e.g.:

   .. code:: console

      chemsmart sub -s shared gaussian -p test -f test.com -l custom_job_name opt

   will create input file named ``custom_job_name.com`` instead of the default ``test_opt.com``.

**Append String to Input File Name**

-  Users can also simply append a string to the base name of the filename supplied, e.g.:

   .. code:: console

      chemsmart sub -s shared gaussian -p test -f test.com -a append_string ts

   will create input file named ``test_append_string.com`` instead of the default ``test_ts.com``.

**Select the Particular Structure in file**

-  If one has more than one structure in the supplied file for input preparation, one can select the particular
   structure to perform job on by using the ``-i/--index`` option, e.g.:

   .. code:: console

      chemsmart sub -s shared gaussian -p test -f small.db -i 5 -c 0 -m 1 opt

   will take the 5th structure (1-indexed, as in chemsmart) from ase database file, ``small.db``, to create the input
   file for geometry optimization.

.. Warning::

   1-indexed numbers are used, instead of 0-indexed numbers in Python language, since most visualization softwares for
   moleculare are 1-indexed.

**Get molecule from pubchem directly**

-  Users can also get the molecule structure directly from PubChem using name, smiles, cid and conformer information,
   e.g.:

   .. code:: console

      chemsmart sub -s shared gaussian -p test -P 962 -c 0 -m 1 -l water opt

   will create input file named ``water.com`` for optimization calculation of water.

.. list-table:: Molecular Properties Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-c, --charge``
      -  int
      -  Charge of the molecule

   -  -  ``-m, --multiplicity``
      -  int
      -  Multiplicity of the molecule

..
   warning::::

   If there is no charge or multiplicity information in input files, users must specify them via ``-c <charge> -m <multiplicity>``.

**Modify Charge and Multiplicity**

-  Users can also modify the charge and multiplicity from the CLI, e.g.:

   Modify the charge in ``test.com`` to charge of +1 in the newly created input file ``test_charge.com`` via:

   .. code:: console

      chemsmart sub -s shared gaussian -p test -f test.com -c 1 -a charge opt

   Modify the multiplicity in ``test.com`` to multiplicity of 3 in the newly created input file
   ``test_multiplicity.com`` via:

   .. code:: console

      chemsmart sub -s shared gaussian -p test -f test.com -m 3 -a multiplicity opt

   Modify the charge to -1 and multiplicity to 2 in the newly created input ``file test_charge_multiplicity.com`` via:

   .. code:: console

      chemsmart sub -s shared gaussian -p test -f test.com -c -1 -m 2 -l test_charge_multiplicity opt

.. tip::

   This can be useful when, e.g., using optimized structure of a neutral closed-shell (charge 0, multiplicity 1) system
   to run a charged radical ion (e.g., charge +1 and multiplicity 2 in radical cation).

.. list-table:: Method and Basis Set Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-x, --functional``
      -  string
      -  New functional to run

   -  -  ``-b, --basis``
      -  string
      -  New basis set to run

   -  -  ``-s, --semiempirical``
      -  string
      -  Semiempirical method to run.

**Modify the Functional and Basis**

-  Users can also modify the functional and basis from the CLI to differ from those in project settings, e.g.:

   Modify the functional to ``b3lyp`` in the newly created input file ``test_functional.com`` via:

   .. code:: console

      chemsmart sub -s shared gaussian -p test -f test.com -x b3lyp -a functional opt

   Modify the basis to ``6-31G*`` in the newly created input file ``test_basis.com`` via:

   .. code:: console

      chemsmart sub -s shared gaussian -p test -f test.com -b "6-31G*" -a basis opt

**Use semiempirical method to run ts**

-  Users can also use semiempirical method to run ts search, e.g.:

   Modify the method to ``pm6`` in the newly created input file ``test_pm6.com`` via:

   .. code:: console

      chemsmart sub -s shared gaussian -p test -f test.com -s pm6 -a pm6 ts

   will use pm6 instead of the functional and basis in ``-p test``.

.. list-table:: Route and Calculation Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-o, --additional-opt-options``
      -  string
      -  Additional opt options

   -  -  ``-r, --additional-route-parameters``
      -  string
      -  Additional route parameters

   -  -  ``-A, --append-additional-info``
      -  string
      -  Additional information to be appended at the end of the input file. E.g., scrf=read

   -  -  ``-C, --custom-solvent``
      -  string
      -  Additional information to be appended at the end of the input file. E.g., scrf=read

   -  -  ``-d, --dieze-tag``
      -  string
      -  Dieze tag for gaussian job. possible options include "n", "p", "t" to get "#n", "#p", "#t", respectively

   -  -  ``--forces/--no-forces``
      -  bool
      -  Whether to calculate forces (default=False)

**Specify Additional Optimization Options**

-  Users can also specify additional optimization options for opt=() in the route, for example,

   .. code:: console

      chemsmart sub -s shared gaussian -p test -f test.com -o maxstep=8,maxsize=12 -a opt_options opt

   will create ``opt=(maxstep=8,maxsize=12)`` as part of the route in the newly created input file
   ``test_opt_options.com``.

**Add in Additional Route Parameters**

-  Users can also add in additional parameters used in the route, e.g.:

   .. code:: console

      chemsmart sub -s shared gaussian -p test -f test.com --r nosymm -a route_params opt

   will add in ``nosymm`` as part of the route in the newly created input file ``test_route_params.com``.

************************************
 SUBCMD for Different Gaussian Jobs
************************************

.. list-table:: Conformational Sampling & Dynamics
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``crest``
      -  for running Gaussian CREST jobs
   -  -  ``traj``
      -  for running Gaussian trajectory jobs

.. list-table:: Structure Optimization
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``opt``
      -  for optimization calculation for Gaussian

.. list-table:: Transition State Search
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``ts``
      -  for transition state calculation for Gaussian
   -  -  ``modred``
      -  for running Gaussian modred jobs
   -  -  ``irc``
      -  for running Gaussian IRC jobs
   -  -  ``scan``
      -  for running Gaussian scan jobs

.. list-table:: Electronic Structure Properties & Analyses
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``sp``
      -  for single point calculation for Gaussian
   -  -  ``nci``
      -  for NCI for Gaussian
   -  -  ``dias``
      -  for Distortion-Interaction/Activation-Strain analysis
   -  -  ``resp``
      -  for RESP for Gaussian
   -  -  ``td``
      -  for time-dependent DFT calculation for Gaussian
   -  -  ``wbi``
      -  for WBI jobs

.. list-table:: Other Jobs
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``com``
      -  for running Gaussian input file as is
   -  -  ``link``
      -  for Gaussian link jobs
   -  -  ``userjob``
      -  for running Gaussian custom jobs

************
 Next Steps
************

For specific calculation types, see the detailed tutorials:

-  Submit Conformational Sampling & Dynamics Jobs
-  Submit Structure Optimization Jobs
-  Submit Transition State Search Jobs
-  Submit Electronic Structure Properties & Analyses Jobs
-  Submit Other Jobs
