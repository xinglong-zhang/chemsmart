######################
 Gaussian CLI Options
######################

This page documents the CLI options available for all Gaussian jobs. Use ``chemsmart sub gaussian --help`` for the
complete list.

*************************
 Basic Command Structure
*************************

.. code:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] <SUBCMD> [SUBCMD_OPTIONS]

******************
 Gaussian Options
******************

Project and File Options
========================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-p, --project``
      -  string
      -  Project settings from ``~/.chemsmart/gaussian/*.yaml``

   -  -  ``-f, --filename``
      -  string
      -  Input file for job preparation

   -  -  ``-l, --label``
      -  string
      -  Custom output filename (without extension)

   -  -  ``-a, --append-label``
      -  string
      -  String to append to the base filename

   -  -  ``-i, --index``
      -  string
      -  Structure index (1-based, default: last structure)

   -  -  ``-t, --title``
      -  string
      -  Gaussian job title

   -  -  ``-P, --pubchem``
      -  string
      -  Query structure from PubChem (name, SMILES, CID)

.. note::

   -  ``-p`` uses the project name without the ``.yaml`` extension.
   -  ``-f`` accepts various formats: ``.xyz``, ``.com``, ``.gjf``, ``.log``, ``.inp``, or ``.out``.

Specifying Output Filenames
---------------------------

Use ``-l`` to set a custom filename:

.. code:: bash

   chemsmart sub gaussian -p test -f test.com -l custom_name opt

This creates ``custom_name.com`` instead of ``test_opt.com``.

Use ``-a`` to append a string to the base filename:

.. code:: bash

   chemsmart sub gaussian -p test -f test.com -a suffix ts

This creates ``test_suffix.com`` instead of ``test_ts.com``.

Selecting Structures
--------------------

Use ``-i`` to select a specific structure from multi-structure files:

.. code:: bash

   chemsmart sub gaussian -p test -f molecules.db -i 5 -c 0 -m 1 opt

This uses the 5th structure (1-indexed) from the ASE database file.

.. warning::

   Chemsmart uses 1-based indexing to match most molecular visualization software.

Using PubChem
-------------

Fetch structures directly from PubChem:

.. code:: bash

   chemsmart sub gaussian -p test -P 962 -c 0 -m 1 -l water opt

This fetches water (CID 962) and creates ``water.com``.

Molecular Properties Options
============================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-c, --charge``
      -  int
      -  Molecular charge

   -  -  ``-m, --multiplicity``
      -  int
      -  Molecular multiplicity

.. warning::

   If the input file lacks charge/multiplicity information, you must specify them with ``-c`` and ``-m``.

Examples:

.. code:: bash

   # Modify charge to +1
   chemsmart sub gaussian -p test -f molecule.com -c 1 opt

   # Modify multiplicity to 3 (triplet)
   chemsmart sub gaussian -p test -f molecule.com -m 3 opt

   # Modify both charge (-1) and multiplicity (2)
   chemsmart sub gaussian -p test -f molecule.com -c -1 -m 2 opt

.. tip::

   This is useful when using an optimized neutral closed-shell structure to run calculations on a charged radical ion.

Method and Basis Set Options
============================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-x, --functional``
      -  string
      -  DFT functional

   -  -  ``-b, --basis``
      -  string
      -  Basis set

   -  -  ``-s, --semiempirical``
      -  string
      -  Semiempirical method (PM6, AM1, etc.)

Examples:

.. code:: bash

   # Use a different functional
   chemsmart sub gaussian -p test -f molecule.com -x b3lyp opt

   # Use a different basis set
   chemsmart sub gaussian -p test -f molecule.com -b "6-31G*" opt

   # Use semiempirical method for TS search
   chemsmart sub gaussian -p test -f molecule.com -s pm6 ts

Route and Calculation Options
=============================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-o, --additional-opt-options``
      -  string
      -  Additional optimization options

   -  -  ``-r, --additional-route-parameters``
      -  string
      -  Additional route parameters

   -  -  ``-A, --append-additional-info``
      -  string
      -  Information appended after coordinates

   -  -  ``-C, --custom-solvent``
      -  string
      -  Custom solvent parameters

   -  -  ``-d, --dieze-tag``
      -  string
      -  Route prefix: "n", "p", or "t" for #n, #p, #t

   -  -  ``--forces/--no-forces``
      -  bool
      -  Calculate forces (default: disabled)

Examples:

.. code:: bash

   # Add optimization options
   chemsmart sub gaussian -p test -f molecule.com -o maxstep=8,maxsize=12 opt

   # Add route parameters
   chemsmart sub gaussian -p test -f molecule.com -r nosymm opt

Solvent Options
===============

Solvent settings can be specified at the Gaussian group level, which means they apply to **any** subcommand (``opt``,
``td``, ``sp``, ``ts``, etc.). This is useful when the project settings define a gas-phase calculation but you want to
add solvation for a particular run without modifying the project file.

.. list-table::
   :header-rows: 1
   :widths: 35 15 50

   -  -  Option
      -  Type
      -  Description

   -  -  ``--remove-solvent/--no-remove-solvent``
      -  bool
      -  Remove solvent from the job, overriding project settings (default: disabled)

   -  -  ``-sm, --solvent-model``
      -  string
      -  Implicit solvent model (e.g. ``smd``, ``cpcm``, ``iefpcm``)

   -  -  ``-si, --solvent-id``
      -  string
      -  Solvent identifier (e.g. ``water``, ``toluene``, ``dichloromethane``)

   -  -  ``-so, --solvent-options``
      -  string
      -  Additional options appended inside the ``scrf=()`` keyword (e.g. ``iterative``)

.. note::

   Both ``-sm`` and ``-si`` must be provided together; specifying only one will raise an error when Gaussian processes
   the input file. ``-so`` is only applied when solvent is active — it is ignored when ``--remove-solvent`` is used.

Examples:

.. code:: bash

   # Gas phase optimization (project default)
   chemsmart sub gaussian -p anomer -f molecule.xyz -c 0 -m 1 -a no_solv opt

   # Solvated optimization with SMD/water (overrides project gas-phase settings)
   chemsmart sub gaussian -p anomer -f molecule.xyz -c 0 -m 1 -sm smd -si water -a solv opt

   # Solvated optimization with additional iterative SCRF option
   chemsmart sub gaussian -p anomer -f molecule.xyz -c 0 -m 1 -sm smd -si water -so iterative -a solv opt

   # Solvated TD-DFT calculation
   chemsmart sub gaussian -p anomer -f molecule.xyz -c 0 -m 1 -sm smd -si water -so iterative -a solv_td td

   # Remove solvent when project settings include one
   chemsmart sub gaussian -p solv_project -f molecule.xyz -c 0 -m 1 --remove-solvent opt

The solvated opt example above produces a route line of the form:

.. code:: text

   # opt freq b3lyp def2svp scrf=(smd,solvent=water,iterative)

and the solvated TD-DFT example gives:

.. code:: text

   # freq cam-b3lyp def2svp scrf=(smd,solvent=water,iterative) TD(singlets,nstates=3,root=1)

***********************
 Available Subcommands
***********************

Conformational Sampling & Dynamics
==================================

.. list-table::
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``crest``
      -  CREST conformational search jobs
   -  -  ``traj``
      -  Trajectory analysis jobs

Structure Optimization
======================

.. list-table::
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``opt``
      -  Geometry optimization

Transition State Search
=======================

.. list-table::
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``ts``
      -  Transition state optimization
   -  -  ``modred``
      -  Modified redundant coordinate optimization
   -  -  ``irc``
      -  Intrinsic reaction coordinate calculations
   -  -  ``scan``
      -  Potential energy surface scanning

Electronic Structure Properties
===============================

.. list-table::
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``sp``
      -  Single point energy calculation
   -  -  ``nci``
      -  Non-covalent interaction analysis
   -  -  ``dias``
      -  Distortion-Interaction/Activation-Strain analysis
   -  -  ``resp``
      -  RESP charge fitting
   -  -  ``td``
      -  Time-dependent DFT calculations
   -  -  ``wbi``
      -  Wiberg Bond Index analysis

Other Jobs
==========

.. list-table::
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``com``
      -  Run Gaussian input file as-is
   -  -  ``link``
      -  Multi-step link jobs
   -  -  ``userjob``
      -  Custom user-defined jobs

************
 Next Steps
************

For detailed information on each job type:

-  :doc:`gaussian-structure-optimization`
-  :doc:`gaussian-transition-state`
-  :doc:`gaussian-conformational-sampling`
-  :doc:`gaussian-qrc`
-  :doc:`gaussian-electronic-structure`
-  :doc:`gaussian-other-jobs`
