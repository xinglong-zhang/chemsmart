#################
 xTB CLI Options
#################

This page documents the CLI options available for all xTB jobs. Use ``chemsmart run xtb --help`` for the complete list.

*************************
 Basic Command Structure
*************************

.. code:: bash

   chemsmart run [OPTIONS] xtb [XTB_OPTIONS] <SUBCMD> [SUBCMD_OPTIONS]

*************
 xTB Options
*************

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
      -  Project settings from ``~/.chemsmart/xtb/*.yaml``

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

   -  -  ``-P, --pubchem``
      -  string
      -  Query structure from PubChem (name, SMILES, CID)

   -  -  ``--ri, --record-index``
      -  int
      -  Select a record from a chemsmart database by its 1-based index

   -  -  ``--rid, --record-id``
      -  string
      -  Select a record from a chemsmart database by its ID

   -  -  ``--sid, --structure-id``
      -  string
      -  Select a structure from a chemsmart database by its ID

.. note::

   -  ``-p`` uses the project name without the ``.yaml`` extension.
   -  ``-f`` accepts various formats: ``.xyz``, ``.com``, ``.gjf``, ``.log``, ``.inp``, ``.out``, or a chemsmart
      database ``.db`` file.

Specifying Output Filenames
---------------------------

Use ``-l`` to set a custom label:

.. code:: bash

   chemsmart run xtb -p test -f water.xyz -l custom_name opt

This creates a job folder/label ``custom_name`` instead of the default ``water_opt``.

Use ``-a`` to append a string to the base filename:

.. code:: bash

   chemsmart run xtb -p test -f water.xyz -a solv opt

This creates ``water_solv`` instead of ``water_opt``.

Selecting Structures
--------------------

Use ``-i`` to select a specific structure from multi-structure files:

.. code:: bash

   chemsmart run xtb -p test -f molecules.xyz -i 5 -c 0 -m 1 opt

This uses the 5th structure (1-indexed) from the XYZ file.

.. warning::

   CHEMSMART uses 1-based indexing to match most molecular visualization software.

Using PubChem
-------------

Fetch structures directly from PubChem:

.. code:: bash

   chemsmart run xtb -p test -P 962 -c 0 -m 1 -l water opt

This fetches water (CID 962) and runs an xTB optimization labeled ``water``.

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

.. note::

   If the input lacks charge/multiplicity, specify them with ``-c`` and ``-m``. CHEMSMART passes charge as ``--chrg``
   and unpaired electrons as ``--uhf`` (``multiplicity - 1``).

Examples:

.. code:: bash

   # Anion doublet
   chemsmart run xtb -p test -f molecule.xyz -c -1 -m 2 opt

   # Triplet
   chemsmart run xtb -p test -f molecule.xyz -c 0 -m 3 sp

Method Options
==============

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-g, --gfn-version``
      -  choice
      -  GFN method: ``gfn0``, ``gfn1``, ``gfn2``, or ``gfnff``

Examples:

.. code:: bash

   # GFN1-xTB optimization
   chemsmart run xtb -p test -f molecule.xyz -g gfn1 opt

   # GFN-FF single point
   chemsmart run xtb -p test -f molecule.xyz -g gfnff sp

Route and Calculation Options
=============================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-r, --additional-route-parameters``
      -  string
      -  Extra xTB CLI flags appended to the generated command

   -  -  ``--grad/--no-grad``
      -  bool
      -  Enable or disable gradient output (``--grad``)

Examples:

.. code:: bash

   # Append extra xTB flags (whitespace-separated)
   chemsmart run xtb -p test -f molecule.xyz -r "--copy --json" sp

   # Write gradient information
   chemsmart run xtb -p test -f molecule.xyz --grad opt

Solvent Options
===============

Solvent options are specified at the xTB group level and apply to all job types.

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-sm, --solvent-model``
      -  string
      -  Implicit solvent model (e.g. ``alpb``, ``gbsa``, ``cosmo``)

   -  -  ``-si, --solvent-id``
      -  string
      -  Solvent name recognized by xTB (e.g. ``water``, ``toluene``)

   -  -  ``--remove-solvent/--no-remove-solvent``
      -  bool
      -  Remove solvent settings inherited from the project YAML

.. important::

   CHEMSMART renders solvent flags only when **both** ``solvent_model`` and ``solvent_id`` are set. Specifying only one
   of them leaves the calculation in the gas phase.

Examples:

.. code:: bash

   # ALPB(water) optimization
   chemsmart run xtb -p test -f molecule.xyz -sm alpb -si water opt

   # Override a solvated project to gas phase
   chemsmart run xtb -p solv_project -f molecule.xyz --remove-solvent opt

Database Input
==============

xTB jobs can take geometries from a chemsmart ``.db`` file using the selectors:

.. code:: bash

   # By record index (last structure of that record by default)
   chemsmart run xtb -p test -f results.db --ri 3 -c 0 -m 1 opt

   # By structure ID
   chemsmart run xtb -p test -f results.db --sid c4d5e6f78a9b -c 0 -m 1 sp

From Other Program Outputs
==========================

xTB can also start from Gaussian/ORCA outputs, or from an existing xTB main ``.out``:

.. code:: bash

   # From Gaussian log
   chemsmart run xtb -p test -f water_opt.log hess

   # From xTB main output (parent folder is resolved automatically)
   chemsmart run xtb -p test -f co2_ohess/co2_ohess.out sp

See :doc:`molecule-input-formats` for the full list of supported geometry sources.

************************
 How Commands Are Built
************************

For a typical GFN2-xTB single point, CHEMSMART writes ``{label}.xyz`` and runs a command similar to:

.. code:: text

   xtb {label}.xyz --gfn 2 --chrg 0 --uhf 0

Job-type flags:

-  ``opt`` adds ``--opt <optimization_level>``
-  ``hess`` adds ``--hess``
-  ``sp`` adds no additional job-type flags

Solvent (when both model and id are set) adds ``--<model> <id>``, for example ``--alpb water``. Additional command-line
arguments provided via ``-r`` are appended at the end.

************
 Next Steps
************

-  :doc:`xtb-structure-optimization` — ``opt``, ``sp``, and ``hess`` workflows
-  :doc:`configuration-project-settings` — ``~/.chemsmart/xtb/*.yaml``
-  :doc:`pymol-visualization` — visualize xTB structures and trajectories
-  :doc:`thermochemistry-analysis` — post-process xTB calculations
-  :doc:`database-assemble` — assemble xTB folders into a chemsmart database
