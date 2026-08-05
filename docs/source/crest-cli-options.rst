###################
 CREST CLI Options
###################

This page documents the CLI options available for all CREST jobs. Use ``chemsmart sub crest --help`` for the complete
list.

*************************
 Basic Command Structure
*************************

.. code:: bash

   chemsmart sub [OPTIONS] crest [CREST_OPTIONS] <SUBCMD> [SUBCMD_OPTIONS]

***************
 CREST Options
***************

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
      -  Project settings from ``~/.chemsmart/crest/*.yaml``

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

   chemsmart sub -s server crest -p test -f molecule.xyz -l custom_name conformers

This creates a job folder labeled ``custom_name`` instead of the default ``molecule_conformers``.

Use ``-a`` to append a string to the base filename:

.. code:: bash

   chemsmart sub -s server crest -p test2 -f molecule.xyz -a solv conformers

This creates ``molecule_solv`` instead of ``molecule_conformers``.

Selecting Structures
--------------------

Use ``-i`` to select a specific structure from multi-structure files:

.. code:: bash

   chemsmart sub -s server crest -p test -f molecules.xyz -i 5 -c 0 -m 1 conformers

This uses the 5th structure (1-indexed) from the XYZ file.

Using PubChem
-------------

Fetch structures directly from PubChem:

.. code:: bash

   chemsmart sub -s server crest -p test -P 356 -c 0 -m 1 -l octane conformers

This fetches octane (CID 356) and runs CREST labeled ``octane``.

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
   chemsmart sub -s server crest -p test -f molecule.xyz -c -1 -m 2 conformers

   # Triplet
   chemsmart sub -s server crest -p test -f molecule.xyz -c 0 -m 3 conformers

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
      -  GFN method: ``gfn1``, ``gfn2``, ``gfnff``, or the ``gfn2//gfnff`` composite protocol

   -  -  ``-O, --optimization-level``
      -  choice
      -  Optimization level: ``crude``, ``sloppy``, ``loose``, ``lax``, ``normal``, ``tight``, ``vtight``, ``extreme``

Examples:

.. code:: bash

   # GFN2-xTB conformational search with very tight optimization level
   chemsmart sub -s server crest -p test -f molecule.xyz -g gfn2 -O vtight conformers

   # GFN-FF method
   chemsmart sub -s server crest -p test -f molecule.xyz -g gfnff conformers

   # GFN2-xTB//GFN-FF composite method
   chemsmart sub -s server crest -p test -f molecule.xyz -g gfn2//gfnff conformers

Route and Calculation Options
=============================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-r, --additional-flags``
      -  string
      -  Extra CREST CLI flags appended to the generated command

   -  -  ``--nci/--no-nci``
      -  bool
      -  Enable or disable non-covalent interaction mode

   -  -  ``--ewin``
      -  float
      -  Energy window in kcal/mol for conformer selection

Examples:

.. code:: bash

   # Append extra CREST flags (whitespace-separated)
   chemsmart sub -s server crest -p test -f molecule.xyz -r "--norotmd --niceprint" conformers

   # Non-covalent interaction mode
   chemsmart sub -s server crest -p test -f molecule.xyz --nci conformers

   # Energy window for conformer selection
   chemsmart sub -s server crest -p test -f molecule.xyz --ewin 6.0 conformers

Solvent Options
===============

Solvent options are specified at the CREST group level and apply to all job types.

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-sm, --solvent-model``
      -  string
      -  Implicit solvent model (e.g. ``alpb``, ``gbsa``)

   -  -  ``-si, --solvent-id``
      -  string
      -  Solvent name recognized by CREST (e.g. ``water``, ``toluene``)

   -  -  ``--remove-solvent/--no-remove-solvent``
      -  bool
      -  Remove solvent settings inherited from the project YAML

.. important::

   CHEMSMART renders solvent flags only when **both** ``solvent_model`` and ``solvent_id`` are set. Specifying only one
   of them leaves the calculation in the gas phase.

Examples:

.. code:: bash

   # ALPB(water) conformational search
   chemsmart sub -s server crest -p test -f molecule.xyz -sm alpb -si water conformers

   # Override a solvated project to gas phase
   chemsmart sub -s server crest -p solv_project -f molecule.xyz --remove-solvent conformers

Database Input
==============

CREST jobs can take geometries from a chemsmart ``.db`` file using the selectors:

.. code:: bash

   # By record index (last structure of that record by default)
   chemsmart sub -s server crest -p test -f results.db --ri 3 -c 0 -m 1 conformers

   # By structure ID
   chemsmart sub -s server crest -p test -f results.db --sid c4d5e6f78a9b -c 0 -m 1 conformers

From Other Program Outputs
==========================

CREST can also start from Gaussian, ORCA, or xTB outputs:

.. code:: bash

   # From Gaussian log
   chemsmart sub -s server crest -p test -f water_opt.log conformers

   # From ORCA output
   chemsmart sub -s server crest -p test -f molecule.out conformers

   # From xTB main output
   chemsmart sub -s server crest -p test -f water_ohess/water_ohess.out conformers

See :doc:`molecule-input-formats` for the full list of supported geometry sources.

************************
 How Commands Are Built
************************

For a typical CREST conformational search, CHEMSMART builds a command like:

.. code:: text

   crest molecule.xyz --gfn2 --chrg 0 --uhf 0

Solvent (when both model and id are set) adds ``--<model> <id>``, for example ``--alpb water``. Additional command-line
arguments provided via ``-r`` are appended at the end.

For constrained searches, CHEMSMART additionally generates a ``constraints.inp`` file containing the specified distance,
angle, and/or dihedral constraints and, if specified, the requested force constant.

************
 Next Steps
************

-  :doc:`crest-conformational-search` — Detailed workflows for free and constrained conformational searches
-  :doc:`configuration-project-settings` — ``~/.chemsmart/crest/*.yaml``
-  :doc:`molecule-input-formats` — Supported input formats
