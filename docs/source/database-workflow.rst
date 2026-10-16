######################
 Workflow Integration
######################

Once a CHEMSMART database has been assembled, its records and structures can be reused directly as inputs for other
CHEMSMART commands without any manual export step. This page describes how database files fit into visualization and new
calculation workflows.

***************************
 Using a Database as Input
***************************

Database files (``.db``) can be used as input for supported CHEMSMART workflows, including Gaussian jobs (see
:doc:`gaussian-cli-options`), ORCA jobs (see :doc:`orca-cli-options`), and PyMOL visualization (see
:doc:`pymol-cli-options`), in the same way as ordinary molecule files.

When ``-f`` points to a ``.db`` file, database selectors can be used to choose the input structure(s):

.. list-table::
   :header-rows: 1
   :widths: 20 28 52

   -  -  Selector
      -  Meaning
      -  Supported workflows

   -  -  ``--ri, --record-index``
      -  1-based record index inside the database.
      -  Visualization, Gaussian, ORCA

   -  -  ``--rid, --record-id``
      -  Record ID, or a unique prefix of at least 12 characters.
      -  Visualization, Gaussian, ORCA

   -  -  ``--sid, --structure-id``
      -  Structure ID, or a unique prefix of at least 12 characters.
      -  Visualization, Gaussian, ORCA

   -  -  ``--mid, --molecule-id``
      -  Molecule ID, or a unique prefix of at least 16 characters.
      -  Visualization only

   -  -  ``-i, --index``
      -  Used together with ``--ri``/``--rid`` to select a specific structure within the selected record.
      -  Visualization, Gaussian, ORCA

By default, when only a record is selected with ``--ri`` or ``--rid``, the *last* structure of that record is used as
the input geometry — typically the optimized geometry. Use ``-i`` to override this default.

.. note::

   -  For Gaussian and ORCA jobs, select database structures with ``--ri``/``--rid``/``--sid`` and, when needed, ``-i``.
   -  The ``--mid`` selector is currently supported by ``chemsmart run mol`` for visualization workflows.

**********
 Examples
**********

The examples below show how to reuse structures from an existing database file for PyMOL visualization and new
Gaussian/ORCA calculations.

**PyMOL visualization**

.. code:: bash

   # Visualize the final structure of a selected record.
   chemsmart run mol -f chemsmart.db --ri 3 visualize

   # Visualize a specific structure within a selected record.
   chemsmart run mol -f chemsmart.db --rid a1b2c3d4e5f6 -i 5 visualize

   # Visualize all structures within a selected record.
   chemsmart run mol -f chemsmart.db --ri 3 -i ':' visualize

   # Visualize all conformers of a molecule by molecule ID.
   chemsmart run mol -f chemsmart.db --mid ABCDEFGHIJKLMN-U visualize

   # Visualize a specific structure by structure ID.
   chemsmart run mol -f chemsmart.db --sid c4d5e6f78a9b visualize

**Gaussian jobs from database structures**

.. code:: bash

   # Submit a Gaussian optimization from the final structure of a selected record.
   chemsmart sub -s server gaussian -p project -f chemsmart.db --ri 3 opt

   # Submit a Gaussian single-point calculation from a specific structure within a record.
   chemsmart sub -s server gaussian -p project -f chemsmart.db --rid a1b2c3d4e5f6 -i 5 sp

   # Submit a Gaussian transition-state optimization from a structure ID.
   chemsmart sub -s server gaussian -p project -f chemsmart.db --sid c4d5e6f78a9b ts

**ORCA jobs from database structures**

.. code:: bash

   # Submit an ORCA optimization from the final structure of a selected record.
   chemsmart sub -s server orca -p project -f chemsmart.db --ri 3 opt

   # Submit an ORCA single-point calculation from a structure ID.
   chemsmart sub -s server orca -p project -f chemsmart.db --sid c4d5e6f78a9b sp

************************
 End-to-End Walkthrough
************************

**Basic workflow: assemble, query, inspect, export, reuse**

.. code:: bash

   # 1. Assemble calculation outputs into a project database.
   chemsmart run database assemble -d ./calculations -o project.db

   # 2. Find Gaussian single-point records with small FMO gaps.
   chemsmart run database query -f project.db \
       -q "program = 'gaussian' AND jobtype = 'sp' AND fmo_gap < 5"

   # 3. Inspect a specific record in detail using the record index from query output.
   chemsmart run database inspect -f project.db --ri 12

   # 4. Export selected scalar properties as a CSV table for plotting / ML.
   chemsmart run database export -f project.db \
       -o project.csv \
       -k total_energy,homo_energy,lumo_energy,fmo_gap

   # 5. Visualize or reuse a selected database structure.
   chemsmart run mol -f project.db --ri 12 visualize
   chemsmart sub -s server gaussian -p project -f project.db --ri 12 sp

This pattern — assemble, query, inspect, export, reuse — captures the most common interaction with CHEMSMART databases:
calculation outputs are collected into a database, inspected or filtered, exported when needed, and reused as inputs for
visualization or new calculations.

**Multi-program workflow: Gaussian → ORCA**

.. code:: bash

   # 1. Assemble existing Gaussian outputs.
   chemsmart run database assemble -d ./gaussian_results -p gaussian -o gaussian.db

   # 2. Select Gaussian optimization and transition-state records.
   chemsmart run database query -f gaussian.db \
       -q "jobtype = 'opt' OR jobtype = 'ts'"

   # 3. Submit new ORCA jobs from selected Gaussian geometries.
   chemsmart sub -s server orca -p project -f gaussian.db --ri 3 opt
   chemsmart sub -s server orca -p project -f gaussian.db --ri 5 ts

   # 4. Assemble the completed ORCA outputs.
   chemsmart run database assemble -d ./orca_results -p orca -o orca.db

   # 5. Export comparable energy tables from both databases.
   chemsmart run database export -f gaussian.db -o gaussian_energies.csv -k total_energy
   chemsmart run database export -f orca.db -o orca_energies.csv -k total_energy

**Conformer screening workflow:**

.. code:: bash

   # 1. Optimize geometries from a multi-structure CREST conformer file.
   chemsmart sub -s server gaussian -p conformer -f crest_conformers.xyz -c 0 -m 1 crest -j opt

   # 2. Assemble the completed conformer optimization outputs.
   chemsmart run database assemble -d ./conformer_opt -o conformers.db

   # 3. List unique molecules and their conformer counts.
   chemsmart run database query -f conformers.db -t molecules

   # 4. Visualize all conformers of one molecule, sorted by increasing energy when available.
   chemsmart run mol -f conformers.db --mid molecule_id_prefix visualize

   # 5. Export all conformers of the molecule in XYZ format
   #    for downstream analysis, e.g. MLIP training.
   chemsmart run database export -f conformers.db --mid molecule_id_prefix -o conformers.xyz
