######################
 Inspecting Databases
######################

The ``inspect`` subcommand provides human-readable summaries and detailed views of a chemsmart database. Without
selectors, it prints a database-level overview. With selectors, it drills down to a single record, structure, or
molecule.

*******
 Usage
*******

.. code:: text

   chemsmart run database inspect -f database.db
                                  [--ri record_index | --rid record_id |
                                   --mid molecule_id | --sid structure_id]
                                  [--si structure_index]

*********
 Options
*********

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   -  -  Option
      -  Type
      -  Description

   -  -  ``-f, --file``
      -  string
      -  Path to the input database file (``.db``). Required.

   -  -  ``--ri, --record-index``
      -  int
      -  1-based index of the record to inspect.

   -  -  ``--rid, --record-id``
      -  string
      -  Record ID, or unique prefix of at least 12 characters, to inspect.

   -  -  ``--si, --structure-index``
      -  int
      -  1-based index of a structure within the selected record. Requires ``--ri`` or ``--rid``.

   -  -  ``--mid, --molecule-id``
      -  string
      -  Molecule ID, or unique prefix of at least 16 characters, to inspect.

   -  -  ``--sid, --structure-id``
      -  string
      -  Structure ID, or unique prefix of at least 12 characters, to inspect.

.. note::

   -  The selectors ``--ri``, ``--rid``, ``--mid`` and ``--sid`` are mutually exclusive.
   -  Use ``--si`` only together with ``--ri`` or ``--rid``.

***************
 Display Modes
***************

.. list-table::
   :header-rows: 1
   :widths: 30 70

   -  -  Selectors
      -  Output

   -  -  *(none)*
      -  Database-level overview: total counts of records, molecules, and structures, plus breakdowns by program and job
         type.

   -  -  ``--ri`` / ``--rid``
      -  Detailed view of a single record: meta, results, thermochemistry, provenance, and a summary table of all
         structures in the record.

   -  -  ``--ri``/``--rid`` + ``--si``
      -  Detailed view of one structure within a record, including available per-structure properties such as energy,
         forces, and Cartesian coordinates.

   -  -  ``--mid``
      -  Detailed view of a molecule, plus a list of records that reference it.

   -  -  ``--sid``
      -  Detailed view of a structure, plus a list of records that reference it.

**********
 Examples
**********

**Database-level overview:**

.. code:: bash

   chemsmart run database inspect -f chemsmart.db

**Detailed view of a record by index:**

.. code:: bash

   chemsmart run database inspect -f chemsmart.db --ri 3

**Detailed view of a record by ID prefix:**

.. code:: bash

   chemsmart run database inspect -f chemsmart.db --rid a1b2c3d4e5f6

**Detailed view of a single structure inside a record:**

.. code:: bash

   chemsmart run database inspect -f chemsmart.db --ri 3 --si 1

**Detailed view of a molecule:**

.. code:: bash

   chemsmart run database inspect -f chemsmart.db --mid ABCDEFGHIJKLMN-U

**Detailed view of a stand-alone structure:**

.. code:: bash

   chemsmart run database inspect -f chemsmart.db --sid c4d5e6f78a9b
