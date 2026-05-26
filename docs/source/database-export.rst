#####################
 Exporting Databases
#####################

The ``export`` subcommand writes data from a chemsmart database to a portable file format. The output format is inferred
from the extension of the output file specified with ``-o/--output``.

*******
 Usage
*******

.. code:: text

   chemsmart run database export -f database.db -o outfile.{json|csv|xyz|extxyz}
                                 [--ri record_index | --rid record_id | --mid molecule_id | --sid structure_id]
                                 [--si structure_index]
                                 [-k keys]
                                 [-x method/basis]

****************
 Output Formats
****************

.. list-table::
   :header-rows: 1
   :widths: 12 28 60

   -  -  Extension
      -  Scope
      -  Description

   -  -  ``.json``
      -  Whole database
      -  Full structured dump of all records, including meta, results, provenance, and all stored structures with their
         available per-structure properties. Suitable for archival and programmatic post-processing.

   -  -  ``.csv``
      -  Whole database
      -  Tabular dump of scalar properties. Default columns are ``record_index``, ``record_id``, and
         ``chemical_formula``. Additional columns can be added with ``-k/--keys``.

   -  -  ``.xyz``
      -  Selected structures
      -  Cartesian coordinates of one or more selected structures. The selection is controlled by
         ``--ri``/``--rid``/``--si``/``--sid``/``--mid``.

   -  -  ``.extxyz``
      -  Selected structures
      -  Extended XYZ format. Each frame includes a machine-readable ``key=value`` header with the per-frame energy and
         per-atom forces. The selection is controlled by ``--ri``/``--rid``/``--si``/``--sid``/``--mid``.

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

   -  -  ``-o, --output``
      -  string
      -  Output file path. Format is inferred from the extension. Required.

   -  -  ``--ri, --record-index``
      -  int
      -  1-based index of the record to export. Only valid for XYZ/extended XYZ export.

   -  -  ``--rid, --record-id``
      -  string
      -  Record ID, or unique prefix of at least 12 characters, to export. Only valid for XYZ/extended XYZ export.

   -  -  ``--si, --structure-index``
      -  string
      -  1-based index of a structure within the selected record. Accepts an integer or a slice, such as ``:`` for the
         full trajectory. Defaults to the last structure when omitted with ``--ri``/``--rid``.

   -  -  ``--mid, --molecule-id``
      -  string
      -  Molecule ID, or unique prefix of at least 16 characters, to export. Only valid for XYZ/extended XYZ export.

   -  -  ``--sid, --structure-id``
      -  string
      -  Structure ID, or unique prefix of at least 12 characters, to export. Only valid for XYZ/extended XYZ export.

   -  -  ``-k, --keys``
      -  string
      -  Comma-separated list of additional scalar columns for CSV export.

   -  -  ``-x, --method-basis``
      -  string
      -  Filter XYZ/extended XYZ export by a specific ``method/basis`` pair. Only valid with ``--sid`` or ``--mid``.

.. note::

   -  ``--ri``, ``--rid``, ``--sid``, and ``--mid`` are mutually exclusive and only valid for XYZ/extended XYZ export.
   -  Use ``--si`` only together with ``--ri`` or ``--rid``.
   -  Use ``-x``/``--method-basis`` only together with ``--sid`` or ``--mid``.
   -  JSON and CSV exports always cover the entire database.

**Supported extra CSV keys** (``-k``):

``program``, ``method``, ``basis``, ``charge``, ``multiplicity``, ``smiles``, ``total_energy``, ``homo_energy``,
``lumo_energy``, ``fmo_gap``, ``zero_point_energy``, ``enthalpy``, ``entropy``, ``gibbs_free_energy``.

.. note::

   **Units:** ``total_energy``, ``zero_point_energy``, ``enthalpy``, and ``gibbs_free_energy`` are in Hartree (Eh).
   ``homo_energy``, ``lumo_energy``, and ``fmo_gap`` are in eV. ``entropy`` is Hartree/K.

**********
 Examples
**********

**Full structured dump as JSON:**

.. code:: bash

   chemsmart run database export -f my.db -o data.json

**CSV table for ML training with extra scalar columns:**

.. code:: bash

   chemsmart run database export -f my.db -o training.csv -k total_energy,homo_energy,fmo_gap

**Final geometry of a record (last structure) as XYZ:**

.. code:: bash

   chemsmart run database export -f my.db --rid a1b2c3d4e5f6 -o final.xyz

**Specific snapshot (5th structure) of record 3 as XYZ:**

.. code:: bash

   chemsmart run database export -f my.db --ri 3 --si 5 -o step5.xyz

**Full optimization trajectory of record 3 as a multi-frame XYZ file:**

.. code:: bash

   chemsmart run database export -f my.db --ri 3 --si ':' -o traj.xyz

**Single structure by structure ID as a single-frame XYZ file:**

.. code:: bash

   chemsmart run database export -f my.db --sid c4d5e6f78a9b -o struct.xyz

**All conformers of one molecule as a multi-frame XYZ file, sorted by increasing energy when available:**

.. code:: bash

   chemsmart run database export -f my.db --mid ABCDEFGHIJKLMN-U -o conformers.xyz

**Single structure exported as extended XYZ:**

.. code:: bash

   chemsmart run database export -f my.db --sid 0df6b2ea4bdc -o struct.extxyz

**All conformers of a molecule as extended XYZ, filtered to a specific level of theory:**

.. code:: bash

   chemsmart run database export -f my.db --mid BLQJIBCZHWBKSL-U -x 'MN15/def2tzvp' -o conformers.extxyz

.. note::

   For ``.extxyz`` export, each exported frame contains an ``energy`` field in **Hartree** and a ``forces`` field in
   **Hartree/Bohr**. These units are written explicitly in the header. Structures that do not have both energy and
   forces at the chosen level of theory are skipped, and a warning summarizing the skipped structures is logged.

   When ``-x``/``--method-basis`` is omitted for extended XYZ export with ``--mid``, chemsmart automatically picks the
   ``(method, basis)`` pair that covers the most structures with both energy and forces stored. With ``--sid``,
   chemsmart picks an available ``(method, basis)`` pair with both energy and forces for that structure. With
   ``--ri``/``--rid``, the selected record's own ``method/basis`` is used.

.. tip::

   Use ``chemsmart run database inspect -f my.db`` to verify which levels of theory have forces stored before attempting
   an extended XYZ export.
