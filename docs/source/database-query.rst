####################
 Querying Databases
####################

The ``query`` subcommand searches a chemsmart database using a compact expression language. The query target can be
switched between records, molecules, or structures.

*******
 Usage
*******

.. code:: text

   chemsmart run database query -f database.db
                                [-t records|molecules|structures]
                                [-q "query expression"]
                                [-l limit]

*********
 Options
*********

.. list-table::
   :header-rows: 1
   :widths: 18 12 70

   -  -  Option
      -  Type
      -  Description

   -  -  ``-f, --file``
      -  string
      -  Path to the input database file (``.db``). Required.

   -  -  ``-t, --target``
      -  string
      -  Query target: ``records`` for calculation records, ``molecules`` for unique chemical species, or ``structures``
         for unique 3D structures. Default: ``records``.

   -  -  ``-q, --query``
      -  string
      -  Query expression. If omitted, all entries for the chosen target are listed.

   -  -  ``-l, --limit``
      -  int
      -  Maximum number of results to display. Must be a positive integer.

***************************
 Query Expression Language
***************************

A query expression is a logical combination of conditions of the form ``FIELD OPERATOR VALUE``. String values must be
quoted with single quotes. Conditions can be combined with the case-insensitive ``AND`` and ``OR`` keywords.

**Supported operators:**

.. list-table::
   :header-rows: 1
   :widths: 12 88

   -  -  Operator
      -  Meaning
   -  -  ``<``
      -  Less than. Numeric fields only.
   -  -  ``<=``
      -  Less than or equal to. Numeric fields only.
   -  -  ``>``
      -  Greater than. Numeric fields only.
   -  -  ``>=``
      -  Greater than or equal to. Numeric fields only.
   -  -  ``=``
      -  Equal to. Numeric or string fields.
   -  -  ``!=``
      -  Not equal to. Numeric or string fields.
   -  -  ``~``
      -  Case-insensitive substring match. String fields only.

**Queryable fields by target:**

.. list-table::
   :header-rows: 1
   :widths: 18 82

   -  -  Target
      -  Fields

   -  -  ``records``

      -  ``program``, ``method``, ``basis``, ``jobtype``, ``solvent_on``, ``solvent_model``, ``total_energy``,
         ``homo_energy``, ``lumo_energy``, ``fmo_gap``, ``zero_point_energy``, ``enthalpy``, ``entropy``,
         ``gibbs_free_energy``, ``source_file``

   -  -  ``molecules``
      -  ``chemical_formula``, ``smiles``, ``inchi``, ``number_of_atoms``, ``mass``

   -  -  ``structures``
      -  ``chemical_formula``, ``number_of_atoms``, ``charge``, ``multiplicity``

.. note::

   **Units for numeric fields:**

   -  ``total_energy``, ``zero_point_energy``, ``enthalpy``, ``gibbs_free_energy``: Hartree (Eh).
   -  ``homo_energy``, ``lumo_energy``, ``fmo_gap``: eV.
   -  ``entropy``: Hartree/K.
   -  ``mass``: amu.

**********
 Examples
**********

**List all records:**

.. code:: bash

   chemsmart run database query -f chemsmart.db

**Filter records by solvent model:**

.. code:: bash

   chemsmart run database query -f chemsmart.db -q "solvent_on = 1 AND solvent_model = 'smd'"

**Filter records by source file name:**

.. code:: bash

   chemsmart run database query -f chemsmart.db -q "source_file ~ 'benzene'"

**List all unique molecules:**

.. code:: bash

   chemsmart run database query -f chemsmart.db -t molecules

**Filter molecules by molecular mass:**

.. code:: bash

   chemsmart run database query -f chemsmart.db -t molecules -q "mass > 100"

**List the first 10 unique structures in the database:**

.. code:: bash

   chemsmart run database query -f chemsmart.db -t structures -l 10

**Filter structures by charge and multiplicity:**

.. code:: bash

   chemsmart run database query -f chemsmart.db -t structures -q "charge = 0 AND multiplicity = 1"
