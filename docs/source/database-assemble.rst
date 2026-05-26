####################
 Assembling Records
####################

The ``assemble`` subcommand scans a project directory for supported quantum chemistry output files, parses the extracted
calculation data, and stores the results as records in a CHEMSMART database.

Each record is uniquely identified by its ``record_id``. If a newly parsed output file generates a ``record_id`` that
already exists in the database, the existing record is updated with the newly parsed data rather than duplicated. This
allows ``assemble`` to be safely re-run on the same project directory after additional calculations have finished or
existing output files have been updated.

*******
 Usage
*******

.. code:: text

   chemsmart run database assemble [-d path/to/directory] [-p gaussian|orca]
                                   [-i index] [-o outfile.db] [--include-failed]

*********
 Options
*********

.. list-table::
   :header-rows: 1
   :widths: 18 12 70

   -  -  Option
      -  Type
      -  Description

   -  -  ``-d, --directory``
      -  string
      -  Root directory to scan recursively for supported output files. Defaults to the current working directory.

   -  -  ``-p, --program``
      -  string
      -  Restrict parsing to output files from a specific program. Supported values are ``gaussian`` and ``orca``. If
         omitted, outputs from all supported programs are scanned.

   -  -  ``-i, --index``
      -  string
      -  1-based index or slice used to select structures from files containing multiple structures. Defaults to ``:``
         (all structures).

   -  -  ``-o, --output``
      -  string
      -  Output database filename. The ``.db`` extension is appended automatically if it is not provided. Defaults to
         ``database.db``.

   -  -  ``--include-failed``
      -  flag
      -  Include parseable partial data from failed calculations, which are skipped by default and marked by termination
         status in the database.

**********
 Examples
**********

**Assemble all output files under the current directory:**

.. code:: bash

   chemsmart run database assemble

**Assemble only Gaussian outputs from a project directory:**

.. code:: bash

   chemsmart run database assemble -d results/ -p gaussian -o gaussian_results.db

**Assemble only the final structure from each file:**

.. code:: bash

   chemsmart run database assemble -d results/ -i -1 -o final_geoms.db
