####################
 Assembling Records
####################

The ``assemble`` subcommand scans a project directory for supported quantum chemistry output files and calculation
folders, parses the available calculation data, and stores the results as records in a CHEMSMART database.

Each record is uniquely identified by its ``record_id``. If a newly parsed calculation generates a ``record_id`` that
already exists in the database, the existing record is updated with the newly parsed data rather than duplicated. This
allows ``assemble`` to be safely re-run on the same project directory after additional calculations have finished or
existing outputs have been updated.

*******
 Usage
*******

.. code:: text

   chemsmart run database assemble [-d path/to/directory] [-p gaussian|orca|xtb] [-f dataset.npz]
                                   [--dataset revqm9|qm9] [--method method] [--basis basis]
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
      -  Root directory to scan recursively for supported output files and folders. Defaults to the current working
         directory.

   -  -  ``-p, --program``
      -  string
      -  Restrict parsing to output files/folders from a specific program. Supported values are ``gaussian``, ``orca``,
         and ``xtb``. If omitted, outputs from all supported programs are scanned.

   -  -  ``-f, --file``
      -  string
      -  Import an NPZ dataset. Cannot be combined with ``--directory``, ``--filetype``, ``--program``, or
         ``--include-failed``.

   -  -  ``--dataset``
      -  string
      -  With ``--file``, declare the NPZ dataset and its units: ``revqm9`` or ``qm9``.

   -  -  ``--method``
      -  string
      -  With ``--file``, the calculation method to record. Defaults to the method of the declared dataset.

   -  -  ``--basis``
      -  string
      -  With ``--file``, the basis set to record. Defaults to the basis set of the declared dataset.

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

**Assemble all supported calculations under the current directory:**

.. code:: bash

   chemsmart run database assemble

**Assemble only Gaussian outputs from a project directory:**

.. code:: bash

   chemsmart run database assemble -d results/ -p gaussian -o gaussian_results.db

**Assemble only xTB calculation outputs from a project directory:**

For xTB calculations, ``assemble`` recursively scans the root directory to arbitrary depth. Each recognized xTB
calculation directory is treated as one record.

Example directory layout:

.. code::

   xtb_results/
   ├── molecule1_ohess/
   ├── project_a/
   │   ├── molecule2_hess/
   │   └── nested/
   │       └── molecule3_opt/
   └── project_b/
       └── molecule4_sp/

.. code:: bash

   chemsmart run database assemble -d xtb_results/ -p xtb -o xtb_results.db

Every xTB calculation directory must contain exactly one xTB main output together with the required auxiliary files.
Directories containing multiple xTB main outputs are ambiguous: CHEMSMART logs an error, skips that directory, and
continues assembling other valid calculations.

**Import a supported NPZ dataset:**

.. code:: bash

   chemsmart run database assemble -f revQM9.npz --dataset revqm9 -o dataset.db

Each selected structure becomes a separate record. For qm9_eV, orbital energies are stored in eV, while the
atom-reference-adjusted U0 is not stored as total energy. Per-molecule dataset properties are preserved as original
values with units where available, accessible through ``Database.get_record()`` in
``record["molecules"][0]["dataset_properties"]``. revQM9 orbital and atomization energies use Hartree, and dipole
moments use Debye. The dataset-level ``meta`` and unmatched ``atomization_test`` array are not imported. The supported
QM9 and revQM9 layouts are treated as neutral closed-shell molecules (charge 0, multiplicity 1).

Use ``--dataset revqm9`` or ``--dataset qm9`` to declare the dataset and its units. The default method and basis are
aPBE0/cc-pVTZ for revQM9 and B3LYP/6-31G(2df,p) for QM9; ``--method`` and ``--basis`` override these values. Without
``--dataset``, the layout is detected and energies are still recorded, assuming the units of that layout, and CHEMSMART
logs a warning. The source, the units of the stored dataset properties, and the default method and basis are then left
unset. The method and basis can still be supplied explicitly. Dataset declarations must match the file layout. These
options apply only to database import; using an NPZ as a calculation input does not require them and expects coordinates
in Angstrom.

Calculated values are stored in the matching database columns, in each column's unit:

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   -  -  Database field
      -  QM9
      -  revQM9

   -  -  ``total_energy``
      -  Not stored
      -  ``energies``

   -  -  ``homo_energy``, ``lumo_energy``
      -  ``homo``, ``lumo``
      -  ``homo``, ``lumo``

   -  -  ``fmo_gap``
      -  ``gap``
      -  ``lumo`` minus ``homo``

   -  -  ``zero_point_energy``
      -  ``zpve``
      -  Not in the dataset

   -  -  ``rotational_constants_in_Hz``
      -  ``A``, ``B``, ``C``
      -  Not in the dataset

   -  -  ``dipole_moment_magnitude``
      -  ``mu``
      -  ``dipole``

QM9's ``U0``, ``U``, ``H``, and ``G`` are atom-reference-adjusted, so they are not comparable with the absolute
``internal_energy``, ``enthalpy``, and ``gibbs_free_energy`` columns and are not stored there. These, and values with no
matching column (``alpha``, ``r2``, ``Cv``, and ``atomization``), are available in ``dataset_properties``.

Reimporting a record with a declared dataset updates that record and retains its record index. If the existing record
has a declared dataset, reimporting the same structure without ``--dataset`` skips it to avoid discarding known units
and metadata. Specify the dataset again when updating the method or basis.

A structure that cannot be assembled is skipped and listed in a warning, and the rest of the import continues. Progress
is kept if an import stops early, and running it again updates the existing records.

**Assemble only the final structure from each calculation source:**

.. code:: bash

   chemsmart run database assemble -d results/ -i -1 -o final_geoms.db
