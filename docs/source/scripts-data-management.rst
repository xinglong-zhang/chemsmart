Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to
contact our team if you have questions or feedback.

#################
 File Management
#################

**************************
 File Organization Script
**************************

The ``file_organizer.py`` script organizes computational chemistry output files based on an Excel spreadsheet. It
creates folders, renames files, and moves them into their corresponding directories.

USAGE
=====

.. code:: console

   file_organizer.py [-d path/to/directory] [-f excel_file] [-t filetype]
                     [-n sheet_name] [-c columns] [-s skip_row(s)]
                     [-r organize_row(s)] [--keep-default-na|--no-keep-default-na]

OPTIONS
=======

.. list-table::
   :header-rows: 1
   :widths: 20 10 70

   -  -  Option
      -  Type
      -  Description

   -  -  ``-d, --directory``
      -  string
      -  |  Specifies the directory containing files to be organized.
         |  Defaults to the current working directory if not provided.

   -  -  ``-f, --filename``
      -  string
      -  |  **Required**. Specifies the Excel file containing metadata for organizing files.

   -  -  ``-t, --type``
      -  string
      -  |  Specifies the extension of files to organize (e.g., ``log``, ``out``, ``xyz``).
         |  Default is ``log``.

   -  -  ``-n, --name``
      -  string
      -  |  **Required**. Specifies the sheet name in the Excel file to read metadata from.

   -  -  ``-c, --cols``

      -  string

      -  |  Specifies a comma-separated list of Excel columns or column ranges (e.g. ``A:C``
         |  or ``A,C,E``) to use for organizing files.
         |  The standard structure consists of three columns: the target folder, the new
         |  file name, and the old file name.
         |  Default is ``B:D``.

   -  -  ``-s, --skip``
      -  int
      -  |  Number of rows in the Excel sheet to skip at the start when reading metadata.

   -  -  ``-r, --row``
      -  int
      -  |  Number of rows to read from the Excel sheet for organizing.
         |  If not provided, all rows except skipped rows are organized.

   -  -  ``--keep-default-na / --no-keep-default-na``

      -  bool

      -  |  Specifies whether to include the default *NaN* values (e.g., `#N/A`, `NULL`) when
         |  reading the Excel file.
         |  By default, these values are treated as *NaN* and ignored.

EXAMPLE
=======

In the ``co2`` worksheet of `test.xlsx
<https://github.com/xinglong-zhang/chemsmart/blob/main/tests/data/IOTests/test.xlsx>`_, rows 32–43 contain the
computational data of twelve conformers of *dimer_sm_epoxide*. After applying the energy corrections, the relative
energies are listed in column Y.

   .. image:: _static/file_organizer_example.png
      :width: 700
      :align: center

To rename these conformers in ascending order of the corrected energies — that is, to rename the files listed in column
**D** (*original runtime name*) to those in column **C** (*final file name*), and move them into the subdirectories
specified in column **B** (*folder name*) — run:

.. code:: console

   file_organizer.py -f test.xlsx -n co2 -c B:D -s 2 -r 45

This command skips the first **2** rows of the sheet and processes up to the first **45** rows.

If the target folder ``final_logs`` does not exist, it will be created automatically in the current working directory.
The script then locates the file ``dimer_sm_epoxide_c5.log`` in the current directory, and saves it to the
``final_logs`` folder with the new name ``dimer_sm_epoxide.log``, while keeping the original file unchanged. The
remaining conformers are processed in the same way.

************************
 File Conversion Script
************************

The ``file_converter.py`` script converts structure files from one format to another. It supports various file types and
allows batch processing of multiple files.

USAGE
=====

.. code:: console

   file_converter.py [-d path/to/directory] [-t filetype] [-f filename] [-o output_type] [-i]

OPTIONS
=======

.. list-table::
   :header-rows: 1
   :widths: 20 10 70

   -  -  Option
      -  Type
      -  Description

   -  -  ``-d, --directory``
      -  string
      -  |  Specifies the directory containing files to be converted.
         |  This option is mutually exclusive with ``-f``.

   -  -  ``-t, --type``

      -  string

      -  |  Specifies the file type in the directory. The currently supported file types
         |  include ``log``, ``com``, ``gjf``, ``out``, ``inp``, ``xyz``, and ``sdf``.
         |  This option is only used together with ``-d``.

   -  -  ``-f, --filename``
      -  string
      -  |  Specifies the file name(s) to be converted.
         |  This option is mutually exclusive with ``-d``.

   -  -  ``-o, --output-filetype``

      -  string

      -  |  Specifies the desired output file type. The currently supported output file
         |  types are ``xyz`` and ``com``.
         |  Default is ``xyz``.

   -  -  ``-i, --include-intermediate-structures / --no-include-intermediate-structures``

      -  bool

      -  |  Specifies whether to include intermediate structures found in the input
         |  files during conversion.
         |  By default, only the final structure is converted.

EXAMPLE
=======

**Single file conversion**

This example converts the Gaussian log file `co2.log
<https://github.com/xinglong-zhang/chemsmart/blob/main/tests/data/GaussianTests/outputs/co2.log>`_ to an XYZ file.

.. code:: console

   Entering Gaussian System, Link 0=/jet/home/xzhangj/programs/g16/g16
   Input=/jet/home/xzhangj/scratch/co2/co2.com
   Output=/jet/home/xzhangj/scratch/co2/co2.log
   Initial command:
   /jet/home/xzhangj/programs/g16/l1.exe "/ocean/projects/che220035p/xzhangj/scratch/co2/Gau-59
   913.inp" -scrdir="/ocean/projects/che220035p/xzhangj/scratch/co2/"
   Entering Link 1 = /jet/home/xzhangj/programs/g16/l1.exe PID=     59915.

   Copyright (c) 1988-2017, Gaussian, Inc.  All Rights Reserved.

   ...

                           Standard orientation:
   ---------------------------------------------------------------------
   Center     Atomic      Atomic             Coordinates (Angstroms)
   Number     Number       Type             X           Y           Z
   ---------------------------------------------------------------------
        1          8           0        0.000000    0.000000    1.163062
        2          8           0        0.000000    0.000000   -1.163062
        3          6           0        0.000000    0.000000    0.000000
   ---------------------------------------------------------------------
   Rotational constants (GHZ):           0.0000000          11.6788340          11.6788340

   ...

   Job cpu time:       0 days  0 hours  4 minutes 31.3 seconds.
   Elapsed time:       0 days  0 hours  0 minutes  4.9 seconds.
   File lengths (MBytes):  RWF=      6 Int=      0 D2E=      0 Chk=      1 Scr=      1
   Normal termination of Gaussian 16 at Wed Nov 15 21:13:40 2023.

To convert the file:

.. code:: console

   file_converter.py -f co2.log

The output file ``co2.xyz`` will be generated in the same directory:

.. code:: console

   3
   co2.xyz    Empirical formula: CO2    Energy(Hartree): -188.444680
   O        0.0000000000    0.0000000000    1.1630620000
   O        0.0000000000    0.0000000000   -1.1630620000
   C        0.0000000000    0.0000000000    0.0000000000

**Batch conversion for all files of a specified type**

In addition to single-file conversion, batch mode can be used to process multiple files of a given type.

.. code:: console

   file_converter.py -d . -t log -o com -i

This command converts all Gaussian log files (``.log``) in the current directory to Gaussian input files (``.com``),
including all intermediate structures.
