Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to
contact our team if you have questions or feedback.

##########################
 File and Data Management
##########################

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

In the ``co2`` worksheet of ``test.xlsx``, rows 32–43 contain the computational data of twelve conformers of
*dimer_sm_epoxide*. After applying the energy corrections, the relative energies are listed in column Y.

   .. image:: _static/file_organizer_example.png
      :width: 700
      :align: center

To rename these conformers in ascending order of the corrected energies — that is, to rename the files listed in column
**D** (*original runtime name*) to those in column **C** (*final file name*), and move them into the subdirectories
specified in column **B** (*folder name*) — run:

.. code:: console

   file_organizer.py -f test.xlsx -n co2 -c B:D -s 2 -r 45

This command skips the first **2** rows of the sheet and processes up to the first **45** rows. If the target folder
``final_logs`` does not exist, it will be created automatically in the current working directory. The script then
locates the file ``dimer_sm_epoxide_c5.log`` in the current directory, and saves it to the ``final_logs`` folder with
the new name ``dimer_sm_epoxide.log``, while keeping the original file unchanged. The remaining conformers are processed
in the same way.
