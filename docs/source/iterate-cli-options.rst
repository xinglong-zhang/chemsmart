#####################
 Iterate CLI Options
#####################

This page documents the CLI options available for the ``iterate`` command.

*************************
 Basic Command Structure
*************************

.. code:: bash

   chemsmart run iterate -f <CONFIG_FILE> [OPTIONS]

********************************
 Role of the Configuration File
********************************

The configuration file (specified by ``-f``) acts as a **manifest** or **catalog** rather than a standard structure
input. Instead of containing raw molecular data, it lists the **file paths** for skeleton and substituent molecules,
along with the **link indices** that define their connection points. This structure allows ``iterate`` to systematically
combine multiple building blocks defined in external files.

You can generate a template configuration file using the ``-g`` flag:

.. code:: bash

   chemsmart run iterate -g my_template.toml

For a detailed explanation of each section within the configuration file, please refer to the :doc:`Structure Generation
<iterate-structure-generation>` page.

*************
 CLI Options
*************

Job Control Options
===================

.. list-table::
   :header-rows: 1
   :widths: 35 15 50

   -  -  Option
      -  Type
      -  Description

   -  -  ``-f, --filename``
      -  string
      -  Path to the TOML configuration file (required)

   -  -  ``-np, --nprocs``
      -  int
      -  Number of parallel processes to use (default: 1)

   -  -  ``-t, --timeout``
      -  int
      -  Timeout in seconds for each molecule to be generated (default: 120)

   -  -  ``-m, --method``
      -  choice
      -  Optimization method for substituent positioning (default: ``lagrange_multipliers``)

   -  -  ``-g, --generate-template``
      -  string
      -  Generate a template configuration file (in TOML format) and exit. Defaults to ``iterate_template.toml`` if not
         specified.

Output Options
==============

.. list-table::
   :header-rows: 1
   :widths: 35 15 50

   -  -  Option
      -  Type
      -  Description

   -  -  ``--separate-outputs`` / ``--no-separate-outputs``
      -  flag
      -  Whether to save each structure as a separate XYZ file or merge into one file (default:
         ``--no-separate-outputs``)

   -  -  ``-o, --outputfile``
      -  string
      -  Output filename for merged output (default: ``iterate_out``). Only valid with ``--no-separate-outputs``.

   -  -  ``-d, --directory``
      -  string
      -  Directory to save separate output files. Only valid with ``--separate-outputs``.

Algorithm Parameters
====================

.. list-table::
   :header-rows: 1
   :widths: 35 15 50

   -  -  Option
      -  Type
      -  Description

   -  -  ``-s, --sphere-direction-samples-number``
      -  int
      -  Number of points to sample on the unit sphere for orientation optimization (default: 96)

   -  -  ``-a, --axial-rotations-sample-number``
      -  int
      -  Number of axial rotations to sample per sphere point (default: 6)

.. note::

   Increasing the sampling numbers for ``-s`` and ``-a`` will increase the time required for molecule generation.
