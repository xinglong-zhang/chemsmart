###################
 PyMOL CLI Options
###################

This page documents the CLI options for molecular visualization and
analysis using PyMOL. Use ``chemsmart run mol --help`` for the complete
list.

*************************
 Basic Command Structure
*************************

.. code:: bash

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] <SUBCMD> [SUBCMD_OPTIONS]

*************
 MOL Options
*************

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-f, --filenames``
      -  string
      -  Input file(s) for visualization

   -  -  ``-d, --directory``
      -  string
      -  Directory for align jobs

   -  -  ``-t, --filetype``
      -  string
      -  File type filter for directory processing

   -  -  ``-l, --label``
      -  string
      -  Custom output filename

   -  -  ``-a, --append-label``
      -  string
      -  String to append to filename

   -  -  ``-i, --index``
      -  string
      -  Structure index, 1-based (default: -1, last)

   -  -  ``--pubchem``
      -  string
      -  Query structure from PubChem

   -  -  ``-o, --overwrite``
      -  bool
      -  Overwrite existing files (default: disabled)

***********************
 Available Subcommands
***********************

Basic Visualization
===================

.. list-table::
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``visualize``
      -  Create PyMOL visualization and save as PSE file
   -  -  ``movie``
      -  Generate rotating molecule movie

Reaction Analysis
=================

.. list-table::
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``irc``
      -  Generate IRC movie and trajectory visualization

Electronic Structure Analysis
=============================

.. list-table::
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``mo``
      -  Generate molecular orbital visualizations
   -  -  ``spin``
      -  Generate spin density visualizations

Interaction Analysis
====================

.. list-table::
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``nci``
      -  Generate non-covalent interaction plots

************
 Next Steps
************

For detailed information on each visualization type:

-  :doc:`pymol-visualization`
-  :doc:`pymol-reaction-analysis`
-  :doc:`pymol-electronic-structure`
-  :doc:`pymol-interaction-analysis`
