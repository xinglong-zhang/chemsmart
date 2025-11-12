Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to
contact our team if you have questions or feedback.

##############################
 Run Basic Visualization Jobs
##############################

ChemSmart provides powerful molecular visualization capabilities using PyMOL for creating high-quality molecular
graphics, movies, and interactive visualizations.

********************
 Visualization Jobs
********************

Create static PyMOL visualizations and interactive session files.

.. code:: console

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] visualize [SUBCMD_OPTIONS]

Visualization-Specific OPTIONS
==============================

.. list-table:: Visualization Job Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-f, --file``
      -  string
      -  PyMOL file script or style. If not specified, defaults to use zhang_group_pymol_style.py (default=None)

   -  -  ``-s, --style``
      -  string
      -  PyMOL render style. Options: pymol, cylview (default=None)

   -  -  ``-t, --trace/--no-trace``
      -  bool
      -  PyMOL options to ray trace or not (default=True)

   -  -  ``-v, --vdw``
      -  bool
      -  Add Van der Waals surface (default=False)

   -  -  ``-q, --quiet/--no-quiet``
      -  bool
      -  Run PyMOL in quiet mode (default=False)

   -  -  ``--command-line-only/--no-command-line-only``
      -  bool
      -  Run PyMOL in command line only (default=True)

   -  -  ``-c, --coordinates <string>``
      -  string
      -  List of coordinates (bonds, angles and dihedrals) for labelling. 1-indexed (default=None)

Visualization Basic Usage
=========================

**Basic molecular visualization**

   .. code:: console

      chemsmart run mol -f molecule.xyz visualize

**Quiet mode visualization**

   .. code:: console

      chemsmart run mol -f calculation.log visualize -q

**Visualization with coordinate labeling**

   .. code:: console

      chemsmart run mol -f structure.xyz visualize -c [[1,2,3]]

**Visualization using custom style or script file**

   .. code:: console

      chemsmart run mol -f molecule.log visualize -f custom_style.py

Visualization Examples
======================

************
 Movie Jobs
************

Generate rotating movie animations of molecular structures.

.. code:: console

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] movie [SUBCMD_OPTIONS]

Movie-Specific OPTIONS
======================

Movie jobs inherit all options from visualization jobs and use the same parameters.

Movie Basic Usage
=================

**Basic rotating movie**

   .. code:: console

      chemsmart run mol -f molecule.xyz movie

Movie Examples
==============

************
 Align Jobs
************

Align multiple molecular structures for comparison.

.. code:: console

   chemsmart run [OPTIONS] mol align [SUBCMD_OPTIONS]

Align-Specific OPTIONS
======================

Inherit all options from visualization jobs and use the same parameters.

Align Basic Usage
=================

**align files**

   .. code:: console

      chemsmart run mol -f mol1.xyz -f mol2.gjf -f mol3.log -i 1 align

This command will align the subsequent molecules to the first molecule sequentially (align mol2 and mol3 to mol1).

**align all files in same type**

   .. code:: console

      chemsmart run mol -t xyz -l xyz_alignment align

.. note::

   When using the ``-i n`` option, user should ensure that every input file contains the n-th molecule.
