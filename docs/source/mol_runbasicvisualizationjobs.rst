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

   -  -  ``--file``
      -  string
      -  PyMOL file script or style. If not specified, defaults to use zhang_group_pymol_style.py (default=None)

   -  -  ``-s, --style``
      -  string
      -  PyMOL render style. Options: pymol, cylview (default=None)

   -  -  ``--trace/--no-trace``
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

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] align [SUBCMD_OPTIONS]

Align-Specific OPTIONS
======================

.. list-table:: Visualization Job Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-f, --filenames``
      -  string
      -  Filenames from which new Gaussian inputs are prepared (multiple=True).

   -  -  ``-t, --filetype``
      -  string
      -  Input file pattern, e.g. '.log','.xyz','.gjf' (default=None). Only for align jobs.

Can use ``-t, --filetype``, also inherit all options from visualization jobs and use the same parameters.

Align Basic Usage
=================

**align files**

   .. code:: console

      chemsmart run mol align -f molecule1.xyz -f molecule2.gjf -f molecule3.log

**align all files in same type**

   .. code:: console

      chemsmart run mol align -t .xyz
