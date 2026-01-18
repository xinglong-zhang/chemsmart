#############################
 Basic Visualization (PyMOL)
#############################

This page covers molecular visualization capabilities using PyMOL for creating high-quality graphics and interactive
session files.

********************
 Visualization Jobs
********************

Create static PyMOL visualizations and interactive PSE session files.

.. code:: bash

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] visualize [SUBCMD_OPTIONS]

Visualization Options
=====================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-f, --file``
      -  string
      -  PyMOL style script (default: zhang_group_pymol_style.py)

   -  -  ``-s, --style``
      -  string
      -  Render style: pymol or cylview

   -  -  ``-t, --trace/--no-trace``
      -  bool
      -  Ray trace rendering (default: enabled)

   -  -  ``-v, --vdw``
      -  bool
      -  Add Van der Waals surface (default: disabled)

   -  -  ``-q, --quiet/--no-quiet``
      -  bool
      -  Quiet mode (default: disabled)

   -  -  ``--command-line-only/--no-command-line-only``
      -  bool
      -  Run without GUI (default: enabled)

   -  -  ``-c, --coordinates``
      -  string
      -  Coordinates for labeling (1-indexed)

Basic Usage
===========

Standard visualization:

.. code:: bash

   chemsmart run mol -f molecule.xyz visualize

Quiet mode:

.. code:: bash

   chemsmart run mol -f calculation.log visualize -q

With coordinate labeling:

.. code:: bash

   chemsmart run mol -f structure.xyz visualize -c [[1,2,3]]

Custom style:

.. code:: bash

   chemsmart run mol -f molecule.log visualize -f custom_style.py

************
 Movie Jobs
************

Generate rotating movie animations.

.. code:: bash

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] movie [SUBCMD_OPTIONS]

Movie jobs inherit all visualization options.

Basic Usage
===========

.. code:: bash

   chemsmart run mol -f molecule.xyz movie

**********************
 Hybrid Visualization
**********************

Create visualizations with different groups displayed in different styles.

.. code:: bash

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] visualize --hybrid [SUBCMD_OPTIONS]

Hybrid Options
==============

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   -  -  Option
      -  Type
      -  Description

   -  -  ``-H, --hybrid``
      -  bool
      -  Enable hybrid visualization mode

   -  -  ``-G, --group``
      -  string
      -  Atom indices for a group (repeatable)

   -  -  ``-C, --color``
      -  string
      -  Color for each group (repeatable)

   -  -  ``-SC, --surface-color``
      -  string
      -  Surface color (default: grey)

   -  -  ``-ST, --surface-transparency``
      -  string
      -  Surface transparency (default: 0.7)

   -  -  ``-NC, --new-color-carbon``
      -  string
      -  Carbon atom color (RGB list)

   -  -  ``-NN, --new-color-nitrogen``
      -  string
      -  Nitrogen atom color (RGB list)

   -  -  ``-NO, --new-color-oxygen``
      -  string
      -  Oxygen atom color (RGB list)

Basic Usage
===========

Basic hybrid visualization:

.. code:: bash

   chemsmart run mol -f molecule.xyz visualize --hybrid -G '1,2,3'

Custom colors:

.. code:: bash

   chemsmart run mol -f molecule.log visualize --hybrid -G '1,2,3' -G '4,5,6' -C cbay -C cbak

Custom background settings:

.. code:: bash

   chemsmart run mol -f structure.xyz visualize --hybrid -G '1,2,3' -ST 0.8 -NC '[0.8, 0.8, 0.9]'

Example
-------

.. code:: bash

   chemsmart run mol -f molecule.xyz visualize --hybrid -G '417-418,422-424' -G '336,397-412'

.. image:: _static/B_in_R.png
   :width: 60%
   :align: center

************
 Align Jobs
************

Align multiple molecular structures for comparison (Alignment reference is the first structure).

.. code:: bash

   chemsmart run [OPTIONS] mol align [SUBCMD_OPTIONS]

Basic Usage
===========

Align multiple files:

.. code:: bash

   chemsmart run mol -f mol1.xyz -f mol2.gjf -f mol3.log -i 1 align

Align all files of the same type:

.. code:: bash

   chemsmart run mol -t xyz -l xyz_alignment align

.. note::

   When using ``-i n``, ensure every input file contains the nth structure.

Align multiple structures in one file:

.. code:: bash

   chemsmart run mol -f conformers.xyz -i 1,3-6,-1 align

.. note::

   If there is no additional index, align all structures in the file by default.
