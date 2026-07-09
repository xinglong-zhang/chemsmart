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
      -  Render style: ``pymol`` or ``cylview``; ``visualize`` also accepts ``glossy``, ``comic``, or ``hybrid`` (see
         below)

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

Batch visualization from a directory (by file type):

.. code:: bash

   chemsmart run mol -d /path/to/outputs -t log visualize

This creates a single PyMOL session visualizing the last structure of every ``.log`` file in the directory.

Batch visualization from a directory (by program):

.. code:: bash

   chemsmart run mol -d /path/to/outputs -p gaussian visualize

This creates a single PyMOL session visualizing the last structure of all Gaussian output files in the directory.

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

.. _glossy-visualization:

**********************
 Glossy Visualization
**********************

Create semi-metallic publication-style figures for metal complexes using ``-s glossy`` on the ``visualize`` subcommand.
Like hybrid mode, glossy mode selects a dedicated job type and accepts additional options that are only valid when the
glossy style is active.

.. code:: bash

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] visualize -s glossy [SUBCMD_OPTIONS]

Glossy Options
==============

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   -  -  Option
      -  Type
      -  Description

   -  -  ``-s glossy``
      -  string
      -  Enable glossy semi-metallic visualization mode

   -  -  ``--style-background``
      -  string
      -  Background: ``white`` (publication) or ``dark`` (slides). Default: ``white``

.. note::

   ``--style-background`` can only be used with ``-s glossy``. Hybrid mode (``-H/--hybrid``) and glossy mode cannot be
   combined.

Basic Usage
===========

Publication figure (white background):

.. code:: bash

   chemsmart run mol -f mn_complex.xyz visualize -s glossy --style-background white

Slide / presentation (dark background):

.. code:: bash

   chemsmart run mol -f mn_complex.xyz visualize -s glossy --style-background dark

Manual PyMOL usage
==================

When editing a saved ``.pse`` session:

.. code:: text

   run glossy_metal_style.py
   metallic_poster_render all
   metallic_poster_render all, elem Mn, None, 2.6, N+O+S+P+H, dark

For Mn coordination complexes and fac/mer comparisons, pass additional arguments to ``metallic_poster_render`` directly
in PyMOL. See ``glossy_metal_style.py`` for the full command signature.

.. _comic-visualization:

*********************
 Comic Visualization
*********************

Create flat, illustrated comic figures for metal complexes using ``-s comic`` (alias ``-s hybrid``) on the ``visualize``
subcommand. Like glossy and group-hybrid modes, this selects a dedicated job type.

.. code:: bash

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] visualize -s comic [SUBCMD_OPTIONS]

The style renders thick sticks and scaled spheres with black ray-traced outlines, centered element labels on the metal
and donor atoms, and an orthoscopic (flat) camera view.

Basic Usage
===========

White background (default in the PyMOL script):

.. code:: bash

   chemsmart run mol -f mn_complex.xyz visualize -s comic --style-background white

Dark background:

.. code:: bash

   chemsmart run mol -f mn_complex.xyz visualize -s comic --style-background dark

Manual PyMOL usage
==================

.. code:: text

   run comic_style.py
   render_comic_metallic_labeled_final all
   render_comic_metallic_labeled_final all, white
   ray 1200, 1200

.. note::

   ``--style-background``, ``-H/--hybrid``, and special ``-s`` styles are mutually exclusive where noted. Only one
   visualization mode can be active.

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

   chemsmart run mol -d . -t xyz -l xyz_alignment align

.. note::

   When using ``-i n``, ensure every input file contains the nth structure. Here ``-t/--filetype`` is extension-based:
   it filters filenames such as ``.xyz`` or ``.log`` and does not infer a Gaussian or ORCA program identity.

Align multiple structures in one file:

.. code:: bash

   chemsmart run mol -f conformers.xyz -i 1,3-6,-1 align

.. note::

   If there is no additional index, align all structures in the file by default.
