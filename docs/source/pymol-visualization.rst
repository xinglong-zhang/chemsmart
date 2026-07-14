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
      -  Render style: ``pymol`` or ``cylview``; ``visualize`` also accepts ``glossy``, ``comic``, ``soft-cartoon``, and
         scientific styles (see below)

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
Like other derived styles, glossy uses :class:`~chemsmart.jobs.mol.visualize.PyMOLScientificStyleVisualizationJob`.

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

.. note::

   Group hybrid mode (``-H/--hybrid``) and derived ``-s`` styles cannot be combined. Ray-traced PNG export uses a
   transparent background.

Basic Usage
===========

.. code:: bash

   chemsmart run mol -f complex.xyz visualize -s glossy

Example
=======

.. figure:: _static/pymol_styles/style_glossy.png
   :alt: Glossy semi-metallic visualization example
   :align: center
   :width: 70%

   Glossy semi-metallic style applied to ``complex.xyz`` via
   ``chemsmart run mol -f complex.xyz visualize -s glossy``.

Manual PyMOL usage
==================

When editing a saved ``.pse`` session:

.. code:: text

   run scientific_styles.py
   metallic_poster_render all
   metallic_poster_render all, elem Mn, None, 2.6, N+O+S+P+H, dark

For Mn coordination complexes and fac/mer comparisons, pass additional arguments to ``metallic_poster_render`` directly
in PyMOL. See ``scientific_styles.py`` for the full command signature.

.. _comic-visualization:

*********************
 Comic Visualization
*********************

Create flat, illustrated comic figures for metal complexes using ``-s comic`` on the ``visualize`` subcommand. Like
other derived styles, comic uses :class:`~chemsmart.jobs.mol.visualize.PyMOLScientificStyleVisualizationJob`.

.. code:: bash

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] visualize -s comic [SUBCMD_OPTIONS]

The style renders thick sticks and scaled spheres with black ray-traced outlines, centered element labels on the metal
and donor atoms, and an orthoscopic (flat) camera view. Ray-traced PNG export uses a transparent background for
compositing into manuscripts or slides.

Basic Usage
===========

.. code:: bash

   chemsmart run mol -f complex.xyz visualize -s comic

Highlight specific metal-ligand bonds in coordination-core styling with ``-c`` (distance labels are suppressed for comic
by default):

.. code:: bash

   chemsmart run mol -f complex.xyz visualize -s comic \\
       -c '[[1,2],[1,5],[1,36],[1,3],[1,15],[1,8]]'

Example
=======

.. figure:: _static/pymol_styles/style_comic.png
   :alt: Comic visualization example
   :align: center
   :width: 70%

   Comic style applied to ``complex.xyz`` with Mn coordination-bond highlighting and a transparent background via
   ``chemsmart run mol -f complex.xyz visualize -s comic -c '[[1,2],[1,5],[1,36],[1,3],[1,15],[1,8]]'``.

Manual PyMOL usage
==================

.. code:: text

   run scientific_styles.py
   render_comic_metallic_labeled_final all
   ray 1200, 1200
   png comic.png, dpi=300

.. _soft-cartoon-visualization:

****************************
 Soft Cartoon Visualization
****************************

Create soft premium ball-and-stick figures for metal complexes using ``-s soft-cartoon`` on the ``visualize``
subcommand. Like glossy and comic modes, this selects a dedicated job type.

.. code:: bash

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] visualize -s soft-cartoon [SUBCMD_OPTIONS]

The style renders scaled sticks and spheres with soft metallic shading, light ray shadows, and a perspective camera
view. Ray-traced PNG export uses a transparent background for compositing.

Basic Usage
===========

.. code:: bash

   chemsmart run mol -f complex.xyz visualize -s soft-cartoon

Example
=======

.. figure:: _static/pymol_styles/style_soft_cartoon.png
   :alt: Soft cartoon visualization example
   :align: center
   :width: 70%

   Soft cartoon ball-and-stick style applied to ``complex.xyz`` via
   ``chemsmart run mol -f complex.xyz visualize -s soft-cartoon``.

Manual PyMOL usage
==================

.. code:: text

   run scientific_styles.py
   render_soft_cartoon all
   render_soft_cartoon all, white
   ray 1200, 1200

.. _scientific-visualization:

**************************
 Scientific Visualization
**************************

CHEMSMART applies ``scientific_styles.py`` for ``visualize -s`` choices including ``glossy``, ``comic``,
``soft-cartoon``, ``editorial-minimal``, ``soft-ceramic``, and other scientific styles. These are routed through
:class:`~chemsmart.jobs.mol.visualize.PyMOLScientificStyleVisualizationJob`.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   -  -  ``-s`` value
      -  Description
   -  -  ``glossy``
      -  Semi-metallic poster rendering with transparent PNG export
   -  -  ``comic``
      -  Comic ball-and-stick rendering with black outlines and labels
   -  -  ``soft-cartoon``
      -  Soft premium ball-and-stick rendering with light metallic shading
   -  -  ``editorial-minimal``
      -  Matte minimal ball-and-stick style for main-text mechanistic figures
   -  -  ``soft-ceramic``
      -  Soft ceramic / studio ball-and-stick style for coordination complexes
   -  -  ``neon-coordination-core``
      -  Neon coordination-core emphasis using radius-ratio shells and a transparent background
   -  -  ``matte-clay``
      -  Matte clay ball-and-stick using radius-ratio coordination spheres and core hydrogens
   -  -  ``xray-wire``
      -  Monochrome wireframe style for SI structure verification
   -  -  ``steric-surface``
      -  Transparent molecular surface for steric / space-filling views
   -  -  ``quasi-chemdraw-bold``
      -  Bold ChemDraw-like style using radius-ratio cores and flat illustrative shading
   -  -  ``labeled-coordination-core``
      -  Coordination core with explicit element labels

Basic Usage
===========

.. code:: bash

   chemsmart run mol -f complex.xyz visualize -s editorial-minimal
   chemsmart run mol -f complex.xyz visualize -s soft-ceramic
   chemsmart run mol -f complex.xyz visualize -s labeled-coordination-core

Examples
========

Editorial minimal
-----------------

.. figure:: _static/pymol_styles/style_editorial_minimal.png
   :alt: Editorial minimal visualization example
   :align: center
   :width: 70%

   ``chemsmart run mol -f complex.xyz visualize -s editorial-minimal``

Soft ceramic
------------

.. figure:: _static/pymol_styles/style_soft_ceramic.png
   :alt: Soft ceramic visualization example
   :align: center
   :width: 70%

   ``chemsmart run mol -f complex.xyz visualize -s soft-ceramic``

Neon coordination core
----------------------

.. figure:: _static/pymol_styles/style_neon_coordination_core.png
   :alt: Neon coordination core visualization example
   :align: center
   :width: 70%

   ``chemsmart run mol -f complex.xyz visualize -s neon-coordination-core``

Matte clay
----------

.. figure:: _static/pymol_styles/style_matte_clay.png
   :alt: Matte clay visualization example
   :align: center
   :width: 70%

   ``chemsmart run mol -f complex.xyz visualize -s matte-clay``

X-ray wire
----------

.. figure:: _static/pymol_styles/style_xray_wire.png
   :alt: X-ray wire visualization example
   :align: center
   :width: 70%

   ``chemsmart run mol -f complex.xyz visualize -s xray-wire``

Steric surface
--------------

.. figure:: _static/pymol_styles/style_steric_surface.png
   :alt: Steric surface visualization example
   :align: center
   :width: 70%

   ``chemsmart run mol -f complex.xyz visualize -s steric-surface``

Quasi-ChemDraw bold
-------------------

.. figure:: _static/pymol_styles/style_quasi_chemdraw_bold.png
   :alt: Quasi-ChemDraw bold visualization example
   :align: center
   :width: 70%

   ``chemsmart run mol -f complex.xyz visualize -s quasi-chemdraw-bold``

Labeled coordination core
-------------------------

.. figure:: _static/pymol_styles/style_labeled_coordination_core.png
   :alt: Labeled coordination core visualization example
   :align: center
   :width: 70%

   ``chemsmart run mol -f complex.xyz visualize -s labeled-coordination-core``

Manual PyMOL usage
==================

.. code:: text

   run scientific_styles.py
   render_editorial_minimal all
   render_soft_ceramic all
   render_labeled_coordination_core all
   ray 1200, 1200

.. note::

   ``-H/--hybrid`` and derived ``-s`` styles are mutually exclusive. Only one visualization mode can be active.

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
