#######################################
 Electronic Structure Analysis (PyMOL)
#######################################

This page covers electronic structure visualization using PyMOL, including molecular orbitals and spin density plots.

.. Warning::

   MO job and Spin job require both ``.log`` and ``.chk`` files in the same folder. If ``-l/--label`` is provided,
   CHEMSMART still processes files from the source filename basename (for example, ``output.log`` ->
   ``output.chk``/``output.fchk``), while the final spin output/session filename follows the custom label.

*****************************
 Molecular Orbital (MO) Jobs
*****************************

Generate molecular orbital visualizations for frontier orbitals and other electronic states.

.. code:: bash

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] mo [SUBCMD_OPTIONS]

.. note::

   For MO jobs, user must provide one of the following [SUBCMD_OPTIONS]: ``-n, --number``, ``-h, --homo`` or ``-l,
   --lumo``.

MO Options
==========

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-n, --number``
      -  int
      -  Specific MO number to visualize

   -  -  ``-h, --homo``
      -  bool
      -  Plot HOMO (default: disabled)

   -  -  ``-l, --lumo``
      -  bool
      -  Plot LUMO (default: disabled)

   -  -  ``-sw, --swap``
      -  bool
      -  Swap positive and negative orbital phase colors (default: disabled)

   -  -  ``-cp, --color-positive``
      -  str
      -  Color for the positive orbital phase isosurface, e.g. ``-cp '[0,1,0]'`` or ``-cp blue`` (default: blue)

   -  -  ``-cn, --color-negative``
      -  str
      -  Color for the negative orbital phase isosurface, e.g. ``-cn '[1,0,0]'`` or ``-cn red`` (default: red)

   -  -  ``-i, --isosurface-value``
      -  float
      -  Set isosurface value to be used in PyMOL .pml file (default: 0.05).

   -  -  ``-tv, --transparency-value``
      -  float
      -  Set transparency value to be used in PyMOL .pml file. Value range: 0.0 - 1.0; 0.0 = fully opaque; 1.0 = fully
         transparent (default: 0.2)

   -  -  ``-sq, --surface-quality``
      -  int
      -  Set surface quality in PyMOL .pml file. Controls the quality of molecular surfaces. value range: 0 (Low
         quality) - 4 (Ultra quality) (default: 3)

   -  -  ``-a, --antialias-value``
      -  int
      -  Set antialias value in PyMOL .pml file. Controls smoothing of edges. value range: 0 (Off, jagged edges) - 4
         (Ultra quality anti-aliasing) (default: 3)

   -  -  ``-m, --ray-trace-mode``

      -  int

      -  Set ray trace mode in PyMOL .pml file. Controls quality of ray-traced images. value range: 0 (standard
         photorealistic render), 1 (outlines around objects, like cell-shading), 2 (no shading, wireframe-like
         appearance), 3 (for figures on dark backgrounds) (default: 1)

.. note::

   MO jobs inherit all visualization options including styling, ray tracing, and surface rendering. Users can further
   modify the *.pml file* after the *.pse file* and *.pml file* have been generated and then reapply the updated
   settings to the PyMOL session.

   Orbital phase colors are for visualization only. Because an orbital's overall sign is arbitrary, swapping both phase
   colors with ``-sw/--swap`` is scientifically valid. When comparing separate orbitals with swapped colors, describe
   them as phase-aligned for visualization rather than implying that a specific color has an absolute physical phase.

Basic Usage
===========

HOMO visualization:

.. code:: bash

   chemsmart run mol -f molecule.log mo -h

LUMO visualization:

.. code:: bash

   chemsmart run mol -f molecule.log mo -l

Specific orbital:

.. code:: bash

   chemsmart run mol -f molecule.log mo -n 5 -m 2

Invert phase colors:

.. code:: bash

   chemsmart run mol -f molecule.log mo -h -sw

Custom phase colors:

.. code:: bash

   chemsmart run mol -f molecule.log mo -h -cp '[0,1,0]' -cn '[1,0,0]'

*******************
 Spin Density Jobs
*******************

Generate spin density visualizations for open-shell systems.

.. code:: bash

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] spin [SUBCMD_OPTIONS]

Spin Options
============

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-i, --isosurface-value``
      -  float
      -  Set isosurface value to be used in PyMOL .pml file (default: 0.05).

   -  -  ``-tv, --transparency-value``
      -  float
      -  Set transparency value to be used in PyMOL .pml file. Value range: 0.0 - 1.0; 0.0 = fully opaque; 1.0 = fully
         transparent (default: 0.2)

   -  -  ``-sq, --surface-quality``
      -  int
      -  Set surface quality in PyMOL .pml file. Controls the quality of molecular surfaces. value range: 0 (Low
         quality) - 4 (Ultra quality) (default: 3)

   -  -  ``-a, --antialias-value``
      -  int
      -  Set antialias value in PyMOL .pml file. Controls smoothing of edges. value range: 0 (Off, jagged edges) - 4
         (Ultra quality anti-aliasing) (default: 3)

   -  -  ``-m, --ray-trace-mode``

      -  int

      -  Set ray trace mode in PyMOL .pml file. Controls quality of ray-traced images. value range: 0 (standard
         photorealistic render), 1 (outlines around objects, like cell-shading), 2 (no shading, wireframe-like
         appearance), 3 (for figures on dark backgrounds) (default: 1)

.. note::

   Spin jobs inherit all visualization options including styling, ray tracing, and surface rendering. Users can further
   modify the *.pml file* after the *.pse file* and *.pml file* have been generated and then reapply the updated
   settings to the PyMOL session.

Basic Usage
===========

Standard spin density:

.. code:: bash

   chemsmart run mol -f radical.log spin

With custom output label while still processing source files:

.. code:: bash

   chemsmart run mol -f output.log -l new_name_new_spin_isovalue spin -i 0.1

This command processes ``output.log`` (and related ``output.chk``/``output.fchk``) and writes the spin session as
``new_name_new_spin_isovalue.pse``.
