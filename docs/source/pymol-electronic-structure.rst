#######################################
 Electronic Structure Analysis (PyMOL)
#######################################

This page covers electronic structure visualization using PyMOL, including molecular orbitals and spin density plots.

*************
 PML Options
*************

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

   The ``.pml options`` can be added directly to the end of the MO or spin job command. Users can also modify the *.pml
   file* after the *.pse file* and *.pml file* have been generated and then reapply the updated settings to the PyMOL
   session.

Molecular Orbital (MO) Jobs
===========================

Generate molecular orbital visualizations for frontier orbitals and other electronic states.

.. code:: bash

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] mo [SUBCMD_OPTIONS]

************
 MO Options
************

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

.. note::

   MO jobs inherit all visualization options including styling, ray tracing, and surface rendering.

*************
 Basic Usage
*************

Standard MO visualization:

.. code:: bash

   chemsmart run mol -f calculation.log mo

HOMO visualization:

.. code:: bash

   chemsmart run mol -f molecule.log mo -h

LUMO visualization:

.. code:: bash

   chemsmart run mol -f molecule.log mo -l

Specific orbital:

.. code:: bash

   chemsmart run mol -f molecule.log mo -n 5 -m 2

Spin Density Jobs
=================

Generate spin density visualizations for open-shell systems.

.. code:: bash

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] spin

.. note::

   Requires both ``.log`` and ``.chk`` files in the same folder. If ``-l/--label`` is provided, CHEMSMART still
   processes files from the source filename basename (for example, ``output.log`` -> ``output.chk``/``output.fchk``),
   while the final spin output/session filename follows the custom label.

Spin density jobs inherit all visualization options.

*************
 Basic Usage
*************

Standard spin density:

.. code:: bash

   chemsmart run mol -f radical.log spin

With ray tracing:

.. code:: bash

   chemsmart run mol -f radical.log spin -t

With custom output label while still processing source files:

.. code:: bash

   chemsmart run mol -f output.log -l new_name_new_spin_isovalue spin -i 0.1

This command processes ``output.log`` (and related ``output.chk``/``output.fchk``) and writes the spin session as
``new_name_new_spin_isovalue.pse``.
