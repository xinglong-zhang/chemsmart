###############################
 Interaction Analysis (PyMOL)
###############################

This page covers non-covalent interaction visualization using PyMOL.

**********
 NCI Jobs
**********

Generate non-covalent interaction (NCI) plots for analyzing intermolecular
and intramolecular interactions.

.. code-block:: bash

   chemsmart run [OPTIONS] mol [MOL_OPTIONS] nci [SUBCMD_OPTIONS]

NCI Options
===========

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-i, --isosurface``
     - float
     - Isosurface value (default: 0.5)
   * - ``-r, --color-range``
     - float
     - Color ramp range (default: 1.0)
   * - ``-b, --binary``
     - bool
     - Two-color plot (default: disabled)
   * - ``--intermediate``
     - bool
     - Intermediate range colors (default: disabled)

.. note::

   NCI jobs inherit all visualization options.

Basic Usage
===========

Standard NCI plot:

.. code-block:: bash

   chemsmart run mol -f complex.xyz nci

Custom color range:

.. code-block:: bash

   chemsmart run mol -f complex.xyz nci -r 1.5

Binary color plot:

.. code-block:: bash

   chemsmart run mol -f hydrogen_bond.xyz nci -b

High-quality with ray tracing:

.. code-block:: bash

   chemsmart run mol -f complex.xyz nci -i 0.4 -t -s pymol
