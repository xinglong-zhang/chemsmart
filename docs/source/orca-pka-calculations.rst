.. _orca-pka-calculations:

#######################
 ORCA pKa Calculations
#######################

This section describes how to run pKa calculations using ORCA. The command structure is designed to be very similar to
the Gaussian pKa calculations for a consistent user experience.

.. contents:: Table of Contents
   :local:
   :depth: 2

*************
 Quick Start
*************

To run a pKa calculation for an acid molecule from an XYZ file with ORCA:

.. code:: bash

   chemsmart run orca -p my_project -f acid.xyz -c 0 -m 1 pka -pi 10

Where:

-  ``-p my_project``: Project settings file (defines functional, basis set, etc. for ORCA)
-  ``-f acid.xyz``: Input geometry file (XYZ, LOG, or other formats)
-  ``-c 0 -m 1``: Charge and multiplicity of the protonated acid (HA)
-  ``-pi 10``: 1-based index of the proton to remove

Proton Exchange Cycle with Reference Acid
=========================================

.. code:: bash

   chemsmart run orca -p my_project -f acid.xyz -c 0 -m 1 pka \
       -pi 10 \
       -t "proton exchange" \
       -r water.xyz \
       -rpi 1 \
       -rc 0 \
       -rm 1 \
       -rp 14.0

This will run the full set of calculations for HA, A⁻, HB, and B⁻ using ORCA.

************************************
 Batch Processing with Input Tables
************************************

For high-throughput pKa calculations, you can use a table-driven approach to submit multiple jobs at once.

Using an Input Table
====================

Provide the ``--input-table`` flag and pass a ``.csv`` or ``.txt`` file to the main ``-f`` option.

.. code:: bash

   chemsmart run orca -p my_project -f pka_input_table.csv -i

The input table should contain the following columns:

-  ``filepath``: Path to the input geometry file for the acid (HA).
-  ``proton_index``: 1-based index of the proton to remove.
-  ``charge``: Charge of the HA molecule.
-  ``multiplicity``: Multiplicity of the HA molecule.

Example ``pka_input_table.csv``:

.. code:: text

   filepath,proton_index,charge,multiplicity
   /path/to/acid1.xyz,12,0,1
   /path/to/acid2.xyz,8,-1,2
   /path/to/acid3.cdxml,,0,1

.. note::

   When using a ChemDraw file (``.cdxml``) in the table, you can leave the ``proton_index`` blank to enable automatic
   proton detection based on atom coloring.

************
 Parameters
************

The parameters for ORCA pKa calculations are identical to the Gaussian pKa parameters. Please refer to the main
:ref:`pka-calculations` page for a full list of options. The main difference is that the default solvent model for ORCA
is ``CPCM``.

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Short
      -  Long
      -  Description

   -  -  ``-sm``
      -  ``--solvent-model``
      -  Solvation model for solution phase SP. Default: ``CPCM``.

   -  -  ``-si``
      -  ``--solvent-id``
      -  Solvent identifier. Default: ``water``.
