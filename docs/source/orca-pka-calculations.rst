.. _orca-pka-calculations:

#######################
 ORCA pKa Calculations
#######################

This section describes how to run pKa calculations using ORCA. The command structure mirrors the Gaussian pKa
calculations for a consistent user experience.

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
-  ``-f acid.xyz``: Input geometry file (XYZ, ORCA output, or other formats)
-  ``-c 0 -m 1``: Charge and multiplicity of the protonated acid (HA)
-  ``-pi 10``: 1-based index of the proton to remove

This will:

#. Optimize HA in gas phase (opt + freq)
#. Optimize A- in gas phase (opt + freq)
#. Run single-point on optimized HA in solution (CPCM/water)
#. Run single-point on optimized A- in solution (CPCM/water)

Proton Exchange Cycle with Reference Acid
=========================================

For more accurate results, use a reference acid with known pKa:

.. code:: bash

   chemsmart run orca -p my_project -f acid.xyz -c 0 -m 1 pka \
       -pi 10 \
       -s "proton exchange" \
       -r reference.xyz \
       -rpi 1 \
       -rc 0 \
       -rm 1 \
       -rp 14.0

This will run the full set of calculations for HA, A-, HRef, and Ref- using ORCA.

************************************
 Batch Processing with Input Tables
************************************

For high-throughput pKa calculations, use a table-driven approach to submit multiple jobs at once.

Using a CSV Input Table
=======================

Pass a ``.csv`` (or whitespace-delimited ``.txt``) file via the ``-f`` option and add the ``pka batch`` subcommand.

.. code:: bash

   chemsmart run orca -p my_project -f pka_input_table.csv pka \
       -s "proton exchange" -r ref.xyz -rpi 5 -rc 0 -rm 1 batch

Or with the direct cycle (no reference acid required):

.. code:: bash

   chemsmart run orca -p my_project -f pka_input_table.csv pka batch

.. note::

   When a ``.csv`` file is detected, ``orca`` defers molecule loading to the ``batch`` subcommand. The ``-c`` / ``-m``
   options on the parent ``orca`` command are **not** required; charge and multiplicity are read from the table rows
   instead.

Table Format
------------

The input table must contain the following four columns (comma or whitespace delimited):

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

For the proton exchange cycle, the reference acid options (``-r``, ``-rpi``, ``-rc``, ``-rm``) are still required on the
``pka`` group even in batch mode.

Batch Execution Policy (Serial vs Non-Serial)
---------------------------------------------

Batch execution policy is controlled by the top-level run/sub flag:

.. code:: bash

   chemsmart run --run-in-serial orca -p my_project -f molecules.csv pka \
       -s "proton exchange" -r ref.xyz -rpi 5 -rc 0 -rm 1 batch

Use ``--run-in-serial`` to force one-by-one execution of table entries. Use the default ``--run-in-parallel`` for
non-serial execution.

Computing pKa from Existing Output Files
========================================

If you already have completed ORCA output files, compute pKa directly using the backend-independent command:

.. code:: bash

   chemsmart run pka analyze \
       -ha acid_opt.out \
       -a conjugate_base_opt.out \
       -hr reference_acid_opt.out \
       -r reference_base_opt.out \
       -has acid_sp_cpcm.out \
       -as conjugate_base_sp_cpcm.out \
       -hrs reference_acid_sp_cpcm.out \
       -rs reference_base_sp_cpcm.out \
       -rp 6.75 \
       -T 298.15

The program automatically detects ORCA output files. See :ref:`pka-calculations` for the full list of ``chemsmart run
pka`` options.

************
 Parameters
************

The parameters for ORCA pKa calculations mirror those for Gaussian. The key difference is the default solvation model
(``CPCM`` instead of ``SMD``).

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Short
      -  Long
      -  Description

   -  -  ``-pi``
      -  ``--proton-index``
      -  **Required** (single-molecule mode). 1-based index of the proton to remove.

   -  -  ``-s``
      -  ``--scheme``
      -  Thermodynamic cycle: ``"direct"`` or ``"proton exchange"`` (default).

   -  -  ``-cc``
      -  ``--conjugate-base-charge``
      -  Charge of A-. Defaults to ``charge - 1``.

   -  -  ``-cm``
      -  ``--conjugate-base-multiplicity``
      -  Multiplicity of A-. Defaults to same as HA.

   -  -  ``-sm``
      -  ``--solvent-model``
      -  Solvation model for solution phase SP. Default: ``CPCM``.

   -  -  ``-si``
      -  ``--solvent-id``
      -  Solvent identifier. Default: ``water``.

   -  -  ``-r``
      -  ``--reference``
      -  Path to geometry file for reference acid (HRef). Required for proton exchange cycle.

   -  -  ``-rpi``
      -  ``--reference-proton-index``
      -  1-based index of proton to remove from HRef.

   -  -  ``-rc``
      -  ``--reference-charge``
      -  Charge of HRef.

   -  -  ``-rm``
      -  ``--reference-multiplicity``
      -  Multiplicity of HRef.

   -  -  ``-rcc``
      -  ``--reference-conjugate-base-charge``
      -  Charge of Ref-. Defaults to ``reference_charge - 1``.

   -  -  ``-rcm``
      -  ``--reference-conjugate-base-multiplicity``
      -  Multiplicity of Ref-. Defaults to ``reference_multiplicity``.

   -  -  ``-rp``
      -  ``--reference-pka``
      -  Experimental pKa of HRef. Required for output file parsing mode.

   -  -  ``-T``
      -  ``--temperature``
      -  Temperature in Kelvin. Default: ``298.15`` K.

   -  -  ``-conc``
      -  ``--concentration``
      -  Concentration in mol/L. Default: ``1.0`` mol/L.

   -  -  ``-csg``
      -  ``--cutoff-entropy-grimme``
      -  Cutoff frequency (cm-1) for quasi-RRHO entropy. Default: ``100.0``.

   -  -  ``-ch``
      -  ``--cutoff-enthalpy``
      -  Cutoff frequency (cm-1) for Head-Gordon enthalpy. Default: ``100.0``.

   -  -  ``-dG``
      -  ``--delta-g-proton``
      -  Absolute free energy of H+ in kcal/mol for direct cycle. Default: ``-265.9``.

**********
 Examples
**********

Example 1: Simple ORCA pKa with Direct Cycle
============================================

.. code:: bash

   chemsmart run orca -p orca_m062x -f phenol.xyz -c 0 -m 1 pka \
       -pi 13 \
       -s direct \
       -T 298.15

Example 2: pKa with Proton Exchange Cycle
=========================================

.. code:: bash

   chemsmart run orca -p orca_m062x -f benzoic_acid.xyz -c 0 -m 1 pka \
       -pi 15 \
       -s "proton exchange" \
       -r acetic_acid.xyz \
       -rpi 10 \
       -rc 0 \
       -rm 1 \
       -T 310.15 \
       -sm CPCM \
       -si water

Example 3: Batch Submission from CSV
====================================

.. code:: bash

   chemsmart run orca -p orca_m062x -f pka_scale.csv pka \
       -s "proton exchange" \
       -r collidine.xyz \
       -rpi 21 \
       -rc 1 \
       -rm 1 \
       -rp 6.75 \
       batch

In this example, ``charge`` and ``multiplicity`` for each molecule are read from the CSV file rather than from the CLI.
The parent ``orca`` command does not require ``-c`` / ``-m``.

Example 4: Extract pKa from Completed ORCA Calculations
=======================================================

.. code:: bash

   chemsmart run pka analyze \
       -ha phenol_opt.out \
       -a phenolate_opt.out \
       -hr collidine_H_opt.out \
       -r collidine_opt.out \
       -has phenol_sp_cpcm.out \
       -as phenolate_sp_cpcm.out \
       -hrs collidine_H_sp_cpcm.out \
       -rs collidine_sp_cpcm.out \
       -rp 6.75 \
       -T 298.15 \
       -csg 100 \
       -ch 100

**********
 See Also
**********

-  :ref:`pka-calculations`
-  :ref:`gaussian-pka-calculations`
-  :doc:`thermochemistry-analysis`
