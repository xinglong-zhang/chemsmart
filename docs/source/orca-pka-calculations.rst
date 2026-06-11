.. _orca-pka-calculations:

#######################
 ORCA pKa Calculations
#######################

This page covers **ORCA pKa job submission**. The command structure mirrors Gaussian pKa for a consistent experience.

.. note::

   Output-file analysis is backend-independent. After calculations finish, use ``chemsmart run pka analyze`` or
   ``chemsmart run pka batch-analyze``. See :ref:`pka-calculations` for the full analysis workflow, table formats, and
   thermochemistry options.

.. contents:: Table of Contents
   :local:
   :depth: 2

*************
 Quick Start
*************

**Proton exchange (default)**

Unless ``-s direct`` is given, CHEMSMART uses the **proton exchange** scheme and expects a reference acid.

.. code:: bash

   chemsmart run orca -p my_project -f acid.xyz -c 0 -m 1 pka \
       -pi 10 \
       -r ref_acid.xyz \
       -rpi 21 \
       -rc 1 \
       -rm 1

Where:

-  ``-p my_project``: ORCA project settings
-  ``-f acid.xyz``: Input geometry
-  ``-c 0 -m 1``: Charge and multiplicity of HA
-  ``-pi 10``: 1-based proton index
-  ``-r`` / ``-rpi`` / ``-rc`` / ``-rm``: Reference acid HRef

This runs gas-phase opt+freq and CPCM/water solvent single-points for HA, A⁻, HRef, and Ref⁻.

**Direct cycle**

.. code:: bash

   chemsmart run orca -p my_project -f acid.xyz -c 0 -m 1 pka \
       -pi 10 \
       -s direct

**ChemDraw CDXML / CDX input**

ChemDraw ``.cdxml`` and ``.cdx`` files are supported. Colour the acidic proton in ChemDraw; CHEMSMART reads the drawing
and detects it automatically, so ``-pi`` can be omitted for single-fragment inputs.

.. code:: bash

   chemsmart run orca -p my_project -f phenol.cdxml -c 0 -m 1 pka \
       -r ref_acid.xyz -rpi 21 -rc 1 -rm 1

   chemsmart run orca -p my_project -f phenol.cdxml -c 0 -m 1 pka \
       -r ref_acid.cdxml -rc 1 -rm 1

See :ref:`pka-calculations` for multi-molecule CDXML workflows and colour-code options.

************************
 Job Output File Naming
************************

ORCA batch submission appends ``_pka`` to the input stem when forming the job label. For input ``acid1.xyz`` the label
is ``acid1_pka`` and typical outputs are:

.. code:: text

   acid1_pka_HA_opt.out
   acid1_pka_A_opt.out
   acid1_pka_HA_sp.out
   acid1_pka_A_sp.out
   acid1_pka_HRef_opt.out    # proton exchange
   acid1_pka_Ref_opt.out
   acid1_pka_HRef_sp.out
   acid1_pka_Ref_sp.out

These names align with the ``batch-analyze`` autodiscovery convention ``<basename>_pka_*`` when ``basename`` is
``acid1``.

************************************
 Batch Processing with Input Tables
************************************

Proton exchange (default — reference acid required):

.. code:: bash

   chemsmart run orca -p my_project -f pka_input_table.csv pka \
       -r ref_acid.xyz -rpi 5 -rc 0 -rm 1 batch

Direct cycle:

.. code:: bash

   chemsmart run orca -p my_project -f pka_input_table.csv pka -s direct batch

.. note::

   When ``-f`` is a submission table, the parent ``orca`` command does not require ``-c`` / ``-m``.

.. note::

   For proton exchange with multiple batch rows, only the first row uses the reference acid; subsequent rows switch to
   the direct cycle (same behaviour as Gaussian batch).

Table Format
============

Required columns: ``filepath``, ``proton_index``, ``charge``, ``multiplicity``. See :ref:`gaussian-pka-calculations` for
CSV examples. For single-molecule ``.cdxml`` / ``.cdx`` rows, ``proton_index`` may be left blank for coloured-proton
auto-detection.

ChemDraw CDXML / CDX batch input
================================

Pass a multi-molecule CDXML file directly as ``-f`` with ``pka batch``. Each ChemDraw fragment becomes one pKa job with
per-fragment coloured-proton detection. Labels are ``<basename>_frag<N>_pka`` (outputs such as
``acids_frag1_pka_HA_opt.out``).

.. code:: bash

   chemsmart run orca -p my_project -f acids.cdxml -c 0 -m 1 pka \
       -r ref_acid.xyz -rpi 21 -rc 1 -rm 1 batch

   chemsmart run orca -p my_project -f acids.cdxml -c 0 -m 1 pka -s direct batch

******************************************
 Computing pKa from Existing Output Files
******************************************

.. code:: bash

   chemsmart run pka analyze \
       -ha acid1_pka_HA_opt.out \
       -hr ref_acid_pka_HRef_opt.out \
       -rp 6.75 \
       -T 298.15

ORCA ``.out`` and Gaussian ``.log`` files can be combined in ``batch-analyze``. See :ref:`pka-calculations`.

************
 Parameters
************

ORCA pKa options mirror Gaussian submission options. The main default difference is the solvent model (**CPCM** instead
of **SMD**).

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Short
      -  Long
      -  Description

   -  -  ``-pi``
      -  ``--proton-index``
      -  **Required** in single-molecule mode (unless CDXML auto-detection applies).

   -  -  ``-cc``
      -  ``--color-code``
      -  CDXML colour-table index for the target proton.

   -  -  ``-s``
      -  ``--scheme``
      -  ``direct`` or ``proton exchange`` (default).

   -  -
      -  ``--conjugate-base-charge``
      -  Charge of A⁻. Defaults to ``charge - 1``.

   -  -
      -  ``--conjugate-base-multiplicity``
      -  Multiplicity of A⁻.

   -  -  ``-sm``
      -  ``--solvent-model``
      -  Solvation model for solvent SP. Default: ``CPCM``.

   -  -  ``-si``
      -  ``--solvent-id``
      -  Solvent identifier. Default: ``water``.

   -  -  ``-r``
      -  ``--reference``
      -  Reference acid geometry (proton exchange).

   -  -  ``-rpi``
      -  ``--reference-proton-index``
      -  Proton index on HRef.

   -  -  ``-rcc``
      -  ``--reference-color-code``
      -  CDXML colour index for the reference proton.

   -  -  ``-rc`` / ``-rm``
      -  ``--reference-charge`` / ``--reference-multiplicity``
      -  Required with ``-r`` for proton exchange batch/submit.

   -  -
      -  ``--reference-conjugate-base-charge``
      -  Charge of Ref⁻.

   -  -
      -  ``--reference-conjugate-base-multiplicity``
      -  Multiplicity of Ref⁻.

   -  -  ``-rp``
      -  ``--reference-pka``
      -  Experimental pKa of HRef (required for ``chemsmart run pka analyze`` only).

   -  -  ``-T``
      -  ``--temperature``
      -  Default: ``298.15`` K.

   -  -  ``-c``
      -  ``--concentration``
      -  Default: ``1.0`` mol/L.

   -  -  ``-P``
      -  ``--pressure``
      -  Default: ``1.0`` atm.

   -  -  ``-csg`` / ``-cst`` / ``-ch``
      -  ``--cutoff-entropy-grimme``, etc.
      -  Entropy and enthalpy cutoffs (cm⁻¹). See :ref:`pka-calculations`.

   -  -  ``-dG``
      -  ``--delta-g-proton``
      -  Default ``-265.9`` kcal/mol for submission; explicit for direct-cycle analysis.

**********
 Examples
**********

Example 1: Direct Cycle Submission
==================================

.. code:: bash

   chemsmart run orca -p orca_m062x -f phenol.xyz -c 0 -m 1 pka \
       -pi 13 \
       -s direct \
       -T 298.15

Example 2: Proton Exchange Submission
=====================================

.. code:: bash

   chemsmart run orca -p orca_m062x -f benzoic_acid.xyz -c 0 -m 1 pka \
       -pi 15 \
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
       -r ref_acid.xyz \
       -rpi 21 \
       -rc 1 \
       -rm 1 \
       batch

Example 4: Analyze Completed ORCA Outputs
=========================================

.. code:: bash

   chemsmart run pka analyze \
       -ha phenol_pka_HA_opt.out \
       -hr ref_acid_pka_HRef_opt.out \
       -rp 6.75 \
       -T 298.15 -c 1.0 -csg 100 -ch 100

Example 5: Mixed Gaussian/ORCA Batch Analysis
=============================================

.. code:: bash

   chemsmart run pka batch-analyze -o pka_output.csv

With ``-p auto`` (default), ORCA target ``.out`` files and Gaussian reference ``.log`` files in the same table are
supported. See :ref:`pka-calculations`.

**********
 See Also
**********

-  :ref:`pka-calculations`
-  :ref:`gaussian-pka-calculations`
-  :doc:`thermochemistry-analysis`
