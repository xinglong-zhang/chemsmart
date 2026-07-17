.. _gaussian-pka-calculations:

###########################
 Gaussian pKa Calculations
###########################

This page covers **Gaussian pKa job submission** — generating input files and running the dual-level HA / A⁻ /
(optional) HRef / Ref⁻ workflow.

.. note::

   Output-file analysis is backend-independent. After calculations finish, use ``chemsmart run pka analyze`` or
   ``chemsmart run pka batch-analyze``. See :ref:`pka-calculations` for the full analysis workflow, table formats, and
   thermochemistry options.

.. contents:: Table of Contents
   :local:
   :depth: 2

********
 Theory
********

For thermodynamic cycles, the dual-level approach, and analysis methodology, see :ref:`pka-calculations`.

*************
 Quick Start
*************

**Proton exchange (default)**

Unless ``-s direct`` is given, CHEMSMART uses the **proton exchange** scheme and expects a reference acid. Provide
``-r``, ``-rpi``, ``-rc``, and ``-rm`` (or a coloured-proton CDXML reference with auto-detected ``-rpi``).

.. code:: bash

   chemsmart run gaussian -p my_project -f acid.xyz -c 0 -m 1 pka \
       -pi 10 \
       -r ref_acid.xyz \
       -rpi 21 \
       -rc 1 \
       -rm 1

Where:

-  ``-p my_project``: Project settings (functional, basis set, etc.)
-  ``-f acid.xyz``: Input geometry (XYZ, LOG, COM, CDXML, …)
-  ``-c 0 -m 1``: Charge and multiplicity of the protonated acid (HA)
-  ``-pi 10``: 1-based index of the proton to remove
-  ``-r`` / ``-rpi`` / ``-rc`` / ``-rm``: Reference acid HRef and its deprotonation settings

This runs gas-phase opt+freq and solvent single-points for HA, A⁻, HRef, and Ref⁻ (default solvent: SMD/water from
project or CLI).

**Direct cycle**

Use ``-s direct`` when you do **not** want a reference acid. Only HA and A⁻ calculations are submitted.

.. code:: bash

   chemsmart run gaussian -p my_project -f acid.xyz -c 0 -m 1 pka \
       -pi 10 \
       -s direct

**ChemDraw CDXML / CDX input**

Structures drawn in ChemDraw can be submitted directly. Mark the acidic proton with a **distinct atom colour**;
CHEMSMART reads the drawing and auto-detects it, so ``-pi`` is optional for single-fragment files.

.. code:: bash

   # Single fragment — proton index auto-detected from colour
   chemsmart run gaussian -p my_project -f phenol.cdxml -c 0 -m 1 pka \
       -r ref_acid.xyz -rpi 21 -rc 1 -rm 1

   # Reference acid can also be CDXML (omit -rpi when uniquely coloured)
   chemsmart run gaussian -p my_project -f phenol.cdxml -c 0 -m 1 pka \
       -r ref_acid.cdxml -rc 1 -rm 1

See :ref:`pka-calculations` for multi-molecule CDXML batch submission and the ``-cc`` / ``-rcc`` colour options.

************************
 Job Output File Naming
************************

Sub-job labels determine output filenames. For a job with label ``acid1`` (the default when submitting ``acid1.xyz``):

.. code:: text

   acid1_HA_opt.log      # HA gas-phase opt+freq
   acid1_A_opt.log       # A- gas-phase opt+freq
   acid1_HA_sp.log       # HA solvent single-point
   acid1_A_sp.log        # A- solvent single-point
   acid1_HRef_opt.log    # HRef gas-phase (proton exchange)
   acid1_Ref_opt.log     # Ref- gas-phase (proton exchange)
   acid1_HRef_sp.log     # HRef solvent SP (proton exchange)
   acid1_Ref_sp.log      # Ref- solvent SP (proton exchange)

When building a ``batch-analyze`` output table, either list these paths explicitly or use a ``basename`` and suffix
convention documented in :ref:`pka-calculations` (the ``_pka_*`` autodiscovery pattern matches ORCA-labelled outputs;
Gaussian-labelled outputs are typically given as explicit paths).

************************************
 Batch Processing with Input Tables
************************************

Pass a ``.csv`` or whitespace-delimited ``.txt`` file via ``-f`` and invoke ``pka batch`` (or omit the subcommand —
batch is selected automatically when ``-f`` points to a submission table).

Proton exchange (default — reference acid required on the ``pka`` group):

.. code:: bash

   chemsmart run gaussian -p my_project -f pka_input_table.csv pka \
       -r ref_acid.xyz -rpi 5 -rc 0 -rm 1 batch

Direct cycle (``-s direct`` must be set explicitly):

.. code:: bash

   chemsmart run gaussian -p my_project -f pka_input_table.csv pka -s direct batch

.. note::

   When ``-f`` is a submission table, the parent ``gaussian`` command does not require ``-c`` / ``-m``; charge and
   multiplicity are read from each table row.

.. note::

   In batch mode, each table row is an independent pKa job. Locally (``chemsmart run``), rows run one after another. On
   a cluster, use ``chemsmart sub`` so each row can run as its own array task. See :ref:`pka-hpc-batch-submission`.

Table Format
============

Required columns (comma or whitespace delimited):

-  ``filepath``: Path to the input geometry for HA
-  ``proton_index``: 1-based index of the proton to remove
-  ``charge``: Charge of HA
-  ``multiplicity``: Multiplicity of HA

Example ``pka_input_table.csv``:

.. code:: text

   filepath,proton_index,charge,multiplicity
   /path/to/acid1.xyz,12,0,1
   /path/to/acid2.xyz,8,-1,2
   /path/to/acid3.cdxml,,0,1

.. note::

   For single-molecule ``.cdxml`` / ``.cdx`` table rows, leave ``proton_index`` blank to auto-detect the coloured
   proton. Multi-fragment CDXML expansion applies only when the CDXML file is passed directly as ``-f`` (see **ChemDraw
   CDXML / CDX batch input** below).

.. note::

   For proton exchange with batch input (the default scheme), ``-r``, ``-rpi``, ``-rc``, and ``-rm`` are required on the
   ``pka`` group. When the table has multiple rows, only the **first** row uses the reference acid; later rows switch to
   the direct cycle automatically.

ChemDraw CDXML / CDX batch input
================================

When ``-f`` is a ``.cdxml`` or ``.cdx`` file (not a CSV table), ``pka batch`` reads **every ChemDraw fragment** in the
file, auto-detects the coloured proton in each fragment independently, and submits one pKa job per molecule. Labels
follow ``<basename>_frag<N>_pka`` (e.g. ``acids_frag1_pka_HA_opt.log``).

.. code:: bash

   # Multi-molecule CDXML — proton exchange (default)
   chemsmart run gaussian -p my_project -f acids.cdxml -c 0 -m 1 pka \
       -r ref_acid.xyz -rpi 21 -rc 1 -rm 1 batch

   # Multi-molecule CDXML — direct cycle
   chemsmart run gaussian -p my_project -f acids.cdxml -c 0 -m 1 pka -s direct batch

   # Optional: pin the colour table index when auto-detection is ambiguous
   chemsmart run gaussian -p my_project -f acids.cdxml -c 0 -m 1 pka -cc 4 -s direct batch

A CSV table may list single-molecule ``.cdxml`` paths per row (see Table Format above); blank ``proton_index`` triggers
coloured-proton auto-detection for that row. ``charge`` and ``multiplicity`` still come from the table columns for CSV
rows; for multi-fragment CDXML passed directly as ``-f``, see :ref:`pka-calculations` (Charge and multiplicity).

On clusters, submit the CDXML batch with ``chemsmart sub ... pka batch`` (one array task per fragment). See
:ref:`pka-hpc-batch-submission`.

******************************************
 Computing pKa from Existing Output Files
******************************************

Use the backend-independent analysis command (not ``gaussian pka``):

.. code:: bash

   chemsmart run pka analyze \
       -ha acid1_HA_opt.log \
       -hr ref_acid_HRef_opt.log \
       -rp 6.75 \
       -T 298.15 -c 1.0 -csg 100 -ch 100

See :ref:`pka-calculations` for auto-discovery rules, direct-cycle syntax, and ``batch-analyze``.

************
 Parameters
************

Core Options
============

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Short
      -  Long
      -  Description

   -  -  ``-pi``
      -  ``--proton-index``
      -  **Required** in single-molecule mode (unless CDXML auto-detection applies). 1-based proton index.

   -  -  ``-cc``
      -  ``--color-code``
      -  CDXML colour-table index for the proton to remove. Auto-detected when uniquely coloured.

   -  -  ``-s``
      -  ``--scheme``
      -  ``direct`` or ``proton exchange`` (default).

   -  -
      -  ``--conjugate-base-charge``
      -  Charge of A⁻. Defaults to ``charge - 1``.

   -  -
      -  ``--conjugate-base-multiplicity``
      -  Multiplicity of A⁻. Defaults to HA multiplicity.

Reference Acid Options (Proton Exchange)
========================================

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Short
      -  Long
      -  Description

   -  -  ``-r``
      -  ``--reference``
      -  Geometry file for HRef.

   -  -  ``-rpi``
      -  ``--reference-proton-index``
      -  1-based proton index on HRef. Required with ``-r`` (unless CDXML auto-detection applies).

   -  -  ``-rcc``
      -  ``--reference-color-code``
      -  CDXML colour index for the reference proton.

   -  -  ``-rc``
      -  ``--reference-charge``
      -  Charge of HRef. Required with ``-r``.

   -  -  ``-rm``
      -  ``--reference-multiplicity``
      -  Multiplicity of HRef. Required with ``-r``.

   -  -
      -  ``--reference-conjugate-base-charge``
      -  Charge of Ref⁻. Defaults to ``reference_charge - 1``.

   -  -
      -  ``--reference-conjugate-base-multiplicity``
      -  Multiplicity of Ref⁻. Defaults to reference multiplicity.

   -  -  ``-rp``
      -  ``--reference-pka``
      -  Experimental pKa of HRef. Required for ``chemsmart run pka analyze`` (not for job submission).

Solvent Options
===============

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Short
      -  Long
      -  Description

   -  -  ``-sm``
      -  ``--solvent-model``
      -  Solvation model for solvent SP. Default: ``SMD`` (or project setting).

   -  -  ``-si``
      -  ``--solvent-id``
      -  Solvent identifier. Default: ``water``.

Thermochemistry Options
=======================

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Short
      -  Long
      -  Description

   -  -  ``-T``
      -  ``--temperature``
      -  Temperature in Kelvin. Default: ``298.15`` K.

   -  -  ``-c``
      -  ``--concentration``
      -  Concentration in mol/L. Default: ``1.0`` mol/L.

   -  -  ``-P``
      -  ``--pressure``
      -  Pressure in atm. Default: ``1.0`` atm.

   -  -  ``-csg``
      -  ``--cutoff-entropy-grimme``
      -  Grimme quasi-RRHO entropy cutoff (cm⁻¹). Default: ``100.0``.

   -  -  ``-cst``
      -  ``--cutoff-entropy-truhlar``
      -  Truhlar quasi-RRHO entropy cutoff (cm⁻¹). Mutually exclusive with ``-csg``.

   -  -  ``-ch``
      -  ``--cutoff-enthalpy``
      -  Head-Gordon enthalpy cutoff (cm⁻¹). Default: ``100.0``.

   -  -  ``-dG``
      -  ``--delta-g-proton``
      -  :math:`\Delta G^{\circ}(\text{H}^{+})_{\text{aq}}` in kcal/mol for direct-cycle **submission**. Default:
         ``-265.9``. For output analysis, pass explicitly with ``-s direct`` (see :ref:`pka-calculations`).

**********
 Examples
**********

Example 1: Direct Cycle Submission
==================================

.. code:: bash

   chemsmart run gaussian -p b3lyp_project -f phenol.xyz -c 0 -m 1 pka \
       -pi 13 \
       -s direct \
       -T 298.15

Example 2: Proton Exchange with Custom Temperature
==================================================

.. code:: bash

   chemsmart run gaussian -p m062x_project -f benzoic_acid.xyz -c 0 -m 1 pka \
       -pi 15 \
       -r acetic_acid.xyz \
       -rpi 10 \
       -rc 0 \
       -rm 1 \
       -T 310.15 \
       -sm SMD \
       -si water

Example 3: Batch Submission from CSV
====================================

.. code:: bash

   chemsmart run gaussian -p m062x_project -f pka_scale.csv pka \
       -r ref_acid.xyz \
       -rpi 21 \
       -rc 1 \
       -rm 1 \
       batch

Example 4: Analyze Completed Gaussian Outputs
=============================================

.. code:: bash

   chemsmart run pka analyze \
       -ha phenol_HA_opt.log \
       -hr ref_acid_HRef_opt.log \
       -rp 6.75 \
       -T 298.15 -c 1.0 -csg 100 -ch 100

Example 5: Direct-Cycle Analysis
================================

.. code:: bash

   chemsmart run pka -s direct -dG -265.9 analyze \
       -ha phenol_HA_opt.log \
       -T 298.15 -csg 100 -ch 100

**********
 See Also
**********

-  :ref:`pka-calculations`
-  :ref:`orca-pka-calculations`
-  :doc:`thermochemistry-analysis`
