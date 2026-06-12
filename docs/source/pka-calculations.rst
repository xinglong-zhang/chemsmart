.. _pka-calculations:

##################
 pKa Calculations
##################

CHEMSMART provides pKa workflows in two separate stages:

#. **Job submission** — generate and run Gaussian or ORCA calculations for HA, A⁻, and (optionally) a reference acid.
   See :ref:`gaussian-pka-calculations` and :ref:`orca-pka-calculations`.

#. **Output analysis** — compute pKa values from completed output files using the backend-independent command
   ``chemsmart run pka``. Analysis is **program-agnostic**: the same workflow reads Gaussian ``.log`` and ORCA ``.out``
   files and extracts the energies and thermal corrections needed for the pKa cycle.

.. toctree::
   :maxdepth: 2
   :caption: pKa Calculations

   gaussian-pka-calculations
   orca-pka-calculations

.. contents:: Table of Contents
   :local:
   :depth: 2

************************
 Execution Architecture
************************

**Job submission**

-  ``chemsmart run/sub gaussian ... pka [submit|batch]`` — prepare and run Gaussian pKa calculations.
-  ``chemsmart run/sub orca ... pka [submit|batch]`` — prepare and run ORCA pKa calculations.
-  A single structure yields one job; batch input (CSV table or multi-molecule CDXML) can produce multiple jobs in one
   invocation.
-  When ``pka`` is invoked without an explicit subcommand, a submission table triggers ``batch``; otherwise ``submit``
   runs.

**Output analysis**

-  ``chemsmart run pka analyze`` — single-system analysis from up to eight output files.
-  ``chemsmart run pka batch-analyze`` — table-driven batch analysis.
-  Both commands use the same pKa analysis workflow. Gaussian ``.log`` and ORCA ``.out`` files can be mixed in the same
   batch table; each file is read and interpreted on its own.

Theory
======

The pKa of an acid HA in aqueous solution is defined by the equilibrium:

.. math::

   \text{HA}_{(\text{aq})} \rightleftharpoons \text{A}^{-}_{(\text{aq})} + \text{H}^{+}_{(\text{aq})}

The pKa is related to the standard Gibbs free energy change by:

.. math::

   \text{p}K_{\text{a}} = \frac{\Delta G^{\circ}_{\text{aq}}}{2.303 \cdot R \cdot T}

where :math:`R` is the gas constant and :math:`T` is the temperature.

**********************
 Thermodynamic Cycles
**********************

CHEMSMART supports two thermodynamic cycles for pKa calculations:

**1. Proton Exchange (Isodesmic) Cycle** (Default, Recommended)

This is the default when ``-s`` is omitted for both job submission and output analysis. A reference acid HRef with known
experimental pKa is required to cancel systematic errors:

.. math::

   \text{HA} + \text{Ref}^{-} \rightarrow \text{A}^{-} + \text{HRef}

The pKa is computed as:

.. math::

   \text{p}K_{\text{a}}(\text{HA}) = \text{p}K_{\text{a}}(\text{HRef}) + \frac{\Delta G_{\text{soln}}}{2.303 \cdot R \cdot T}

where:

.. math::

   \Delta G_{\text{soln}} = \left[ G(\text{A}^{-})_{\text{soln}} + G(\text{HRef})_{\text{soln}} \right] - \left[ G(\text{HA})_{\text{soln}} + G(\text{Ref}^{-})_{\text{soln}} \right]

**2. Direct Cycle**

Uses the absolute free energy of a proton in water:

.. math::

   \text{p}K_{\text{a}} = \frac{G(\text{A}^{-})_{\text{aq}} - G(\text{HA})_{\text{aq}} + \Delta G^{\circ}(\text{H}^{+})_{\text{aq}}}{2.303 \cdot R \cdot T}

Default value: :math:`\Delta G^{\circ}(\text{H}^{+})_{\text{aq}} = -265.9` kcal/mol (Tissandier et al., 1998).

*********************
 Dual-Level Approach
*********************

CHEMSMART implements a dual-level approach for accurate solvation free energies:

#. **Thermal corrections** (:math:`G_{\text{corr}}`) from gas-phase frequency calculations using quasi-harmonic Gibbs
   free energy:

   .. math::

      G_{\text{corr}} = G_{\text{qh}}(T) - E_{\text{gas}}

#. **Solvent energies** (:math:`E_{\text{solv}}`) from high-level single-point calculations in implicit solvent (e.g.,
   SMD or CPCM).

#. **Total free energy in solution**:

   .. math::

      G_{\text{soln}} = E_{\text{solv}} + G_{\text{corr}}

.. note::

   All internal energies are stored in Hartree (au). :math:`\Delta G_{\text{soln}}` (proton exchange) and :math:`\Delta
   G_{\text{diss}}` (direct dissociation) are converted to kcal/mol for the pKa formula (1 Hartree = 627.5094740631
   kcal/mol).

**********************************
 Job Submission (Gaussian / ORCA)
**********************************

Job submission is backend-specific. Use the dedicated pages for full examples and parameter tables:

-  :ref:`gaussian-pka-calculations`
-  :ref:`orca-pka-calculations`

**Commands**

The default scheme is **proton exchange**, which requires a reference acid (``-r``, ``-rpi``, ``-rc``, ``-rm``). Use
``-s direct`` only when you want the direct dissociation cycle without a reference acid.

.. code:: bash

   # Proton exchange (default) — reference acid required
   chemsmart run gaussian -p my_project -f acid.xyz -c 0 -m 1 pka \
       -pi 10 -r ref_acid.xyz -rpi 21 -rc 1 -rm 1

   chemsmart run orca -p my_project -f acid.xyz -c 0 -m 1 pka \
       -pi 10 -r ref_acid.xyz -rpi 21 -rc 1 -rm 1

   # Direct cycle — no reference acid; must set -s direct explicitly
   chemsmart run gaussian -p my_project -f acid.xyz -c 0 -m 1 pka -pi 10 -s direct
   chemsmart run orca -p my_project -f acid.xyz -c 0 -m 1 pka -pi 10 -s direct

   # Batch submission (proton exchange requires reference options on the pka group)
   chemsmart run gaussian -p my_project -f pka_input.csv pka \
       -r ref_acid.xyz -rpi 21 -rc 1 -rm 1 batch

   # Batch submission (direct cycle)
   chemsmart run gaussian -p my_project -f pka_input.csv pka -s direct batch

**Submission input table** (``pka batch``)

Comma- or whitespace-delimited table with columns ``filepath``, ``proton_index``, ``charge``, ``multiplicity``.

***************************************
 ChemDraw CDXML / CDX Input (pKa Jobs)
***************************************

pKa job submission can read structures directly from ChemDraw ``.cdxml`` and ``.cdx`` files. CHEMSMART reads atom
colours in the drawing to identify the **acidic proton** to remove. Colour the proton (or the ``H`` in a functional
group such as –OH) in ChemDraw with a distinct colour; CHEMSMART auto-detects it so ``-pi`` is often unnecessary.

**Single-molecule submit**

When ``-f`` points to one CDXML structure with a single fragment, omit ``-pi`` if the coloured proton is unique:

.. code:: bash

   chemsmart run gaussian -p my_project -f phenol.cdxml -c 0 -m 1 pka \
       -r ref_acid.xyz -rpi 21 -rc 1 -rm 1

   chemsmart run gaussian -p my_project -f phenol.cdxml -c 0 -m 1 pka -s direct

Use ``-cc`` / ``--color-code`` when several hydrogens share similar styling and you need to select a specific ChemDraw
colour-table index. The reference acid may also be a CDXML file; in that case ``-rpi`` can be omitted when the reference
proton is uniquely coloured (or use ``-rcc`` / ``--reference-color-code``).

**Multi-molecule CDXML (one job per fragment)**

A single ``.cdxml`` / ``.cdx`` file may contain **multiple molecules** (multiple ChemDraw fragments). CHEMSMART performs
**per-fragment** coloured-proton detection and creates **one pKa job per fragment**.

Pass the file with ``pka batch`` (or ``pka submit`` for a single-fragment file):

.. code:: bash

   chemsmart run gaussian -p my_project -f acids.cdxml -c 0 -m 1 pka \
       -r ref_acid.xyz -rpi 21 -rc 1 -rm 1 batch

   chemsmart run orca -p my_project -f acids.cdxml -c 0 -m 1 pka -s direct batch

Job labels are derived from the filename, e.g. ``acids_frag1_pka`` (Gaussian) or ``acids_frag1_pka`` (ORCA). Charge and
multiplicity come from the parent ``-c`` / ``-m`` options and apply to every fragment unless overridden in a CSV table.

**CDXML paths inside a CSV batch table**

You can mix XYZ and CDXML inputs in the same submission table. Each row becomes one pKa job. For a **single-molecule**
``.cdxml`` / ``.cdx`` row, leave ``proton_index`` blank to auto-detect the coloured proton; an explicit value overrides
the detection (same behaviour as single-file CDXML submit). Non-CDXML rows still require ``proton_index``.

.. code:: text

   filepath,proton_index,charge,multiplicity
   /path/to/acid1.xyz,12,0,1
   /path/to/acid2.cdxml,,0,1

Multi-molecule CDXML files cannot be expanded from a table row. Pass them directly as ``-f`` with ``pka batch`` (see
above) to create one job per ChemDraw fragment.

.. note::

   If ``-f`` is a CDXML file (not a CSV table), CHEMSMART routes to coloured-proton batch expansion automatically. For
   general CDXML structure handling outside pKa, see :doc:`chemdraw-organometallic`.

**Proton and reference options for CDXML**

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Short
      -  Long
      -  Description

   -  -  ``-pi``
      -  ``--proton-index``
      -  Optional for CDXML when a uniquely coloured proton is present. Required for XYZ/LOG/COM inputs.

   -  -  ``-cc``
      -  ``--color-code``
      -  ChemDraw colour-table index for the target acidic proton (``.cdxml`` / ``.cdx`` only).

   -  -  ``-rpi``
      -  ``--reference-proton-index``
      -  Optional when ``-r`` is a CDXML file with a uniquely coloured reference proton.

   -  -  ``-rcc``
      -  ``--reference-color-code``
      -  ChemDraw colour-table index for the reference proton (``.cdxml`` / ``.cdx`` reference only).

**Job output file naming**

Each pKa job creates gas-phase opt+freq and solvent single-point sub-jobs. Output filenames follow the sub-job label:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   -  -  Sub-job label suffix
      -  Species
   -  -  ``_HA_opt``
      -  Target acid HA (gas-phase opt+freq)
   -  -  ``_A_opt``
      -  Target conjugate base A⁻ (gas-phase opt+freq)
   -  -  ``_HA_sp``
      -  HA solvent single-point
   -  -  ``_A_sp``
      -  A⁻ solvent single-point
   -  -  ``_HRef_opt`` / ``_Ref_opt``
      -  Reference acid / conjugate base (proton exchange only)
   -  -  ``_HRef_sp`` / ``_Ref_sp``
      -  Reference solvent single-points (proton exchange only)

Gaussian batch jobs use the input stem as the job label (e.g. ``acid1_HA_opt.log``). ORCA batch jobs append ``_pka`` to
the stem (e.g. ``acid1_pka_HA_opt.out``). The output-analysis autodiscovery convention below is aligned with the
``{basename}_pka_*`` pattern used by ORCA submission and by typical batch output tables.

*****************************************
 Output Analysis (``chemsmart run pka``)
*****************************************

All post-processing lives under ``chemsmart run pka``. No Gaussian or ORCA backend is invoked during analysis.

Thermochemistry extraction
==========================

For each output file, analysis:

#. Open the file and detect the program (Gaussian or ORCA).
#. Reads the gas-phase SCF energy and quasi-harmonic Gibbs free energy (for opt+freq outputs).
#. Reads the solvent-phase SCF energy (for single-point outputs).
#. Raises a clear error if a required quantity cannot be extracted.

Computing pKa from Output Files (``analyze``)
=============================================

If you have completed output files for a single acid, compute pKa with ``analyze``.

**Proton exchange (default)**

Only ``-ha`` and ``-hr`` are strictly required; the remaining six companion files are auto-discovered when they follow
the naming convention below. ``-rp`` / ``--reference-pka`` is required.

.. code:: bash

   chemsmart run pka analyze \
       -ha acid1_pka_HA_opt.log \
       -hr ref_acid_pka_HRef_opt.log \
       -rp 6.75 \
       -T 333.15 -c 1.0 -csg 100 -ch 100

Provide all eight files explicitly when auto-discovery is not appropriate:

.. code:: bash

   chemsmart run pka analyze \
       -ha acid1_pka_HA_opt.log \
       -a acid1_pka_A_opt.log \
       -hr ref_acid_pka_HRef_opt.log \
       -r ref_acid_pka_Ref_opt.log \
       -has acid1_pka_HA_sp.log \
       -as acid1_pka_A_sp.log \
       -hrs ref_acid_pka_HRef_sp.log \
       -rs ref_acid_pka_Ref_sp.log \
       -rp 6.75 \
       -T 298.15

**Direct dissociation**

Four output files are required (HA, A⁻, and their solvent single-points). Both ``-s direct`` and ``-dG`` must be
specified on the ``pka`` group **before** the ``analyze`` subcommand:

.. code:: bash

   chemsmart run pka -s direct -dG -265.9 analyze \
       -ha acid1_pka_HA_opt.log \
       -T 298.15

Only ``-ha`` is strictly required; ``-a``, ``-has``, and ``-as`` are auto-discovered from the target-acid suffix
convention when omitted.

File autodetection (``analyze``)
================================

When companion paths are omitted, CHEMSMART derives them from the HA and HRef gas-phase files using the same suffix
patterns as ``batch-analyze``:

**From the HA gas-phase file (``-ha``)**

-  ``<basename>_pka_A_opt.<ext>`` — conjugate base gas-phase
-  ``<basename>_pka_HA_sp.<ext>`` — HA solvent single-point
-  ``<basename>_pka_A_sp.<ext>`` — conjugate base solvent SP

**From the HRef gas-phase file (``-hr``)**

-  ``<basename>_pka_Ref_opt.<ext>`` — reference conjugate base
-  ``<basename>_pka_HRef_sp.<ext>`` — reference acid solvent SP
-  ``<basename>_pka_Ref_sp.<ext>`` — reference conjugate base solvent SP

Alternative suffixes (``_HRef_opt``, ``_pka_cb``, etc.) are also recognised. The file extension (``.log`` or ``.out``)
is chosen from the detected program. Override any auto-discovered path with the corresponding flag.

If a required file is missing or cannot be parsed, analysis stops with a clear error (missing paths or missing
thermochemistry data).

Batch Processing of Output Files (``batch-analyze``)
====================================================

Parse a table of pre-computed output file paths to calculate pKa values in batch.

**Proton exchange (default)**

.. code:: bash

   chemsmart run pka -T 333.15 -c 1.0 -csg 100 -ch 100 batch-analyze \
       -o pka_output_table.csv \
       -O results.dat

**Direct dissociation**

Both ``-s direct`` and ``-dG`` are required:

.. code:: bash

   chemsmart run pka -s direct -dG -265.9 batch-analyze \
       -o pka_output_table_direct.csv \
       -O results_direct.dat

The formatted batch summary table is printed to stdout. When ``-O`` / ``--output-results`` is given, the same formatted
report is written to that file (not a wide CSV of input columns).

Output table format
-------------------

The output table (``-o`` / ``--output-table``) must contain at least a ``basename`` column. Other file paths may be
omitted and are auto-discovered from ``basename`` when blank.

**Required column**

-  ``basename``: Unique identifier for each acid. Used for file auto-discovery.

**Target-acid columns (both schemes)**

Accepted header aliases include ``ha_opt``, ``a_opt``, ``ha_solv``, ``a_solv``, etc. (see column list below).

When blank, CHEMSMART searches for ``<basename><suffix>.<ext>`` in the current working directory. Suffixes are tried in
order; both ``.log`` and ``.out`` are tested.

.. list-table::
   :header-rows: 1
   :widths: 20 50 30

   -  -  Column
      -  Description
      -  Auto-discovery suffixes (first match wins)

   -  -  ``ha_gas``
      -  HA gas-phase opt+freq output
      -  ``_pka_HA_opt``, ``_pka_HA``, ``_pka``

   -  -  ``a_gas``
      -  A⁻ gas-phase opt+freq output
      -  ``_pka_A_opt``, ``_pka_A``, ``_pka_cb``

   -  -  ``ha_sp``
      -  HA solvent single-point output
      -  ``_pka_HA_sp``, ``_pka_sp``

   -  -  ``a_sp``
      -  A⁻ solvent single-point output
      -  ``_pka_A_sp``, ``_pka_cb_sp``

**Reference-acid columns (proton exchange only)**

Ignored for direct dissociation. Blank reference columns inherit values from the previous row.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   -  -  Column
      -  Description
   -  -  ``href_gas``
      -  HRef gas-phase opt+freq output (not auto-discovered from ``basename``; provide explicitly or inherit)
   -  -  ``ref_gas``
      -  Ref⁻ gas-phase opt+freq output
   -  -  ``href_sp``
      -  HRef solvent single-point output
   -  -  ``ref_sp``
      -  Ref⁻ solvent single-point output
   -  -  ``pka_ref``
      -  Experimental pKa of the reference acid

**Example** ``pka_output_table.csv``:

.. code:: text

   basename,ha_gas,a_gas,ha_sp,a_sp,href_gas,ref_gas,href_sp,ref_sp,pka_ref
   phenol,,,,,ref_acid_pka_HRef_opt.log,ref_acid_pka_Ref_opt.log,ref_acid_pka_HRef_sp.log,ref_acid_pka_Ref_sp.log,6.75
   benzoic_acid,,,,,,,,,6.75

``batch-analyze`` options
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Short
      -  Long
      -  Description

   -  -  ``-o``
      -  ``--output-table``
      -  **Required.** Path to the output-file table.

   -  -  ``-O``
      -  ``--output-results``
      -  Optional path for the formatted results report. Stdout always receives the summary table.

   -  -  ``-p``
      -  ``--program``
      -  Require every populated output path to match ``gaussian`` or ``orca``. Default: ``auto`` (per-file detection;
         supports mixed tables).

Mixed Gaussian / ORCA tables
----------------------------

With ``-p auto`` (default), the program behind each output file is detected automatically. A batch table may contain
Gaussian ``.log`` targets and ORCA ``.out`` reference files in the same run. Use ``-p gaussian`` or ``-p orca`` only
when you want to **validate** that all populated paths belong to one backend.

*************************
 Analysis Scheme Options
*************************

These options apply to ``chemsmart run pka`` (``analyze`` and ``batch-analyze``). They are separate from submission
options on ``chemsmart run/sub gaussian ... pka`` and ``chemsmart run/sub orca ... pka``.

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Short
      -  Long
      -  Description

   -  -  ``-s``
      -  ``--scheme``
      -  Thermodynamic cycle: ``direct`` or ``proton exchange``. Default: ``proton exchange``.

   -  -  ``-dG``
      -  ``--delta-g-proton``
      -  :math:`G_{\text{soln}}(\text{H}^{+})` in kcal/mol for the direct cycle. **Required** when ``-s direct`` is used
         for analysis. Has no default during analysis.

.. note::

   For job submission, ``-dG`` defaults to ``-265.9`` kcal/mol. For output analysis you must pass ``-dG`` explicitly
   whenever ``-s direct`` is used.

*********************
 Output File Options
*********************

Used by ``analyze`` (not ``batch-analyze``, which reads paths from the table).

**Gas-phase optimization + frequency files**

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Short
      -  Long
      -  Description

   -  -  ``-ha``
      -  ``--ha``
      -  HA gas-phase opt+freq output.

   -  -  ``-a``
      -  ``--a``
      -  A⁻ gas-phase opt+freq output.

   -  -  ``-hr``
      -  ``--href``
      -  HRef gas-phase opt+freq output.

   -  -  ``-r``
      -  ``--ref``
      -  Ref⁻ gas-phase opt+freq output.

**Solvent single-point files**

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Short
      -  Long
      -  Description

   -  -  ``-has``
      -  ``--ha-solv``
      -  HA solvent single-point output.

   -  -  ``-as``
      -  ``--a-solv``
      -  A⁻ solvent single-point output.

   -  -  ``-hrs``
      -  ``--href-solv``
      -  HRef solvent single-point output.

   -  -  ``-rs``
      -  ``--ref-solv``
      -  Ref⁻ solvent single-point output.

*************************
 Thermochemistry Options
*************************

Shared by ``analyze`` and ``batch-analyze``.

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
      -  Cutoff frequency (cm⁻¹) for entropy using Grimme's quasi-RRHO. Default: ``100.0``.

   -  -  ``-cst``
      -  ``--cutoff-entropy-truhlar``
      -  Cutoff frequency (cm⁻¹) for entropy using Truhlar's quasi-RRHO. Mutually exclusive with ``-csg``.

   -  -  ``-ch``
      -  ``--cutoff-enthalpy``
      -  Cutoff frequency (cm⁻¹) for enthalpy using Head-Gordon's method. Default: ``100.0``.

   -  -  ``-rp``
      -  ``--reference-pka``
      -  Experimental pKa of HRef. Required for proton exchange analysis.

Output Format
=============

When computing pKa from output files, CHEMSMART prints a detailed summary. The format depends on the analysis scheme.

**Proton exchange**

.. code:: text

   ==============================================================================
   pKa Calculation - Dual-level Proton Exchange Scheme
   ==============================================================================
   Reaction: HA + Ref- -> A- + HRef
   Temperature: 373.15 K

   Method:
     G_corr = qh-G(T) - E_gas  (from gas-phase freq calculation)
     G_soln = E_solv + G_corr  (solution free energy)
     DG_soln = [G(A-)_soln + G(HRef)_soln] - [G(HA)_soln + G(Ref-)_soln]
     pKa = pKa_ref + DG_soln / (RT * ln10)
   ------------------------------------------------------------------------------

   Gas-Phase Electronic Energies (E_gas, au):
     HA:    -345.7419436500
     A-:    -344.9153986020
     HRef:  -365.8436493070
     Ref-:  -365.4561783660

   Thermal Corrections (G_corr = qh-G - E_gas, au):
     HA:    0.0931931305
     A-:    0.0758935969
     HRef:  0.1404528844
     Ref-:  0.1267467582

   Solvent Single-Point Energies (E_solv, au):
     HA:    -346.4882221850
     A-:    -345.8989956310
     HRef:  -366.5974351550
     Ref-:  -366.1368369100

   Solution Free Energies (G_soln = E_solv + G_corr, au):
     HA:    -346.3950290545
     A-:    -345.8231020341
     HRef:  -366.4569822706
     Ref-:  -366.0100901518
   ------------------------------------------------------------------------------

   pKa Calculation:
     DG_soln = 0.1250349015 au
             = 78.4606 kcal/mol
     pKa(HRef)_ref = 6.75

     *** Computed pKa(HA) = 52.70 ***
   ==============================================================================

**Direct dissociation**

.. code:: text

   ==============================================================================
   pKa Calculation - Direct Dissociation Scheme
   ==============================================================================
   Reaction: HA -> A- + H+
   Temperature: 298.15 K

   Method:
     G_corr = qh-G(T) - E_gas  (from gas-phase freq calculation)
     G_soln = E_solv + G_corr  (solution free energy)
     DG_diss = G_soln(A-) + G_soln(H+) - G_soln(HA)
     pKa = DG_diss / (2.303 * R * T)
   ------------------------------------------------------------------------------

   Gas-Phase Electronic Energies (E_gas, au):
     HA:  -345.7419436500
     A-:  -344.9153986020

   Thermal Corrections (G_corr = qh-G - E_gas, au):
     HA:  0.0931931305
     A-:  0.0758935969

   Solvent Single-Point Energies (E_solv, au):
     HA:  -346.4882221850
     A-:  -345.8989956310

   Solution Free Energies (G_soln = E_solv + G_corr, au):
     HA:  -346.3950290545
     A-:  -345.8231020341
   ------------------------------------------------------------------------------

   pKa Calculation:
     G_soln(H+) = -265.9000 kcal/mol
     DG_diss = 0.1250349015 au
             = 78.4606 kcal/mol

     *** Computed pKa(HA) = 52.70 ***
   ==============================================================================

**Batch analyze output**

``batch-analyze`` prints a compact table whose ΔG column header matches the scheme. The same formatted report is written
to ``-O`` when provided:

.. code:: text

   ==============================================================================
   Batch pKa Results (Dual-level Proton Exchange)
   ==============================================================================
   Temperature: 298.15 K
   Pressure: 1.0 atm
   basename                              pKa     DG_soln (kcal/mol)
   ------------------------------------------------------------------------------
   phenol                               10.12              13.4567
   benzoic_acid                          4.20               5.7890
   ==============================================================================

For direct dissociation, the header reads ``Batch pKa Results (Direct Dissociation)`` and the column is labeled
``DG_diss (kcal/mol)``.

References
==========

#. Tissandier, M. D. et al. (1998). *J. Phys. Chem. A*, 102, 7787. (Absolute proton solvation energy)
#. Grimme, S. (2012). *Chem. Eur. J.*, 18, 9955. (Quasi-RRHO method)
#. Marenich, A. V.; Cramer, C. J.; Truhlar, D. G. (2009). *J. Phys. Chem. B*, 113, 6378. (SMD solvation model)

See Also
========

-  :doc:`thermochemistry-analysis`
