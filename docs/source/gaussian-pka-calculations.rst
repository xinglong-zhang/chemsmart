.. _gaussian-pka-calculations:

###########################
 Gaussian pKa Calculations
###########################

This module provides tools for computing acid dissociation constants (pKa) using Gaussian electronic structure
calculations with proper thermodynamic cycles.

.. contents:: Table of Contents
   :local:
   :depth: 2

********
 Theory
********

The pKa of an acid HA in aqueous solution is defined by the equilibrium:

.. math::

   \text{HA}_{(\text{aq})} \rightleftharpoons \text{A}^{-}_{(\text{aq})} + \text{H}^{+}_{(\text{aq})}

The pKa is related to the standard Gibbs free energy change by:

.. math::

   \text{p}K_{\text{a}} = \frac{\Delta G^{\circ}_{\text{aq}}}{2.303 \cdot R \cdot T}

where :math:`R` is the gas constant and :math:`T` is the temperature.

Thermodynamic Cycles
====================

CHEMSMART supports two thermodynamic cycles for pKa calculations:

**1. Proton Exchange (Isodesmic) Cycle** (Default, Recommended)

Uses a reference acid HB with known experimental pKa to cancel systematic errors:

.. math::

   \text{HA} + \text{B}^{-} \rightarrow \text{A}^{-} + \text{HB}

The pKa is computed as:

.. math::

   \text{p}K_{\text{a}}(\text{HA}) = \text{p}K_{\text{a}}(\text{HB}) + \frac{\Delta G_{\text{soln}}}{2.303 \cdot R \cdot T}

where:

.. math::

   \Delta G_{\text{soln}} = \left[ G(\text{A}^{-})_{\text{soln}} + G(\text{HB})_{\text{soln}} \right] - \left[ G(\text{HA})_{\text{soln}} + G(\text{B}^{-})_{\text{soln}} \right]

**2. Direct Cycle**

Uses the absolute free energy of a proton in water:

.. math::

   \text{p}K_{\text{a}} = \frac{G(\text{A}^{-})_{\text{aq}} - G(\text{HA})_{\text{aq}} + \Delta G^{\circ}(\text{H}^{+})_{\text{aq}}}{2.303 \cdot R \cdot T}

Default value: :math:`\Delta G^{\circ}(\text{H}^{+})_{\text{aq}} = -265.9` kcal/mol (Tissandier et al., 1998).

Dual-Level Approach
===================

CHEMSMART implements a dual-level approach for accurate solvation free energies:

#. **Thermal corrections** (:math:`G_{\text{corr}}`) from gas-phase frequency calculations using quasi-harmonic Gibbs
   free energy:

   .. math::

      G_{\text{corr}} = G_{\text{qh}}(T) - E_{\text{gas}}

#. **Solvent energies** (:math:`E_{\text{solv}}`) from high-level single-point calculations in implicit solvent (e.g.,
   SMD).

#. **Total free energy in solution**:

   .. math::

      G_{\text{soln}} = E_{\text{solv}} + G_{\text{corr}}

.. note::

   All internal energies are stored in Hartree (au). Only :math:`\Delta G_{\text{soln}}` is converted to kcal/mol for
   the pKa formula (1 Hartree = 627.5094740631 kcal/mol).

*************
 Quick Start
*************

Basic pKa Calculation
=====================

To run a pKa calculation for an acid molecule from an XYZ file:

.. code:: bash

   chemsmart run gaussian -p my_project -f acid.xyz -c 0 -m 1 pka -pi 10

Where:

-  ``-p my_project``: Project settings file (defines functional, basis set, etc.)
-  ``-f acid.xyz``: Input geometry file (XYZ, LOG, or COM format)
-  ``-c 0 -m 1``: Charge and multiplicity of the protonated acid (HA)
-  ``-pi 10``: 1-based index of the proton to remove

Example with Acetic Acid
------------------------

.. code:: bash

   # Run pKa calculation for acetic acid (proton at index 10)
   chemsmart run gaussian -p my_project -f acetic_acid.xyz -c 0 -m 1 pka -pi 10

This will:

#. Optimize HA in gas phase (opt + freq)
#. Optimize A⁻ in gas phase (opt + freq)
#. Run single-point on optimized HA in solution (SMD/water)
#. Run single-point on optimized A⁻ in solution (SMD/water)

Proton Exchange Cycle with Reference Acid
=========================================

For more accurate results, use a reference acid with known pKa:

.. code:: bash

   chemsmart run gaussian -p my_project -f acid.xyz -c 0 -m 1 pka \
       -pi 10 \
       -t "proton exchange" \
       -r water.xyz \
       -rpi 1 \
       -rc 0 \
       -rm 1 \
       -rp 14.0

This additionally runs calculations for the reference acid HB and its conjugate base B⁻.

************************************
 Batch Processing with Input Tables
************************************

For high-throughput pKa calculations, you can use a table-driven approach to submit multiple jobs at once.

Using an Input Table
====================

Provide the ``--input-table`` flag and pass a ``.csv`` or ``.txt`` file to the main ``-f`` option.

.. code:: bash

   chemsmart run gaussian -p my_project -f pka_input_table.csv -i

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

Computing pKa from Existing Output Files
========================================

If you already have completed Gaussian output files, compute pKa directly:

.. code:: bash

   chemsmart run gaussian pka \
       -pi 1 \
       -ha acid_opt.log \
       -a conjugate_base_opt.log \
       -hb reference_acid_opt.log \
       -b reference_base_opt.log \
       -has acid_sp_smd.log \
       -as conjugate_base_sp_smd.log \
       -hbs reference_acid_sp_smd.log \
       -bs reference_base_sp_smd.log \
       -rp 6.75 \
       -T 298.15

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
      -  **Required.** 1-based index of the proton to remove for deprotonation.

   -  -  ``-t``
      -  ``--thermodynamic-cycle``
      -  Thermodynamic cycle type: ``"direct"`` or ``"proton exchange"`` (default).

   -  -  ``-cc``
      -  ``--conjugate-base-charge``
      -  Charge of A⁻. Defaults to ``charge - 1``.

   -  -  ``-cm``
      -  ``--conjugate-base-multiplicity``
      -  Multiplicity of A⁻. Defaults to same as HA.

Reference Acid Options (Proton Exchange Cycle)
==============================================

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Short
      -  Long
      -  Description

   -  -  ``-r``
      -  ``--reference``
      -  Path to geometry file for reference acid (HB).

   -  -  ``-rpi``
      -  ``--reference-proton-index``
      -  1-based index of proton to remove from HB. Required with ``-r``.

   -  -  ``-rc``
      -  ``--reference-charge``
      -  Charge of HB. Required with ``-r``.

   -  -  ``-rm``
      -  ``--reference-multiplicity``
      -  Multiplicity of HB. Required with ``-r``.

   -  -  ``-rcc``
      -  ``--reference-conjugate-base-charge``
      -  Charge of B⁻. Defaults to ``reference_charge - 1``.

   -  -  ``-rcm``
      -  ``--reference-conjugate-base-multiplicity``
      -  Multiplicity of B⁻. Defaults to ``reference_multiplicity``.

   -  -  ``-rp``
      -  ``--reference-pka``
      -  Experimental pKa of HB. Required for output file parsing mode.

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
      -  Solvation model for solution phase SP. Default: ``SMD``.

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

   -  -  ``-conc``
      -  ``--concentration``
      -  Concentration in mol/L. Default: ``1.0`` mol/L.

   -  -  ``-csg``
      -  ``--cutoff-entropy-grimme``
      -  Cutoff frequency (cm⁻¹) for entropy using Grimme's quasi-RRHO. Default: ``100.0``.

   -  -  ``-ch``
      -  ``--cutoff-enthalpy``
      -  Cutoff frequency (cm⁻¹) for enthalpy using Head-Gordon's method. Default: ``100.0``.

   -  -  ``-dG``
      -  ``--delta-g-proton``
      -  :math:`\Delta G^{\circ}(\text{H}^{+})_{\text{aq}}` in kcal/mol for direct cycle. Default: ``-265.9``.

Output File Options (Parsing Mode)
==================================

For computing pKa from existing Gaussian output files:

**Gas-Phase Optimization + Frequency Files:**

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Short
      -  Long
      -  Description

   -  -  ``-ha``
      -  ``--ha-output``
      -  HA (protonated acid) gas-phase opt+freq output.

   -  -  ``-a``
      -  ``--a-output``
      -  A⁻ (conjugate base) gas-phase opt+freq output.

   -  -  ``-hb``
      -  ``--hb-output``
      -  HB (reference acid) gas-phase opt+freq output.

   -  -  ``-b``
      -  ``--b-output``
      -  B⁻ (reference conjugate base) gas-phase opt+freq output.

**Solvent Single-Point Files:**

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Short
      -  Long
      -  Description

   -  -  ``-has``
      -  ``--ha-solv-output``
      -  HA solvent single-point output.

   -  -  ``-as``
      -  ``--a-solv-output``
      -  A⁻ solvent single-point output.

   -  -  ``-hbs``
      -  ``--hb-solv-output``
      -  HB solvent single-point output.

   -  -  ``-bs``
      -  ``--b-solv-output``
      -  B⁻ solvent single-point output.

**********
 Examples
**********

Example 1: Simple pKa with Direct Cycle
=======================================

.. code:: bash

   chemsmart run gaussian -p b3lyp_project -f phenol.xyz -c 0 -m 1 pka \
       -pi 13 \
       -t direct \
       -T 298.15

Example 2: pKa with Proton Exchange and Custom Temperature
==========================================================

.. code:: bash

   chemsmart run gaussian -p m062x_project -f benzoic_acid.xyz -c 0 -m 1 pka \
       -pi 15 \
       -t "proton exchange" \
       -r acetic_acid.xyz \
       -rpi 10 \
       -rc 0 \
       -rm 1 \
       -T 310.15 \
       -sm SMD \
       -si water

Example 3: Extract pKa from Completed Calculations
==================================================

.. code:: bash

   chemsmart run gaussian pka \
       -pi 1 \
       -ha phenol_opt.log \
       -a phenolate_opt.log \
       -hb collidine_H_opt.log \
       -b collidine_opt.log \
       -has phenol_sp_smd.log \
       -as phenolate_sp_smd.log \
       -hbs collidine_H_sp_smd.log \
       -bs collidine_sp_smd.log \
       -rp 6.75 \
       -T 298.15 \
       -csg 100 \
       -ch 100

***************
 Output Format
***************

When computing pKa from output files, CHEMSMART prints a detailed summary:

.. code:: text

   ==============================================================================
   pKa Calculation - Dual-level Proton Exchange Scheme
   ==============================================================================
   Reaction: HA + B⁻ → A⁻ + HB
   Temperature: 373.15 K

   Method:
     G_corr = qh-G(T) - E_gas  (from gas-phase freq calculation)
     G_soln = E_solv + G_corr  (solution free energy)
     ΔG_soln = [G(A⁻)_soln + G(HB)_soln] - [G(HA)_soln + G(B⁻)_soln]
     pKa = pKa_ref + ΔG_soln / (RT × ln10)
   ------------------------------------------------------------------------------

   Gas-Phase Electronic Energies (E_gas, au):
     HA:  -345.7419436500
     A⁻:  -344.9153986020
     HB:  -365.8436493070
     B⁻:  -365.4561783660

   Thermal Corrections (G_corr = qh-G - E_gas, au):
     HA:  0.0931931305
     A⁻:  0.0758935969
     HB:  0.1404528844
     B⁻:  0.1267467582

   Solvent Single-Point Energies (E_solv, au):
     HA:  -346.4882221850
     A⁻:  -345.8989956310
     HB:  -366.5974351550
     B⁻:  -366.1368369100

   Solution Free Energies (G_soln = E_solv + G_corr, au):
     HA:  -346.3950290545
     A⁻:  -345.8231020341
     HB:  -366.4569822706
     B⁻:  -366.0100901518
   ------------------------------------------------------------------------------

   pKa Calculation:
     ΔG_soln = 0.1250349015 au
            = 78.4606 kcal/mol
     pKa(HB)_ref = 6.75

     *** Computed pKa(HA) = 52.70 ***
   ==============================================================================

************
 References
************

#. Tissandier, M. D. et al. (1998). *J. Phys. Chem. A*, 102, 7787. (Absolute proton solvation energy)
#. Grimme, S. (2012). *Chem. Eur. J.*, 18, 9955. (Quasi-RRHO method)
#. Marenich, A. V.; Cramer, C. J.; Truhlar, D. G. (2009). *J. Phys. Chem. B*, 113, 6378. (SMD solvation model)

**********
 See Also
**********

-  :doc:`thermochemistry-analysis`
