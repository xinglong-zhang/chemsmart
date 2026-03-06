.. _pka-calculations:

##################
 pKa Calculations
##################

This module provides tools for computing acid dissociation constants (pKa) using quantum chemistry calculations with
proper thermodynamic cycles. Both Gaussian and ORCA are supported.

.. toctree::
   :maxdepth: 2
   :caption: pKa Calculations

   gaussian-pka-calculations
   orca-pka-calculations

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
   SMD or CPCM).

#. **Total free energy in solution**:

   .. math::

      G_{\text{soln}} = E_{\text{solv}} + G_{\text{corr}}

.. note::

   All internal energies are stored in Hartree (au). Only :math:`\Delta G_{\text{soln}}` is converted to kcal/mol for
   the pKa formula (1 Hartree = 627.5094740631 kcal/mol).

**********************************
 Batch Processing of Output Files
**********************************

You can also parse a table of pre-computed output files to calculate pKa values in batch.

.. code:: bash

   chemsmart run pka -O pka_output_table.csv --output-results computed_pka.csv -rp 6.75

The output table (e.g., ``pka_output_table.csv``) must contain columns for the basename and paths to all required output
files.

Required columns:

-  ``basename``: A unique identifier for each acid.
-  ``ha_gas``: Path to the gas-phase optimization output for HA.
-  ``a_gas``: Path to the gas-phase optimization output for A⁻.
-  ``ha_sp``: Path to the solvent single-point output for HA.
-  ``a_sp``: Path to the solvent single-point output for A⁻.
-  ``hb_gas`` (optional): Path to the gas-phase optimization output for the reference acid HB.
-  ``b_gas`` (optional): Path to the gas-phase optimization output for B⁻.
-  ``hb_sp`` (optional): Path to the solvent single-point output for HB.
-  ``b_sp`` (optional): Path to the solvent single-point output for B⁻.
-  ``pka_ref`` (optional): The pKa of the reference acid.

.. note::

   If reference acid columns (``hb_gas``, ``b_gas``, etc.) are left blank for a row, the values from the most recently
   defined reference acid in a previous row will be used. The ``pka_ref`` can also be provided via the ``-rp``
   command-line option.

The results, including the computed pKa for each entry, will be saved to the file specified by ``--output-results``.

*********************************
 Computing pKa from Output Files
*********************************

If you have already completed output files for a single acid, you can compute the pKa directly without resubmitting
jobs. This is the recommended way to get the pKa value after your calculations are done.

.. code:: bash

   chemsmart run pka \
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

The program will automatically detect whether the output files are from Gaussian or ORCA.

Output File Options
===================

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

   -  -  ``-rp``
      -  ``--reference-pka``
      -  Experimental pKa of HB. Required for output file parsing mode.

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

-  :ref:`thermochemistry-analysis`
