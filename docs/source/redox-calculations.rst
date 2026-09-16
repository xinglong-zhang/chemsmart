.. _redox-calculations:

####################
 Redox Calculations
####################

CHEMSMART provides redox workflows in two separate stages:

#. **Job submission** — generate and run Gaussian or ORCA calculations for the oxidized and reduced target and a
   reference couple. See :ref:`gaussian-redox-calculations` and :ref:`orca-redox-calculations`.

#. **Output analysis** — compute the exchange redox potential from completed output files using the backend-independent
   command ``chemsmart run redox analyze``. Analysis is **program-agnostic**: the same workflow reads Gaussian ``.log``
   and ORCA ``.out`` files. The same analysis is also registered under ``chain``; see :doc:`chain-jobs`.

The scheme is **exchange-only**. There is no isolated-electron :math:`G(\mathrm{e}^{-})` term and no bundled Fc/Fc+
geometries.

.. toctree::
   :maxdepth: 2
   :caption: Redox Calculations

   gaussian-redox-calculations
   orca-redox-calculations

.. contents:: Table of Contents
   :local:
   :depth: 2

************************
 Execution Architecture
************************

**Job submission**

-  ``chemsmart run/sub gaussian ... redox`` — prepare and run Gaussian redox calculations.

-  ``chemsmart run/sub orca ... redox`` — prepare and run ORCA redox calculations.

-  ``chemsmart run/sub chain -p combined ... redox --program {gaussian,orca}`` — same jobs, with theory and solvent from
   the chain YAML alias for ``--program``. See :ref:`chain-workflow-subcommands`.

-  Use ``chemsmart run`` for local preparation and execution; use ``chemsmart sub`` on HPC clusters to generate
   scheduler scripts.

-  The oxidized target comes from parent ``-f``. The reduced target uses the same geometry (or ``--red``) with charge
   ``ox − n``. The oxidized reference comes from ``--ref-ox`` (or registry ``ox_file``); the reduced reference uses the
   same geometry (or ``--ref-red``) with charge ``ox − n``.

**Output analysis**

-  ``chemsmart run redox analyze`` — single-system analysis. Only ``--ox-gas`` and ``--ref-ox-gas`` are required; the
   other six outputs are auto-discovered from CHEMSMART redox job labels (same convention as submit). The reference
   couple is inferred from Ref_ox/Ref_red formulas; ``-r`` overrides.

-  Analysis never invokes ``gaussian`` or ``orca`` job submission — only reads completed output files.

********
 Theory
********

The exchange reaction is:

.. math::

   \mathrm{Ox} + \mathrm{Ref_{red}} \rightarrow \mathrm{Red} + \mathrm{Ref_{ox}}

The exchange free-energy change is:

.. math::

   \Delta G_{\mathrm{exchange}} = G(\mathrm{Red}) + G(\mathrm{Ref_{ox}}) - G(\mathrm{Ox}) - G(\mathrm{Ref_{red}})

The target reduction potential on the reference scale is:

.. math::

   E_{\mathrm{target}} = E_{\mathrm{ref}} - \frac{\Delta G_{\mathrm{exchange}}}{n F}

:math:`n` must match the target couple and the reference couple. :math:`F` is Faraday's constant in C/mol. :math:`\Delta
G_{\mathrm{exchange}}` is converted from Hartree to J/mol before dividing by :math:`n F`. The summary reports
:math:`\Delta G` in au, kcal/mol, eV, and J/mol.

*********************
 Dual-Level Approach
*********************

Solution free energies follow the same dual-level construction as pKa:

#. **Thermal corrections** (:math:`G_{\mathrm{corr}}`) from gas-phase frequency calculations using quasi-harmonic Gibbs
   free energy:

   .. math::

      G_{\mathrm{corr}} = G_{\mathrm{qh}}(T) - E_{\mathrm{gas}}

#. **Solvent energies** (:math:`E_{\mathrm{solv}}`) from high-level single-point calculations in implicit solvent.

#. **Total free energy in solution**:

   .. math::

      G_{\mathrm{soln}} = E_{\mathrm{solv}} + G_{\mathrm{corr}}

**********************************
 Job Submission (Gaussian / ORCA)
**********************************

Job submission is backend-specific. Use the dedicated pages for full examples:

-  :ref:`gaussian-redox-calculations`
-  :ref:`orca-redox-calculations`

**Commands**

The built-in ``fc_fc+`` couple has :math:`E_{\mathrm{ref}} = 0.0` V and :math:`n = 1` on the Fc/Fc+ scale. It does not
bundle geometries, so ``--ref-ox`` is required unless another registered couple supplies ``ox_file``. ``--ref-red``
defaults to the same geometry with charge ``ox − n``.

.. code:: bash

   chemsmart run gaussian -p my_project -f ox.xyz -c 1 -m 2 redox \
       --ref-ox ref_ox.xyz

   chemsmart run orca -p my_project -f ox.xyz -c 1 -m 2 redox \
       --ref-ox ref_ox.xyz

   chemsmart run chain -p combined -f ox.xyz -c 1 -m 2 \
       redox --program gaussian --ref-ox ref_ox.xyz

Phases: Opt (Ox, Red) → Ref Opt → SP → Ref SP. Child labels for a job labelled ``mol_redox``:

.. code:: text

   mol_redox_ox_opt / mol_redox_red_opt
   mol_redox_RefOx_opt / mol_redox_RefRed_opt
   mol_redox_ox_sp / mol_redox_red_sp
   mol_redox_RefOx_sp / mol_redox_RefRed_sp

***************************
 Reference Couple Registry
***************************

Core jobs consume a ``RedoxReference`` from the public registry. Register additional couples without editing the redox
job class:

.. code:: python

   from chemsmart.analysis.redox import (
       RedoxReference,
       get_redox_reference,
       list_redox_references,
       register_redox_reference,
   )

   register_redox_reference(
       RedoxReference(
           name="custom_she",
           E_ref_V=0.40,
           n_electrons=1,
           scale="SHE",
           couple_label="Custom/Custom+",
           ox_file="ref_ox.xyz",
           red_file="ref_red.xyz",
           ox_charge=1,
           ox_multiplicity=2,
           red_charge=0,
           red_multiplicity=1,
       )
   )

``-r/--reference`` selects a registry name for **submit** (default ``fc_fc+``). ``-n/--n-electrons`` defaults to the
couple and must match it when given. Analyze does not use the registry: pass ``--e-ref`` (the reference potential in
volts).

*******************************************
 Output Analysis (``chemsmart run redox``)
*******************************************

All post-processing lives under ``chemsmart run redox analyze``. No Gaussian or ORCA backend is invoked during analysis.

Only ``--ox-gas`` and ``--ref-ox-gas`` are required. The remaining six outputs are auto-discovered when they follow
CHEMSMART redox labels (``<basename>_redox_red_opt``, ``<basename>_redox_ox_sp``, ``<basename>_redox_RefRed_opt``, and
so on). Override any path with the corresponding flag. ``--e-ref`` is the reference potential in volts used as
:math:`E_{\mathrm{ref}}`.

.. code:: bash

   chemsmart run redox analyze --e-ref 0.0 \
       --ox-gas mol_redox_ox_opt.log \
       --ref-ox-gas mol_redox_RefOx_opt.log

All eight outputs can still be given explicitly:

.. code:: bash

   chemsmart run redox analyze --e-ref 0.0 \
       --ox-gas mol_redox_ox_opt.log \
       --red-gas mol_redox_red_opt.log \
       --ref-ox-gas mol_redox_RefOx_opt.log \
       --ref-red-gas mol_redox_RefRed_opt.log \
       --ox-solv mol_redox_ox_sp.log \
       --red-solv mol_redox_red_sp.log \
       --ref-ox-solv mol_redox_RefOx_sp.log \
       --ref-red-solv mol_redox_RefRed_sp.log \
       -T 333.15 -c 1.0 -csg 100 -ch 100 \
       -o redox.dat

``-n`` defaults to 1. Example with a non-zero experimental reference potential:

.. code:: bash

   chemsmart run redox analyze --e-ref 0.2 -n 1 --ox-gas ...

Output Format
=============

``chemsmart run redox analyze`` prints a detailed summary to stdout (or writes it with ``-o``). The report lists
gas-phase energies, thermal corrections, solvent single-point energies, solution free energies, and the exchange redox
potential on the scale implied by ``--e-ref``.

.. code:: text

   ==============================================================================
   Redox Potential - Dual-level Exchange Scheme
   ==============================================================================
   Reaction: Ox + Ref_red → Red + Ref_ox
   Reference: fc_fc+ (Fc/Fc+, Fc/Fc+)
   n = 1
   Temperature: 298.15 K

   Method:
     G_corr = qh-G(T) - E_gas  (from gas-phase freq calculation)
     G_soln = E_solv + G_corr  (solution free energy)
     ΔG_exchange = G(Red) + G(Ref_ox) − G(Ox) − G(Ref_red)
     E_target = E_ref − ΔG_exchange / (n F)
   ------------------------------------------------------------------------------

   Gas-Phase Electronic Energies (E_gas, au):
     Ox:      1.0000000000
     Red:     1.1000000000
     Ref_ox:  2.0000000000
     Ref_red: 2.2000000000

   Thermal Corrections (G_corr = qh-G - E_gas, au):
     Ox:      0.0100000000
     Red:     0.0200000000
     Ref_ox:  0.0300000000
     Ref_red: 0.0400000000

   Solvent Single-Point Energies (E_solv, au):
     Ox:      0.9000000000
     Red:     1.0000000000
     Ref_ox:  1.8000000000
     Ref_red: 2.0000000000

   Solution Free Energies (G_soln = E_solv + G_corr, au):
     Ox:      0.9100000000
     Red:     1.0200000000
     Ref_ox:  1.8300000000
     Ref_red: 2.0400000000
   ------------------------------------------------------------------------------

   Redox Potential:
     ΔG_exchange = -0.1000000000 au
                 = -62.7509 kcal/mol
                 = -2.7211 eV
                 = -262549.9452 J/mol
     E_ref = 0.0000 V (Fc/Fc+)

     *** E_target = 2.7211 V (Fc/Fc+) ***
   ==============================================================================

When ``--e-ref`` is supplied without a registry couple during analyze, the reference line reads ``Reference: E_ref =
<value> V`` instead of the registry name. The scale suffix after ``E_ref`` and ``E_target`` is omitted in that case.

*************
 CLI Options
*************

Submit Options
==============

.. list-table::
   :header-rows: 1
   :widths: 30 70

   -  -  Option
      -  Description

   -  -  ``-r, --reference``
      -  Registry name of the reference couple (default ``fc_fc+``).

   -  -  ``-n, --n-electrons``
      -  Electrons transferred. Defaults to the reference couple; must match it when given.

   -  -  ``-rd, --red``
      -  Reduced target geometry. Defaults to the oxidized structure from parent ``-f`` with charge ``ox − n``.

   -  -  ``-rdc, --red-charge`` / ``-rdm, --red-multiplicity``
      -  Charge and multiplicity of the reduced target.

   -  -  ``-ro, --ref-ox`` / ``-rr, --ref-red``
      -  Oxidized reference geometry (required unless the registry provides ``ox_file``). Reduced reference defaults to
         the same structure with charge ``ox − n`` (or use ``--ref-red``).

   -  -  ``-roc, --ref-ox-charge`` / ``-rom, --ref-ox-multiplicity``
      -  Charge and multiplicity of the oxidized reference.

   -  -  ``-rrc, --ref-red-charge`` / ``-rrm, --ref-red-multiplicity``
      -  Charge and multiplicity of the reduced reference.

   -  -  ``-T, --temperature`` / ``-c, --concentration`` / ``-csg`` / ``-ch``
      -  Thermochemistry options stored on the job (same defaults as pKa analysis).

Analyze Options
===============

.. list-table::
   :header-rows: 1
   :widths: 30 70

   -  -  Option
      -  Description
   -  -  ``-er, --e-ref``
      -  Reference reduction potential in volts. Required. Used as :math:`E_{\mathrm{ref}}`.
   -  -  ``-n, --n-electrons``
      -  Electrons transferred (default 1).
   -  -  ``-oxg, --ox-gas`` / ``-rdg, --red-gas``
      -  Gas-phase opt+freq outputs for the target couple.
   -  -  ``-rog, --ref-ox-gas`` / ``-rrg, --ref-red-gas``
      -  Gas-phase opt+freq outputs for the reference couple.
   -  -  ``-oxs, --ox-solv`` / ``-rds, --red-solv``
      -  Solution-phase SP outputs for the target couple.
   -  -  ``-ros, --ref-ox-solv`` / ``-rrs, --ref-red-solv``
      -  Solution-phase SP outputs for the reference couple.
   -  -  ``-T, --temperature`` / ``-c, --concentration`` / ``-csg`` / ``-ch``
      -  Thermochemistry options for quasi-harmonic G (same defaults as pKa analysis).
   -  -  ``-o, --output``
      -  Write the formatted summary to this file instead of printing it.

**********
 See Also
**********

-  :ref:`gaussian-redox-calculations`
-  :ref:`orca-redox-calculations`
-  :ref:`chain-workflow-subcommands`
-  :ref:`pka-calculations`
-  :doc:`thermochemistry-analysis`
