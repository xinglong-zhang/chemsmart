.. _orca-reaction-calculations:

############################
 ORCA Reaction Calculations
############################

This page covers **ORCA reaction job submission**. Path search uses project NEB-TS instead of QST. Submit with
``chemsmart sub/run chain … reaction --program orca``.

.. note::

   Output-file analysis is not part of ``chain … reaction``. After calculations finish, use ``chemsmart run
   thermochemistry`` on the child outputs. See :ref:`reaction-calculations`. A dedicated ``chemsmart run reaction
   analyze`` command is planned as a follow-up.

.. contents:: Table of Contents
   :local:
   :depth: 2

*************
 Quick Start
*************

**Case 1 — TS guess**

``-f`` is the TS structure. Guess (NEB) is skipped.

.. code:: bash

   chemsmart sub chain -p combined -f ts_guess.xyz -c 0 -m 1 \
       reaction --program orca

Where:

-  ``-p combined``: Chain project (ORCA theory from the ``orca:`` YAML alias)
-  ``-f ts_guess.xyz``: TS guess
-  ``-c 0 -m 1``: Charge and multiplicity of the TS

This runs OptTS + freq (``ts_settings()``) and a solvent single-point (``sp_settings()``). Hessian/ScanTS flags from
``orca ts`` are not re-exposed; TS children use project ``ts_settings()``.

**Case 2 — reactant + product (NEB-TS)**

``-f`` is the reactant. ``--product`` selects path search. Reactant and product are optimized first. The Guess phase
reuses ``ORCANEBJob`` with project ``neb_settings()`` (default ``NEB-TS``): optimized reactant start, optimized product
``ending_xyzfile``, optional ``--ts-guess`` as ``intermediate_xyzfile``. Endpoint optimization is handled by the
reaction chain (not ORCA ``preopt_ends``). Then OptTS characterization.

.. code:: bash

   chemsmart sub chain -p combined -f reactant.xyz -c 0 -m 1 \
       reaction --program orca --product product.xyz

   chemsmart sub chain -p combined -f reactant.xyz -c 0 -m 1 \
       reaction --program orca --product product.xyz --ts-guess ts.xyz

``orca ts`` remains the single-structure TS search and ``orca neb`` remains the standalone NEB command. See
:doc:`orca-transition-state`.

************************
 Job Output File Naming
************************

For a job with label ``sn2``:

.. code:: text

   sn2_R_opt.out                 # reactant endpoint opt+freq (case 2)
   sn2_P_opt.out                 # product endpoint opt+freq (case 2)
   sn2_neb.inp / sn2_neb.out     # NEB-TS (case 2; optimized endpoints)
   sn2_neb_P.xyz                 # optimized product geometry written for NEB
   sn2_neb_TS.xyz                # optional TS intermediate
   sn2_TS_opt.out                # TS gas-phase opt+freq
   sn2_TS_sp.out                 # TS solvent single-point
   sn2_R_sp.out / sn2_P_sp.out   # matching SP children

************************************
 Batch Processing with Input Tables
************************************

.. code:: bash

   chemsmart sub chain -p combined -f reactions.csv reaction --program orca batch

.. note::

   When ``-f`` is a submission table, chain ``-c`` / ``-m`` are not required.

Table format is the same as Gaussian. See :ref:`reaction-calculations`.

On HPC clusters, ``chemsmart sub ... reaction batch`` writes one scheduler script per ``reaction_id`` with a
reconstructed ``reaction submit`` command.

************
 Parameters
************

ORCA reaction options match Gaussian submission options. Guess uses project ``neb_settings()`` rather than QST.

.. list-table::
   :header-rows: 1
   :widths: 15 20 65

   -  -  Short
      -  Long
      -  Description

   -  -  ``-r``
      -  ``--reactant``
      -  Reactant geometry file. Repeatable. With ``--product``, parent ``-f`` is the NEB intermediate.

   -  -
      -  ``--product``
      -  Product geometry (NEB ending structure). Repeatable. Presence selects path search when no ``--reactant`` is
         given.

   -  -
      -  ``--ts-guess``
      -  NEB intermediate when ``-f`` is the reactant. Requires ``--product``.

   -  -  ``-S`` / ``-R``
      -  ``--skip-completed`` / ``--no-skip-completed``
      -  Skip completed child jobs (default) or rerun them.

The full standalone NEB CLI (``-j``, ``-e``, ``-i``, image count, …) is not re-exposed on ``reaction``.

**********
 Examples
**********

Example 1: TS Guess Only
========================

.. code:: bash

   chemsmart sub chain -p orca_b3lyp -f ts_guess.xyz -c 0 -m 1 \
       reaction --program orca

Example 2: NEB-TS then OptTS
============================

.. code:: bash

   chemsmart sub chain -p orca_b3lyp -f reactant.xyz -c 0 -m 1 \
       reaction --program orca --product product.xyz

Example 3: NEB with TS Intermediate
===================================

.. code:: bash

   chemsmart sub chain -p orca_b3lyp -f reactant.xyz -c 0 -m 1 \
       reaction --program orca --product product.xyz --ts-guess ts.xyz

Example 4: Batch Submission from CSV
====================================

.. code:: bash

   chemsmart sub chain -p orca_b3lyp -f reactions.csv reaction --program orca batch

Example 5: Thermochemistry on Completed Outputs
===============================================

Until ``chemsmart run reaction analyze`` exists:

.. code:: bash

   chemsmart run thermochemistry -f sn2_TS_opt.out
   chemsmart run thermochemistry -f sn2_R_opt.out
   chemsmart run thermochemistry -f sn2_P_opt.out

**********
 See Also
**********

-  :ref:`reaction-calculations`
-  :ref:`gaussian-reaction-calculations`
-  :doc:`orca-transition-state`
-  :doc:`thermochemistry-analysis`
