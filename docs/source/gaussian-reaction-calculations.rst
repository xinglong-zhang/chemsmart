.. _gaussian-reaction-calculations:

################################
 Gaussian Reaction Calculations
################################

This page covers **Gaussian reaction job submission** — endpoint optimization, QST2/QST3 path search on minimized
reactant and product, TS opt+freq, and solvent single-points. Submit with ``chemsmart sub/run chain … reaction --program
gaussian``.

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

``-f`` is the TS structure. Guess (QST) is skipped. Optional ``--reactant`` without ``--product`` adds reactant optima
to optimize.

.. code:: bash

   chemsmart sub chain -p combined -f ts_guess.xyz -c 0 -m 1 \
       reaction --program gaussian

Where:

-  ``-p combined``: Chain project (Gaussian theory from the ``gaussian:`` YAML alias)
-  ``-f ts_guess.xyz``: TS guess (XYZ, LOG, COM, …)
-  ``-c 0 -m 1``: Charge and multiplicity of the TS

This runs TS opt+freq (``ts_settings()``) and a solvent single-point (``sp_settings()``).

**Case 2 — reactant + product (QST2)**

``-f`` is the reactant. ``--product`` selects path search. Reactant and product are optimized first, then Gaussian
builds a QST2 ``.com`` (two geometry blocks from the minimized endpoints) and runs it as a COM job, then characterizes
the located TS.

.. code:: bash

   chemsmart sub chain -p combined -f reactant.xyz -c 0 -m 1 \
       reaction --program gaussian --product product.xyz

**Case 2 — QST3**

Pass a TS guess as well. Atom order must match across reactant, product, and TS guess.

.. code:: bash

   chemsmart sub chain -p combined -f reactant.xyz -c 0 -m 1 \
       reaction --program gaussian --product product.xyz --ts-guess ts.xyz

   # Equivalent: -f is the QST3 guess when both --reactant and --product are set
   chemsmart sub chain -p combined -f ts.xyz -c 0 -m 1 \
       reaction --program gaussian --reactant reactant.xyz --product product.xyz

``gaussian ts`` remains the single-structure TS search. See :doc:`gaussian-transition-state`.

************************
 Job Output File Naming
************************

Sub-job labels determine output filenames. For a job with label ``sn2`` (the default when submitting ``sn2.xyz``):

.. code:: text

   sn2_R_opt.log                 # reactant endpoint opt+freq (case 2)
   sn2_R1_opt.log / sn2_R2_opt   # extra reactant fragments (case 1)
   sn2_P_opt.log                 # product endpoint opt+freq (case 2)
   sn2_qst.com / sn2_qst.log     # QST2 or QST3 (case 2; optimized endpoints)
   sn2_TS_opt.log                # TS gas-phase opt+freq
   sn2_TS_sp.log                 # TS solvent single-point
   sn2_R_sp.log / sn2_P_sp.log   # matching SP children

************************************
 Batch Processing with Input Tables
************************************

Pass a ``.csv`` or whitespace-delimited ``.txt`` file via ``-f`` and invoke ``reaction batch`` (or omit the subcommand —
batch is selected automatically when ``-f`` points to a submission table).

.. code:: bash

   chemsmart sub chain -p combined -f reactions.csv reaction --program gaussian batch

.. note::

   When ``-f`` is a submission table, chain ``-c`` / ``-m`` are not required; charge and multiplicity are read from each
   table row.

On HPC clusters, use ``chemsmart sub`` instead of ``chemsmart run``; each ``reaction_id`` receives its own scheduler
script with a reconstructed ``reaction submit`` command. See :ref:`reaction-calculations`.

Table format is documented in :ref:`reaction-calculations`.

************
 Parameters
************

Reaction Options
================

.. list-table::
   :header-rows: 1
   :widths: 15 20 65

   -  -  Short
      -  Long
      -  Description

   -  -  ``-r``
      -  ``--reactant``
      -  Reactant geometry file. Repeatable for extra fragments. With ``--product``, parent ``-f`` is the TS guess.

   -  -  ``-p``
      -  ``--product``
      -  Product geometry file. Repeatable. Presence selects path search (case 2) when no ``--reactant`` is given.

   -  -  ``-ts``
      -  ``--ts-guess``
      -  QST3 intermediate when ``-f`` is the reactant. Requires ``--product``.

   -  -  ``-S`` / ``-R``
      -  ``--skip-completed`` / ``--no-skip-completed``
      -  Skip completed child jobs (default) or rerun them.

QST notes
=========

-  QST2 uses two coordinate blocks (reactant, product). QST3 adds a TS-guess block.
-  Structures must have the same number of atoms and the same atom order.
-  The Guess job reuses project TS theory with ``opt=qst2`` or ``opt=qst3``; it does not change the shared Gaussian
   writer or route builder.

**********
 Examples
**********

Example 1: TS Guess Only
========================

.. code:: bash

   chemsmart sub chain -p b3lyp_project -f ts_guess.xyz -c 0 -m 1 \
       reaction --program gaussian

Example 2: QST2 from Reactant and Product
=========================================

.. code:: bash

   chemsmart sub chain -p b3lyp_project -f reactant.xyz -c 0 -m 1 \
       reaction --program gaussian --product product.xyz

Example 3: QST3
===============

.. code:: bash

   chemsmart sub chain -p b3lyp_project -f reactant.xyz -c 0 -m 1 \
       reaction --program gaussian --product product.xyz --ts-guess ts.xyz

Example 4: Batch Submission from CSV
====================================

.. code:: bash

   chemsmart sub chain -p b3lyp_project -f reactions.csv reaction --program gaussian batch

Example 5: Thermochemistry on Completed Outputs
===============================================

Until ``chemsmart run reaction analyze`` exists:

.. code:: bash

   chemsmart run thermochemistry -f sn2_TS_opt.log
   chemsmart run thermochemistry -f sn2_R_opt.log
   chemsmart run thermochemistry -f sn2_P_opt.log

**********
 See Also
**********

-  :ref:`reaction-calculations`
-  :ref:`orca-reaction-calculations`
-  :doc:`gaussian-transition-state`
-  :doc:`thermochemistry-analysis`
