################################
 CREST and Trajectory Workflows
################################

This page covers workflows for submitting Gaussian jobs from CREST conformer ensembles and MD trajectories, with
optional pre-grouping to reduce redundant calculations.

.. note::

   **Behavior of -N parameter:**

   -  **Without -g (no grouping):** ``-N`` selects the N lowest-energy conformers directly.
   -  **With -g (grouping enabled):** ``-N`` targets N groups using complete-linkage distance levels, then uses the
      first member of each group (``group[0]``) as representative.

   This means ``-g irmsd -N 10`` gives 10 *structurally diverse* conformers, while ``-N 10`` alone gives the 10
   *lowest-energy* conformers (which may be similar). Because tied linkage levels are not split, grouping may return
   slightly more than N groups.

****************
 CREST Workflow
****************

The ``crest`` subcommand processes CREST conformer ensembles (multi-structure xyz files), optionally groups them, and
generates Gaussian input files for each unique conformer.

Basic Usage
===========

.. code:: bash

   # Submit all conformers without grouping, using project setting in proj.yaml
   chemsmart sub -s server gaussian -p proj -f crest_conformers.xyz -c 0 -m 1 crest

   # Group conformers first, then submit unique ones
   chemsmart sub -s server gaussian -p proj -f crest_conformers.xyz -c 0 -m 1 crest -g irmsd -N 10

   # Use threshold-based grouping
   chemsmart sub -s server gaussian -p proj -f crest_conformers.xyz -c 1 -m 3 crest -g irmsd -T 0.3

Grouping Options
================

Use ``-g`` to specify grouping strategy before submission:

.. code:: bash

   # iRMSD grouping (recommended for symmetric molecules)
   chemsmart sub -s server gaussian -p proj -f conformers.xyz -c 0 -m 1 crest -g irmsd -N 15

   # Torsion Fingerprint Deviation (TFD) grouping (good for flexible molecules)
   chemsmart sub -s server gaussian -p proj -f conformers.xyz -c 0 -m 1 crest -g tfd -T 0.1

   # Tanimoto similarity
   chemsmart sub -s server gaussian -p proj -f conformers.xyz -c 0 -m 1 crest -g tanimoto -T 0.9

Available strategies: ``rmsd``, ``hrmsd``, ``spyrmsd``, ``irmsd``, ``pymolrmsd``, ``tanimoto``, ``tfd``,
``isomorphism``, ``formula``, ``connectivity``, ``energy``

Representative workflow with grouping:

``grouping -> representative strategy reorders each group -> group[0] is the representative -> unique() takes group[0]
-> CREST/traj uses that representative directly``

Representative strategy (``-r, --representative``):

-  ``lowest`` (default): groups are ordered by energy.
-  ``center`` (matrix-based strategies only): groups are ordered by centrality.
-  ``top3`` (matrix-based strategies only): representative chosen from the three lowest-energy members by centrality;
   representative first, then remaining members in energy order.

Non-matrix strategies (``formula``, ``connectivity``, ``isomorphism``) support only ``lowest`` and reject
``center``/``top3``.

Strategy-specific Options
=========================

.. list-table::
   :header-rows: 1
   :widths: 25 75

   -  -  Option
      -  Description
   -  -  ``--inversion``
      -  For iRMSD: auto/on/off (default: auto)
   -  -  ``-ft, --fingerprint-type``
      -  For tanimoto: rdkit/morgan/maccs/usr/usrcat (default: rdkit)
   -  -  ``--use-weights/--no-use-weights``
      -  For TFD (Torsion Fingerprint Deviation): use torsion weights (default: True)
   -  -  ``--max-dev``
      -  For TFD (Torsion Fingerprint Deviation): equal/spec (default: equal)

*********************
 Trajectory Workflow
*********************

The ``traj`` subcommand processes MD trajectory or optimization trajectory files, selecting structures from the end of
the trajectory.

Basic Usage
===========

.. code:: bash

   # Use last 10% of trajectory (default)
   chemsmart sub -s server gaussian -p proj -f trajectory.xyz -c 0 -m 1 traj

   # Use last 50% of trajectory
   chemsmart sub -s server gaussian -p proj -f trajectory.xyz -c 0 -m 1 traj -x 0.5

   # Group selected structures before submission
   chemsmart sub -s server gaussian -p proj -f trajectory.xyz -c 0 -m 1 traj -x 0.3 -g spyrmsd -N 5

Trajectory Options
==================

.. list-table::
   :header-rows: 1
   :widths: 30 70

   -  -  Option
      -  Description
   -  -  ``-x, --proportion-structures-to-use``
      -  Proportion of structures from end of trajectory (0.0 < x <= 1.0, default: 0.1)

**Note:** Original conformer indices are preserved in the output to track which trajectory frames were selected.

****************
 Common Options
****************

Both ``crest`` and ``traj`` share these grouping options:

-  ``-g, --grouping-strategy``: choose grouping strategy
-  ``-T, --threshold`` or ``-N, --num-groups``: threshold-based grouping or target group count
-  ``-np, --num-procs``: requested process count (strategy-dependent; serial-only strategies warn and fall back to
   ``1``)
-  ``-r, --representative``: representative selection strategy (default ``lowest``)
-  ``-ih, --ignore-hydrogens``: ignore H atoms where supported by the selected strategy

Group XYZ files preserve representative-defined ordering (not always pure energy ordering):

-  ``lowest``: energy order
-  ``center``: centrality order
-  ``top3``: representative first, remaining members in energy order

********************
 Practical Examples
********************

Example 1: CREST Post-processing for DFT
========================================

Process CREST output, group to 15 conformers, and submit for DFT optimization:

.. code:: bash

   chemsmart sub -s server gaussian -p proj -f crest_conformers.xyz -c 0 -m 1 crest -g irmsd -N 15 -np 4

Example 2: Trajectory Analysis
==============================

Select last 30% of MD trajectory, group, and submit for single-point calculations:

.. code:: bash

   chemsmart sub -s server gaussian -p proj -f md_trajectory.xyz -c 0 -m 1 traj -x 0.3 -g tfd -T 0.5 --no-use-weights
