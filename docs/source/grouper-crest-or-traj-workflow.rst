################################
 CREST and Trajectory Workflows
################################

This page covers workflows for submitting Gaussian jobs from CREST conformer ensembles and MD trajectories, with
optional pre-grouping to reduce redundant calculations.

.. note::

   **Behavior of -N parameter:**

   -  **Without -g (no grouping):** ``-N`` selects the N lowest-energy conformers directly.
   -  **With -g (grouping enabled):** ``-N`` first groups conformers into N groups, then selects the lowest-energy
      conformer from each group as representative.

   This means ``-g irmsd -N 10`` gives 10 *structurally diverse* conformers, while ``-N 10`` alone gives the 10
   *lowest-energy* conformers (which may be similar).

****************
 CREST Workflow
****************

The ``crest`` subcommand processes CREST conformer ensembles (multi-structure xyz files), optionally groups them, and
generates Gaussian input files for each unique conformer.

Basic Usage
===========

.. code:: bash

   # Submit all conformers without grouping
   chemsmart sub gaussian -f crest_conformers.xyz crest

   # Group conformers first, then submit unique ones
   chemsmart sub gaussian -f crest_conformers.xyz crest -g irmsd -N 10

   # Use threshold-based grouping
   chemsmart sub gaussian -f crest_conformers.xyz crest -g irmsd -T 0.3

Grouping Options
================

Use ``-g`` to specify grouping strategy before submission:

.. code:: bash

   # iRMSD grouping (recommended for symmetric molecules)
   chemsmart sub gaussian -f conformers.xyz crest -g irmsd -N 15

   # TFD grouping (good for flexible molecules)
   chemsmart sub gaussian -f conformers.xyz crest -g torsion -T 0.1

   # Tanimoto similarity
   chemsmart sub gaussian -f conformers.xyz crest -g tanimoto -T 0.9

Available strategies: ``rmsd``, ``hrmsd``, ``spyrmsd``, ``irmsd``, ``pymolrmsd``, ``tanimoto``, ``torsion``,
``isomorphism``, ``formula``, ``connectivity``, ``energy``

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
      -  For torsion: use torsion weights (default: True)
   -  -  ``--max-dev``
      -  For torsion: equal/spec (default: equal)

*********************
 Trajectory Workflow
*********************

The ``traj`` subcommand processes MD trajectory or optimization trajectory files, selecting structures from the end of
the trajectory.

Basic Usage
===========

.. code:: bash

   # Use last 10% of trajectory (default)
   chemsmart sub gaussian -f trajectory.xyz traj

   # Use last 50% of trajectory
   chemsmart sub gaussian -f trajectory.xyz traj -x 0.5

   # Group selected structures before submission
   chemsmart sub gaussian -f trajectory.xyz traj -x 0.3 -g spyrmsd -N 5

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

********************
 Practical Examples
********************

Example 1: CREST Post-processing for DFT
========================================

Process CREST output, group to 15 conformers, and submit for DFT optimization:

.. code:: bash

   chemsmart sub gaussian -f crest_conformers.xyz crest -g irmsd -N 15 -np 4

Example 2: Trajectory Analysis
==============================

Select last 30% of MD trajectory, group, and submit for single-point calculations:

.. code:: bash

   chemsmart sub gaussian -f md_trajectory.xyz traj -x 0.3 -g tfd -T 0.5 --no-use-weights
