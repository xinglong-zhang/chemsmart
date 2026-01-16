###############################################
 Conformational Sampling & Dynamics (Gaussian)
###############################################

This page covers CREST conformational searches and trajectory analysis using Gaussian.

************
 CREST Jobs
************

CREST (Conformer-Rotamer Ensemble Sampling Tool) performs systematic conformational searches to find low-energy
conformers. Chemsmart combines CREST with Gaussian for subsequent calculations.

.. code:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] crest -j <jobtype> [SUBCMD_OPTIONS]

CREST Options
=============

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-j, --jobtype``
      -  string
      -  Job type: opt, ts, modred, scan, sp (required)

   -  -  ``-N, --num-confs-to-run``
      -  int
      -  Number of conformers to process

.. warning::

   The ``-j, --jobtype`` option is required for CREST jobs.

Basic Usage
===========

Run optimization on CREST conformers:

.. code:: bash

   chemsmart sub gaussian -p project -f crest_conformers.xyz -c 1 -m 1 crest -j opt

Run TS optimization on conformers:

.. code:: bash

   chemsmart sub gaussian -p project -f crest_conformers.xyz -c 0 -m 1 crest -j ts

Run modred optimization with constraints:

.. code:: bash

   chemsmart sub gaussian -p project -f crest_conformers.xyz -c 0 -m 2 crest -j modred -c [1,4]

.. note::

   The ``<input_file>`` is typically the CREST output file named ``crest_conformers.xyz``.

Examples
========

Process 10 lowest-energy conformers:

.. code:: bash

   chemsmart sub gaussian -p project -f crest_conformers.xyz -l structure_from_lowest -c 0 -m 1 crest -j opt -N 10

Output files are named ``structure_from_lowest_opt_c1`` through ``structure_from_lowest_opt_c10``.

.. note::

   If a job terminates before all conformers are processed (e.g., walltime limit), resubmitting will continue from where
   it left off.

*********************
 Trajectory Analysis
*********************

Process molecular dynamics trajectories and extract structures for further analysis:

.. code:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] traj -j <jobtype> [SUBCMD_OPTIONS]

Trajectory Options
==================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-j, --jobtype``
      -  string
      -  Job type: opt, ts, modred, scan, sp (required)

   -  -  ``-N, --num-structures-to-run``
      -  int
      -  Number of unique structures to process

   -  -  ``-x, --proportion-structures-to-use``
      -  float
      -  Proportion of trajectory to use (0.0 < x <= 1.0, default: 0.1)

.. warning::

   The ``-j, --jobtype`` option is required for trajectory analysis.

Basic Usage
===========

Process trajectory structures:

.. code:: bash

   chemsmart sub gaussian -p project -f trajectory.xyz -c 0 -m 1 traj -j opt

Use the last 30% of trajectory:

.. code:: bash

   chemsmart sub gaussian -p project -f md.traj -c 0 -m 1 traj -j opt -x 0.3

*******************************************
 Grouper Options for CREST/Trajectory Jobs
*******************************************

Apply molecular similarity-based grouping to filter structures:

.. code:: bash

   chemsmart sub gaussian [GAUSSIAN_OPTIONS] crest/traj -j <jobtype> -g <strategy> [OPTIONS]

Grouper Options
===============

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-g, --grouping-strategy``
      -  string
      -  Strategy: rmsd, tanimoto, formula, isomorphism, connectivity

   -  -  ``-i, --ignore-hydrogens``
      -  bool
      -  Ignore H atoms (only for rmsd, default: disabled)

   -  -  ``-t, --threshold``
      -  float
      -  Grouping threshold (default: 0.5 for rmsd, 0.9 for tanimoto)

   -  -  ``-p, --num-procs``
      -  int
      -  Number of processors (default: 1)

.. note::

   No grouper is used by default. Specify ``-g`` to enable grouping.

Basic Usage
===========

CREST with RMSD grouping:

.. code:: bash

   chemsmart sub gaussian -p project -f crest_conformers.xyz -c 0 -m 1 crest -j opt -g rmsd -t 1 -p 4

Trajectory with Tanimoto grouping:

.. code:: bash

   chemsmart sub gaussian -p project -f trajectory.xyz -c 0 -m 1 traj -j opt -x 0.5 -g tanimoto

Examples
========

Use RMSD grouper with tight threshold and hydrogen exclusion:

.. code:: bash

   chemsmart run gaussian -p local -f crest_conformers.xyz -l grouped -c 0 -m 1 crest -j opt -g rmsd -t 0.2 -p 4 -i

Output files are named ``grouped_opt_c1.com``, ``grouped_opt_c2.com``, etc.
