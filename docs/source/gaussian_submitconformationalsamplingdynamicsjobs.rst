Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to
contact our team if you have questions or feedback.

################################################
 Submit Conformational Sampling & Dynamics Jobs
################################################

ChemSmart provides powerful tools for conformational sampling and molecular dynamics calculations using Gaussian. This
section covers CREST conformational searches and trajectory analysis workflows.

************
 CREST jobs
************

CREST (Conformer-Rotamer Ensemble Sampling Tool) is used for systematic conformational searches to find low-energy
conformers. In chemsmart, the crest job allow you to combine the crest program with gaussian.

.. code:: console

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] crest -j <jobtype> [SUBCMD_OPTIONS]

CREST-Specific OPTIONS
======================

.. list-table:: CREST Job Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-j, --jobtype``
      -  string
      -  Gaussian job type. Options: ["opt", "ts", "modred", "scan", "sp"]

   -  -  ``-N, --num-confs-to-run``
      -  int
      -  Number of conformers to optimize

.. warning::

   ``-j, --Jobtype`` must be provided for Crest job!

Crest Basic Usage
=================

**Basic CREST-Gaussian jobs**

-  To run opt or modred or ts conformers from CREST run output, do:

   .. code:: console

      chemsmart sub -s shared gaussian -p test-f <input_file> -c 1 -m 1 crest -j opt

   or

   .. code:: console

      chemsmart sub -s shared gaussian -p test -f <input_file> -c 0 -m 2 crest -j modred -c [1,4]

   or

   .. code:: console

      chemsmart sub -s shared gaussian -p test -f <input_file> -c 0 -m 1 crest -j ts

   respectively

.. note::

   Typically, the ``<input_file>`` is a list of all conformers obtained by CREST program and named
   ``crest_conformers.xyz``.

Crest Examples
==============

**Use CREST-opt job for 10 conformers**

-  use crest_conformers.xyz file to run the opt job for 10 conformers from the lowest energy conformer:

   .. code:: console

      chemsmart sub gaussian -p ti -f crest_conformers.xyz -l RR_T4_from_lowest -c 0 -m 1 crest -j opt -N 10

   Output files will be saved as *RR_T4_from_lowest_opt_c1* to *RR_T4_from_lowest_opt_c10*.

.. note::

   If the job terminates before ``<n_conformers_to_opt>`` are all optimized, perhaps due to walltime limit, resubmitting
   the job will continue crest opt job until all ``<n_conformers_to_opt>`` are optimized. Charge and multiplicity need
   to be specified.

*********************
 Trajectory Analysis
*********************

Trajectory analysis allows you to process molecular dynamics trajectories and extract specific structures for further
analysis.

.. code:: console

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] traj -j <jobtype> [SUBCMD_OPTIONS]

Trajectory-Specific OPTIONS
===========================

.. list-table:: Trajectory Job Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-j, --jobtype``
      -  string
      -  Gaussian job type. Options: ["opt", "ts", "modred", "scan", "sp"]

   -  -  ``-N, --num-structures-to-run``
      -  int
      -  Number of structures from the list of unique structures to run the job on

   -  -  ``-x, --proportion-structures-to-use``
      -  float
      -  Proportion of structures from the end of trajectory to use. Values ranges from 0.0 < x <= 1.0. Defaults to 0.1
         (last 10% of structures)

.. warning::

   ``-j, --Jobtype`` must be provided for Trajectory analysis!

Trajectory Basic Usage
======================

**Basic trajectory analysis**

   .. code:: console

      chemsmart sub gaussian -p trajectory -f trajectory.xyz -c 0 -m 1 traj -j opt

**Trajectory analysis with specific proportion of structures**

-  to consider the last 30% of the structures in md.traj trajectory file:

   .. code:: console

      chemsmart sub gaussian -p traj_analysis -f md.traj -c 0 -m 1 traj -j opt -x 0.3

***********************************************
 Additional Grouper Option for crest/traj Jobs
***********************************************

Process the results of the crest/traj task further using multiple molecular similarity-based grouping strategies.

.. code:: console

   chemsmart sub gaussian [GAUSSIAN OPTIONS] crest/traj -j <jobtype> -g <> [SUBCMD_OPTIONS]

Grouper-Specific OPTIONS
========================

.. list-table:: Grouper Job Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-g, --grouping-strategy``
      -  string
      -  Grouping strategy to use for grouping. Options: "rmsd", "tanimoto", "formula", "isomorphism", "connectivity"

   -  -  ``-i, --ignore-hydrogens``
      -  bool
      -  Ignore H atoms in the grouping (Default = False, only for "rmsd")

   -  -  ``-t, --threshold``
      -  float
      -  Threshold value for grouping (Default = 0.5 for "rmsd", 0.9 for "tanimoto" and 0.0 for "connectivity")

   -  -  ``-p, --num-procs``
      -  int
      -  Number of processors to use for grouper (Default=1)

Grouper Basic Usage
===================

**Basic grouping for crest job**

   .. code:: console

      chemsmart sub gaussian -p test -f crest_conformers.xyz -c 0 -m 1 crest -j opt -g rmsd -t 1 -p 4

**basic grouping for traj job**

   .. code:: console

      chemsmart sub gaussian -p traj_test -f trajectory.xyz -c 0 -m 1 traj -j opt -x 0.5 -g tanimoto

Grouper Examples
================

**Use RMSD grouper for crest job**

-  To run crest job in local and save the output files with "grouped" label, tight threshold is used:

   .. code:: console

      chemsmart run gaussian -p local -f crest_conformers.xyz -l grouped -c 0 -m 1 crest -j opt -g rmsd -t 0.2 -p 4 -i

   Output files will be saved as *grouped_opt_c1.com, grouped_opt_c1.log, ..., grouped_opt_cN.com, grouped_opt_cN.log*
