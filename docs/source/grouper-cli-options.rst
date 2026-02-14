#####################
 Grouper CLI Options
#####################

This page documents the CLI options for molecular structure grouping/clustering. Use ``chemsmart run grouper --help``
for the complete list.

*************************
 Basic Command Structure
*************************

.. code:: bash

   chemsmart run/sub [OPTIONS] grouper [GROUPER_OPTIONS] <STRATEGY> [STRATEGY_OPTIONS]

*****************
 Grouper Options
*****************

These options are shared across all grouping strategies and must be placed BEFORE the strategy subcommand.

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-f, --filenames``
      -  string
      -  Input file containing structures to group (e.g. xyz)

   -  -  ``-d, --directory``
      -  string
      -  Directory containing structure files to group

   -  -  ``-t, --filetype``
      -  string
      -  File type filter for directory processing (only support log files currently)

   -  -  ``-l, --label``
      -  string
      -  Custom output label/filename

   -  -  ``-a, --append-label``
      -  string
      -  String to append to auto-generated label

   -  -  ``-T, --threshold``
      -  float
      -  Threshold for grouping (strategy-specific defaults apply)

   -  -  ``-N, --num-groups``
      -  int
      -  Target number of groups (adaptive threshold finding)

   -  -  ``-np, --num-procs``
      -  int
      -  Number of processors for parallel calculation (default: 1)

   -  -  ``-ih, --ignore-hydrogens``
      -  flag
      -  Ignore hydrogen atoms in grouping calculations

   -  -  ``-o, --output-format``
      -  string
      -  Output format for group results (xlsx, csv, txt, default: xyz)

**********************
 Available Strategies
**********************

RMSD-based Strategies
=====================

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Strategy
      -  Default T
      -  Description

   -  -  ``rmsd``
      -  0.5 Å
      -  Simple Kabsch RMSD alignment

   -  -  ``hrmsd``
      -  0.5 Å
      -  Hungarian RMSD (optimal atom mapping)

   -  -  ``spyrmsd``
      -  0.5 Å
      -  spyrmsd library-based RMSD with symmetry handling

   -  -  ``irmsd``
      -  0.125 Å
      -  Invariant RMSD (considers molecular symmetry)

   -  -  ``pymolrmsd``
      -  0.5 Å
      -  PyMOL-based RMSD alignment

Fingerprint-based Strategies
============================

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Strategy
      -  Default T
      -  Description

   -  -  ``tfd``
      -  0.1
      -  Torsion Fingerprint Deviation

   -  -  ``tanimoto``
      -  0.9
      -  Tanimoto similarity with molecular fingerprints

Energy-based Strategies
=======================

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Strategy
      -  Default T
      -  Description

   -  -  ``energy``
      -  1.0 kcal/mol
      -  Energy-based grouping (requires energy data)

Other Strategies
================

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Strategy
      -  Default T
      -  Description

   -  -  ``formula``
      -  N/A
      -  Group by molecular formula

   -  -  ``connectivity``
      -  N/A
      -  Group by molecular connectivity/topology

   -  -  ``isomorphism``
      -  N/A
      -  Graph isomorphism-based grouping

*************
 Input Modes
*************

Single File Mode
================

Load all structures from a single multi-structure file:

.. code:: bash

   chemsmart run grouper -f conformers.xyz rmsd

Directory Mode
==============

Load structures from a directory of files (Only for Gaussian TS/Opt log files currently). Files should follow the naming
convention ``xxx_c1_xxx.log``, ``xxx_c2_xxx.log``, or ``xxx_c1.log``, ``xxx_c2.log``, etc.

.. code:: bash

   chemsmart run grouper -d . -t log rmsd

For TS/Opt log files, the program extracts Gibbs free energy (SCF + thermal correction)

**************
 Output Files
**************

After grouping, the following files are created in ``{label}_group_result/`` folder:

Excel File
==========

``{label}_{strategy}_T{threshold}.xlsx`` or ``{label}_{strategy}_N{num_groups}.xlsx``

Contains:

-  **Parameters**: Grouping parameters (threshold, num_procs, ignore_hydrogens, etc.)
-  **Matrix**: Pairwise distance/similarity matrix with conformer IDs as labels
-  **Groups**: Group assignments with member lists

Group XYZ Files
===============

``{label}_group_1.xyz``, ``{label}_group_2.xyz``, ...

Each file contains:

-  All molecules in that group, sorted by energy (lowest first)
-  Comment line with: Group number, Original_Index, Energy (Hartree)

Example comment line:

.. code:: text

   Group 1 Molecule 1 Original_Index: c3 Energy(Hartree): -126.25755080

************
 Next Steps
************

For detailed information on each grouping strategy:

-  :doc:`grouper-strategies`
-  :doc:`grouper-crest-or-traj-workflow`
