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

   -  -  ``-p, --program``
      -  string
      -  File type filter for directory processing (for ``gaussian`` and ``orca``)

   -  -  ``-t, --filetype``
      -  string
      -  File type filter for directory processing (e.g. log)

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

   -  -  ``-m, --matrix-format``
      -  string
      -  Output format for group results (xlsx, csv, txt, default: xyz)

   -  -  ``-E, --energy-type``
      -  string
      -  Energy type to extract (E, H, G, qhH, qhG, sp_qhG)

   -  -  ``-csg, --cutoff-entropy-grimme``
      -  float
      -  Cutoff frequency for entropy (cm⁻¹) using Grimme's qRRHO method (Default: 100)

   -  -  ``-cst, --cutoff-entropy-truhlar``
      -  float
      -  Cutoff frequency for entropy (cm⁻¹) using Truhlar's qRRHO method

   -  -  ``-ch, --cutoff-enthalpy``
      -  float
      -  Cutoff frequency for enthalpy (cm⁻¹), Head-Gordon qRRHO method (Default: 100)

   -  -  ``-c, --concentration``
      -  float
      -  Concentration in mol/L (Default: 1.0)

   -  -  ``-P, --pressure``
      -  float
      -  Pressure in atm (Default: 1.0)

   -  -  ``-temp, --temperature``
      -  float
      -  Temperature in Kelvin (Default: 298.15)

   -  -  ``--alpha``
      -  int
      -  Interpolator exponent used in the qRRHO approximation (Default: 4)

   -  -  ``--weighted/--no-weighted``
      -  flag
      -  Use natural abundance weighted masses vs most abundant masses (Default: weighted)

   -  -  ``--energy-units``
      -  string
      -  Units of energetic values (hartree, eV, kcal/mol, kJ/mol. Default: hartree)

   -  -  ``--check-imaginary-frequencies``
      -  flag
      -  Whether to check for imaginary frequencies (Default: True)

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

Load structures from a directory of Gaussian or ORCA output files. The program extracts the last structure from each
file along with the Gibbs free energy.

**Supported file types:**

-  ``gaussian``: Gaussian output files (``.log``, ``.out``)
-  ``orca``: ORCA output files (``.out``)

**File naming:**

Files with conformer pattern (``xxx_c1_xxx.log``, ``xxx_c2.log``, etc.) will use ``c1``, ``c2``, ... as conformer IDs
and be sorted numerically. Files without this pattern will use the filename (without extension) as the conformer ID and
be sorted alphabetically after the numbered conformers.

.. code:: bash

   chemsmart run grouper -d . -p gaussian rmsd
   chemsmart run grouper -d . -p orca rmsd
   chemsmart run grouper -d . -t log rmsd
   chemsmart run grouper -d . -t out rmsd

**Validation:**

-  Files must have normal termination
-  Frequency calculation is required (for Gibbs energy extraction)
-  For optimization jobs: no imaginary frequencies allowed
-  For TS jobs: exactly one imaginary frequency required

**Energy extraction:**

For Gaussian and ORCA output files, the program can extract various types of energies directly from the output using the
``-E`` or ``--energy-type`` option.

Supported energy types:

-  ``E``: SCF energy (default)
-  ``H``: Enthalpy
-  ``G``: Gibbs free energy
-  ``qhH``: Quasi-harmonic enthalpy (requires frequency calculation)
-  ``qhG``: Quasi-harmonic Gibbs free energy (requires frequency calculation)
-  ``sp_qhG``: Single-point corrected quasi-harmonic Gibbs free energy (qhG - Egas + Esolv)

**Thermochemistry corrections:**

When using ``qhH``, ``qhG``, or ``sp_qhG``, chemsmart leverages the internal thermochemistry module to compute corrected
thermodynamic values. Additional thermochemistry CLI options (e.g., temperature, concentration, frequency cutoffs,
entropy method) become available.

For ``sp_qhG``, the program will automatically search for a matching single-point file in a subfolder named 'sp'
(containing "sp" in its name) to extract the energy in solvation or using better basis and compute the corrected free
energy.

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
