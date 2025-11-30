##################
 ORCA CLI Options
##################

This page documents the CLI options available for all ORCA jobs. Use
``chemsmart sub orca --help`` for the complete list.

*************************
 Basic Command Structure
*************************

.. code:: bash

   chemsmart sub [OPTIONS] orca [ORCA_OPTIONS] <SUBCMD> [SUBCMD_OPTIONS]

.. tip::

   ORCA options are largely similar to Gaussian options with some
   ORCA-specific parameters.

**************
 ORCA Options
**************

Project and File Options
========================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-p, --project``
      -  string
      -  Project settings from ``~/.chemsmart/orca/*.yaml``

   -  -  ``-f, --filename``
      -  string
      -  Input file for job preparation

   -  -  ``-l, --label``
      -  string
      -  Custom output filename (without extension)

   -  -  ``-a, --append-label``
      -  string
      -  String to append to the base filename

   -  -  ``-t, --title``
      -  string
      -  ORCA job title

   -  -  ``-i, --index``
      -  string
      -  Structure index (1-based, default: last structure)

   -  -  ``-P, --pubchem``
      -  string
      -  Query structure from PubChem

.. note::

   -  ``-p`` uses the project name without the ``.yaml`` extension.
   -  ``-f`` accepts various formats: ``.xyz``, ``.com``, ``.gjf``,
      ``.log``, ``.inp``, or ``.out``.

Molecular Properties Options
============================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-c, --charge``
      -  int
      -  Molecular charge

   -  -  ``-m, --multiplicity``
      -  int
      -  Molecular multiplicity

Method and Basis Set Options
============================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-A, --ab-initio``
      -  string
      -  Ab initio method (e.g., DLPNO-CCSD(T))

   -  -  ``-x, --functional``
      -  string
      -  DFT functional

   -  -  ``-D, --dispersion``
      -  string
      -  Dispersion correction

   -  -  ``-b, --basis``
      -  string
      -  Basis set

   -  -  ``-a, --aux-basis``
      -  string
      -  Auxiliary basis set

   -  -  ``-e, --extrapolation-basis``
      -  string
      -  Extrapolation basis set

SCF and Grid Options
====================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-d, --defgrid``
      -  choice
      -  Grid: defgrid1, defgrid2, defgrid3 (default: defgrid2)

   -  -  ``--scf-tol``
      -  choice
      -  SCF tolerance: NormalSCF, LooseSCF, TightSCF, etc.

   -  -  ``--scf-algorithm``
      -  choice
      -  SCF algorithm: GDIIS, DIIS, SOSCF, AutoTRAH

   -  -  ``--scf-maxiter``
      -  int
      -  Maximum SCF iterations

   -  -  ``--scf-convergence``
      -  float
      -  SCF convergence criterion

Property Calculation Options
============================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``--dipole/--no-dipole``
      -  bool
      -  Dipole moment calculation

   -  -  ``--quadrupole/--no-quadrupole``
      -  bool
      -  Quadrupole moment calculation

   -  -  ``--forces/--no-forces``
      -  bool
      -  Forces calculation (default: disabled)

MDCI Options
============

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``--mdci-cutoff``
      -  choice
      -  MDCI cutoff: loose, normal, tight

   -  -  ``--mdci-density``
      -  choice
      -  MDCI density: none, unrelaxed, relaxed

Additional Options
==================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-r, --additional-route-parameters``
      -  string
      -  Additional route parameters

***********************
 Available Subcommands
***********************

Structure Optimization
======================

.. list-table::
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``opt``
      -  Geometry optimization
   -  -  ``sp``
      -  Single point calculation

Transition State Search
=======================

.. list-table::
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``ts``
      -  Transition state optimization
   -  -  ``modred``
      -  Modified redundant coordinate optimization
   -  -  ``irc``
      -  Intrinsic reaction coordinate calculations
   -  -  ``scan``
      -  Coordinate scanning

Direct Input
============

.. list-table::
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``inp``
      -  Run ORCA input file as-is

************
 Next Steps
************

For detailed information on each job type:

-  :doc:`orca-structure-optimization`
-  :doc:`orca-transition-state`
-  :doc:`orca-direct-input`
