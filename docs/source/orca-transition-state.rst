################################
 Transition State Search (ORCA)
################################

This page covers transition state optimization, NEB calculations, IRC calculations, and coordinate scanning using ORCA.

*************************
 Transition State Search
*************************

Locate saddle points on the potential energy surface.

.. code:: bash

   chemsmart sub [OPTIONS] orca [ORCA_OPTIONS] ts [SUBCMD_OPTIONS]

TS Options
==========

.. list-table:: Hessian Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-i, --inhess/--no-inhess``
      -  flag
      -  Read Hessian file (default: disabled)

   -  -  ``-f, --inhess-filename``
      -  string
      -  Hessian filename

   -  -  ``--numhess/--no-numhess``
      -  flag
      -  Use numerical Hessian (default: disabled)

   -  -  ``-h, --hybrid-hess/--no-hybrid-hess``
      -  flag
      -  Use hybrid Hessian (default: disabled)

   -  -  ``-a, --hybrid-hess-atoms``
      -  string
      -  Atoms for hybrid Hessian (0-indexed)

   -  -  ``-s, --recalc-hess``
      -  int
      -  Hessian recalculation frequency (default: 5)

.. list-table:: Search Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-t, --trust-radius``
      -  float
      -  Trust radius for optimization

   -  -  ``-ts, --tssearch-type``
      -  string
      -  Search type: optts or scants (default: optts)

   -  -  ``-fs, --full-scan/--no-full-scan``
      -  flag
      -  Perform full scan (default: disabled)

Basic Usage
===========

Standard TS search:

.. code:: bash

   chemsmart sub orca -p project -f ts_guess.xyz ts

TS search with existing Hessian:

.. code:: bash

   chemsmart sub orca -p project -f molecule.xyz ts -i -f hessian.hess

TS search with numerical Hessian:

.. code:: bash

   chemsmart sub orca -p project -f molecule.xyz ts --numhess

ScanTS mode with full scan:

.. code:: bash

   chemsmart sub orca -p project -f molecule.xyz ts -ts scants -fs

***************************
 Nudged Elastic Band (NEB)
***************************

Find reaction pathways and transition states by optimizing a chain of molecular structures connecting reactant and
product geometries.

.. code:: bash

   chemsmart sub [OPTIONS] orca [ORCA_OPTIONS] neb [SUBCMD_OPTIONS]

NEB Options
===========

.. list-table:: Job Type Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-j, --job-type``
      -  choice
      -  NEB calculation type: NEB, NEB-CI, NEB-TS, FAST-NEB-TS, TIGHT-NEB-TS, LOOSE-NEB, ZOOM-NEB, ZOOM-NEB-CI,
         ZOOM-NEB-TS, NEB-IDPP

.. list-table:: Structure Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-e, --ending-xyzfile``
      -  string
      -  Filename of ending geometry (product)

   -  -  ``-i, --intermediate-xyzfile``
      -  string
      -  Filename of intermediate geometry (TS guess)

   -  -  ``-r, --restarting-xyzfile``
      -  string
      -  Filename of geometry for restarting calculation

.. list-table:: Optimization Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-o, --pre-optimization/--no-pre-optimization``
      -  flag
      -  Whether to optimize input geometries (default: disabled)

   -  -  ``-A, --semiempirical``
      -  choice
      -  Semiempirical method: XTB0, XTB1, XTB2

Basic Usage
===========

Standard NEB calculation:

.. code:: bash

   chemsmart sub orca -p project -f reactant.xyz -c 0 -m 1 neb -j NEB-TS -e product.xyz -A XTB2

NEB with climbing image:

.. code:: bash

   chemsmart sub orca -p project -f reactant.xyz -c 0 -m 1 neb -j NEB-CI -e product.xyz -A XTB1

NEB with intermediate guess:

.. code:: bash

   chemsmart sub orca -p project -f reactant.xyz -c 0 -m 1 neb -j NEB-CI -e product.xyz -i ts_guess.xyz

Restart from previous calculation:

.. code:: bash

   chemsmart sub orca -p project -f reactant.xyz -c 0 -m 1 neb -j NEB -r restart.allxyz

NEB with geometry pre-optimization:

.. code:: bash

   chemsmart sub orca -p project -f reactant.xyz -c 0 -m 1 neb -j NEB-TS -e product.xyz -o -A XTB2

Job Types
=========

.. list-table::
   :header-rows: 1
   :widths: 20 80

   -  -  Job Type
      -  Description
   -  -  ``NEB``
      -  Standard nudged elastic band calculation
   -  -  ``NEB-CI``
      -  NEB with climbing image to find exact saddle point
   -  -  ``NEB-TS``
      -  NEB calculation optimized for transition state search
   -  -  ``FAST-NEB-TS``
      -  Fast convergence NEB for TS search with loose criteria
   -  -  ``TIGHT-NEB-TS``
      -  Tight convergence NEB for accurate TS search
   -  -  ``LOOSE-NEB``
      -  NEB calculation with loose convergence criteria
   -  -  ``ZOOM-NEB``
      -  Variable density NEB with focus on important regions
   -  -  ``ZOOM-NEB-CI``
      -  ZOOM-NEB with climbing image
   -  -  ``ZOOM-NEB-TS``
      -  ZOOM-NEB optimized for TS search
   -  -  ``NEB-IDPP``
      -  NEB with image-dependent pair potential interpolation

Requirements
============

-  Starting geometry (``-f`` option): reactant structure
-  Ending geometry (``-e`` option): product structure
-  At minimum, NEB requires reactant and product geometries
-  For TS searches, use NEB-TS or NEB-CI job types
-  Semiempirical methods (XTB) are recommended for initial exploration

*************
 Modred Jobs
*************

Constrained geometry optimization using modified redundant coordinates.

.. code:: bash

   chemsmart sub [OPTIONS] orca [ORCA_OPTIONS] modred [SUBCMD_OPTIONS]

*************************************
 Intrinsic Reaction Coordinate (IRC)
*************************************

Trace the minimum energy pathway from a transition state.

.. code:: bash

   chemsmart sub [OPTIONS] orca [ORCA_OPTIONS] irc [SUBCMD_OPTIONS]

IRC Options
===========

.. list-table:: General Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``--maxiter``
      -  int
      -  Maximum iterations

   -  -  ``-p, --printlevel``
      -  int
      -  Print level

   -  -  ``-d, --direction``
      -  choice
      -  Direction: both, forward, backward, down

.. list-table:: Hessian Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-i, --inithess``
      -  choice
      -  Initial Hessian: read, calc_anfreq, calc_numfreq

   -  -  ``-f, --hess-filename``
      -  string
      -  Initial Hessian filename

   -  -  ``-m, --hessmode``
      -  int
      -  Hessian mode for initial displacement (default: 0)

.. list-table:: Displacement Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``--init-displ``
      -  choice
      -  Initial displacement: DE or length

   -  -  ``--scale-init-displ``
      -  float
      -  Step size from TS (default: 0.1 a.u.)

   -  -  ``--de-init-displ``
      -  float
      -  Energy difference for displacement (default: 2 mEh)

.. list-table:: Convergence Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``--tolrmsg``
      -  float
      -  RMS gradient tolerance (default: 5e-4)

   -  -  ``--tolmaxg``
      -  float
      -  Max gradient tolerance (default: 2e-3)

   -  -  ``-M, --monitor-internals/--no-monitor-internals``
      -  flag
      -  Monitor internal coordinates

   -  -  ``-I, --internal-modred``
      -  string
      -  Internal coordinates to monitor

Basic Usage
===========

Standard IRC calculation:

.. code:: bash

   chemsmart sub orca -p project -f ts.xyz irc

IRC in both directions:

.. code:: bash

   chemsmart sub orca -p project -f ts.xyz irc -d both

IRC with existing Hessian:

.. code:: bash

   chemsmart sub orca -p project -f ts.xyz irc -i read -f hessian.hess

*********************
 Coordinate Scanning
*********************

Explore potential energy surfaces by varying specific coordinates.

.. code:: bash

   chemsmart sub [OPTIONS] orca [ORCA_OPTIONS] scan [SUBCMD_OPTIONS]

Scan Options
============

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-j, --jobtype``
      -  string
      -  Job type: opt, ts, modred, scan, sp

   -  -  ``-c, --coordinates``
      -  string
      -  Coordinates to scan (1-indexed)

   -  -  ``-x, --dist-start``
      -  string
      -  Starting distance/angle

   -  -  ``-y, --dist-end``
      -  string
      -  Ending distance/angle

   -  -  ``-n, --num-steps``
      -  string
      -  Number of scan steps

   -  -  ``-cc, --constrained-coordinates``
      -  string
      -  Coordinates to keep fixed

Basic Usage
===========

Scan a bond distance:

.. code:: bash

   chemsmart sub orca -p project -f molecule.xyz scan -j scan -c [[1,2]] -x 1.0 -y 3.0 -n 20

Scan with constraints:

.. code:: bash

   chemsmart sub orca -p project -f molecule.xyz scan -j scan -c [[1,2]] -x 1.0 -y 3.0 -n 20 -cc [[5,8]]

Multi-coordinate scan:

.. code:: bash

   chemsmart sub orca -p project -f molecule.xyz scan -j scan -c [[1,2],[3,4,5]] -x [1.5,70.0] -y [3.5,85.0] -n [10,15]
