################################
 Transition State Search (ORCA)
################################

This page covers transition state optimization, IRC calculations, and coordinate scanning using ORCA.

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
