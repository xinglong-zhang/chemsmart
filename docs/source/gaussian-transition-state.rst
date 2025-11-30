#####################################
 Transition State Search (Gaussian)
#####################################

This page covers transition state optimization, IRC calculations, and
potential energy surface scanning using Gaussian.

*******************************
 Transition State Optimization
*******************************

Run transition state optimization to find saddle points on the potential
energy surface:

.. code-block:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] ts [SUBCMD_OPTIONS]

TS Options
==========

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-f, --freeze-atoms``
     - string
     - Atom indices to freeze (1-based indexing)

Basic Usage
===========

Standard TS optimization:

.. code-block:: bash

   chemsmart sub gaussian -p project -f ts_guess.xyz -c 0 -m 1 ts

TS optimization with frozen atoms:

.. code-block:: bash

   chemsmart sub gaussian -p project -f ts_guess.xyz ts -f 2,5,8

Using semiempirical method for pre-screening:

.. code-block:: bash

   chemsmart sub gaussian -p project -f ts_guess.gjf -c 0 -m 2 -s PM6 ts

*********************************************
 Modified Redundant Coordinate (Modred) Jobs
*********************************************

Run constrained geometry optimization using modified redundant coordinates:

.. code-block:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] modred [SUBCMD_OPTIONS]

Modred Options
==============

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-c, --coordinates``
     - string
     - Coordinate indices to constrain (1-based indexing)

Basic Usage
===========

Constrain a bond between atoms 4 and 17:

.. code-block:: bash

   chemsmart sub gaussian -p project -f input.com modred -c [[4,17]]

Constrain multiple bonds:

.. code-block:: bash

   chemsmart sub gaussian -p project -f input.gjf -c 0 -m 2 modred -c [[85,100],[100,101],[101,89]]

**************************************************
 Intrinsic Reaction Coordinate (IRC) Calculations
**************************************************

Run IRC calculations to follow the reaction path from a transition state:

.. code-block:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] irc [SUBCMD_OPTIONS]

IRC Options
===========

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-fl, --flat-irc/--no-flat-irc``
     - bool
     - Run flat IRC (default: disabled)
   * - ``-pt, --predictor``
     - string
     - Predictor type: LQA, HPC, EulerPC, DVV, Euler
   * - ``-rc, --recorrect``
     - string
     - Recorrection step: Never, Always, Test
   * - ``-rs, --recalc-step``
     - int
     - Hessian recalculation frequency (default: 6)
   * - ``-mp, --maxpoints``
     - int
     - Points along reaction path (default: 512)
   * - ``-mc, --maxcycles``
     - int
     - Maximum IRC steps (default: 128)
   * - ``-ss, --stepsize``
     - int
     - Step size in 0.01 Bohr units (default: 20)
   * - ``-d, --direction``
     - string
     - Direction: forward, reverse, or both (default: both)

Basic Usage
===========

Standard IRC calculation:

.. code-block:: bash

   chemsmart sub gaussian -p project -f ts.xyz irc

***********************************
 Potential Energy Surface Scanning
***********************************

Run coordinate scanning to explore potential energy surfaces:

.. code-block:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] scan [SUBCMD_OPTIONS]

.. note::

   Scanning coordinates, step size, and number of steps are all required.

Scan Options
============

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description
   * - ``-c, --coordinates``
     - string
     - Coordinates to scan (1-based indexing)
   * - ``-s, --step-size``
     - float
     - Step size for scanning
   * - ``-n, --num-steps``
     - int
     - Number of scan steps
   * - ``-cc, --constrained-coordinates``
     - string
     - Coordinates to keep fixed (1-based indexing)

Basic Usage
===========

Scan a bond distance:

.. code-block:: bash

   chemsmart sub gaussian -p project -f molecule.xyz scan -c [[2,3]] -s 0.1 -n 10

This scans the bond between atoms 2 and 3 for 10 steps of 0.1 Å each.

Multi-coordinate scan:

.. code-block:: bash

   chemsmart sub gaussian -p project -f molecule.xyz scan -c [[2,3],[5,6]] -s [0.1,0.15] -n [10,12]

Scan with constraints:

.. code-block:: bash

   # Keep bond 5-8 fixed while scanning bond 2-3
   chemsmart sub gaussian -p project -f molecule.xyz scan -c [[2,3]] -s 0.1 -n 10 -cc [[5,8]]

   # Multiple constraints
   chemsmart sub gaussian -p project -f molecule.xyz scan -c [[2,3]] -s 0.1 -n 10 -cc [[5,8],[1,4,6]]
