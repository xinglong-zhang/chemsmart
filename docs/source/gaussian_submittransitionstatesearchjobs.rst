Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to
contact our team if you have questions or feedback.

#####################################
 Submit Transition State Search Jobs
#####################################

ChemSmart provides comprehensive transition state search capabilities using Gaussian. This section covers transition
state optimization, IRC calculations, and potential energy surface scanning.

*******************************
 Transition State Optimization
*******************************

Run transition state optimization to find saddle points on the potential energy surface.

.. code:: console

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] ts [SUBCMD_OPTIONS]

TS-Specific OPTIONS
===================

.. list-table:: TS Frozen Atom Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-f, --freeze-atoms``
      -  string
      -  Indices of atoms to freeze for constrained optimization

TS Basic Usage
==============

**Transition state with specific charge and multiplicity**

   .. code:: console

      chemsmart sub gaussian -p ts_opt -f ts_guess.xyz -c 0 -m 1 ts

**TS optimization with frozen atoms**

   .. code:: console

      chemsmart sub gaussian -p constrained_ts -f ts_guess.xyz ts -f 2,5,8

TS Examples
===========

**Use semiempirical method for pre-TS**

   .. code:: console

      chemsmart sub -s SLURM gaussian -p ti -f radical_opening_ts2.gjf -c 0 -m 2 -s PM6 ts

   Can search ts in pm6 level first.

*********************************************
 Modified Redundant Coordinate (Modred) Jobs
*********************************************

Run modified redundant coordinate optimization for constrained geometry optimization and transition state searches using
the ``modred`` command.

.. code:: console

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] modred [SUBCMD_OPTIONS]

Modred-Specific OPTIONS
=======================

.. list-table:: Modred Job Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-c, --coordinates``
      -  string
      -  List of coordinates to be fixed for modred job (1-based indexing)

Modred Basic Usage
==================

**Defining Coordinate Constraints**

-  to submit a modredundant job with constraints on bond between atom 4 and atom 17 and on bond between atom 9 and atom
   10, do:

   .. code:: console

      chemsmart sub -s shared gaussian -p test -f input.com modred -c [[4,17]]

Modred Examples
===============

**Modred optimization**

-  The structure can be optimized while keeping the bonds involved in the transition state fixed：

   .. code:: console

      chemsmart sub -s SLURM gaussian -p ti -f I_6m_ts_guess3_new.gjf -c 0 -m 2 modred -c [[85,100],[100,101],[101,89],[89,90],[90,88],[88,85]]

**************************************************
 Intrinsic Reaction Coordinate (IRC) Calculations
**************************************************

Run IRC calculations to follow the reaction path from a transition state.

.. code:: console

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] irc [SUBCMD_OPTIONS]

IRC-Specific OPTIONS
====================

.. list-table:: IRC Job Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-fl, --flat-irc/--no-flat-irc``
      -  bool
      -  Whether to run flat IRC or not (default=False)

   -  -  ``-pt, --predictor``
      -  string
      -  Type of predictors used for IRC. Options: LQA, HPC, EulerPC, DVV, Euler (default=none)

   -  -  ``-rc, --recorrect``
      -  string
      -  Recorrection step of HPC and EulerPC IRCs. Options: Never, Always, Test (default=none)

   -  -  ``-rs, --recalc-step``
      -  int
      -  Compute the Hessian analytically every N predictor steps or every N corrector steps if N<0 (default=6)

   -  -  ``-mp, --maxpoints``
      -  int
      -  Number of points along reaction path to examine (default=512)

   -  -  ``-mc, --maxcycles``
      -  int
      -  Maximum number of steps along IRC to run (default=128)

   -  -  ``-ss, --stepsize``
      -  int
      -  Step size along reaction path, in units of 0.01 Bohr (default=20)

   -  -  ``-d, --direction``
      -  string
      -  Only run the forward or reverse IRC. Options: forward, reverse (default=None, run both directions)

IRC Basic Usage
===============

**Basic IRC calculation**:

   .. code:: console

      chemsmart sub -s shared gaussian -p test -f irc.xyz irc

***********************************
 Potential Energy Surface Scanning
***********************************

Run coordinate scanning to explore potential energy surfaces and locate transition states.

.. code:: console

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] scan [SUBCMD_OPTIONS]

.. note::

   Scanning coordinates, step size and number of steps are all required!

Scan-Specific OPTIONS
=====================

.. list-table:: Scan Job Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-c, --coordinates``
      -  string
      -  List of coordinates to be fixed for scan job. 1-indexed (default=None)

   -  -  ``-s, --step-size``
      -  float
      -  Step size of coordinates to scan (default=None)

   -  -  ``-n, --num-steps``
      -  int
      -  Number of steps to scan (default=None)

   -  -  ``-cc, --constrained-coordinates``
      -  string
      -  List of coordinates to be fixed in scan job (1-based indexing)

Scan Basic Usage
================

**Basic coordinate scan**

-  For example, to submit the PES scan job with along bond between atom 2 and atom 3 for 10 steps with 0.1Å increment
   per step:

   .. code:: console

      chemsmart sub gaussian -p pes_scan -f molecule.xyz scan -c [[2,3]] -s 0.1 -n 10

**Multi-coordinate scan**

-  To scan multiple coordinates simultaneously with different step sizes and numbers:

   .. code:: console

      chemsmart sub gaussian -p pes_scan -f molecule.xyz scan -c [[2,3],[5,6]] -s [0.1,0.15] -n [10,12]

   Here, bond 2-3 is scanned with 0.1Å steps for 10 points, and bond 5-6 with 0.15Å steps for 12 points. Note that the
   number of distances/angles/dihedrals being scanned should be the same as the number of steps/step-sizes supplied,
   such that the first degree of freedom being scanned takes the first scanning step and step-size, and so on
   sequentially.

**Coordinate scan with additional constraints**

-  To submit the PES scan job along bond between atom 2 and atom 3 while keeping bond distance between atom 5 and atom 8
   fixed, one may use

   .. code:: console

      chemsmart sub gaussian -p pes_scan -f molecule.xyz scan -c [2,3] -s 0.1 -n 10 -cc [5,8]

   where a single pair of bonding atoms is enclosed by a pair of square brackets, or

   .. code:: console

      chemsmart sub gaussian -p pes_scan -f molecule.xyz scan -c [[2,3]] -s 0.1 -n 10 -cc [[5,8]]

   where the pair of bonding atoms is enclosed as a list of list.

-  You can also specify multiple constraints (distance, angle, dihedral) as follows:

   .. code:: console

      chemsmart sub gaussian -p pes_scan -f molecule.xyz scan -c [[2,3]] -s 0.1 -n 10 -cc [[5,8],[1,4,6],[2,3,5,7]]

**Multi-coordinate scan with constraints**

-  Combine multi-coordinate scanning with constraints:

   .. code:: console

      chemsmart sub gaussian -p pes_scan -f molecule.xyz scan -c [[1,2],[3,4,5]] -s [0.1,5.0] -n [10,18] -cc [[6,7],[8,9,10]]

   Here, bond 1-2 is scanned with 0.1Å steps for 10 points and angle 3-4-5 with 5.0° steps for 18 points, while bond 6-7
   and angle 8-9-10 are kept fixed.
