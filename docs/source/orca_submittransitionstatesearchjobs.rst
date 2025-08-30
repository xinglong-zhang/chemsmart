Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to
contact our team if you have questions or feedback.

#####################################
 Submit Transition State Search Jobs
#####################################

ChemSmart provides comprehensive tools for transition state searches and reaction pathway analysis using ORCA. This
section covers transition state optimization, intrinsic reaction coordinate (IRC) calculations, and coordinate scanning
workflows.

*************************
 Transition State Search
*************************

Transition state (TS) search is used to locate saddle points on the potential energy surface that correspond to
transition states in chemical reactions.

.. code:: console

   chemsmart sub [OPTIONS] orca [ORCA_OPTIONS] ts [SUBCMD_OPTIONS]

TS-Specific OPTIONS
===================

.. list-table:: Hessian Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-i, --inhess/--no-inhess``
      -  flag
      -  Option to read in Hessian file (default=False)

   -  -  ``-f, --inhess-filename <string>``
      -  string
      -  Filename of Hessian file (default=None)

   -  -  ``--numhess/--no-numhess``
      -  flag
      -  Option to use numerical Hessian (default=False)

   -  -  ``-h, --hybrid-hess/--no-hybrid-hess``
      -  flag
      -  Option to use hybrid Hessian (default=False)

   -  -  ``-a, --hybrid-hess-atoms <string>``
      -  string
      -  List of atoms to use for hybrid Hessian. Zero-indexed, e.g. [0, 1, 2, 3] (default=None)

   -  -  ``-s, --recalc-hess <int>``
      -  int
      -  Number of steps to recalculate Hessian (default=5)

.. list-table:: Search Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-t, --trust-radius <float>``
      -  float
      -  Trust radius for TS optimization (default=None)
<<<<<<< HEAD

   -  -  ``-ts, --tssearch-type <string>``
      -  string
      -  Type of TS search to perform. Options: optts, scants (default="optts")

   -  -  ``-fs, --full-scan/--no-full-scan``
      -  flag
      -  Option to perform a full scan (default=False)

**Search Modes:** - **OptTS mode (optts)**: Direct transition state optimization - **ScanTS mode (scants)**: Transition
state search through scanning

Basic Usage
===========
=======

   -  -  ``-ts, --tssearch-type <string>``
      -  string
      -  Type of TS search to perform. Options: optts, scants (default="optts")

   -  -  ``-fs, --full-scan/--no-full-scan``
      -  flag
      -  Option to perform a full scan (default=False)

**Search Modes:** - **OptTS mode (optts)**: Direct transition state optimization - **ScanTS mode (scants)**: Transition
state search through scanning

TS Basic Usage
==============
>>>>>>> 8a3e0a1 (Docs setup (#263))

**Basic transition state search**:

   .. code:: console

      chemsmart sub orca -p project_name -f input.xyz ts

**TS search with existing Hessian**:

   .. code:: console

      chemsmart sub orca -p ts_hess -f molecule.xyz ts -i -f hessian.hess

**TS search with numerical Hessian**:

   .. code:: console

      chemsmart sub orca -p ts_numhess -f molecule.xyz ts --numhess

**ScanTS mode with full scan**:

   .. code:: console

      chemsmart sub orca -p scan_ts -f molecule.xyz ts -ts scants -fs

*************
 Modred jobs
*************

Modred-Specific OPTIONS
=======================

Jobtype mission

<<<<<<< HEAD
Basic Usage
===========
=======
Modred Basic Usage
==================
>>>>>>> 8a3e0a1 (Docs setup (#263))

Whats the Difference TODO

*************************************
 Intrinsic Reaction Coordinate (IRC)
*************************************

IRC calculations trace the minimum energy pathway from a transition state to reactants and products.

Run "chemsmart sub [GENERAL_OPTIONS] orca [ORCA_OPTIONS] irc [OPTIONS]" to perform IRC calculations.

IRC-Specific OPTIONS
====================

.. list-table:: General IRC Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``--maxiter <int>``
      -  int
      -  Maximum number of iterations (default=None)

   -  -  ``-p, --printlevel <int>``
      -  int
      -  Print level (default=None)

   -  -  ``-d, --direction <choice>``
      -  choice
      -  IRC direction. Options: both, forward, backward, down (default=None)

.. list-table:: Hessian and Initial Settings
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-i, --inithess <choice>``
      -  choice
      -  Initial Hessian. Options: read, calc_anfreq, calc_numfreq (default=None)

   -  -  ``-f, --hess-filename <string>``
      -  string
      -  Filename of initial Hessian (default=None)

   -  -  ``-m, --hessmode <int>``
      -  int
      -  Hessian mode used for the initial displacement. Default 0 (default=None)

.. list-table:: Displacement Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``--init-displ <choice>``
      -  choice
      -  Initial displacement. Options: DE, length. DE for energy difference, length for step size (default=None)

   -  -  ``--scale-init-displ <float>``
      -  float
      -  Step size for initial displacement from TS. Default 0.1 a.u. (default=None)

   -  -  ``--de-init-displ <float>``
      -  float
      -  Energy difference for initial displacement based on provided Hessian. Default: 2 mEh (default=None)

   -  -  ``--scale-displ-sd <float>``
      -  float
      -  Scaling factor for scaling the 1st SD step. Default to 0.15 (default=None)

   -  -  ``--adapt-scale-displ/--no-adapt-scale-displ``
      -  flag
      -  Modify Scale_Displ_SD when the step size becomes smaller or larger (default=False)

.. list-table:: Steepest Descent Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``--sd-parabolicfit/--no-sd-parabolicfit``
      -  flag
      -  Do a parabolic fit for finding an optimal SD step length (default=False)

   -  -  ``--interpolate-only/--no-interpolate-only``
      -  flag
      -  Only allow interpolation for parabolic fit, not extrapolation (default=False)

   -  -  ``--do-sd-corr/--no-do-sd-corr``
      -  flag
      -  Do SD correction to 1st step (default=False)

   -  -  ``--scale-displ-sd-corr``
      -  float
      -  Scaling factor for scaling the correction step to the SD step (default=None)

   -  -  ``--sd-corr-parabolicfit/--no-sd-corr-parabolicfit``
      -  flag
      -  Do a parabolic fit for finding an optimal correction step length (default=False)

.. list-table:: Convergence and Monitoring
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description
<<<<<<< HEAD

   -  -  ``--tolrmsg``
      -  float
      -  Tolerance for RMS gradient (a.u.). Default 5.e-4 (default=None)

   -  -  ``--tolmaxg``
      -  float
      -  Tolerance for maximum gradient (a.u.). Default 2.e-3 (default=None)

   -  -  ``-M, --monitor-internals/--no-monitor-internals``
      -  flag
      -  Monitor internals to print out up to three internal coordinates (default=False)

   -  -  ``-I, --internal-modred``
      -  string
      -  Internal modred. Up to three internal coordinates can be defined and values printed (default=None)

   -  -  ``--follow-coordtype``
      -  string
      -  Follow coordinate type. Default cartesian. The only option (default=None)

Basic Usage
===========
=======

   -  -  ``--tolrmsg``
      -  float
      -  Tolerance for RMS gradient (a.u.). Default 5.e-4 (default=None)

   -  -  ``--tolmaxg``
      -  float
      -  Tolerance for maximum gradient (a.u.). Default 2.e-3 (default=None)

   -  -  ``-M, --monitor-internals/--no-monitor-internals``
      -  flag
      -  Monitor internals to print out up to three internal coordinates (default=False)

   -  -  ``-I, --internal-modred``
      -  string
      -  Internal modred. Up to three internal coordinates can be defined and values printed (default=None)

   -  -  ``--follow-coordtype``
      -  string
      -  Follow coordinate type. Default cartesian. The only option (default=None)

IRC Basic Usage
===============
>>>>>>> 8a3e0a1 (Docs setup (#263))

**Basic IRC calculation**:

   .. code:: console

      chemsmart sub orca -p project_name -f ts_structure.xyz irc

**IRC in both directions**:

   .. code:: console

      chemsmart sub orca -p irc_both -f ts.xyz irc -d both

**IRC with existing Hessian**:

   .. code:: console

      chemsmart sub orca -p irc_hess -f ts.xyz irc -i read -f hessian.hess

**IRC with monitoring internal coordinates**:

   .. code:: console

      chemsmart sub orca -p irc_monitor -f ts.xyz irc -M -I [[1,2,3,4],[2,3,4,5]]

*********************
 Coordinate Scanning
*********************

Coordinate scanning performs a systematic exploration of the potential energy surface by varying specific coordinates.

.. code:: console

   chemsmart sub [OPTIONS] orca [ORCA_OPTIONS] scan [SUBCMD_OPTIONS]

Scan-Specific OPTIONS
=====================

.. list-table:: Scan Job Options (Required)
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description
<<<<<<< HEAD

   -  -  ``-j, --jobtype``
      -  string
      -  ORCA job type. Options: opt, ts, modred, scan, sp (default=None)

   -  -  ``-c, --coordinates``
      -  string
      -  List of coordinates to be fixed for modred or scan job. 1-indexed (default=None)

   -  -  ``-x, --dist-start``
      -  string
      -  Starting distance to scan, in Angstroms (default=None)

   -  -  ``-y, --dist-end``
      -  string
      -  Ending distance to scan, in Angstroms (default=None)

   -  -  ``-n, --num-steps``
      -  string
      -  Number of steps for coordinate scanning (default=None)

Basic Usage
===========
=======

   -  -  ``-j, --jobtype``
      -  string
      -  ORCA job type. Options: opt, ts, modred, scan, sp (default=None)

   -  -  ``-c, --coordinates``
      -  string
      -  List of coordinates to be fixed for modred or scan job. 1-indexed (default=None)

   -  -  ``-x, --dist-start``
      -  string
      -  Starting distance to scan, in Angstroms (default=None)

   -  -  ``-y, --dist-end``
      -  string
      -  Ending distance to scan, in Angstroms (default=None)

   -  -  ``-n, --num-steps``
      -  string
      -  Number of steps for coordinate scanning (default=None)

Scan Basic Usage
================
>>>>>>> 8a3e0a1 (Docs setup (#263))

**Basic distance scan**:

   .. code:: console

      chemsmart sub orca -p scan_job -f molecule.xyz scan -j scan -c [[1,2]] -x 1.0 -y 3.0 -n 20

**Bond optimization with constrained distance**:

   .. code:: console

      chemsmart sub orca -p modred_opt -f molecule.xyz scan -j modred -c [[1,2]]
