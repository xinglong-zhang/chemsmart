######################################
 Quick Reaction Coordinate (QRC) Jobs
######################################

The Quick Reaction Coordinate (QRC) method is a lightweight alternative to full Intrinsic Reaction Coordinate (IRC)
calculations for linking a transition state (TS) to its adjacent minima.

**************
 Introduction
**************

The QRC method was introduced by Jonathan M. Goodman and M. A. Silva in *Tetrahedron Letters*, 2003, 44, 8233.

The core idea:

#. Start from a TS with a completed frequency calculation.
#. Displace along the imaginary mode (negative eigenvalue direction) in both + and − directions by a small amplitude.
#. Optimize each displaced geometry to find the adjacent minima.
#. Plot energy vs. displacement for a quick reaction profile.

Chemsmart automates this process, taking a TS frequency calculation output and directly submitting QRC jobs.

.. code:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] qrc [SUBCMD_OPTIONS]

*************
 QRC Options
*************

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-j, --jobtype``
      -  string
      -  Job type: opt, ts, modred, scan, sp, irc (default: opt)

   -  -  ``-c, --coordinates``
      -  string
      -  Coordinates to freeze/scan (1-indexed, for modred/scan)

   -  -  ``-s, --step-size``
      -  float
      -  Scan step size in Å or radians (for scan)

   -  -  ``-n, --num-steps``
      -  int
      -  Number of scan steps (for scan)

   -  -  ``-m, --mode-idx``
      -  int
      -  Vibrational mode index, 1-indexed (default: 1)

   -  -  ``-a, --amp``
      -  float
      -  Displacement amplitude in Å (default: 0.5)

   -  -  ``-N, --nframes``
      -  int
      -  Number of frames for vibrational movie

   -  -  ``-p, --phase``
      -  float
      -  Phase angle in radians (default: π/2)

   -  -  ``--normalize/--no-normalize``
      -  flag
      -  Normalize eigenvector to 1.0 Å max displacement

   -  -  ``--return-xyz/--no-return-xyz``
      -  flag
      -  Output multi-frame XYZ with ``--nframes``

*******************
 Workflow Overview
*******************

#. **Read TS geometry and modes**

   Supply a Gaussian output (``ts.log``) or ORCA output (``ts.out``) containing an optimized TS and frequency
   calculation.

#. **Select eigenmode**

   By default, mode 1 is used (the imaginary mode). Use ``-m`` to select a different mode.

#. **Generate displaced guesses**

   Two geometries are created:

   .. math::

      R(\pm) = R_\mathrm{TS} \pm \mathrm{amp} \times v_\mathrm{mode}

   If ``--normalize`` is used, the eigenvector is rescaled so the largest single-atom displacement is 1.0 Å.

#. **Launch downstream jobs**

   Each displaced structure is used for the specified job type:

   -  ``opt``: Optimize to find adjacent minima
   -  ``ts``: TS optimization (for removing extra imaginary frequencies)
   -  ``modred``/``scan``: Constrained optimizations or scans
   -  ``sp``: Single-point energy calculations
   -  ``irc``: Full IRC calculations

#. **(Optional) Generate movie**

   With ``--nframes``, sample multiple points along the vibrational mode. Use ``--return-xyz`` to output a multi-frame
   XYZ file.

*************
 Basic Usage
*************

Standard QRC to find reactant and product minima:

.. code:: bash

   chemsmart sub gaussian -p project -f ts.log qrc -j opt

Remove an extra imaginary frequency from a TS:

.. code:: bash

   chemsmart sub gaussian -p project -f ts.log qrc -m 2 -a 1.2 -j ts

Generate a vibrational movie:

.. code:: bash

   chemsmart sub gaussian -p project -f ts.log qrc -N 10 --return-xyz

***************************
 Extended QRC Applications
***************************

**Removing spurious imaginary frequencies from geometry optimizations:**

If an optimization produces an unwanted imaginary frequency, use QRC to displace along that mode and re-optimize:

.. code:: bash

   chemsmart sub gaussian -p project -f opt_with_imag.log qrc -m 1 -a 0.8 -j opt

**Removing second imaginary frequency from TS:**

If a TS has two imaginary frequencies, displace along the second mode:

.. code:: bash

   chemsmart sub gaussian -p project -f ts.log qrc -m 2 -a 1.2 -j ts

*********
 Summary
*********

Chemsmart QRC automates the workflow of:

-  Extracting vibrational modes from TS calculations
-  Displacing the TS geometry by ± amplitude along a mode
-  Generating Gaussian inputs for downstream jobs
-  Submitting those jobs automatically

This provides:

-  Fast connectivity checks between TS and minima
-  Approximate reaction profiles without full IRC
-  Tools for removing spurious imaginary frequencies
-  Ready-to-visualize mode animations
