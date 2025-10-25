Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to
contact our team if you have questions or feedback.

#############################################
 Submit Quick Reaction Coordinate (QRC) Jobs
#############################################

**************
 Introduction
**************
The quick reaction coordinate (QRC) method is a lightweight alternative to the full intrinsic
reaction coordinate (IRC) calculations for linking a located transition structure (TS) to its adjacent 
minima (reactant- and product-like structures). It was introduced by Jonathan M. Goodman and M. A. Silva 
in Tetrahedron Letters, 2003, 44, 8233. 

The core idea is that one can start from a TS with a frequency calculation done and displace along the imaginary mode 
(the Hessian eigenvector with a negative eigenvalue) in the + and - directions by a small amplitude to make two perturbed geometries, which can be optimized and fall into the basin of the adjacent minimum. The energy vs a distance measure (often mass-weighted RMS displacement from the TS) can be plotted to get a quick reaction profile and to confirm connectivity.

Chemsmart is able to use a located TS with frequencies calculation file (such as Gaussian or ORCA output) and submit the
QRC job directly, so that two optimization jobs (or other job types) can be created.


.. code:: console

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] qrc [SUBCMD_OPTIONS]

QRC-Specific OPTIONS
====================

.. list-table:: QRC Job Options
   :header-rows: 1
   :widths: 30 15 55

   * - Option
     - Type
     - Description

   * - ``-j, --jobtype``
     - string
     - Gaussian job type. Options: ``["opt", "ts", "modred", "scan", "sp", "irc"]``.
       ``opt`` (default) will run standard geometry optimizations on the displaced
       structures. ``ts`` will attempt transition-state optimizations from the
       displaced guesses. ``modred`` and ``scan`` generate constrained jobs
       (see below). ``sp`` runs single-point energy calculations. ``irc``
       launches a full IRC.

   * - ``-c, --coordinates``
     - string / list
     - List of internal coordinate indices (1-indexed) to freeze or scan.
       Only used for ``jobtype=modred`` or ``jobtype=scan``.
       Example: ``-c [1,2]``.

   * - ``-s, --step-size``
     - float
     - Step size for each scan increment (in Å or radians depending on the
       coordinate). Only meaningful for ``jobtype=scan``.

   * - ``-n, --num-steps``
     - integer
     - Number of scan steps to perform for ``jobtype=scan``. Combined
       with ``--step-size`` this defines the grid.

   * - ``-m, --mode-idx``
     - integer
     - Vibrational mode index (1-indexed) to displace along. Defaults to ``1``
       (the first mode, typically the imaginary mode of the TS).

   * - ``-a, --amp``
     - float
     - Displacement amplitude in Å applied along the chosen vibrational mode
       in both the + and − directions. Defaults to ``0.5`` Å. Larger values may
       drive the structure more clearly downhill but can lead to unrealistic
       geometries.

   * - ``-N, --nframes``
     - integer
     - Number of structures to generate *along* the chosen vibrational mode
       (a "vibrational movie"). If omitted, Chemsmart only makes the two
       ± displaced endpoints for QRC jobs. If provided, you can also use
       ``--return-xyz`` to collect these frames into a multi-frame XYZ.

   * - ``-p, --phase``
     - float (radians)
     - Phase angle for sampling along the vibrational mode. The default is
       ``pi/2`` (so that ``sin(phase) = 1``), which gives maximum displacement
       for the generated frame(s). Advanced use only.

   * - ``--normalize`` / ``--no-normalize``
     - flag
     - If ``--normalize`` is set, Chemsmart rescales the raw (un-weighted)
       vibrational eigenvector so that the atom with the largest magnitude
       displacement moves by exactly 1.0 Å. After that rescaling, ``--amp``
       becomes the literal maximum per-atom displacement. Default:
       ``--no-normalize``.

   * - ``-S, --skip-completed`` / ``-R, --no-skip-completed``
     - flag
     - By default Chemsmart will try to avoid resubmitting jobs it believes
       have already successfully completed in the project directory. Use
       ``-R`` / ``--no-skip-completed`` to force reruns.

   * - ``--return-xyz`` / ``--no-return-xyz``
     - flag
     - If ``--return-xyz`` is given together with ``--nframes``, Chemsmart
       will also emit a concatenated multi-frame XYZ trajectory of the sampled
       structures for quick visualization. Default: ``--no-return-xyz``.


Workflow and Theory
===================

At a high level, a QRC job in Chemsmart automates a standard
"displace-the-TS-along-its-imaginary-mode" workflow:

1. **Read TS geometry and modes**

   You supply a Gaussian output (e.g. ``ts.log``) that already contains:
   (i) an optimized transition structure (first-order saddle point) and
   (ii) a frequency calculation.

2. **Select an eigenmode**

   Chemsmart extracts the vibrational eigenvectors from that output file.
   By default, it chooses mode index ``1``, which in a well-characterized TS
   corresponds to the single imaginary mode (negative curvature direction).
   You can override this using ``-m/--mode-idx``.

3. **Generate displaced guesses**

   Two new geometries are created by taking the TS Cartesian coordinates
   ``R_TS`` and displacing them along the chosen eigenvector ``v_mode``
   by ± ``amp``:

   .. math::

      R(\pm) = R_\mathrm{TS} \; \pm \; \mathrm{amp} \times v_\mathrm{mode}

   If ``--normalize`` is used, ``v_mode`` is first rescaled so that its largest
   single-atom displacement is 1.0 Å. Otherwise, the raw eigenvector from the
   Gaussian frequencies block is used.

4. **Launch downstream jobs**

   Each displaced structure becomes the starting point for whatever Gaussian
   job type you requested via ``-j/--jobtype``:

   - ``opt``:
     Relax the geometry downhill into the nearest minimum. This typically lands
     you in the "reactant-like" and "product-like" wells on either side of the
     TS.

   - ``ts``:
     Attempt a transition-state (TS) optimization from both + and − guesses.
     This is helpful if you suspect your current TS is not fully converged or
     if you want to test robustness of the TS.

   - ``modred`` / ``scan``:
     Prepare constrained optimizations or relaxed scans along specified
     internal coordinates. Useful for mapping approximate reaction coordinates
     or enforcing specific bond distances.

   - ``sp``:
     Single-point energy only, no geometry change.

   - ``irc``:
     Launch a full IRC calculation from each displaced guess, if you want
     a rigorous intrinsic reaction path rather than the "quick" approximation.

5. **(Optional) Movie / sanity check**

   If you provide ``--nframes``, Chemsmart samples multiple evenly spaced
   points along that vibrational mode (using trigonometric sampling of a phase
   angle; ``--phase`` sets the offset). With ``--return-xyz``, it produces a
   multi-frame XYZ block that you can visualize to confirm that the selected
   mode actually corresponds to bond formation / cleavage, rather than e.g.
   rotation of a remote substituent.

This removes the manual steps of:
copying displaced coordinates,
hand-editing multiple Gaussian input files for the "+" and "−" directions,
and tracking filenames for each guess.


Basic Usage
===========

Goal: starting from an existing Gaussian transition state ``ts.log`` (which
already contains a frequency calculation), generate ± displaced structures
along the first imaginary mode and run standard geometry optimizations.

.. code-block:: console

   chemsmart sub gaussian \
      -p my_project \
      -f ts.log \
      qrc

Explanation:

* ``sub gaussian``:
  Submit Gaussian jobs through Chemsmart.

* ``-p my_project``:
  Name of your Chemsmart project / working directory. Inputs and outputs
  will be organized there.

* ``-f ts.log``:
  Path to the Gaussian log file that holds the TS geometry and its
  vibrational analysis.

* ``qrc``:
  Activates the QRC workflow.

With no further flags:

* ``jobtype`` defaults to ``opt``.
* ``mode-idx`` defaults to ``1``.
* ``amp`` defaults to ``0.5`` Å.

Chemsmart will:

1. Displace the TS coordinates along mode 1 by +0.5 Å and −0.5 Å.
2. Write two Gaussian input files for geometry optimization (``opt``).
3. Submit both jobs under ``my_project``.

After completion, you typically obtain two relaxed structures that
approximate the reactant- and product-like minima connected by your TS.


Choosing a Different Mode or Amplitude
======================================

If you suspect the first imaginary mode mixes irrelevant motion, or you
want to push harder along a particular coordinate, you can override the
mode index and amplitude:

.. code-block:: console

   chemsmart sub gaussian \
      -p my_project \
      -f ts.log \
      qrc \
      -m 2 \
      -a 1.2

Here:

* ``-m 2`` selects vibrational mode #2 (1-indexed).
* ``-a 1.2`` increases the displacement magnitude to 1.2 Å in both + and −
  directions.

Larger amplitudes can help "kick" the structure more clearly into
separate wells, but overly large values can lead to distorted or
chemically unreasonable geometries. Always visualize the displaced
structures if you increase ``--amp`` significantly.


TS Refinement From Both Sides
=============================

Instead of downhill optimizations to minima, you can ask Chemsmart to
re-optimize the TS starting from both displaced guesses. This is useful
to check whether the TS you found is robust:

.. code-block:: console

   chemsmart sub gaussian \
      -p my_project \
      -f ts.log \
      qrc \
      -j ts \
      -m 2 \
      -a 1.5

Notes:

* ``-j ts`` changes the downstream Gaussian job type from ``opt`` to
  ``ts`` (TS optimization).
* ``-m 2`` and ``-a 1.5`` are shown here for illustration; you can
  keep the defaults if you prefer.

If both + and − guesses converge back to essentially the same TS
(geometry, energy, imaginary frequency), that is a strong indication
your TS is a proper first-order saddle and not just a shoulder on a
flatter surface.


Constrained Coordinates: ``modred`` and ``scan``
================================================

For mechanistic mapping, you may want constrained optimizations
(``modred``) or relaxed scans (``scan``) of specific internal
coordinates. This is commonly used to approximate a reaction coordinate
by controlling key bond distances or dihedrals.

``modred`` example (freeze internal coordinates during optimization):

.. code-block:: console

   chemsmart sub gaussian \
      -p my_project \
      -f ts.log \
      qrc \
      -j modred \
      -c [1,2] \
      -m 2 \
      -a 1.5

Explanation:

* ``-j modred`` requests a constrained optimization using Gaussian's
  "modredundant" style input.
* ``-c [1,2]`` is a Python-like list of internal coordinate indices
  (1-indexed) to be constrained/fixed.
* ``-m`` and ``-a`` again control which mode is used for the initial ±
  displacements.

``scan`` example (relaxed scan over one coordinate):

.. code-block:: console

   chemsmart sub gaussian \
      -p my_project \
      -f ts.log \
      qrc \
      -j scan \
      -c [1] \
      -s 0.05 \
      -n 10 \
      -m 1 \
      -a 0.5

Explanation:

* ``-j scan`` enables a coordinate scan job.
* ``-c [1]`` chooses the internal coordinate to scan.
* ``-s 0.05`` sets the scan step size for that coordinate.
* ``-n 10`` requests 10 steps.
* ``-m 1 -a 0.5`` control displacement along the selected vibrational
  mode as before.

A relaxed scan from the + and − displaced structures can generate a
coarse reaction profile suitable for mechanism discussion and SI
figures, at the cost of launching multiple Gaussian jobs in series.


Generating a Vibrational Movie
==============================

You can also generate multiple frames along the vibrational mode
purely for visualization or sanity checking. This does *not* submit
Gaussian jobs for each frame unless you explicitly do so; instead it
samples along the eigenvector and can write a multi-frame XYZ.

.. code-block:: console

   chemsmart sub gaussian \
      -p my_project \
      -f ts.log \
      qrc \
      -m 1 \
      -a 0.5 \
      -N 20 \
      --return-xyz

Here:

* ``-N 20`` tells Chemsmart to generate 20 evenly spaced frames along
  the vibrational mode by scanning a phase angle.
* ``--return-xyz`` instructs Chemsmart to output a concatenated
  multi-frame XYZ trajectory. You can save this to a file such as
  ``mode1_movie.xyz`` and open it in any molecular viewer that supports
  trajectories (VMD, IQmol, Avogadro, etc.).
* ``--phase`` (optional) can shift the starting phase if you want your
  first frame to be at some particular displacement rather than the
  maximum.


Practical Recommendations
=========================

1. Always visualize the displaced guesses
   --------------------------------------

   Before trusting any optimization, inspect the ``+amp`` and ``-amp``
   structures. If you see clearly unphysical bond lengths or severe
   steric clashes, reduce ``--amp``.

2. Use ``--normalize`` for cross-comparisons
   -----------------------------------------

   Vibrational eigenvectors from different levels of theory (different
   basis sets / solvents / functionals) can have different absolute
   scales. With ``--normalize``, the largest per-atom displacement is
   rescaled to 1.0 Å before applying ``--amp``. This makes
   ``--amp`` values more comparable across multiple TS calculations.

3. ``scan`` is slower but gives a figure
   -------------------------------------

   A relaxed scan produces an energy profile along a chosen coordinate,
   which is very convenient for supporting mechanistic claims in a
   manuscript or SI. But note that it can explode into many Gaussian
   jobs, especially if you run scans from both + and − displacements.

4. ``ts`` from ± is a robustness test
   ----------------------------------

   If both the ``+amp`` and ``-amp`` starting guesses re-converge to
   essentially the same TS (same geometry, same imaginary frequency),
   that's strong evidence you have located the correct first-order
   saddle for that step.


Summary
=======

* Chemsmart QRC takes an existing TS + frequency calculation and
  automatically:
  
  - extracts a vibrational mode,
  - displaces the TS geometry by ± ``amp`` along that mode,
  - generates Gaussian inputs for the requested job type,
  - and submits those jobs for you.

* You control:
  
  - which mode (``-m/--mode-idx``),
  - how hard you push (``-a/--amp``),
  - what downstream job to run (``-j/--jobtype``),
  - whether to apply constraints or scans (``-c``, ``-s``, ``-n``),
  - and whether to generate a multi-frame XYZ preview (``-N``,
    ``--return-xyz``).

* This gives you:
  
  - A fast connectivity check between a TS and nearby minima,
  - approximate reaction profiles without a full IRC,
  - constrained scans suitable for mechanistic discussion,
  - and ready-to-visualize mode animations.

Future sections can document:
``Project layout and output files``, i.e. how Chemsmart names the
``+`` and ``-`` jobs in ``-p my_project``, and
``Troubleshooting``, e.g. what to do if Gaussian fails to converge or
if the displaced guess is too distorted.

***************************
 Extension of the QRC Jobs
***************************

Given an output file with extraneous imaginary frequencies, we can use this file as input for direct job submission, aiming to remove the extra imaginary frequencies after submission of the QRC job. 
For example, suppose we have a geometry optimisation output file with one imaginary frequency, which should not be there, we can run QRC to displace along the vibration and submit it for further geometry optimization, with the desired outcome of re-optimized structures reaching a PES minimum without any imaginary frequency.
