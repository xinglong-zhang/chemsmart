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
QRC job directly, so that two optimization jobs can be created and run.


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

   You supply a Gaussian output (e.g. ``ts.log``) or an ORCA output (e.g. ``ts.out``) that already contains:
   (i) an optimized transition structure (first-order saddle point) and
   (ii) a frequency calculation.

2. **Select an eigenmode**

   Chemsmart extracts the vibrational eigenvectors from that output file.
   By default, it chooses mode index ``1``, which in a well-characterized TS
   corresponds to the single imaginary mode (negative curvature direction).
   You can change this using ``-m/--mode-idx``.

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
     Advanced use: if suppose you have a geometry optimization but has one imaginary frequency (a properly optimized structure should have no imaginary frequency) and you want to remove it. You may submit the job by creating initial guess structures by displacement along the imaginary frequency and subjecting the resulting guess structures for geometry optimization, via

    .. code:: console

       chemsmart sub [OPTIONS] gaussian -p <project_settings> -f <erroneous_opt_with_one_im_freq.log> qrc -m 1 -a 0.8 -j opt

This will prepare two guess structures that are obtained from displacement along the vibrational mode 1 (``-m 1``) with an amplitude of ±0.8 and subjecting those structures for geometry optimization (``-j opt``), with the intended results of successfully optimized structure with no imaginary frequency.

   - ``ts``:
     Attempt a transition-state (TS) optimization from both + and − guesses.
     This is helpful if you want to locate a TS but instead obtain a TS output file having e.g., two imaginary frequencies, in which case, you would want to displace along the second imaginary frequency (via ``-m 2`` or ``--mode-idx 2``) to create new input structures to redo the TS optimization, in hopes of removing the second imaginary frequency mode.

.. code-block:: console

   chemsmart sub gaussian -p my_project -f ts.log qrc -m 2 -a 1.2 -j ts

Here:

* ``-m 2`` selects vibrational mode #2 (1-indexed).
* ``-a 1.2`` increases the displacement magnitude to 1.2 Å in both + and −
  directions.


   - ``modred`` / ``scan``:
     Prepare constrained optimizations or relaxed scans along specified
     internal coordinates. Again, useful when there is any imaginary frequency modes that one wishes to remove.

   - ``sp``:
     Single-point energy only, no geometry optimization after displacement.

   - ``irc``:
     Launch a full IRC calculation from each displaced guess, if you want
     a rigorous intrinsic reaction path rather than the "quick" approximation.

5. **(Optional) Movie / sanity check**

   If you provide ``--nframes``, Chemsmart samples multiple evenly spaced
   points along that vibrational mode (using trigonometric sampling of a phase
   angle; ``--phase`` sets the offset). With ``--return-xyz``, it produces a
   multi-frame XYZ block that you can visualize (e.g., via ``chemsmart run mol``) to confirm that the selected
   mode actually corresponds to bond formation / cleavage, rather than e.g.
   rotation of a remote substituent.

This removes the manual steps of:
copying displaced coordinates,
hand-editing multiple Gaussian input files for the "+" and "−" directions,
and tracking filenames for each guess.


***************************
 Extension of the QRC Jobs
***************************

Given an output file with extraneous imaginary frequencies, we can use this file as input for direct job submission, aiming to remove the extra imaginary frequencies after submission of the QRC job. 
For example, suppose we have a geometry optimisation output file with one imaginary frequency, which should not be there, we can run QRC to displace along that vibration to generate new guess structures and submit them for further geometry optimization, with the desired outcome of re-optimized structures reaching a PES minimum without any imaginary frequency.
Similarly, we can do this to remove the second imaginary frequency of a TS-like structure with the correct first imaginary frequency mode, but with an extra second imaginary mode that we intend to remove, by using displacement along the second imaginary mode and creating new input structures for TS re-optimization automatically. 

Summary
=======
 
* Chemsmart QRC takes an existing TS + frequency calculation (or erroneous optimization job with spurious imaginary frequency)  and
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

