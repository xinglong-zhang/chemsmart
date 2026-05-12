#####################
 Other Gaussian Jobs
#####################

This page covers additional Gaussian job types including multi-step link jobs, custom user-defined calculations, and
direct input file execution.

****************************************
 Minimum Energy Cross Point (MECP) Jobs
****************************************

A Minimum Energy Crossing Point (MECP) is a geometry where two potential energy surfaces of different spin multiplicity
are degenerate. It is the spin-forbidden analogue of a transition state and is relevant to intersystem crossing,
spin-state reactivity, and organometallic reaction mechanisms.

ChemSmart provides a **self-contained MECP optimizer** that drives the search entirely in Python, calling Gaussian only
for single-point energies and Cartesian forces at each step. No external MECP code (e.g. MECP2, Harvey's code) is
required.

.. code:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] mecp [SUBCMD_OPTIONS]

Algorithm
=========

At each optimization step the following operations are performed:

#. **State A SP+forces** — Gaussian single-point energy and Cartesian forces for spin state A at the current geometry.

#. **State B SP+forces** — Same for spin state B.

#. **Effective displacement** — The displacement vector is computed using the penalty-function / projected-gradient
   approach:

   -  Let :math:`\Delta E = E_A - E_B` and :math:`\mathbf{g}_\Delta = \nabla E_A - \nabla E_B`.

   -  **Seam correction** (drives geometry toward the crossing seam): :math:`\mathbf{d}_\text{seam} = -\frac{\Delta
      E}{\|\mathbf{g}_\Delta\|^2}\,\mathbf{g}_\Delta`

   -  **Seam-tangent gradient** (component of :math:`\nabla E_A` perpendicular to :math:`\mathbf{g}_\Delta`):
      :math:`\mathbf{g}_\perp = \nabla E_A - \frac{\nabla E_A \cdot
      \mathbf{g}_\Delta}{\|\mathbf{g}_\Delta\|^2}\,\mathbf{g}_\Delta`

   -  **Downhill step** (minimizes energy along the seam): :math:`\mathbf{d}_\text{seam\text{-}min} =
      -\alpha\,\mathbf{g}_\perp` where :math:`\alpha` is ``step_size`` (Bohr²/Hartree).

   -  **Total displacement**: :math:`\mathbf{d} = \mathbf{d}_\text{seam} + \mathbf{d}_\text{seam\text{-}min}`

#. **Trust radius** — Each atom's displacement vector is independently scaled so that its Cartesian norm does not exceed
   ``trust_radius`` (Bohr).

#. **Convergence check** — See the convergence criteria table below.

#. **Geometry update** — :math:`\mathbf{r}_{n+1} = \mathbf{r}_n + \mathbf{d}`.

MECP Options
============

.. list-table::
   :header-rows: 1
   :widths: 30 10 20 40

   -  -  Option
      -  Type
      -  Default
      -  Description

   -  -  ``--multiplicity-a``
      -  int
      -  1 (singlet)
      -  Spin multiplicity of state A. Falls back to ``-m`` if not set.

   -  -  ``--multiplicity-b``
      -  int
      -  multiplicity-a + 2
      -  Spin multiplicity of state B.

   -  -  ``--charge-a``
      -  int
      -  0
      -  Charge of state A. Falls back to ``-c`` if not set.

   -  -  ``--charge-b``
      -  int
      -  charge-a
      -  Charge of state B.

   -  -  ``--title-a``
      -  string
      -  ``"First"``
      -  Title prefix for state A Gaussian jobs.

   -  -  ``--title-b``
      -  string
      -  ``"Second"``
      -  Title prefix for state B Gaussian jobs.

   -  -  ``--max-steps``
      -  int
      -  500
      -  Maximum number of MECP optimization steps.

   -  -  ``--step-size``
      -  float
      -  0.1 Bohr²/Hartree
      -  Scaling factor :math:`\alpha` for the seam-tangent gradient step.

   -  -  ``--convergence``
      -  string
      -  ``"standard"``
      -  Convergence preset: ``"standard"`` (default, general-purpose) or ``"tight"``
         (publication-quality). See the :ref:`convergence-presets` section below.
         Individual tolerance options override the preset.

   -  -  ``--trust-radius``
      -  float
      -  0.3 Bohr (standard) / 0.1 Bohr (tight)
      -  Maximum per-atom Cartesian displacement magnitude (Bohr) per step. Overrides preset.

   -  -  ``--energy-diff-tol``
      -  float
      -  5.0×10⁻⁵ Ha (standard) / 1.0×10⁻⁵ Ha (tight)
      -  Convergence threshold for :math:`|E_A - E_B|`. Overrides preset.

   -  -  ``--force-max-tol``
      -  float
      -  7.0×10⁻⁴ Ha/Bohr (standard) / 3.0×10⁻⁴ Ha/Bohr (tight)
      -  Convergence threshold for :math:`\max|\mathbf{g}_\perp|`. Overrides preset.

   -  -  ``--force-rms-tol``
      -  float
      -  5.0×10⁻⁴ Ha/Bohr (standard) / 1.0×10⁻⁴ Ha/Bohr (tight)
      -  Convergence threshold for :math:`\text{RMS}(\mathbf{g}_\perp)`. Overrides preset.

   -  -  ``--disp-max-tol``
      -  float
      -  4.0×10⁻³ Bohr (standard) / 2.0×10⁻³ Bohr (tight)
      -  Convergence threshold for :math:`\max|\mathbf{d}|`. Overrides preset.

   -  -  ``--disp-rms-tol``
      -  float
      -  2.5×10⁻³ Bohr (standard) / 1.0×10⁻³ Bohr (tight)
      -  Convergence threshold for :math:`\text{RMS}(\mathbf{d})`. Overrides preset.

   -  -  ``--adaptive-step-size / --no-adaptive-step-size``
      -  bool
      -  True
      -  Enable automatic step size adaptation each iteration (default: enabled).

   -  -  ``--step-size-method``
      -  string
      -  ``"bb"``
      -  Step size adaptation algorithm: ``"bb"`` (Barzilai-Borwein, default) or ``"grow_shrink"`` (merit-based
         grow/shrink).

   -  -  ``--step-size-grow``
      -  float
      -  1.2 (dimensionless)
      -  Grow factor for the ``grow_shrink`` method.

   -  -  ``--step-size-shrink``
      -  float
      -  0.7 (dimensionless)
      -  Shrink factor for the ``grow_shrink`` method. Together with grow=1.2, the product 0.84 per oscillation cycle
         provides mild damping to stabilise convergence.

   -  -  ``--step-size-min``
      -  float
      -  1.0×10⁻⁴ Bohr²/Hartree
      -  Floor for the adaptive step size (both methods).

   -  -  ``--step-size-max``
      -  float
      -  1.0 Bohr²/Hartree
      -  Ceiling for the adaptive step size (both methods).

   -  -  ``--verify-seam-minimum / --no-verify-seam-minimum``
      -  bool
      -  False
      -  After convergence, verify the MECP is a true minimum on the crossing seam via effective Hessian
         analysis. Requires ~4×3N additional Gaussian sub-jobs.  Results are written to
         ``<label>_seam_check.log``. See the :ref:`seam-minimum-verification` section below.

   -  -  ``--hess-step-size``
      -  float
      -  1.0×10⁻³ Bohr
      -  Finite-difference step size used by ``--verify-seam-minimum`` for the numerical Hessian.

.. _convergence-presets:

Convergence Presets
===================

Two built-in convergence presets are available via ``--convergence``:

**Standard** (``--convergence standard``, default)
   General-purpose preset suitable for most MECP searches.

   .. list-table::
      :header-rows: 1
      :widths: 35 20 15 30

      -  -  Quantity
         -  Symbol
         -  Value
         -  Unit

      -  -  Energy difference
         -  :math:`|E_A - E_B|`
         -  5.0×10⁻⁵
         -  Hartree

      -  -  Maximum seam-tangent gradient
         -  :math:`\max|\mathbf{g}_\perp|`
         -  7.0×10⁻⁴
         -  Hartree/Bohr

      -  -  RMS seam-tangent gradient
         -  :math:`\text{RMS}(\mathbf{g}_\perp)`
         -  5.0×10⁻⁴
         -  Hartree/Bohr

      -  -  Maximum displacement
         -  :math:`\max|\mathbf{d}|`
         -  4.0×10⁻³
         -  Bohr

      -  -  RMS displacement
         -  :math:`\text{RMS}(\mathbf{d})`
         -  2.5×10⁻³
         -  Bohr

      -  -  Trust radius
         -  —
         -  0.3
         -  Bohr/atom

**Tight** (``--convergence tight``)
   Publication-quality refinement. Use to confirm and report final MECP geometries.

   .. list-table::
      :header-rows: 1
      :widths: 35 20 15 30

      -  -  Quantity
         -  Symbol
         -  Value
         -  Unit

      -  -  Energy difference
         -  :math:`|E_A - E_B|`
         -  1.0×10⁻⁵
         -  Hartree

      -  -  Maximum seam-tangent gradient
         -  :math:`\max|\mathbf{g}_\perp|`
         -  3.0×10⁻⁴
         -  Hartree/Bohr

      -  -  RMS seam-tangent gradient
         -  :math:`\text{RMS}(\mathbf{g}_\perp)`
         -  1.0×10⁻⁴
         -  Hartree/Bohr

      -  -  Maximum displacement
         -  :math:`\max|\mathbf{d}|`
         -  2.0×10⁻³
         -  Bohr

      -  -  RMS displacement
         -  :math:`\text{RMS}(\mathbf{d})`
         -  1.0×10⁻³
         -  Bohr

      -  -  Trust radius
         -  —
         -  0.1
         -  Bohr/atom

Individual options (``--energy-diff-tol``, ``--force-max-tol``, etc.) always override
the preset values when provided.

Convergence is declared when **all** five criteria are simultaneously satisfied.

Adaptive Step Size
==================

When ``--adaptive-step-size`` is enabled (the default), the step size :math:`\alpha` is updated at the end of each
iteration. Two algorithms are available via ``--step-size-method``.

Barzilai-Borwein (``"bb"``, default)
------------------------------------

The BB2 step size is derived from the secant condition and provides near-quadratic convergence near the MECP:

.. math::

   \alpha_{n+1} = \frac{\|\Delta\mathbf{r}\|^2}{\Delta\mathbf{r} \cdot \Delta\mathbf{g}_\perp}

where :math:`\Delta\mathbf{r} = \mathbf{r}_n - \mathbf{r}_{n-1}` and :math:`\Delta\mathbf{g}_\perp =
\mathbf{g}_{\perp,n} - \mathbf{g}_{\perp,n-1}`.

The denominator is the inner product of the position change and the seam-tangent gradient change; it approximates the
local curvature. If the curvature is non-positive (negative-curvature region or first step), the algorithm falls back to
the initial ``step_size``. The result is clamped to ``[step_size_min, step_size_max]``.

Grow-Shrink (``"grow_shrink"``)
-------------------------------

A dimensionless merit function tracks progress:

.. math::

   M_n = \frac{|\Delta E_n|}{\epsilon_{\Delta E}} + \frac{\text{RMS}(\mathbf{g}_{\perp,n})}{\epsilon_{\text{rms}}}

-  If :math:`M_n < M_{n-1}` (progress): :math:`\alpha_{n+1} = \min(\alpha_n \times \texttt{step\_size\_grow},\,
   \texttt{step\_size\_max})`
-  If :math:`M_n \geq M_{n-1}` (stall/overshoot): :math:`\alpha_{n+1} = \max(\alpha_n \times
   \texttt{step\_size\_shrink},\, \texttt{step\_size\_min})`

The current step size is recorded on every line of ``<label>_report.log``.

Output Files
============

Three output files are produced alongside the Gaussian sub-job input/output files; the
third (``<label>_seam_check.log``) is written only when ``--verify-seam-minimum`` is requested:

``<label>_report.log``
   Step-by-step optimization log. The file header records the run settings; each subsequent line reports one step, using
   a **1-indexed** six-digit zero-padded step counter (``000001`` = first step):

   .. code::

      step=NNNNNN E_A=<Hartree> E_B=<Hartree> dE=<±Hartree>
      pgrad_max=<H/Bohr> pgrad_rms=<H/Bohr>
      disp_max=<Bohr> disp_rms=<Bohr>
      seam_max=<Bohr> seam_rms=<Bohr>
      step_size=<Bohr²/Hartree>

   where:

   -  ``dE`` = :math:`E_A - E_B` — energy difference (drives toward the seam).
   -  ``pgrad_max`` / ``pgrad_rms`` — max and RMS of the seam-tangent gradient
      :math:`\mathbf{g}_\perp` (projection of :math:`\nabla E_A` onto the seam; drives
      geometry to the minimum on the seam).
   -  ``disp_max`` / ``disp_rms`` — max and RMS of the total Cartesian displacement
      :math:`\mathbf{d} = \mathbf{d}_\text{seam} + \mathbf{d}_\text{seam-min}` after trust-radius scaling.
   -  ``seam_max`` / ``seam_rms`` — max and RMS of the seam-correction component
      :math:`\mathbf{d}_\text{seam} = -(\Delta E / \|\mathbf{g}_\Delta\|^2)\,\mathbf{g}_\Delta` that
      moves the geometry toward the crossing surface.

   The final line reads ``Converged at step NNNNNN.`` on successful convergence. The presence of this ``Converged``
   marker is used by ``skip_completed`` to avoid re-running a finished job.

``<label>_traj.xyz``
   Multi-frame XYZ trajectory of the MECP geometry at every optimization step (coordinates in Ångström).

``<label>_seam_check.log``
   Written only when ``--verify-seam-minimum`` is requested. Reports the eigenvalues of the effective projected
   Hessian :math:`H_\text{eff}` (translations, rotations, and gradient-difference direction removed) and whether
   the MECP is a true minimum on the seam. See the :ref:`seam-minimum-verification` section below.

.. _seam-minimum-verification:

Seam Minimum Verification
==========================

A converged MECP may be a crossing point anywhere on the crossing seam, not necessarily the **minimum energy**
point on it.  To confirm that the MECP is a true minimum on the seam (analogous to verifying a transition state has
exactly one imaginary frequency), ChemSmart implements an effective Hessian analysis.

Theory
------

At the MECP the crossing seam is a :math:`(3N-1)`-dimensional hypersurface.  Motions along the gradient-difference
direction :math:`\mathbf{g}_\Delta = \nabla E_A - \nabla E_B` take the molecule off the seam.  The remaining
:math:`3N-1` directions span the seam tangent space; after further removal of translations (3) and rotations (up to 3)
there are :math:`3N - 7` (or :math:`3N - 6` for linear molecules) internal seam degrees of freedom.

ChemSmart constructs a projector that removes these constrained directions:

.. math::

   P = I - \sum_i |\mathbf{v}_i\rangle\langle\mathbf{v}_i|

where :math:`\{\mathbf{v}_i\}` is an orthonormal set spanning translations, rotations, and
:math:`\hat{\mathbf{g}}_\Delta`.  The **effective Hessian** is

.. math::

   H_\text{eff} = P\,\bar{H}\,P, \qquad \bar{H} = \tfrac{1}{2}(H_A + H_B)

where :math:`H_A` and :math:`H_B` are the numerical Hessians of the two states, averaged to give a balanced
description.  If all non-zero eigenvalues of :math:`H_\text{eff}` are positive, the point is confirmed as a seam
minimum; any negative eigenvalue indicates a lower-energy MECP elsewhere on the seam.

.. note::

   A **standard Gaussian frequency analysis** at the MECP geometry is **not sufficient** for this check: it does not
   project out the gradient-difference direction, so it will always show one near-zero or spurious mode whose sign
   is ambiguous.  The effective Hessian analysis described here is the correct diagnostic (cf. ORCA manual, §9.40).

Usage
-----

Add ``--verify-seam-minimum`` to the MECP command after the optimization converges:

.. code:: bash

   chemsmart sub gaussian -p project -f structure.log -c 0 -m 1 mecp \
       --convergence tight --verify-seam-minimum

The verification requires **4 × 3N** additional Gaussian sub-jobs (2 displaced geometries × 2 spin states
× 3N Cartesian coordinates), each labelled ``<label>_step900000_A`` etc.  For a 10-atom molecule this is 120
additional Gaussian calculations.  The finite-difference step size (default 1×10⁻³ Bohr) can be adjusted with
``--hess-step-size``.

Results are written to ``<label>_seam_check.log``:

.. code::

   CHEMSMART MECP seam-minimum verification
   label=... hess_step=1.00e-03 Bohr step_prefix=900000
   energy_diff=+1.234567e-06 Hartree n_projected=7
   n_negative_eigenvalues=0  MECP MINIMUM

   Eigenvalues of H_eff (Hartree/Bohr^2):
     mode    1: +1.234567e-03
     mode    2: +2.345678e-03
     ...

Basic Usage
===========

Singlet/triplet MECP from a Gaussian output structure (charge and multiplicity inferred from file):

.. code:: bash

   chemsmart sub gaussian -p project -f structure.log mecp

Set spin states explicitly (singlet ↔ triplet):

.. code:: bash

   chemsmart sub gaussian -p project -f structure.log -c 0 -m 1 mecp --multiplicity-a 1 --multiplicity-b 3

Doublet/quartet MECP for an open-shell cation:

.. code:: bash

   chemsmart sub gaussian -p project -f radical.log -c 1 -m 2 mecp --multiplicity-a 2 --multiplicity-b 4

Use tight convergence (publication quality):

.. code:: bash

   chemsmart sub gaussian -p project -f structure.log -c 0 -m 1 mecp --convergence tight

Use tight convergence and then verify the geometry is a true seam minimum:

.. code:: bash

   chemsmart sub gaussian -p project -f structure.log -c 0 -m 1 mecp \
       --convergence tight --verify-seam-minimum

Override individual thresholds (tight preset + custom energy threshold):

.. code:: bash

   chemsmart sub gaussian -p project -f structure.log -c 0 -m 1 mecp \
       --convergence tight --energy-diff-tol 5.0e-6

Use a fixed step size (disable adaptive scaling):

.. code:: bash

   chemsmart sub gaussian -p project -f structure.log -c 0 -m 1 mecp \
       --no-adaptive-step-size --step-size 0.05

Use the grow/shrink adaptive method instead of the default Barzilai-Borwein:

.. code:: bash

   chemsmart sub gaussian -p project -f structure.log -c 0 -m 1 mecp \
       --step-size-method grow_shrink --step-size-grow 1.2 --step-size-shrink 0.6

.. note::

   Each MECP step generates two Gaussian sub-jobs named ``<label>_step<NNNNNN>_A`` and ``<label>_step<NNNNNN>_B``
   (single-point energy + forces), where ``<NNNNNN>`` is the **1-indexed** six-digit zero-padded step number (e.g.
   ``000001`` for the first step, ``000010`` for the tenth). Steps are numbered starting from 1 up to ``max_steps`` (up
   to 999 999). These sub-jobs are always re-run (``skip_completed=False``), while the outer MECP job itself honours
   ``skip_completed`` via the ``Converged`` marker in the report file.

***********
 Link Jobs
***********

Run multi-step Gaussian calculations with linked job steps.

.. code:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] link [SUBCMD_OPTIONS]

Link Options
============

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-j, --jobtype``
      -  string
      -  Job type: opt, ts, modred, scan, sp, irc

   -  -  ``-st, --stable``
      -  string
      -  Stability test options (default: opt)

   -  -  ``-g, --guess``
      -  string
      -  Guess options (default: mix). Separate multiple options with a comma, e.g. ``mix,always``.

   -  -  ``--route``
      -  string
      -  Route for link section

Basic Usage
===========

Link job with optimization:

.. code:: bash

   chemsmart sub gaussian -p project -f molecule.xyz link -j opt

Link job with single point:

.. code:: bash

   chemsmart sub gaussian -p project -f molecule.xyz -c 0 -m 1 -r scf=qc link -j sp -so iterative

Examples
========

Optimization of singlet open-shell structure:

.. code:: bash

   chemsmart sub -s SLURM gaussian -p project -f dimer.gjf -c 0 -m 1 link -j opt

This creates a multi-step workflow:

.. code:: text

   # um062x def2svp stable=opt guess=mix
   ...
   # opt freq um062x def2svp geom=check guess=read
   ...
   #N Geom=AllCheck Guess=TCheck SCRF=Check GenChk UM062X/def2SVP Freq

To use multiple guess options, separate them with a comma:

.. code:: bash

   chemsmart sub -s SLURM gaussian -p project -f dimer.gjf -c 0 -m 1 link -j opt -g mix,always

This sets ``guess=(mix,always)`` in the route string:

.. code:: text

   # um062x def2svp stable=opt guess=(mix,always)
   ...
   # opt freq um062x def2svp geom=check guess=read
   ...
   #N Geom=AllCheck Guess=TCheck SCRF=Check GenChk UM062X/def2SVP Freq

******************
 Custom User Jobs
******************

Run custom calculations not built into Chemsmart.

.. code:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] userjob [SUBCMD_OPTIONS]

Custom Job Options
==================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-r, --route``
      -  string
      -  User-defined route (required)

   -  -  ``-a, --append-info``
      -  string
      -  Information to append after coordinates

Basic Usage
===========

.. code:: bash

   chemsmart sub gaussian -p project -f molecule.com -l custom_job userjob -r 'opt freq b3lyp/6-31g*' -a 'B 1 2 F'

*****************************
 Direct Input File Execution
*****************************

Run a pre-prepared Gaussian input file without modifications.

.. code:: bash

   chemsmart sub [OPTIONS] gaussian [GAUSSIAN_OPTIONS] com

Basic Usage
===========

Run a ``.com`` file:

.. code:: bash

   chemsmart sub gaussian -p project -f input_file.com com

Run a ``.gjf`` file:

.. code:: bash

   chemsmart sub gaussian -p project -f input_file.gjf com

Modify charge and multiplicity:

.. code:: bash

   chemsmart sub gaussian -p project -f input_file.com -c 1 -m 2 com
