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

CHEMSMART provides a **self-contained MECP optimizer** that drives the search entirely in Python, calling Gaussian only
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

   -  -  ``--multiplicity1``
      -  int
      -  1 (singlet)
      -  Spin multiplicity of state A. Falls back to ``-m`` if not set.

   -  -  ``--multiplicity2``
      -  int
      -  multiplicity1 + 2
      -  Spin multiplicity of state B.

   -  -  ``--charge1``
      -  int
      -  0
      -  Charge of state A. Falls back to ``-c`` if not set.

   -  -  ``--charge2``
      -  int
      -  charge1
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
      -  Convergence preset: ``"standard"`` (default, general-purpose) or ``"tight"`` (publication-quality). See the
         :ref:`convergence-presets` section below. Individual tolerance options override the preset.

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
      -  ``"harvey"``
      -  MECP optimizer: ``"harvey"`` (inverse-BFGS, default), ``"bb"`` (Barzilai-Borwein), or ``"grow_shrink"``
         (merit-based grow/shrink).

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
      -  Floor for the adaptive step size (``bb`` and ``grow_shrink``).

   -  -  ``--step-size-max``
      -  float
      -  1.0 Bohr²/Hartree
      -  Ceiling for the adaptive step size (``bb`` and ``grow_shrink``).

   -  -  ``--restart / --no-restart``
      -  bool
      -  True
      -  Resume an interrupted MECP optimization from ``<label>_state.npz``.

   -  -  ``--mecp-numfreq / --no-mecp-numfreq``
      -  bool

      -  False
      -  Perform the seam-minimum verification and additionally write mass-weighted projected MECP frequencies and
         normal modes to ``<label>_mecp_freq.log``.

   -  -  ``--hess-step-size``
      -  float
      -  1.0×10⁻³ Bohr
      -  Finite-difference step size used by ``--mecp-numfreq`` for the numerical Hessian.

   -  -  ``--follow-seam-imaginary-mode / --no-follow-seam-imaginary-mode``
      -  bool
      -  False
      -  If the projected Hessian contains a significant imaginary mode, follow it in both directions using constrained
         progress planes, then release and retain the lower verified seam minimum. Implies ``--mecp-numfreq``.

   -  -  ``--seam-mode-displacement``
      -  float
      -  0.05 Å
      -  Progress increment along the tracked projected mode at each constrained macro step.

   -  -  ``--seam-mode-max-steps``
      -  int
      -  30
      -  Maximum constrained mode-following macro steps in each positive/negative direction.

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

Individual options (``--energy-diff-tol``, ``--force-max-tol``, etc.) always override the preset values when provided.

Convergence is declared when **all** five criteria are simultaneously satisfied.

Adaptive Step Size
==================

When ``--adaptive-step-size`` is enabled (the default), the step size :math:`\alpha` is updated at the end of each
iteration. The available algorithms are selected via ``--step-size-method``. For the default ``harvey`` optimizer, the
inverse Hessian controls the step; ``--adaptive-step-size`` and the scalar ``step-size-*`` controls apply only to ``bb``
and ``grow_shrink``.

Barzilai-Borwein (``"bb"``)
---------------------------

The BB step size is derived from the secant condition and accelerates convergence near the MECP:

.. math::

   \alpha_{n+1} = \frac{\|\Delta\mathbf{r}\|^2}{\Delta\mathbf{r} \cdot \Delta\mathbf{g}_\perp}

where :math:`\Delta\mathbf{r} = \mathbf{r}_n - \mathbf{r}_{n-1}` and :math:`\Delta\mathbf{g}_\perp =
\mathbf{g}_{\perp,n} - \mathbf{g}_{\perp,n-1}`.

Before applying the secant formula, the position change is projected onto the current seam tangent so that it is paired
consistently with the seam-tangent gradient change. Curvature is accepted only when its normalized magnitude is reliably
positive. An unreliable pair damps the current step by ``step_size_shrink`` instead of resetting it. Even a valid BB
estimate is limited to between 0.5 and 2 times the current step before the configured ``[step_size_min, step_size_max]``
bounds are applied. These safeguards prevent small secant denominators from causing abrupt jumps to ``step_size_max``.

Grow-Shrink (``"grow_shrink"``)
-------------------------------

A dimensionless merit function tracks progress:

.. math::

   M_n = \frac{|\Delta E_n|}{\epsilon_{\Delta E}} + \frac{\text{RMS}(\mathbf{g}_{\perp,n})}{\epsilon_{\text{rms}}}

The update uses relative merit progress rather than reacting to every numerical change:

-  Improvement greater than 10%: grow by ``step_size_grow``.
-  Change between a 2% regression and a 10% improvement: keep the step.
-  Regression greater than 2%: shrink by ``step_size_shrink``.

The dead band prevents small SCF and gradient fluctuations from making the step size oscillate.

The current step size is recorded on every line of ``<label>_report.log``.

Harvey inverse-BFGS (``"harvey"``, default)
-------------------------------------------

This method follows the inverse-BFGS update used by easyMECP: it builds a full inverse Hessian from successive
effective-gradient and Cartesian-displacement pairs, including negative-curvature updates, and limits the largest
Cartesian component using Harvey's ``STPMX`` rule. Numerically singular updates are skipped. Ill-conditioned inverse
Hessians and non-descent directions are reset to the configured diagonal initial inverse Hessian. The corresponding
``bfgs_status=RESET_*`` reason is recorded in the report.

Restarting interrupted calculations
===================================

MECP optimization state is written atomically after every completed step to ``<label>_state.npz``. It contains the next
geometry, inverse Hessian, previous effective gradient, and adaptive-step history. Re-running the same job resumes from
that state by default. The saved atom sequence and optimizer method must match the new invocation; otherwise CHEMSMART
stops with an explicit error. Use ``--no-restart`` to deliberately start from the supplied input geometry. The state
file is removed after successful convergence.

MECP force calculations always include ``nosymm`` so that Gaussian Cartesian forces remain aligned with the optimizer
coordinate frame. An explicit conflicting ``symmetry`` route option is rejected.

Like easyMECP, CHEMSMART maintains one rolling checkpoint per state (``<label>_A.chk`` and ``<label>_B.chk``) instead of
one checkpoint per iteration. It does not add ``guess=read`` automatically: checkpoint orbitals are read only when that
option is explicitly present in the Gaussian route. If the first step requests ``guess=read`` but its checkpoint does
not yet exist, ``read`` is removed for that first step and retained thereafter. Iteration ``.com`` and ``.log`` files
and the two rolling checkpoints are kept in ``<label>_optimization``, along with the report and trajectory.

When Gaussian scratch storage is enabled, initial MECP scratch jobs are grouped below ``<label>_optimization``
and frequency-analysis jobs below ``<label>_numfreq``, instead of being created directly in the scratch root.

Seam-minimum verification uses a separate pair of temporary rolling checkpoints. They are deleted after a successful
verification, so the final ``<label>_A.chk`` and ``<label>_B.chk`` continue to represent the converged MECP geometry.
Temporary checkpoints are retained if verification fails.

Output Files
============

The main job directory contains the final report and, when requested, the thermochemistry-ready frequency file.
Intermediate files are grouped by function: ``<label>_optimization`` for the initial MECP search,
``<label>_numfreq`` for numerical Hessian analysis, and ``<label>_seam_follow`` for optional imaginary-mode following.
Only requested functions create their directories:

``<label>_report.log``
   Located in ``<label>_optimization``. Step-by-step optimization log. The file header records the run settings; each subsequent line reports one step, using
   a **1-indexed** step counter (``1`` = first step):

   .. code::

      step=N E_A=<Hartree> E_B=<Hartree> dE=<±Hartree>
      pgrad_max=<H/Bohr> pgrad_rms=<H/Bohr>
      disp_max=<Bohr> disp_rms=<Bohr>
      seam_max=<Bohr> seam_rms=<Bohr>
      step_size=<Bohr²/Hartree>

   where:

   -  ``dE`` = :math:`E_A - E_B` — energy difference (drives toward the seam).
   -  ``pgrad_max`` / ``pgrad_rms`` — max and RMS of the seam-tangent gradient :math:`\mathbf{g}_\perp` (projection of
      :math:`\nabla E_A` onto the seam; drives geometry to the minimum on the seam).
   -  ``disp_max`` / ``disp_rms`` — max and RMS of the total Cartesian displacement :math:`\mathbf{d} =
      \mathbf{d}_\text{seam} + \mathbf{d}_\text{seam-min}` after trust-radius scaling.
   -  ``seam_max`` / ``seam_rms`` — max and RMS of the seam-correction component :math:`\mathbf{d}_\text{seam} =
      -(\Delta E / \|\mathbf{g}_\Delta\|^2)\,\mathbf{g}_\Delta` that moves the geometry toward the crossing surface.

   After optimization, the report distinguishes the initial crossing from the final verified MECP. If seam-mode
   following was needed, it names the selected branch, reports its number of macro steps, and records the final energy,
   energy gap, and number of
   significant imaginary modes. If seam verification fails, it records the error instead. Only a successful job ends
   with ``Converged at step N.``; this marker is used by ``skip_completed`` to avoid re-running a finished job.

``<label>_final_report.log``
   Concise report for the final selected structure. It contains both state energies, MECP energy, final values and
   thresholds for the energy gap, projected gradients, and displacements, plus the final Cartesian geometry. With
   ``--mecp-numfreq`` it also lists the final projected frequencies and seam-minimum status. It distinguishes the
   initial optimization step count from the number of constrained seam-following macro steps and, when applicable,
   the final branch optimization step count. Intermediate structures and per-step diagnostics remain in the other
   reports.

``<label>_traj.xyz``
   Located in ``<label>_optimization``. Multi-frame XYZ trajectory of the initial MECP optimization
   (coordinates in Ångström). Its A/B Gaussian sub-jobs are in the same directory.

``<label>_seam_check.log``
   Located in ``<label>_numfreq`` when ``--mecp-numfreq`` is requested, alongside the numerical-Hessian Gaussian
   sub-jobs. Reports the eigenvalues of the effective projected Hessian
   :math:`H_\text{eff}` (translations, rotations, and gradient-difference direction removed) and whether the MECP is a
   true minimum on the seam. See the :ref:`seam-minimum-verification` section below.

``<label>_mecp_freq.log``
   Written only by ``--mecp-numfreq``. Contains the converged geometry, both state energies, MECP energy, atomic
   masses, mass-weighted projected frequencies, and Cartesian normal modes. It can be passed directly to CHEMSMART
   thermochemistry.

.. _seam-minimum-verification:

Seam Minimum Verification
=========================

A converged MECP may be a crossing point anywhere on the crossing seam, not necessarily the **minimum energy** point on
it. To confirm that the MECP is a true minimum on the seam (analogous to verifying a transition state has exactly one
imaginary frequency), CHEMSMART implements an effective Hessian analysis.

Theory
------

At the MECP the crossing seam is a :math:`(3N-1)`-dimensional hypersurface. Motions along the gradient-difference
direction :math:`\mathbf{g}_\Delta = \nabla E_A - \nabla E_B` take the molecule off the seam. The remaining :math:`3N-1`
directions span the seam tangent space; after further removal of translations (3) and rotations (up to 3) there are
:math:`3N - 7` (or :math:`3N - 6` for linear molecules) internal seam degrees of freedom.

CHEMSMART constructs a projector that removes these constrained directions:

.. math::

   P = I - \sum_i |\mathbf{v}_i\rangle\langle\mathbf{v}_i|

where :math:`\{\mathbf{v}_i\}` is an orthonormal set spanning translations, rotations, and
:math:`\hat{\mathbf{g}}_\Delta`. The **effective Lagrangian Hessian** is

.. math::

   H_\text{eff} = P\,[(1-\lambda)H_A + \lambda H_B]P

where :math:`H_A` and :math:`H_B` are the numerical Hessians of the two states and :math:`\lambda` is the constrained
MECP Lagrange multiplier. For frequencies, the Hessian and projection vectors are transformed to mass-weighted
coordinates before diagonalisation. If all projected eigenvalues are positive, the point is confirmed as a seam
minimum; any negative eigenvalue indicates a lower-energy MECP elsewhere on the seam.

.. note::

   A **standard Gaussian frequency analysis** at the MECP geometry is **not sufficient** for this check: it does not
   project out the gradient-difference direction, so it will always show one near-zero or spurious mode whose sign is
   ambiguous. The effective Hessian analysis described here is the correct diagnostic (cf. ORCA manual, §9.40).

Usage
-----

Use ``--mecp-numfreq`` to check whether the MECP is a seam minimum and write
projected frequencies and modes:

.. code:: bash

   chemsmart sub gaussian -p project -f structure.log -c 0 -m 1 mecp \
       --convergence tight --mecp-numfreq

To escape a stationary point that is a saddle on the crossing seam, add
``--follow-seam-imaginary-mode``. In each direction CHEMSMART fixes a progress plane normal to the tracked negative
mode, restores :math:`E_A=E_B`, and optimizes every remaining seam coordinate. It then recomputes the effective Hessian,
identifies the continuation of the mode by maximum absolute overlap, and advances the progress plane again. The progress
constraint is released only after no significant negative projected mode remains; a final unconstrained MECP optimization
and frequency check must still confirm the result:

.. code:: bash

   chemsmart sub gaussian -p project -f structure.xyz -c 0 -m 1 mecp \
       --convergence tight --follow-seam-imaginary-mode \
       --seam-mode-displacement 0.05 --seam-mode-max-steps 30

All displaced structures, branch reports, trajectories, frequency logs, and Gaussian sub-jobs are collected under
``<label>_seam_follow/``. Branch labels end in ``_seam_follow_plus`` and
``_seam_follow_minus``. The final projected-frequency log remains in the main calculation directory; the final seam
check is in ``<label>_numfreq`` and the selection record ``<label>_seam_follow.log`` is in
``<label>_seam_follow``.
This is constrained iterative
seam-mode following; it does not invoke Gaussian IRC/QRC.

Branch optimizations use at least the tight convergence thresholds. If neither direction loses its negative mode within
``--seam-mode-max-steps``, or if the released optimization returns to a saddle, the calculation stops and retains both
branch results. The user can then adjust the progress increment or maximum number of macro steps.

Each Hessian evaluation requires **4 × 3N** Gaussian sub-jobs (2 displaced geometries × 2 spin states × 3N Cartesian
coordinates). During seam-mode following, labels identify the macro step, coordinate, displacement sign, and state,
such as ``<label>_seam_follow_plus_macro07_check_coord18_plus_B.log``. Constrained optimization sub-jobs use names such as
``<label>_seam_follow_plus_macro07_inner012_A.log``. Standalone seam checks retain the original ``check_stepN`` naming.
For a 10-atom molecule this is 120
Gaussian calculations per macro step. Iterative following performs this analysis after every constrained macro step in
both directions, plus the initial and final checks, so users should choose ``--seam-mode-max-steps`` conservatively.
The finite-difference step size (default 1×10⁻³ Bohr) can be adjusted with ``--hess-step-size``.

Results are written to ``<label>_numfreq/<label>_seam_check.log``:

.. code::

   CHEMSMART MECP seam-minimum verification
   label=... hess_step=1.00e-03 Bohr check_step_start=1
   energy_diff=+1.234567e-06 Hartree n_projected=7
   n_negative_eigenvalues=0  MECP MINIMUM

   Eigenvalues of H_eff (Hartree/Bohr^2):
     mode    1: +1.234567e-03
     mode    2: +2.345678e-03
     ...

The corresponding ``3N-7`` projected frequencies and normal modes are written to ``<label>_mecp_freq.log``. For a
linear molecule the number of modes is ``3N-6``. Thermochemistry can then be calculated with, for example:

.. code:: bash

   chemsmart run thermochemistry -f <label>_mecp_freq.log -T 298.15

The MECP frequency file uses a conservative electronic degeneracy of 1 by default. If a different statistical weight
is required by the chosen non-adiabatic rate theory, specify it with ``--electronic-degeneracy``.

Basic Usage
===========

Singlet/triplet MECP from a Gaussian output structure (charge and multiplicity inferred from file):

.. code:: bash

   chemsmart sub gaussian -p project -f structure.log mecp

Set spin states explicitly (singlet ↔ triplet):

.. code:: bash

   chemsmart sub gaussian -p project -f structure.log -c 0 -m 1 mecp --multiplicity1 1 --multiplicity2 3

Doublet/quartet MECP for an open-shell cation:

.. code:: bash

   chemsmart sub gaussian -p project -f radical.log -c 1 -m 2 mecp --multiplicity1 2 --multiplicity2 4

Use tight convergence (publication quality):

.. code:: bash

   chemsmart sub gaussian -p project -f structure.log -c 0 -m 1 mecp --convergence tight

Use tight convergence and verify the geometry is a true seam minimum:

.. code:: bash

   chemsmart sub gaussian -p project -f structure.log -c 0 -m 1 mecp \
       --convergence tight --mecp-numfreq

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

   Each MECP step generates two Gaussian sub-jobs in ``<label>_optimization``, named ``<label>_step<N>_A`` and
   ``<label>_step<N>_B`` (single-point energy + forces), where ``<N>`` is the **1-indexed** step number (for example,
   ``step1`` and ``step10``). These sub-jobs are always re-run (``skip_completed=False``), while the outer MECP job
   itself honours ``skip_completed`` via the ``Converged`` marker in the report file.

***********
 Link Jobs
***********

Run multi-step Gaussian calculations with linked job steps. Current Link job only supports openshell singlet broken
symmetry guess as first job, followed by link job to run different job types (opt, modred, ts, scan, sp, etc).

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

To use multiple guess options, separate them with a comma:

.. code:: bash

   chemsmart sub -s SLURM gaussian -p project -f dimer.gjf -c 0 -m 1 link -j opt -g mix,always

This sets ``guess=(mix,always)`` in the route string:

.. code:: text

   # um062x def2svp stable=opt guess=(mix,always)
   ...
   # opt freq um062x def2svp geom=check guess=read
   ...

******************
 Custom User Jobs
******************

Run custom calculations not built into CHEMSMART.

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
