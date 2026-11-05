################
 ORCA MECP Jobs
################

Run minimum energy crossing-point calculations with ORCA using the following command structure:

.. code:: bash

   chemsmart sub orca -p project -f structure.xyz -c 0 mecp -m1 1 -m2 3 [MODE] [OPTIONS]

ORCA MECP calculations use the same charge and level of theory for both spin states. Specify their multiplicities with
``-m1/--multiplicity1`` and ``-m2/--multiplicity2``. The default ``mecp opt`` performs native ``SurfCrossOpt``
optimization; ``mecp numfreq`` retains ORCA's native ``SurfCrossNumFreq`` verification after optimization.

Use ``anafreq`` for two fixed-geometry analytic state Hessians, or ``followseam`` to follow a negative projected seam
mode in both directions:

.. code:: bash

   chemsmart sub orca -p project -f structure.xyz -c 0 mecp -m1 1 -m2 3 anafreq
   chemsmart sub orca -p project -f crossing.xyz -c 0 mecp -m1 1 -m2 3 anafreq --only
   chemsmart sub orca -p project -f structure.xyz -c 0 mecp -m1 1 -m2 3 followseam --hessian-method analytic

Without ``--only``, initial optimization and final unconstrained branch relaxation use native ``SurfCrossOpt`` with
``TightOpt``. The shared CHEMSMART driver controls the intervening constrained macro steps and selects the lower
verified seam minimum. Every initial, macro-step and final Hessian check uses the selected method. ``followseam``
defaults to numerical Hessians; ``--hess-step-size`` (Bohr, default 0.001) controls the Cartesian central differences.
``--seam-mode-displacement`` (Angstrom, default 0.05) and ``--seam-mode-max-steps`` (default 30) control following;
``--max-steps`` (default 200) limits each inner constrained optimization. ``--maxiter`` controls native ORCA
optimization iterations separately.

CHEMSMART seam analysis defaults to the ``standard`` preset, as for Gaussian MECP. Its maximum and RMS projected
gradient limits are ``7e-4`` and ``5e-4`` Hartree/Bohr. These settings control shared analysis and constrained seam
following, not native ORCA optimization. Native ``TightOpt`` is a separate ORCA setting.

When ``followseam`` creates each branch, it takes the smaller of each current threshold and the corresponding backend
``tight`` preset value, preserving stricter user settings. This applies before constrained macro steps, branch
relaxation and verification. ORCA branch defaults are ``3e-4`` and ``1e-4`` Hartree/Bohr for maximum and RMS projected
gradients, ``2e-3`` and ``1e-3`` Bohr for maximum and RMS displacement, and ``0.1`` Bohr for the trust radius.

Both ORCA presets retain an absolute state energy-gap limit of ``1e-4`` Hartree, matching native ``SurfCrossOpt``
convergence. This applies both after native optimization and to ``anafreq --only`` and ``numfreq --only``; no extra
tolerance option is needed for a structure optimized with the native default gap. Projected gradient checks still apply.
To require a stricter gap, use, for example, ``anafreq --only --energy-diff-tol 1e-5`` and supply a structure optimized
to that accuracy. This option changes the shared analysis threshold, not native optimization convergence, and does not
infer custom convergence settings from XYZ files. A failed stationarity check reports the measured energy gap and
projected gradients with their limits before computing Hessians. Gaussian convergence presets are unchanged.

``numfreq --only`` also uses the shared driver, computing numerical Hessians at the supplied geometry. Both
frequency-only modes require a stationary crossing and do not optimize it. They record whether the supplied crossing is
a minimum; optimization workflows fail if no verified minimum is obtained.

Analytic checks run two reference ``EnGrad`` calculations followed by two ``AnFreq`` calculations reading the
corresponding reference orbitals. Raw Cartesian Hessians are read from ORCA ``.hess`` files in Hartree/Bohr²; the
effective Hessian, seam projection, mass weighting and following algorithms are shared with Gaussian. Missing or invalid
files and unsupported analytic calculations fail explicitly without a numerical fallback.

State calculations are stored in ``<label>_analytic_frequency/`` or ``<label>_numerical_frequency/``, constrained
optimization jobs in ``<label>_optimization/``, and branch records in ``<label>_seam_follow/``. The selected final
geometry, projected frequencies and modes are written to ``<label>_mecp_freq.log``; the workflow status and selected
branch are in ``<label>_final_report.log``. ORCA's original output and native quality report describe the initial native
optimization. Read the final workflow report when following changes the geometry. Raw ``.hess``, ``.engrad`` and
``.gbw`` files are retained for inspection.

The shared ORCA workflow supports ordinary and broken-symmetry PES2 spin states with structured settings and
unconstrained molecules. Supply ``--broken-sym NA,NB`` (or ``-bs NA,NB``) with ``--multiplicity2 abs(NA-NB)+1``. For
example, ``-m1 3 -m2 1 -bs 1,1`` describes a triplet and a BS singlet. The first standalone B-state calculation starts
with multiplicity ``NA+NB+1`` and ``%scf BrokenSym NA,NB``; ORCA constructs the BS wavefunction before evaluating its
gradient. Subsequent B-state gradients and analytic Hessians use unrestricted SCF with the target low-spin multiplicity
and read that state's own GBW file. A-state orbitals remain separate. The native optimization continues to use the
existing ``%mecp brokenSym`` setting. This follows the `ORCA BS procedure
<https://www.faccts.de/docs/orca/6.1/manual/contents/spectroscopyproperties/magnx.html>`_.

Explicit PES2 orbital guesses, CASSCF and frozen-atom/modred constraints remain available through native ``opt/numfreq``
but are rejected by the shared analysis/following workflow until their adapters are available. Analytic Hessians require
support from the chosen ORCA method and solvent model; consult the `ORCA frequency documentation
<https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/frequencies.html>`_.

*****************************
 State conventions and units
*****************************

State A means state 1 (``--multiplicity1/-m1``); state B means state 2 (``--multiplicity2/-m2``). Native ORCA MECP
requires one shared charge for both states.

Use ``--hessian-method numerical|analytic`` to select the Hessian source. ``--hess-step-size`` is in Bohr;
``--seam-mode-displacement/--displacement`` is in Angstrom, and ``--seam-mode-max-steps/--max-follow-steps`` limits
seam-following macro steps. ``numfreq --only`` and ``anafreq --only`` analyze a supplied stationary crossing without
optimization.

Energies are in Hartree, Cartesian geometry in Angstrom, projected gradients in Hartree/Bohr, optimization displacements
in Bohr, Cartesian Hessians in Hartree/Bohr², and projected frequencies in cm⁻¹.

*************************
 Output file conventions
*************************

``<label>_final_report.log``
   Final optimization, analysis or seam-following results, including frequency-only jobs.

``<label>_mecp_freq.log``
   Selected geometry, masses, projected frequencies and computed Cartesian normal modes, in the format accepted by
   CHEMSMART thermochemistry and molecule readers.

Native ORCA ``opt/numfreq`` retains ``<label>_report.log`` with its quality assessment, canonical energy fields and
units. ``numfreq`` without ``--only`` reuses native ``SurfCrossNumFreq`` results without additional state or Hessian
calculations. ``numfreq --only`` instead uses CHEMSMART numerical Lagrangian-Hessian analysis.

``<label>_PES1_mecp_freq.log`` and ``<label>_PES2_mecp_freq.log`` without additional state or Hessian calculations.
These files retain the separate native effective-Hessian definitions; the two frequency sets are never averaged. ORCA
``numfreq --only`` retains CHEMSMART's numerical Lagrangian-Hessian analysis of the supplied structure. Gaussian obtains
numerical state Hessians through the CHEMSMART driver.

ORCA ``anafreq`` runs fixed-geometry ``AnFreq`` for each state and reads both ``.hess`` files. CHEMSMART constructs the
Lagrangian Hessian using the state gradients, then projects translations, rotations and the gradient-difference
direction to obtain seam frequencies.

Energies are in Hartree, Cartesian geometry in Angstrom, projected gradients in Hartree/Bohr, optimization displacements
in Bohr, Cartesian Hessians in Hartree/Bohr², and projected frequencies in cm⁻¹.

**************
 Common files
**************

``<label>_final_report.log``
   One schema for Gaussian and ORCA selected MECP results, including native ORCA ``opt/numfreq`` and CHEMSMART
   analysis/following workflows. Native ORCA also retains its separate ``<label>_report.log`` quality summary.

``<label>_mecp_freq.log``
   Selected geometry, masses, projected frequencies and Cartesian normal modes. Both engines use the existing format
   accepted by CHEMSMART thermochemistry and molecule readers.

``<label>_PES1_mecp_freq.log`` and ``<label>_PES2_mecp_freq.log``
   Native ORCA projected frequencies with the projected zero-mode prefix removed. Each file identifies
   ``frequency_source=ORCA_SurfCrossNumFreq`` and ``frequency_surface=PES1|PES2``. Thermochemistry uses its
   ``electronic_energy`` (the corresponding state's energy); ``mecp_energy`` retains the mean for reference. Only
   complete, normally terminated and converged native results are exported. No normal-mode vectors or Lagrange
   multiplier are fabricated. The electronic statistical weight defaults to one, as for other MECP inputs.

   Analyze the surfaces separately, for example:

   .. code::

      chemsmart run thermochemistry -T 298.15 -f crossing_PES1_mecp_freq.log
      chemsmart run thermochemistry -T 298.15 -f crossing_PES2_mecp_freq.log

   These yield separate Gibbs energies, not a unique Lagrangian-seam Gibbs energy.

``<label>_numerical_frequency/`` or ``<label>_analytic_frequency/``
   Numerical or analytic state calculations and ``<label>_seam_check.log``. Reference jobs use
   ``hessian_reference_A/B``; analytic jobs use ``hessian_analytic_A/B``. Cartesian differences use
   ``hessian_a1_x_plus_A`` and analogous atom/axis/sign/state names. Atom numbering starts at 1.

``<label>_seam_follow/``
   Follow summary ``<label>_seam_follow.log``, displaced XYZ structures, and
   ``<label>_seam_follow_plus/minus_record.log`` branch records. Macro-step jobs include prefixes such as ``macro07``;
   inner optimization files use ``macro07_opt_step012_A/B``. Branch labels prefix their own calculation directories.

ORCA uses ``.inp/.out/.gbw/.engrad/.hess`` files. ``state_output_A/B`` report fields identify original engine outputs
using paths relative to the report.

*********************
 Final report fields
*********************

-  ``label``, ``backend``, ``mode``, and ``hessian_method``.
-  ``optimization=CONVERGED|SKIPPED`` and ``initial_optimization_steps``.
-  ``seam_follow_macro_steps``, ``seam_follow_selected_branch`` and, for a followed result,
   ``final_branch_optimization_steps``.
-  ``state_output_A/B`` when available.
-  ``energy_A``, ``energy_B``, signed ``energy_difference`` (A minus B), and ``mecp_energy`` (their arithmetic mean).
-  Convergence values, thresholds, units and PASS/FAIL assessments. Unavailable native metrics use ``value=NA`` and
   ``status=NA``, never a fabricated zero or PASS.
-  ``seam_minimum=PASS|FAIL|NOT_CHECKED`` and, when checked, ``is_minimum``, ``significant_imaginary_modes`` and
   ``frequency_file``.
-  Final geometry and any computed projected frequencies.
-  ``status=COMPLETED`` after the requested workflow finishes. For frequency-only analysis this indicates completion of
   the analysis; inspect ``seam_minimum`` to distinguish a seam minimum from a saddle.

Native fields including ``state_1_energy/state_2_energy`` and the PASSED/WARNING/FAILED assessment describe the native
quality report.

Code API: ``final_report_file`` identifies the selected final result, and ``native_quality_report_file`` identifies the
native ORCA quality summary. The native quality report describes the initial ORCA calculation, which may differ from the
selected result after seam following. ``verify_seam_minimum()`` returns the canonical ``state_output_files`` mapping.

**Python results API**

``job.results`` returns the selected completed MECP result, or ``None`` if not available. Both backends expose
``energy_a``, ``energy_b``, ``energy_difference`` (Hartree), ``symbols``, ``positions_angstrom``, ``is_minimum``,
``projected_frequencies`` (cm⁻¹), ``source``, ``quality_status``, and ``minimum_check_method`` (``projected_hessian``,
``orca_native_frequency``, or ``None`` when no minimum check was completed). Values absent from the underlying output
are ``None``. Shared results are read from the final report, so values have the report's printed precision. For native
ORCA calculations, ``job.output`` is a dictionary with the same keys as Gaussian. ``native_output`` is the ORCA output
parser, ``native_summary`` is the parsed native MECP result. ``energy_gap_history`` preserves the energy gap records
without assigning them to optimization steps. Explicitly reported ORCA optimization cycle count is available as
``native_summary.optimization_steps``; when absent it is ``None``. ORCA native ``optimization_steps`` is empty because
separate A/B outputs are not available for each cycle. Its native frequency tables are not interchangeable with shared
projected frequencies.

Native ORCA reports record ``analysis_source=orca_native`` and ``quality_status``; their native minimum assessment is
not a CHEMSMART projected-Hessian check.
