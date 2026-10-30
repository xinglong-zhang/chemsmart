##################
 ORCA CLI Options
##################

This page documents the CLI options available for all ORCA jobs. Use ``chemsmart sub orca --help`` for the complete
list.

*************************
 Basic Command Structure
*************************

.. code:: bash

   chemsmart sub [OPTIONS] orca [ORCA_OPTIONS] <SUBCMD> [SUBCMD_OPTIONS]

.. tip::

   ORCA options are largely similar to Gaussian options with some ORCA-specific parameters.

**************
 ORCA Options
**************

Project and File Options
========================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-p, --project``
      -  string
      -  Project settings from ``~/.chemsmart/orca/*.yaml``

   -  -  ``-f, --filename``
      -  string
      -  Input file for job preparation

   -  -  ``-l, --label``
      -  string
      -  Custom output filename (without extension)

   -  -  ``-a, --append-label``
      -  string
      -  String to append to the base filename

   -  -  ``-t, --title``
      -  string
      -  ORCA job title

   -  -  ``-i, --index``
      -  string
      -  Structure index (1-based, default: last structure)

   -  -  ``-P, --pubchem``
      -  string
      -  Query structure from PubChem

   -  -  ``--ri, --record-index``
      -  int
      -  Select a record from a CHEMSMART database by its 1-based index

   -  -  ``--rid, --record-id``
      -  string
      -  Select a record from a CHEMSMART database by its ID

   -  -  ``--sid, --structure-id``
      -  string
      -  Select a structure from a CHEMSMART database by its ID

.. note::

   -  ``-p`` uses the project name without the ``.yaml`` extension.
   -  ``-f`` accepts various formats: ``.xyz``, ``.com``, ``.gjf``, ``.log``, ``.inp``, ``.out``, or a CHEMSMART
      database ``.db`` file.

Molecular Properties Options
============================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-c, --charge``
      -  int
      -  Molecular charge

   -  -  ``-m, --multiplicity``
      -  int
      -  Molecular multiplicity

Method and Basis Set Options
============================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-A, --ab-initio``
      -  string
      -  Ab initio method (e.g., DLPNO-CCSD(T))

   -  -  ``-x, --functional``
      -  string
      -  DFT functional

   -  -  ``-D, --dispersion``
      -  string
      -  Dispersion correction

   -  -  ``-b, --basis``
      -  string
      -  Basis set

   -  -  ``-B, --aux-basis``
      -  string
      -  Auxiliary basis set

   -  -  ``-e, --extrapolation-basis``
      -  string
      -  Extrapolation basis set

SCF and Grid Options
====================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-d, --defgrid``
      -  choice
      -  Grid: defgrid1, defgrid2, defgrid3 (default: defgrid2)

   -  -  ``--scf-tol``
      -  choice
      -  SCF tolerance: NormalSCF, LooseSCF, TightSCF, etc.

   -  -  ``--scf-algorithm``
      -  choice
      -  SCF algorithm: GDIIS, DIIS, SOSCF, AutoTRAH

   -  -  ``--scf-maxiter``
      -  int
      -  Maximum SCF iterations

   -  -  ``--scf-convergence``
      -  float
      -  SCF convergence criterion

Property Calculation Options
============================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``--dipole/--no-dipole``
      -  bool
      -  Dipole moment calculation

   -  -  ``--quadrupole/--no-quadrupole``
      -  bool
      -  Quadrupole moment calculation

   -  -  ``--forces/--no-forces``
      -  bool
      -  Forces calculation (default: disabled)

MDCI Options
============

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``--mdci-cutoff``
      -  choice
      -  MDCI cutoff: loose, normal, tight

   -  -  ``--mdci-density``
      -  choice
      -  MDCI density: none, unrelaxed, relaxed

Additional Options
==================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-r, --additional-route-parameters``
      -  string
      -  Additional route parameters

Solvent Options
===============

Solvent settings can be specified at the ORCA group level, which means they apply to **any** subcommand (``sp``,
``opt``, ``ts``, ``irc``, ``scan``, etc.). This is useful when the project settings define a gas-phase calculation but
you want to add solvation for a particular run without modifying the project file.

They can also be specified at the **subcommand level** to override the group-level settings for a single calculation.

.. list-table::
   :header-rows: 1
   :widths: 35 15 50

   -  -  Option
      -  Type
      -  Description

   -  -  ``--remove-solvent/--no-remove-solvent``
      -  bool
      -  Remove solvent from the job, overriding project settings (default: disabled)

   -  -  ``-sm, --solvent-model``
      -  string
      -  Implicit solvent model: ``cpcm`` (CPCM with CPCM epsilon), ``cpcmc`` (CPCM with COSMO epsilon; replaces the
         legacy COSMO model removed in ORCA 4.0), ``smd`` (Minnesota SMD), or ``cosmors`` (openCOSMO-RS interface)

   -  -  ``-si, --solvent-id``
      -  string
      -  Named solvent identifier (e.g. ``water``, ``toluene``, ``cyclohexane``). Omit when using a fully custom
         dielectric via ``-so``

   -  -  ``-so, --solvent-options``
      -  string
      -  Additional parameters for the model's solvent block (see tables below), newline-separated for multiple options

   -  -  ``-sf, --solventfilename``

      -  path

      -  Path to a solvent file for the ``cosmors`` model. Any file format is accepted — it does **not** have to be a
         ``.cosmorsxyz`` file. If the path points to a Gaussian output file (e.g. ``basename.log``) or an ORCA output
         file (e.g. ``basename.out``), CHEMSMART automatically converts it to ``basename.cosmorsxyz`` (via
         ``Molecule.write_cosmorsxyz()``) before use. The ``.cosmorsxyz`` file is then copied to the running directory
         (scratch or job folder) and its basename (without the ``.cosmorsxyz`` extension) is written as
         ``solventfilename "name"`` inside the ``%cosmors`` block.

.. note::

   -  For **CPCM** with a named solvent, ``CPCM(solvent_id)`` is written in the route line.

   -  For **CPCMC** (CPCM + COSMO epsilon, replaces old ``COSMO`` keyword removed in ORCA 4.0), ``CPCMC(solvent_id)`` is
      written in the route line.

   -  For **SMD**, ``SMD(solvent_id)`` is written in the route line (canonical ORCA 6.0 form). Any additional options
      (e.g. ``SurfaceType``, SMD descriptors) go into the ``%cpcm`` block via ``-so``.

   -  For **openCOSMO-RS** (``cosmors``), ``COSMORS(solvent_id)`` is written in the route line and a ``%cosmors`` block
      is added with any ``-so`` parameters.

      .. warning::

         **ORCA 6.1 duplicate-keyword guard (openCOSMO-RS only):** ORCA raises an ``INPUT ERROR`` if
         ``COSMORS(solvent_id)`` is on the route line *and* ``solvent "solvent_id"`` also appears in the ``%cosmors``
         block. When ``-si`` / ``solvent_id`` is set, CHEMSMART automatically filters out any ``solvent "..."`` lines
         from the ``%cosmors`` block to prevent this error. Note that ``solventfilename "..."`` is a **different**
         keyword (it specifies the path to a ``.cosmorsxyz`` file) and is **not** filtered.

   -  For a **custom dielectric** (no named solvent), the bare keyword (``CPCM``, ``CPCMC``, ``SMD``, or ``COSMORS``) is
      written in the route line and the dielectric constants go into the corresponding block via ``-so`` (or
      ``custom_solvent`` in the project YAML).

   -  ``-so`` is only applied when a solvent model is active — it is ignored when ``--remove-solvent`` is used.

Supported ``%cpcm`` block options (via ``-so``, for ``cpcm``, ``cpcmc``, and ``smd`` models):

.. list-table::
   :header-rows: 1
   :widths: 30 70

   -  -  Option
      -  Description

   -  -  ``Epsilon <value>``
      -  Static dielectric constant (used for custom/non-named solvents, e.g. ``Epsilon 78.36``)

   -  -  ``Refrac <value>``
      -  Refractive index (e.g. ``Refrac 1.33``)

   -  -  ``SurfaceType <type>``
      -  Cavity surface type: ``gepol_ses``, ``gepol_sas``, ``vdw_gaussian`` (default since ORCA 5), or
         ``gepol_ses_gaussian``

   -  -  ``Rsolv <value>``
      -  Solvent probe radius in Ångström (e.g. ``Rsolv 1.30``)

   -  -  ``MaxIter <n>``
      -  Maximum iterations (e.g. ``MaxIter 100``)

   -  -  ``Tolerance <value>``
      -  Convergence tolerance

   -  -  ``soln``, ``soln25``, ``sola``, ``solb``, ``solg``, ``solc``, ``solh``
      -  SMD solvent descriptors (refractive index, H-bond acidity/basicity, surface tension, aromaticity,
         halogenicity); used only with the ``smd`` model

Supported ``%cosmors`` block options (via ``-so``, for the ``cosmors`` model):

.. note::

   The ORCA ``%cosmors`` keyword for temperature is ``temp`` (lowercase), as listed in the ORCA 6.0 manual.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   -  -  Option
      -  Description

   -  -  ``temp <value>``
      -  Reference temperature in Kelvin (e.g. ``temp 298.15``)

   -  -  ``aeff <value>``
      -  Effective contact area between surface segments in Å² (default: ``5.925``)

   -  -  ``lnalpha <value>``
      -  Logarithm of the misfit prefactor (default: ``0.202``)

   -  -  ``lnchb <value>``
      -  Hydrogen bond (HB) strength parameter (default: ``0.166``)

   -  -  ``chbt <value>``
      -  Parameter for temperature dependence of HB (default: ``1.50``)

   -  -  ``sigmahb <value>``
      -  HB threshold parameter in e/Å² (default: ``9.61e-3``)

   -  -  ``rav <value>``
      -  Radius to average ideal screening charges in Å (default: ``0.50``)

   -  -  ``fcorr <value>``
      -  Parameter adjusted from dielectric screening energies (default: ``2.40``)

   -  -  ``ravcorr <value>``
      -  Additional radius for misfit energy calculation in Å (default: ``1.00``)

   -  -  ``astd <value>``
      -  Standard surface area normalization factor in Å² (default: ``41.624``)

   -  -  ``zcoord <value>``
      -  Coordination number (default: ``10.0``)

   -  -  ``dgsolv_eta <value>``
      -  Offset for solvation energy calculation (default: ``-4.4480``)

   -  -  ``dgsolv_omegaring <value>``
      -  Correction for solvation energy of molecules with rings (default: ``0.2630``)

   -  -  ``dftfunc "name"``
      -  DFT functional for COSMO-RS sub-calculations (default: ``"BP86"``)

   -  -  ``dftbas "name"``
      -  Basis set for COSMO-RS sub-calculations (default: ``"def2-TZVPD"``)

   -  -  ``solvent "name"``
      -  Solvent from the internal COSMO-RS database (e.g. ``solvent "water"``)

   -  -  ``solventfilename "name"``
      -  Name of the ``.cosmorsxyz`` solvent file to read (prefer the ``-sf`` CLI option, which also handles
         auto-conversion from ``.log``/``.out`` files)

   -  -  ``orbs_vac true|false``
      -  Reuse gas-phase orbitals for the conductor calculation (default: ``false``)

Examples:

.. code:: bash

   # CPCM with a named solvent (group-level, applies to all subcommands)
   chemsmart sub orca -p myproject -f molecule.xyz -c 0 -m 1 -sm cpcm -si water sp

   # CPCMC (CPCM + COSMO epsilon) with a named solvent
   chemsmart sub orca -p myproject -f molecule.xyz -c 0 -m 1 -sm cpcmc -si water sp

   # SMD with a named solvent (route: ! SMD(water))
   chemsmart sub orca -p myproject -f molecule.xyz -c 0 -m 1 -sm smd -si water opt

   # SMD with a surface-type option (goes into %cpcm block)
   chemsmart sub orca -p myproject -f molecule.xyz -c 0 -m 1 -sm smd -si water -so 'SurfaceType gepol_ses' opt

   # CPCMC with custom dielectric (no named solvent; replaces old COSMO usage)
   chemsmart sub orca -p myproject -f molecule.xyz -c 0 -m 1 -sm cpcmc -so $'Epsilon 16.7\nRefrac 1.275' sp

   # openCOSMO-RS with a named solvent and temperature (route: ! COSMORS(water))
   chemsmart sub orca -p myproject -f molecule.xyz -c 0 -m 1 -sm cosmors -si water -so 'temp 298.15' sp

   # openCOSMO-RS with a named solvent and a custom .cosmorsxyz file
   chemsmart sub orca -p myproject -f molecule.xyz -c 0 -m 1 -sm cosmors -si water -sf /path/to/water.cosmorsxyz sp

   # openCOSMO-RS with a named solvent and a Gaussian/ORCA output file (auto-converted to .cosmorsxyz)
   chemsmart sub orca -p myproject -f molecule.xyz -c 0 -m 1 -sm cosmors -si water -sf /path/to/water.log sp
   chemsmart sub orca -p myproject -f molecule.xyz -c 0 -m 1 -sm cosmors -si water -sf /path/to/water.out sp

   # Custom dielectric (no named solvent): remove project solvent first, then set custom Epsilon/Refrac
   chemsmart sub orca -p myproject -f molecule.xyz -c 0 -m 1 --remove-solvent sp -sm cpcm -so $'Epsilon 16.7\nRefrac 1.275'

   # Subcommand-level override (overrides group-level solvent)
   chemsmart sub orca -p myproject -f molecule.xyz -sm smd -si water sp -sm smd -si toluene

   # Remove solvent defined in project settings
   chemsmart sub orca -p solv_project -f molecule.xyz -c 0 -m 1 --remove-solvent sp

The SMD example produces:

.. code:: text

   ! SMD(water) B3LYP def2-SVP ...

The SMD + SurfaceType example produces:

.. code:: text

   ! SMD(water) B3LYP def2-SVP ...
   %cpcm
     SurfaceType gepol_ses
   end

The openCOSMO-RS example produces:

.. code:: text

   ! COSMORS(water) B3LYP def2-SVP ...
   %cosmors
     temp 298.15
   end

The openCOSMO-RS with custom solvent file example produces (regardless of whether ``-sf`` points to a ``.cosmorsxyz``,
``.log``, or ``.out`` file — non-``.cosmorsxyz`` files are auto-converted first):

.. code:: text

   ! COSMORS(water) B3LYP def2-SVP ...
   %cosmors
     solventfilename "water"
   end

The custom-dielectric CPCM example produces:

.. code:: text

   ! CPCM B3LYP def2-SVP ...
   %cpcm
     Epsilon 16.7
     Refrac 1.275
   end

***********************
 Available Subcommands
***********************

Structure Optimization
======================

.. list-table::
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``opt``
      -  Geometry optimization
   -  -  ``sp``
      -  Single point calculation

Transition State Search
=======================

.. list-table::
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``ts``
      -  Transition state optimization
   -  -  ``modred``
      -  Modified redundant coordinate optimization
   -  -  ``irc``
      -  Intrinsic reaction coordinate calculations
   -  -  ``scan``
      -  Coordinate scanning
   -  -  ``neb``
      -  Nudged Elastic Band calculations
   -  -  ``qrc``
      -  Quick reaction coordinate calculations

Direct Input
============

.. list-table::
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``inp``
      -  Run ORCA input file as-is

Excited States
==============

.. list-table::
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description

   -  -  ``td``
      -  TDDFT / TDA excited states. Vertical by default; supports optional excited-state ``Opt``/``Freq`` via ``-r``,
         plus SOC coupling and NTO analysis.

TD-specific options (``chemsmart sub orca ... td ...``):

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-n, --nroots``
      -  int
      -  Number of excited states to solve for (``NRoots``, default 3).

   -  -  ``--tda/--no-tda``
      -  bool
      -  Enable/disable the Tamm–Dancoff approximation. Default ``--no-tda`` writes ``TDA false`` so a bare ``td``
         request runs full TDDFT.

   -  -  ``--triplets/--no-triplets``
      -  bool
      -  Solve for triplet excitations. Implicitly enabled by ``--dosoc`` when not set.

   -  -  ``--dosoc/--no-dosoc``
      -  bool
      -  Write ``DoSOC true`` inside ``%tddft``. Does **not** by itself emit a ``%rel`` block.

   -  -  ``--soc-type``
      -  int
      -  Emit a minimal ``%rel`` block with ``SOCType <value>``. If omitted, no ``%rel`` is written even when
         ``--dosoc`` is enabled.

   -  -  ``--printlevel``
      -  int
      -  ``PrintLevel`` inside ``%tddft``.

   -  -  ``--cpcmeq/--no-cpcmeq``
      -  bool
      -  ``CPCMEQ`` inside ``%tddft``. Explicit ``false`` is preserved.

   -  -  ``--nto/--no-nto``
      -  bool
      -  ``DoNTO`` inside ``%tddft``. Passing ``--nto-states`` or ``--nto-thresh`` implicitly enables it.

   -  -  ``--nto-states``
      -  string
      -  Comma / range list of positive 1-based state indices, e.g. ``"1,2,3"`` or ``"1-3"``.

   -  -  ``--nto-thresh``
      -  float
      -  ``NTOThresh``.

   -  -  ``--td-maxiter``
      -  int
      -  ``MaxIter`` inside ``%tddft``. Independent of the SCF ``MaxIter``.

   -  -  ``--td-maxdim``
      -  int
      -  ``MaxDim`` (Davidson subspace).

   -  -  ``--td-etol``
      -  float
      -  ``ETol`` (energy convergence).

   -  -  ``--td-rtol``
      -  float
      -  ``RTol`` (residual convergence).

   -  -  ``--tprint``
      -  float
      -  ``TPrint`` transition print threshold.

   -  -  ``--root``
      -  int
      -  Target excited-state index (``IRoot``). Defaults to ``1`` when an excited-state task is requested (via YAML or
         ``-r``); must satisfy ``1 <= root <= nroots``.

   -  -  ``--root-mult``
      -  singlet | triplet
      -  Target-state multiplicity (``IRootMult``). Defaults to ``singlet`` for excited-state tasks. ``triplet`` implies
         ``--triplets`` and is distinct from ``--triplets`` alone (which only enables solving for triplet excitations).

   -  -  ``--follow-root/--no-follow-root``
      -  bool
      -  ``FollowIRoot`` inside ``%tddft``. Only meaningful when an excited-state ``Opt`` / ``Freq`` / ``NumFreq`` task
         is requested via ``additional_route_parameters`` (YAML or ``-r``).

Unset options are omitted from the input file; ORCA falls back to its own defaults. Only ``NRoots`` and ``TDA`` are
written by default. When neither the project YAML nor the CLI requests ``Opt`` / ``Freq`` / ``NumFreq`` through
``additional_route_parameters``, ``orca td`` runs a vertical calculation: the CLI clears Opt/Freq/NumFreq that would
otherwise be inherited from ``ORCAJobSettings`` defaults or from parsing a ``.log`` / ``.inp`` file, so the route line
carries no task keyword.

Excited-state optimization and frequencies are opt-in via ``additional_route_parameters``, set either on the project
YAML (``td:`` section) or on the ``orca`` group via ``-r`` / ``--additional-route-parameters``. Any ``Opt`` / ``Freq`` /
``NumFreq`` token found there is consumed as a task request; remaining tokens (e.g. ``TightSCF``) are appended to the
``!`` route line with duplicates deduplicated. ``-r`` **replaces** the YAML value in full — there is no per-token merge,
and only that field is replaced (other structured settings such as functional or solvent are unaffected). Use ``-r ''``
to force a vertical TD when the YAML would otherwise request a task; use ``-r TightSCF`` to keep the extras while
dropping any YAML-provided Opt/Freq. ``Opt`` / ``Freq`` / ``NumFreq`` act on the target excited state selected by
``--root`` / ``--root-mult`` (not on the ground state). ORCA's ``DoSOC`` computes singlet–triplet couplings and is not
equivalent to a SOC-based gradient method; ``Opt`` combined with ``--dosoc`` is refused because CHEMSMART does not drive
``SOCGrad`` — use ``chemsmart run orca inp`` for that workflow. SF-TDA, ESD, XAS, TD-specific restarts and spectrum
post-processing are also intentionally out of scope for this subcommand; use ``orca inp``.

Minimal vertical TDDFT (``td --nroots 10``):

.. code::

   %tddft
     NRoots 10
     TDA false
   end

TD + SOC (``td --nroots 10 --dosoc --printlevel 3``):

.. code::

   %tddft
     NRoots 10
     TDA false
     Triplets true
     DoSOC true
     PrintLevel 3
   end

Adding ``--soc-type 3`` additionally writes:

.. code::

   %rel
     SOCType 3
   end

Excited-state optimization on the first singlet root (``-r Opt td --nroots 5 --root 1``) produces a route line beginning
with ``Opt`` and adds ``IRoot 1`` / ``IRootMult singlet`` inside ``%tddft``:

.. code::

   ! Opt CAM-B3LYP def2-SVP defgrid3 SMD(water)
   %tddft
     NRoots 5
     TDA false
     IRoot 1
     IRootMult singlet
   end

Combining ``-r 'Opt Freq'`` / ``-r 'Opt NumFreq'`` writes both keywords onto the route line, and any user-supplied
extras such as ``-r 'TightSCF'`` are preserved (with duplicate task keywords deduplicated).

NTO (``td --nto-states 1,2,3 --nto-thresh 1e-4``) writes ``DoNTO true``, ``NTOStates 1,2,3`` and ``NTOThresh 0.0001``
inside the same single ``%tddft`` block.

.. note::

   SOC support here covers closed-shell (S=0) reference states with singlet–triplet couplings. Combinations with
   UHF/UKS, double-hybrids, or ECPs may need a manual ORCA input via ``orca inp``; CHEMSMART does not auto-inject
   ``ForceECP`` or reshape the reference multiplicity.

QM/MM
=====

``qmmm`` is nested under a parent job type. Use ``<JOBTYPE> qmmm`` where ``<JOBTYPE>`` is ``opt``, ``ts``, ``sp``,
``scan``, ``modred``, ``qrc``, or ``neb``. See :doc:`orca-multiscale-calculations`.

************
 Next Steps
************

For detailed information on each job type:

-  :doc:`orca-structure-optimization`
-  :doc:`orca-transition-state`
-  :doc:`orca-direct-input`
-  :doc:`orca-multiscale-calculations`
