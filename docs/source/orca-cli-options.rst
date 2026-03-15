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

.. note::

   -  ``-p`` uses the project name without the ``.yaml`` extension.
   -  ``-f`` accepts various formats: ``.xyz``, ``.com``, ``.gjf``, ``.log``, ``.inp``, or ``.out``.

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

   -  -  ``-a, --aux-basis``
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
      -  Implicit solvent model: ``cpcm`` (CPCM with CPCM epsilon), ``cpcmc`` (CPCM with COSMO epsilon;
         replaces the legacy COSMO model removed in ORCA 4.0), ``smd`` (Minnesota SMD), or ``cosmors``
         (openCOSMO-RS interface)

   -  -  ``-si, --solvent-id``
      -  string
      -  Named solvent identifier (e.g. ``water``, ``toluene``, ``cyclohexane``). Omit when using a fully custom
         dielectric via ``-so``

   -  -  ``-so, --solvent-options``
      -  string
      -  Additional parameters for the model's solvent block (see tables below), newline-separated for multiple options

.. note::

   -  For **CPCM** with a named solvent, ``CPCM(solvent_id)`` is written in the route line.
   -  For **CPCMC** (CPCM + COSMO epsilon, replaces old ``COSMO`` keyword removed in ORCA 4.0),
      ``CPCMC(solvent_id)`` is written in the route line.
   -  For **SMD**, ``SMD(solvent_id)`` is written in the route line (canonical ORCA 6.0 form).
      Any additional options (e.g. ``SurfaceType``, SMD descriptors) go into the ``%cpcm`` block via ``-so``.
   -  For **openCOSMO-RS** (``cosmors``), ``COSMORS(solvent_id)`` is written in the route line and a
      ``%cosmors`` block is added with any ``-so`` parameters.
   -  For a **custom dielectric** (no named solvent), the bare keyword (``CPCM``, ``CPCMC``, ``SMD``, or
      ``COSMORS``) is written in the route line and the dielectric constants go into the corresponding block
      via ``-so`` (or ``custom_solvent`` in the project YAML).
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
      -  Cavity surface type: ``gepol_ses``, ``gepol_sas``, ``vdw_gaussian`` (default since ORCA 5),
         or ``gepol_ses_gaussian``
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

.. list-table::
   :header-rows: 1
   :widths: 30 70

   -  -  Option
      -  Description
   -  -  ``Temperature <value>``
      -  Reference temperature in Kelvin (e.g. ``Temperature 298.15``)
   -  -  ``dftfunc <name>``
      -  DFT functional for COSMO-RS sub-calculations (default: ``BP86``)
   -  -  ``dftbas <name>``
      -  Basis set for COSMO-RS sub-calculations (default: ``def2-TZVPD``)

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
   chemsmart sub orca -p myproject -f molecule.xyz -c 0 -m 1 -sm cosmors -si water -so 'Temperature 298.15' sp

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
     Temperature 298.15
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

Direct Input
============

.. list-table::
   :header-rows: 1
   :widths: 15 85

   -  -  Subcommand
      -  Description
   -  -  ``inp``
      -  Run ORCA input file as-is

************
 Next Steps
************

For detailed information on each job type:

-  :doc:`orca-structure-optimization`
-  :doc:`orca-transition-state`
-  :doc:`orca-direct-input`
