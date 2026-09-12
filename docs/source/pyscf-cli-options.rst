###################
 PySCF CLI Options
###################

This page documents the CLI options available for all PySCF jobs. Use ``chemsmart sub pyscf --help`` for the complete
list.

PySCF is a Python library, not an external binary, so chemsmart runs it by generating a standalone driver script and
executing it with a chosen Python interpreter. GPU4PySCF is an execution engine of the same program, selected with
``--gpu``; it is not a separate program.

*************************
 Basic Command Structure
*************************

.. code:: bash

   chemsmart run [RUN_OPTIONS] pyscf [PYSCF_OPTIONS] <SUBCMD> [SUBCMD_OPTIONS]
   chemsmart sub [SUB_OPTIONS] pyscf [PYSCF_OPTIONS] <SUBCMD> [SUBCMD_OPTIONS]

.. note::

   A project is required. PySCF has no defensible default functional or basis set, so chemsmart will not invent one.
   This is unlike xTB, whose GFN2 default is a complete method.

***************
 PySCF Options
***************

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
      -  Project settings from ``~/.chemsmart/pyscf/*.yaml`` (required)

   -  -  ``-f, --filename``
      -  string
      -  Input file for job preparation

   -  -  ``-l, --label``
      -  string
      -  Custom output filename (without extension)

   -  -  ``-a, --append-label``
      -  string
      -  String to append to the base filename

   -  -  ``-i, --index``
      -  string
      -  Structure index (1-based, default: last structure)

   -  -  ``-t, --title``
      -  string
      -  PySCF job title

   -  -  ``-P, --pubchem``
      -  string
      -  Query structure from PubChem (name, SMILES, CID)

.. note::

   -  ``-p`` uses the project name without the ``.yaml`` extension.
   -  Exactly one of ``-f`` and ``-P`` must be given.
   -  ``-l`` and ``-a`` cannot be used together.

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
      -  Spin multiplicity, given as ``2S+1``

.. note::

   PySCF's own ``spin`` keyword is ``2S``, not ``2S+1``. Always give the multiplicity here and let chemsmart convert it,
   so that the value in the command matches the value in the project YAML and in every other program group.

Method and Basis Set Options
============================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-x, --functional``
      -  string
      -  DFT functional, resolved through libxc (e.g. ``b3lyp``, ``m062x``)

   -  -  ``-A, --ab-initio``
      -  string
      -  Ab initio method; use ``hf`` for a Hartree-Fock reference

   -  -  ``-b, --basis``
      -  string
      -  Basis set in hyphenated PySCF spelling (e.g. ``def2-svp``)

   -  -  ``-ab, --aux-basis``
      -  string
      -  Auxiliary basis for density fitting; omit to let PySCF choose

   -  -  ``-d, --dispersion``
      -  string
      -  Dispersion correction (e.g. ``d3bj``); requires ``pyscf-dispersion``

Give either ``-x`` or ``-A``, never both. Basis names follow PySCF and libcint spelling, which is hyphenated:
``def2-svp``, not Gaussian's ``def2svp``. Functional names go through libxc, which accepts ``m062x`` as well as
``M06-2X``.

.. warning::

   PySCF's ``b3lyp`` is the Gaussian VWN3 definition. ORCA's ``B3LYP`` is the VWN5 variant, which PySCF spells
   ``b3lyp5``. For water in def2-SVP the two differ by roughly 23 kcal/mol, so a PySCF ``b3lyp`` energy must not be
   compared with an ORCA ``B3LYP`` energy. chemsmart logs a warning when it resolves ``b3lyp``, but it cannot know which
   convention you meant.

Double-hybrid functionals are rejected, since their perturbative correlation term would not be applied. ``-A`` accepts
only ``hf``; post-Hartree-Fock correlation methods are not part of this integration.

Examples:

.. code:: bash

   # Use a different functional
   chemsmart sub pyscf -p test -f molecule.xyz -x m062x opt

   # Use a Hartree-Fock reference
   chemsmart sub pyscf -p test -f molecule.xyz -A hf -b def2-svp sp

   # Add a dispersion correction
   chemsmart sub pyscf -p test -f molecule.xyz -x b3lyp -d d3bj opt

SCF and Grid Options
====================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``--defgrid``
      -  string
      -  DFT integration grid: ``defgrid1``, ``defgrid2``, or ``defgrid3``

   -  -  ``--scf-tol``
      -  float
      -  SCF convergence tolerance (``mf.conv_tol``)

   -  -  ``--scf-maxiter``
      -  int
      -  Maximum SCF cycles (``mf.max_cycle``)

   -  -  ``--density-fit/--no-density-fit``
      -  bool
      -  Use resolution-of-the-identity density fitting

The three grids map to ``(radial, angular)`` point counts of ``(50, 194)``, ``(75, 302)``, and ``(99, 590)``.
``--defgrid`` is DFT-only and cannot be combined with ``-A``. ``-ab`` cannot be combined with ``--no-density-fit``.

Optimization Options
====================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``--opt-solver``
      -  string
      -  Optimization backend: ``geometric``, ``berny``, or ``ase``

   -  -  ``--opt-maxsteps``
      -  int
      -  Maximum geometry optimization steps

The three backends do not share convergence criteria, so the one that was used is recorded in the results file. The
selected backend must be importable from the interpreter that runs the job.

Engine Options
==============

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``--gpu/--no-gpu``
      -  bool
      -  Run on GPU4PySCF or on CPU

When neither flag is given, the engine follows ``SERVER.NUM_GPUS`` from the server YAML. GPU and CPU results are not
bit-identical, so the resolved engine is recorded in the results file. See :doc:`pyscf-gpu-acceleration`.

Solvent Options
===============

.. list-table::
   :header-rows: 1
   :widths: 35 15 50

   -  -  Option
      -  Type
      -  Description

   -  -  ``-sm, --solvent-model``
      -  string
      -  ``pcm``, ``iefpcm``, ``cpcm``, ``cosmo``, ``ssvpe``, or ``smd``

   -  -  ``-si, --solvent-id``
      -  string
      -  Solvent name, which must exist in ``pyscf.solvent.smd.solvent_db``

   -  -  ``--remove-solvent/--no-remove-solvent``
      -  bool
      -  Drop the solvent model from the project settings (default: disabled)

.. note::

   A solvent model always needs a solvent identifier. chemsmart resolves the dielectric explicitly, because PySCF's PCM
   otherwise falls back to water regardless of which solvent was requested.

Examples:

.. code:: bash

   # SMD water on a gas-phase project
   chemsmart sub pyscf -p test -f molecule.xyz -sm smd -si water sp

   # C-PCM toluene
   chemsmart sub pyscf -p test -f molecule.xyz -sm cpcm -si toluene sp

   # Drop the solvent the project specifies
   chemsmart sub pyscf -p test -f molecule.xyz --remove-solvent sp

***********************
 Project YAML Contract
***********************

A PySCF project uses two sections. ``gas`` feeds the ``opt`` and ``hess`` jobs; ``solv`` feeds ``sp``. That split
matches the usual workflow: optimize and take frequencies at one level, then refine the single point in solvent.

.. code:: yaml

   gas:
     functional: b3lyp
     basis: def2-svp
     freq: True
     density_fit: True
     aux_basis: null
     defgrid: defgrid2
     scf_tol: 1.0e-9
     scf_maxiter: 100
     opt_solver: geometric
     opt_maxsteps: 100
   solv:
     functional: b3lyp
     basis: def2-tzvp
     freq: False
     density_fit: True
     defgrid: defgrid2
     scf_tol: 1.0e-9
     solvent_model: smd
     solvent_id: water

A key that does not exist on the PySCF settings raises an error at load time and lists every valid key, so a typo stops
the run instead of quietly applying a default. Settings inherited from the Gaussian and ORCA schemas are refused rather
than ignored, including ``gen_genecp_file``, ``heavy_elements``, ``numfreq``, ``additional_route_parameters``,
``route_to_be_written``, ``input_string``, and ``custom_solvent``.

**********************
 Server Configuration
**********************

Readiness belongs to the interpreter that runs the calculation, not to the process running chemsmart. The ``PYSCF``
block in the server YAML names that interpreter's ``bin`` directory:

.. code:: yaml

   PYSCF:
     EXEFOLDER: ~/miniconda3/envs/pyscf-gpu/bin
     LOCAL_RUN: True
     SCRATCH: False

Omit ``EXEFOLDER`` when PySCF is installed in the same environment as chemsmart; the running interpreter is used. The
block itself is optional. ``SCRATCH`` is off because a PySCF job writes a handful of small files, so scratch storage
adds nothing.

***********************
 Available Subcommands
***********************

.. list-table::
   :header-rows: 1
   :widths: 20 80

   -  -  Subcommand
      -  Description
   -  -  ``sp``
      -  Single point energy, using the ``solv`` project section
   -  -  ``opt``
      -  Geometry optimization, using the ``gas`` project section
   -  -  ``hess``
      -  Analytic Hessian and harmonic frequencies, using the ``gas`` section

Every subcommand also accepts ``-S/--skip-completed`` and ``-R/--no-skip-completed``. Use ``-R`` to rerun a job that has
already completed.

************
 Next Steps
************

For detailed information on each job type:

-  :doc:`pyscf-structure-optimization`
-  :doc:`pyscf-frequency-calculations`
-  :doc:`pyscf-gpu-acceleration`
-  :doc:`pyscf-results`
