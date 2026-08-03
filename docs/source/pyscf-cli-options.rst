###################
 PySCF CLI Options
###################

ChemSmart exposes PySCF 2.14.0 through the same ``run`` and ``sub`` command
families as its executable-backed programs.  The maintained surface is
deliberately limited to ``sp``, ``opt``, and ``hess``.  GPU4PySCF 1.8.0 is an
execution engine of the PySCF program; it is not a separate program.

*************************
 Basic Command Structure
*************************

.. code:: bash

   chemsmart run [RUN_OPTIONS] pyscf -p PROJECT -f GEOMETRY [PYSCF_OPTIONS] sp
   chemsmart run [RUN_OPTIONS] pyscf -p PROJECT -f GEOMETRY [PYSCF_OPTIONS] opt
   chemsmart run [RUN_OPTIONS] pyscf -p PROJECT -f GEOMETRY [PYSCF_OPTIONS] hess

The same program and leaf commands are available below ``chemsmart sub``.
PySCF requires a validated project YAML; ChemSmart does not invent a default
method or basis.

.. warning::

   ChemSmart generates a standalone Python script and a structured HDF5
   result.  The generated script is an execution artifact, not a supported
   user-editing interface.  Change the project YAML or CLI options and let
   ChemSmart regenerate it.

**********************
 Project YAML Contract
**********************

A PySCF project has exactly three stage-specific sections.  A stage does not
inherit scientific settings from another stage.

.. code:: yaml

   sp:
     ab_initio: hf
     functional: null
     basis: def2-svp
     density_fit: false
     freq: false
   opt:
     ab_initio: hf
     functional: null
     basis: def2-svp
     density_fit: false
     opt_solver: geometric
     opt_maxsteps: 100
     freq: false
   hess:
     ab_initio: hf
     functional: null
     basis: def2-svp
     density_fit: false
     freq: true

Use either ``functional`` or ``ab_initio`` in a stage, never both.  Unknown
keys and inherited Gaussian/ORCA-only settings are rejected.  In particular,
native route text, ``modred``, semiempirical settings, arbitrary mixed-basis
text, and unsupported forces are not silently ignored.

*********************
 Program-Level Options
*********************

.. list-table::
   :header-rows: 1
   :widths: 32 18 50

   -  -  Option
      -  Value
      -  Meaning
   -  -  ``-p, --project``
      -  name or YAML path
      -  Required stage-specific project settings.
   -  -  ``-f, --filename``
      -  molecular artifact
      -  Geometry source; ``--index`` remains 1-based.
   -  -  ``-c, --charge``
      -  integer
      -  Molecular charge override.
   -  -  ``-m, --multiplicity``
      -  positive integer
      -  Spin multiplicity ``2S+1``.  ChemSmart converts it to PySCF ``spin=2S``.
   -  -  ``-A, --ab-initio``
      -  method
      -  Hartree--Fock reference; use ``hf`` for the maintained v1 surface.
   -  -  ``-x, --functional``
      -  libxc functional
      -  DFT functional.  Program-specific definitions remain scientifically distinct.
   -  -  ``-b, --basis``
      -  basis name
      -  PySCF basis spelling, for example ``def2-svp``.
   -  -  ``-ab, --aux-basis``
      -  basis name
      -  Auxiliary basis used only with density fitting.
   -  -  ``--density-fit/--no-density-fit``
      -  boolean
      -  Enable or disable density fitting.
   -  -  ``--scf-tol``
      -  float
      -  SCF convergence tolerance.
   -  -  ``--scf-maxiter``
      -  positive integer
      -  Maximum SCF cycles.
   -  -  ``--defgrid``
      -  defgrid1/2/3
      -  ChemSmart's documented PySCF grid mapping; it is not an ORCA-grid equivalence claim.
   -  -  ``--opt-solver``
      -  geometric/berny/ase
      -  Geometry optimizer; the selected dependency must exist in the compute interpreter.
   -  -  ``--opt-maxsteps``
      -  positive integer
      -  Geometry-optimization step ceiling.
   -  -  ``-sm, --solvent-model``
      -  PCM-family or SMD
      -  Implicit-solvent implementation.
   -  -  ``-si, --solvent-id``
      -  solvent name
      -  Solvent identity resolved by the compute environment.
   -  -  ``--gpu/--no-gpu``
      -  boolean
      -  Explicit GPU4PySCF request or CPU selection.  Missing GPU support never falls back to CPU.

**********************
 Server Configuration
**********************

PySCF is a Python library, so readiness is bound to the exact compute
interpreter rather than to the controller process.  A server block may point
``EXEFOLDER`` at the ``bin`` directory that owns the required Python.

.. code:: yaml

   PYSCF:
     EXEFOLDER: /path/to/pyscf-environment/bin
     LOCAL_RUN: true
     SCRATCH: false

Before execution, ChemSmart records the interpreter and required dependency
versions.  A GPU request additionally requires matching GPU4PySCF, CuPy, CUDA,
cuTENSOR, driver, and device observations.  Merely declaring ``NUM_GPUS`` does
not establish GPU readiness.

************************
 Results and Completion
************************

``LABEL.h5`` is the machine-readable result contract.  The human-readable
``LABEL.out`` is retained as evidence but is not the authority for completion.
The HDF5 record binds requested and applied settings, geometry, charge,
multiplicity, engine, environment, convergence, properties, and artifact
hashes.  A process exit code of zero is insufficient when preflight,
provenance, convergence, or required-property validation is red.

The ``hess`` leaf uses the supplied geometry without optimizing it.  In a
multi-stage workflow, bind it to the exact optimized-geometry artifact from a
validated ``opt`` node rather than reusing the initial geometry.

*********************
 Unsupported Requests
*********************

The initial integration does not offer transition-state search, IRC, TDDFT,
scan, QMMM/ONIOM, NEB, post-HF correlation, double hybrids, arbitrary mixed
basis/ECP input, or unsupported constraints.  These requests must block; they
must not be rewritten as a superficially similar PySCF calculation.
