###################
 PySCF CLI Options
###################

ChemSmart exposes PySCF 2.14.0 through the same ``run`` and ``sub`` command families as its executable-backed programs.
The executable CPU surface is ``sp``, ``opt``, ``hess`` and ``td``: ground-state single points, optimisations and
Hessians on a Hartree--Fock or DFT reference; TDA/TDDFT vertical excitations on a closed-shell reference (singlet or
triplet manifold) or on an open-shell reference (the one ``unrestricted`` manifold), gas phase or with an implicit
solvent; optimisation on an excited root of that manifold; and MP2, CCSD and CCSD(T) as ``ab_initio`` methods on a
Hartree--Fock reference. GPU4PySCF 1.8.0 is an execution engine of the PySCF program; it is not a separate program.
GPU4PySCF configuration and safe preview are available, but this release does not claim a qualified Agent GPU run.

*************************
 Basic Command Structure
*************************

.. code:: bash

   chemsmart run [RUN_OPTIONS] pyscf -p PROJECT -f GEOMETRY [PYSCF_OPTIONS] sp
   chemsmart run [RUN_OPTIONS] pyscf -p PROJECT -f GEOMETRY [PYSCF_OPTIONS] opt
   chemsmart run [RUN_OPTIONS] pyscf -p PROJECT -f GEOMETRY [PYSCF_OPTIONS] hess
   chemsmart run [RUN_OPTIONS] pyscf -p PROJECT -f GEOMETRY [PYSCF_OPTIONS] td

The same program and leaf commands are available below ``chemsmart sub``. PySCF requires a validated project YAML;
ChemSmart does not invent a default method or basis.

.. warning::

   ChemSmart generates a standalone Python script and a structured HDF5 result. The generated script is an execution
   artifact, not a supported user-editing interface. Change the project YAML or CLI options and let ChemSmart regenerate
   it.

***********************
 Project YAML Contract
***********************

A PySCF project uses stage-specific ``sp``, ``opt``, ``hess`` and ``td`` sections. A stage does not inherit scientific
settings from another stage.

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
   td:
     functional: b3lyp
     basis: def2-svp
     density_fit: true
     response_method: tda
     state_manifold: singlet
     nstates: 10
     td_max_cycle: 100
     freq: false

An ``opt`` section carrying ``excited_state_root`` optimises on that root of the manifold it also names, with the same
response fields as ``td``; an ``sp`` or ``opt`` section naming a correlated ``ab_initio`` method may carry
``frozen_core`` and, for coupled cluster, ``cc_max_cycle``:

.. code:: yaml

   opt:
     functional: b3lyp
     basis: def2-svp
     response_method: tda
     state_manifold: singlet
     nstates: 3
     excited_state_root: 1
     opt_solver: geometric
     opt_maxsteps: 100
   sp:
     ab_initio: ccsd(t)
     basis: cc-pvtz
     frozen_core: auto
     cc_max_cycle: 100

Roots are ascending ordinals within the requested manifold at the geometry the stage ran on, never state labels: a
``td`` result reports how many roots were requested and how many PySCF's positive-eigenvalue filter kept, an
excited-root ``opt`` follows its root by index at every step and records, at the reached geometry, that root's gap to
the ground state and to its neighbouring root, and the reported numbers carry no threshold. A correlated stage records
its reference and correlation energies separately (``correlation_energy`` is the final method's whole correlation,
triples included) and the number of orbitals it left uncorrelated; PySCF correlates every electron unless
``frozen_core`` says otherwise, where ORCA and Gaussian freeze the core by default, and ``auto`` applies PySCF's own
chemical-core rule. An excited-root optimisation is gas phase only, a correlated method takes no density fitting or
implicit solvent in this release, and ``ccsd(t)`` optimises nothing; each of these is refused when the project is
validated, naming its route.

Use either ``functional`` or ``ab_initio`` in a stage, never both. Unknown keys and inherited Gaussian/ORCA-only
settings are rejected. In particular, native route text, ``modred``, semiempirical settings, arbitrary mixed-basis text,
and unsupported forces are not silently ignored.

Program-Level Options
=====================

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
      -  Spin multiplicity ``2S+1``. ChemSmart converts it to PySCF ``spin=2S``.

   -  -  ``-A, --ab-initio``
      -  hf, mp2, ccsd, ccsd(t)
      -  Hartree--Fock, or a correlated method on a Hartree--Fock reference; never together with ``--functional``.

   -  -  ``--frozen-core``
      -  non-negative integer or ``auto``
      -  Orbitals a correlated method leaves uncorrelated. Omitted, PySCF correlates every electron; ``auto`` applies
         PySCF's chemical-core rule. The applied count is recorded in the result.

   -  -  ``--cc-max-cycle``
      -  positive integer
      -  Coupled-cluster amplitude iteration cap (PySCF default 50); the repair control behind unconverged amplitudes.

   -  -  ``--nstates``
      -  positive integer
      -  Roots computed per manifold, for ``td`` or an excited-root ``opt``.

   -  -  ``--response-method``
      -  tda/tddft
      -  Tamm--Dancoff approximation or the full response.

   -  -  ``--state-manifold``
      -  singlet/triplet/unrestricted
      -  Singlet or triplet excitations of a closed-shell reference; ``unrestricted`` is the one manifold of an
         open-shell reference.

   -  -  ``--excited-root``
      -  positive integer, at most ``nstates``
      -  ``opt`` only: optimise on that root of the manifold, by index.

   -  -  ``--td-max-cycle``
      -  positive integer
      -  Davidson iteration cap for the response solver (PySCF default 100); the repair control behind an unconverged
         root.

   -  -  ``-x, --functional``
      -  libxc functional
      -  DFT functional. Program-specific definitions remain scientifically distinct.

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
      -  Explicit GPU4PySCF request or CPU selection. Missing GPU support never falls back to CPU.

**********************
 Server Configuration
**********************

PySCF is a Python library, so readiness is bound to the exact compute interpreter rather than to the controller process.
A server block may point ``EXEFOLDER`` at the ``bin`` directory that owns the required Python.

.. code:: yaml

   PYSCF:
     EXEFOLDER: /path/to/pyscf-environment/bin
     LOCAL_RUN: true
     SCRATCH: false

Before execution, ChemSmart records the interpreter and required dependency versions. A GPU request additionally
requires matching GPU4PySCF, CuPy, CUDA, cuTENSOR, driver, and device observations. Merely declaring ``NUM_GPUS`` does
not establish GPU readiness, and a green preview is not a claim that the GPU path has been release-qualified.

************************
 Results and Completion
************************

``LABEL.h5`` is the machine-readable result contract. The human-readable ``LABEL.out`` is retained as evidence but is
not the authority for completion. The HDF5 record binds requested and applied settings, geometry, charge, multiplicity,
engine, environment, convergence, properties, and artifact hashes. A process exit code of zero is insufficient when
preflight, provenance, convergence, or required-property validation is red.

The artifact carries two structures and every quantity belongs to the second. ``spec/positions`` is the geometry the
calculation was handed; ``results/positions`` is where it ended. For an optimisation the driver re-converges the SCF on
the final geometry, from the optimiser's own last density, before any energy, orbital, dipole, population, spin
expectation or frequency is read, so a PySCF result has one structure and it is the final one. For ``sp`` and ``hess``
the two coincide by construction and the validator enforces it. An optimisation that stops on its step limit still
writes the last geometry the optimiser evaluated, never the input; its receipt records the failure, and the structure
stays readable. The typed analysis layer serves the supplied structure as ``supplied_positions`` and the final one as
``positions`` (and, for ``opt``, ``reached_positions`` and ``converged``).

The ``hess`` leaf uses the supplied geometry without optimizing it. In a multi-stage workflow, bind it to the exact
optimized-geometry artifact from a validated ``opt`` node rather than reusing the initial geometry. A Hessian's
frequencies are projected free of translations and rotations, so a spectrum with no imaginary mode does not by itself
prove the geometry is a stationary point; the driver records the gradient at the Hessian geometry beside the frequencies
(``results/forces``), the mass table behind the frequencies (isotope-averaged) and the tolerance behind the detected
point group, so a free energy derived from the result states its conventions. A Hessian that includes a D3 or D4
dispersion correction, or an SMD cavity term, is analytic except those blocks, which PySCF evaluates by finite
differences.

For an open-shell reference the driver also records per-atom Mulliken spin populations; a closed-shell result marks the
property as not applicable rather than reporting zeros.

A ``td`` result records its excitation energies (hartree, ascending within the manifold), oscillator strengths,
transition dipole moments (debye), the manifold multiplicity of each root for a closed-shell reference, and each root's
convergence, beside the count of roots requested and obtained and the iteration cap applied; a solvated spectrum records
the static dielectric applied to the reference and the response dielectric PySCF actually used, which is 1.78 for every
solvent in this build. A run with an unconverged root or an unconverged amplitude set records the failure in its receipt
and stays inspectable; the typed analysis layer serves no quantity from it. An excited-root optimisation re-converges
the SCF and re-evaluates the spectrum at the reached geometry, so its excitation energies, like every other quantity,
belong to the final structure, and its ``energy`` is the followed root's total there while ``scf_energy`` is the
reference's. A correlated result's ``energy`` is the correlated total and ``reference_energy``, ``correlation_energy``,
``ccsd_correlation_energy`` and ``triples_correction`` are the components PySCF returned. The dipole moment,
populations, orbital energies and spin expectation on an excited-root or correlated result belong to the SCF reference,
and the typed analysis layer says so beside each value.

**********************
 Unsupported Requests
**********************

The executable integration does not offer transition-state search, IRC, scan, QMMM/ONIOM, NEB, double hybrids, arbitrary
mixed basis/ECP input, unsupported constraints, excited-state or correlated Hessians, CCSD(T) gradients (present
upstream, not audited through this driver), EOM-CCSD, CASSCF, density fitting or implicit solvent with a correlated
method, or a solvated excited-state gradient. These requests must block; they must not be rewritten as a superficially
similar PySCF calculation. PySCF 2.14 has no analytic Hessian for any ROHF reference, which its ``scf.HF`` selects for
every one-electron system, nor for an open-shell reference under a non-local-correlation functional; both are refused at
preflight rather than inside the engine. Only the ``geometric`` optimiser is installed in this host's compute
environment; ``berny`` and ``ase`` are refused by the environment probe when absent.
