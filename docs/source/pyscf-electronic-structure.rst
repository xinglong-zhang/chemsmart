######################################
 Electronic Structure Methods (PySCF)
######################################

This page covers the correlated and multiconfigurational methods available for PySCF jobs, and the stability analysis of
the SCF reference.

The method is chosen with ``-A`` or the ``ab_initio`` key of the project YAML, as for Hartree-Fock. For the complete
option reference, see :doc:`pyscf-cli-options`.

********************
 Correlated Methods
********************

Run MP2, CCSD or CCSD(T) on the Hartree-Fock reference, which is RHF for a closed shell and UHF otherwise.

.. code:: bash

   chemsmart sub [OPTIONS] pyscf [PYSCF_OPTIONS] -A METHOD sp

``METHOD`` is ``mp2``, ``ccsd`` or ``ccsd(t)``. All three give single-point energies. ``mp2`` and ``ccsd`` also give
forces (``--forces``) and geometry optimizations on their own analytic gradient; ``ccsd(t)`` gives energies only.

Correlation Options
===================

These keys go in a ``correlation`` block of the project section:

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   -  -  Key
      -  Type
      -  Description

   -  -  ``frozen_core``
      -  bool
      -  Freeze the core orbitals (default: ``true``). ``false`` correlates every electron.

   -  -  ``max_cycle``
      -  int
      -  Maximum CCSD iterations (default: PySCF's)

   -  -  ``conv_tol``
      -  float
      -  CCSD energy convergence threshold (default: PySCF's). Not used by ``mp2``, which is not iterative.

The frozen core is PySCF's chemical core, which matches ORCA's default for H–Zn and for the main-group elements Ga–Xe.
Electrons replaced by an effective core potential count as frozen, and the number of frozen orbitals is recorded.
Gaussian's default ``FC`` freezes more orbitals for Li, Be, Sc–Zn and, with Pople basis sets, Ga–Kr, so set the same
frozen orbitals in Gaussian before comparing correlation energies for those elements.

Density fitting must be off, because PySCF would otherwise run the density-fitted variant of the method: pass
``--no-density-fit`` or set ``density_fit: false``. Dispersion corrections, solvent models and the GPU engine cannot be
combined with these methods. PySCF has no analytic Hessian for them, so ``hess``, ``freq: True`` and every other job
that needs a Hessian are refused; set ``freq: false`` in the ``gas`` section to optimize. Optimizations use the
``geometric`` optimizer.

The reported energy is the correlated total, as in Gaussian's output. The Hartree-Fock energy, the correlation energy
and, for CCSD(T), the CCSD part and the triples correction are stored beside it. If the CCSD iterations do not converge,
no correlated energy is reported and the job fails. Mulliken charges and the dipole moment are not reported, since only
the reference density is available.

Basic Usage
===========

Standard MP2 single point:

.. code:: bash

   chemsmart sub pyscf -p project -f molecule.xyz -A mp2 -b cc-pvdz --no-density-fit sp

CCSD(T) single point, with the quotes the shell needs:

.. code:: bash

   chemsmart sub pyscf -p project -f molecule.xyz -A "ccsd(t)" -b cc-pvtz --no-density-fit sp

CCSD geometry optimization, with ``freq: false`` in the ``gas`` section:

.. code:: bash

   chemsmart sub pyscf -p project -f molecule.xyz -A ccsd -b cc-pvdz --no-density-fit opt

The same calculation set in the project YAML, correlating every electron:

.. code:: yaml

   solv:
     ab_initio: ccsd(t)
     basis: cc-pvtz
     density_fit: false
     correlation:
       frozen_core: false
       conv_tol: 1.0e-8

******************************
 Multiconfigurational Methods
******************************

Run CASCI or CASSCF single points, optionally with strongly contracted NEVPT2, in an active space described in the
project YAML.

.. code:: bash

   chemsmart sub [OPTIONS] pyscf [PYSCF_OPTIONS] -A METHOD sp

``METHOD`` is ``casci`` or ``casscf``. The reference is RHF for a closed shell and ROHF otherwise.

Active Space Options
====================

These keys go in an ``active_space`` block of the project section:

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   -  -  Key
      -  Type
      -  Description

   -  -  ``ncas``
      -  int
      -  Number of active orbitals

   -  -  ``nelecas``
      -  int or list
      -  Number of active electrons, or ``[n_alpha, n_beta]``

   -  -  ``orbitals``
      -  block
      -  How the active orbitals are chosen (required; see below)

   -  -  ``nroots``
      -  int
      -  Number of roots (default: 1). More than one averages the orbitals over the roots.

   -  -  ``weights``
      -  list
      -  State-averaging weights (default: equal weights)

   -  -  ``nevpt2``
      -  bool
      -  Add the strongly contracted NEVPT2 correction to every root (default: ``false``)

   -  -  ``max_cycle_macro``
      -  int
      -  Maximum CASSCF iterations (default: PySCF's)

   -  -  ``conv_tol``
      -  float
      -  CASSCF energy convergence threshold (default: PySCF's)

No orbitals are chosen implicitly. ``orbitals.selection`` takes one of three values:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   -  -  Selection
      -  Active orbitals

   -  -  ``mo_list``
      -  The reference orbitals listed in ``mo_list``, numbered from ``index_base``, which must be given: ``1`` if the
         first orbital is 1, as in Gaussian, or ``0`` if it is 0, as in ORCA

   -  -  ``avas``

      -  Orbitals projected onto the atomic orbitals in ``ao_labels``, for example ``['Fe 3d', 'C 2pz']``, with an
         optional ``threshold`` (default: 0.2). AVAS sets the number of active orbitals and electrons, so ``ncas`` and
         ``nelecas`` are left out.

   -  -  ``frontier``
      -  The ``ncas`` orbitals around the highest occupied orbital

All roots have the multiplicity of the molecule. A spin penalty keeps roots of other multiplicities out, and a root
whose ⟨S²⟩ still differs fails the job. With ``nevpt2: true``, each root receives the NEVPT2 correction computed from
its own CI vector in the final CASSCF orbitals, including state-averaged ones.

The reported energy is that of the lowest root: the CASSCF or CASCI energy, plus its NEVPT2 correction when ``nevpt2``
is true. Every root's energy, NEVPT2 correction, ⟨S²⟩, natural-orbital occupations and largest CI coefficients are
stored beside it. If CASSCF does not converge, no energy is reported and the job fails.

Density fitting, dispersion corrections, solvent models and the GPU engine cannot be combined with these methods, and
gradients, optimizations and frequencies are not available. No point group is imposed, so NEVPT2 corrections for
molecules with degenerate orbitals, such as N₂, can differ slightly between runs; a tighter ``conv_tol`` reduces the
spread.

Basic Usage
===========

CASSCF(6,6) with the frontier orbitals, in the ``solv`` section:

.. code:: yaml

   solv:
     ab_initio: casscf
     basis: cc-pvdz
     density_fit: false
     active_space:
       ncas: 6
       nelecas: 6
       orbitals:
         selection: frontier

.. code:: bash

   chemsmart sub pyscf -p project -f molecule.xyz sp

Two-root state-averaged CASSCF on listed orbitals, with NEVPT2:

.. code:: yaml

   solv:
     ab_initio: casscf
     basis: cc-pvdz
     density_fit: false
     active_space:
       ncas: 6
       nelecas: 6
       orbitals:
         selection: mo_list
         mo_list: [5, 6, 7, 8, 9, 10]
         index_base: 1
       nroots: 2
       nevpt2: true

************************
 SCF Stability Analysis
************************

Check whether the converged SCF reference is a minimum in orbital-rotation space.

.. code:: bash

   chemsmart sub [OPTIONS] pyscf [PYSCF_OPTIONS] --scf-stability SUBCMD

A converged SCF need not be a minimum. ``--scf-stability``, or ``scf_stability: true`` in the project YAML, asks PySCF
three questions about the final reference: whether it is stable within its own family (internal), towards complex
orbitals (real to complex), and towards the next less constrained family (external: RHF/RKS to UHF/UKS, or UHF/UKS to
GHF/GKS). Each answer is recorded with its lowest orbital-rotation Hessian eigenvalue, which is comparable only within
the same question.

The analysis only observes. An instability is reported but not followed, and the job's convergence and termination are
unchanged. If you want the lower solution, start a new calculation from it, for example with a different multiplicity or
a broken-symmetry guess.

Basic Usage
===========

Single point with a stability analysis:

.. code:: bash

   chemsmart sub pyscf -p project -f molecule.xyz --scf-stability sp

The answers are stored under ``results/scf_stability``. See :doc:`pyscf-results`.
