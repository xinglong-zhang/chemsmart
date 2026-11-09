######################################
 Electronic Structure Methods (PySCF)
######################################

This page covers the correlated methods available for PySCF jobs, and the stability analysis of the SCF reference.

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
``--no-density-fit`` or set ``density_fit: false``. Dispersion corrections and solvent models cannot be combined with
these methods. PySCF has no analytic Hessian for them, so ``hess``, ``freq: True`` and every other job that needs a
Hessian are refused; set ``freq: false`` in the ``gas`` section to optimize. Optimizations use the ``geometric``
optimizer.

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
