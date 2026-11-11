######################################
 Electronic Structure Methods (PySCF)
######################################

This page covers the correlated methods available for PySCF jobs.

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
