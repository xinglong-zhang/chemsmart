#################################
 Transition State Search (PySCF)
#################################

This page covers transition state searches using PySCF and the geomeTRIC optimizer.

The search reads the ``gas`` section of the project YAML, like ``opt``, and needs ``opt_solver: geometric``. It runs
with Hartree-Fock and DFT. MP2 and coupled cluster have no analytic Hessian in PySCF and are refused. For the complete
option reference, see :doc:`pyscf-cli-options`.

*************************
 Transition State Search
*************************

Find a first-order saddle point from a guess structure that lies near it.

.. code:: bash

   chemsmart sub [OPTIONS] pyscf [PYSCF_OPTIONS] ts

The search uses geomeTRIC's partitioned rational function optimization (P-RFO), starting from the analytic Hessian of
the guess.

A converged search has found a stationary point, not necessarily a transition state: from a poor guess it can fall into
a minimum. With ``freq: True`` in the ``gas`` section, the Hessian at the result counts the imaginary modes, and exactly
one is required. Any other number fails the job. The count is recorded, and the structure and its frequencies are kept
for inspection. Without ``freq`` the number of imaginary modes is not known and is recorded as null.

Every negative frequency counts as an imaginary mode, with no threshold. This is the same rule that CHEMSMART's
thermochemistry applies, so a second-order saddle point is never reported as a transition state. A spurious mode of a
few cm⁻¹ from a coarse DFT grid also fails the job; if the recorded frequencies show this, rerun with a finer grid such
as ``--defgrid defgrid3``.

When the job succeeds with ``freq: True``, thermochemistry treats the imaginary mode as the reaction coordinate and
leaves it out of the vibrational partition function.

Basic Usage
===========

Standard transition state search:

.. code:: bash

   chemsmart sub pyscf -p test -f guess.xyz ts

Results are written to ``guess_ts_gas_phase.h5``.

Transition state search with a finer DFT grid:

.. code:: bash

   chemsmart sub pyscf -p test -f guess.xyz --defgrid defgrid3 ts
