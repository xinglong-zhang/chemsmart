#################################
 Transition State Search (PySCF)
#################################

This page covers transition state searches and intrinsic reaction coordinate (IRC) calculations using PySCF and the
geomeTRIC optimizer.

Both jobs read the ``gas`` section of the project YAML, like ``opt``, and both need ``opt_solver: geometric``. They run
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

*******************************
 Intrinsic Reaction Coordinate
*******************************

Follow the reaction path down from a transition state to the minima on either side.

.. code:: bash

   chemsmart sub [OPTIONS] pyscf [PYSCF_OPTIONS] irc [SUBCMD_OPTIONS]

IRC Options
===========

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-d, --direction``
      -  choice
      -  Branches to follow: ``both`` (default), ``forward`` or ``backward``. The project key ``irc_direction`` sets the
         same choice.

The job first computes the analytic Hessian of the starting structure and stops unless it has exactly one imaginary
mode, so a minimum or a higher-order saddle point is never followed as a transition state. Start from the result of a
``ts`` job. ``freq`` does not apply.

``forward`` and ``backward`` are the two signs of the imaginary mode, which the eigensolver chooses, so they do not say
which side is the reactant. Each branch is followed until it reaches a minimum or ``opt_maxsteps`` steps, and each
records whether it converged. A branch that fails does not discard the other.

The path is stored in the results file, from the end of the backward branch, through the starting structure, to the end
of the forward branch. The job's own energy, forces and geometry are those of the starting structure. See
:doc:`pyscf-results`.

Basic Usage
===========

Follow both branches from a transition state:

.. code:: bash

   chemsmart sub pyscf -p test -f guess_ts_gas_phase.out irc

Results are written to ``guess_ts_gas_phase_irc_gas_phase.h5``; pass ``-l`` to choose a shorter label.

Follow one branch:

.. code:: bash

   chemsmart sub pyscf -p test -f guess_ts_gas_phase.out irc -d forward

The label gains ``f`` for ``forward`` or ``r`` for ``backward``, as Gaussian IRC labels do, so the results go to
``guess_ts_gas_phase_ircf_gas_phase.h5`` and the two branches never share a results file.
