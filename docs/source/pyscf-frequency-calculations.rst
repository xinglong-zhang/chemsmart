################################
 Frequency Calculations (PySCF)
################################

This page covers Hessian and harmonic frequency calculations using PySCF.

Frequencies are controlled by ``freq`` in the ``gas`` section of the project YAML. They can be taken as part of an
optimization, or as a separate job on a geometry you already have.

**********************************
 Frequencies with an Optimization
**********************************

With ``freq: True`` in the ``gas`` section, an ``opt`` job runs the Hessian immediately after the geometry converges, in
the same process:

.. code:: yaml

   gas:
     functional: b3lyp
     basis: def2-svp
     freq: True
     density_fit: True
     defgrid: defgrid2
     scf_tol: 1.0e-9
     opt_solver: geometric

One command then gives both the optimized structure and its frequencies:

.. code:: bash

   chemsmart sub pyscf -p test -f molecule.xyz opt

This is the cheaper route, since the converged wavefunction is still in memory when the Hessian starts. Turn ``freq``
off when you only want the geometry, for example while screening conformers.

********************
 Standalone Hessian
********************

Compute the Hessian and harmonic frequencies on a geometry that is already optimized.

.. code:: bash

   chemsmart sub [OPTIONS] pyscf [PYSCF_OPTIONS] hess

The ``hess`` job reads the same ``gas`` section as ``opt``.

Basic Usage
===========

Frequency calculation on a structure you already trust:

.. code:: bash

   chemsmart sub pyscf -p test -f optimized.xyz hess

Frequencies at a tighter grid than the project specifies:

.. code:: bash

   chemsmart sub pyscf -p test -f optimized.xyz --defgrid defgrid3 hess

Frequencies are reported in cm⁻¹, with imaginary modes kept as negative values, so a minimum has none and a transition
state has exactly one. The results file also stores the Hessian, the normal modes, the reduced masses, and the force
constants.

.. warning::

   ``hess`` uses the geometry it is given and does not optimize it first. A Hessian evaluated away from a stationary
   point is not a meaningful vibrational analysis. Run ``hess`` on a converged structure.

*****************************
 Using an Optimized Geometry
*****************************

To take frequencies at a different level of theory from the optimization, point ``hess`` at the completed optimization.

Optimize first:

.. code:: bash

   chemsmart sub pyscf -p test -f molecule.xyz opt

Then run the frequency job on that output:

.. code:: bash

   chemsmart sub pyscf -p test -f molecule_opt_gas_phase.out -x m062x hess

Give the ``.out`` file, not the ``.h5``. chemsmart recognizes a PySCF log and reads the geometry from the structured
results file stored beside it, so both files need to stay in the same directory. See :doc:`pyscf-results`.

The ``thermochemistry`` command does not yet accept PySCF output; enthalpies, entropies, and Gibbs free energies
currently come from Gaussian, ORCA, or xTB results. See :doc:`thermochemistry-analysis`.
