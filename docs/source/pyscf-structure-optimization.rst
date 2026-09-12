################################
 Structure Optimization (PySCF)
################################

This page covers geometry optimization and single point calculations using PySCF.

A PySCF project splits its settings into two sections. ``gas`` supplies the ``opt`` and ``hess`` jobs, and ``solv``
supplies ``sp``. The examples below use the packaged ``test`` project; substitute your own project name. For the
complete option reference, see :doc:`pyscf-cli-options`.

***********************
 Geometry Optimization
***********************

Find the minimum energy structure by adjusting atomic positions.

.. code:: bash

   chemsmart sub [OPTIONS] pyscf [PYSCF_OPTIONS] opt

The ``gas`` section of the project YAML supplies the method, basis set, and optimizer:

.. code:: yaml

   gas:
     functional: b3lyp
     basis: def2-svp
     freq: True
     density_fit: True
     defgrid: defgrid2
     scf_tol: 1.0e-9
     opt_solver: geometric
     opt_maxsteps: 100

Optimization Options
====================

Optimizer settings are given at the PySCF group level, before the ``opt`` subcommand:

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

Basic Usage
===========

Standard geometry optimization:

.. code:: bash

   chemsmart sub pyscf -p test -f molecule.xyz opt

Results are written to ``molecule_opt_gas_phase.h5``, alongside the PySCF log ``molecule_opt_gas_phase.out``.

Change the optimizer for a single run without editing the project:

.. code:: bash

   chemsmart sub pyscf -p test -f molecule.xyz --opt-solver berny opt

Raise the step ceiling for a shallow or awkward potential energy surface:

.. code:: bash

   chemsmart sub pyscf -p test -f molecule.xyz --opt-maxsteps 300 opt

With ``freq: True`` in the ``gas`` section, the optimization is followed by a Hessian in the same process, so one
command gives both the optimized geometry and its frequencies. See :doc:`pyscf-frequency-calculations`.

***************************
 Single Point Calculations
***************************

Compute energy and properties at a fixed geometry.

.. code:: bash

   chemsmart sub [OPTIONS] pyscf [PYSCF_OPTIONS] sp

Single points read the ``solv`` section, which is normally where the larger basis set and the solvent model live:

.. code:: yaml

   solv:
     functional: b3lyp
     basis: def2-tzvp
     freq: False
     density_fit: True
     defgrid: defgrid2
     scf_tol: 1.0e-9
     solvent_model: smd
     solvent_id: water

Basic Usage
===========

Standard single point calculation:

.. code:: bash

   chemsmart sub pyscf -p test -f molecule.xyz sp

Because the ``solv`` section requests SMD water, the results file is named ``molecule_sp_smd_water.h5``. A section with
no solvent produces ``_gas_phase`` instead.

Single point on one structure of a multi-structure file, using 1-based indexing:

.. code:: bash

   chemsmart sub pyscf -p test -f trajectory.xyz -i 3 sp

Single point on a structure fetched from PubChem, naming the job explicitly:

.. code:: bash

   chemsmart sub pyscf -p test -P 962 -l water sp

This fetches water (CID 962) and writes ``water_sp_smd_water.h5``.

*****************************************
 Changing the Method, Basis, and Solvent
*****************************************

Scientific settings normally live in the project YAML, where the reasoning behind them stays visible and reusable. The
command line overrides them for a single run. Every option below is given at the PySCF group level, before the
subcommand.

Method and Basis Set
====================

Use ``-x`` for a DFT functional or ``-A`` for a Hartree-Fock reference, never both. Basis sets take the hyphenated PySCF
spelling, for example ``def2-svp``, not Gaussian's ``def2svp``.

Optimize with a different functional and basis set:

.. code:: bash

   chemsmart sub pyscf -p test -f molecule.xyz -x m062x -b def2-tzvp opt

Run a Hartree-Fock single point instead:

.. code:: bash

   chemsmart sub pyscf -p test -f molecule.xyz -A hf -b def2-svp sp

.. warning::

   PySCF's ``b3lyp`` follows the Gaussian VWN3 definition, while ORCA's ``B3LYP`` is the VWN5 variant that PySCF calls
   ``b3lyp5``. The two differ by roughly 23 kcal/mol for water in def2-SVP, so energies from the two programs are not
   comparable unless the same variant was used in both.

``--defgrid`` selects the DFT integration grid, so it cannot be combined with ``-A``. An auxiliary basis only applies
when density fitting is on, so ``-ab`` cannot be combined with ``--no-density-fit``.

Solvent
=======

``-sm`` and ``-si`` are given together, so the applied dielectric is always explicit. The available models are ``pcm``,
``iefpcm``, ``cpcm``, ``cosmo``, ``ssvpe``, and ``smd``.

Add SMD water to a gas phase optimization:

.. code:: bash

   chemsmart sub pyscf -p test -f molecule.xyz -sm smd -si water opt

Drop the solvent that the project specifies:

.. code:: bash

   chemsmart sub pyscf -p test -f molecule.xyz --remove-solvent sp

*************************
 Charge and Multiplicity
*************************

By default the charge and multiplicity come from the molecular source. Use ``-c`` and ``-m`` to override them:

.. code:: bash

   chemsmart sub pyscf -p test -f molecule.xyz -c 1 -m 2 sp

Give the multiplicity as ``2S+1``. PySCF's own ``spin`` keyword is ``2S``, and chemsmart converts between the two, so a
multiplicity of 2 becomes ``spin=1`` in the generated script.
