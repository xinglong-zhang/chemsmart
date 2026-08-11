################
 xTB CLI Options
################

ChemSmart's bounded xTB 6.7.1 execution surface supports CPU ``sp``, ``opt``,
and ``hess`` jobs under both ``run`` and ``sub``.  The execution surface is
narrower than the xTB output parser and rejects unsupported native features.

*************************
 Basic Command Structure
*************************

.. code:: bash

   chemsmart run [RUN_OPTIONS] xtb -f GEOMETRY [XTB_OPTIONS] sp
   chemsmart run [RUN_OPTIONS] xtb -f GEOMETRY [XTB_OPTIONS] opt
   chemsmart run [RUN_OPTIONS] xtb -f GEOMETRY [XTB_OPTIONS] hess

A project is optional because GFN2 is a complete default method.  If a project
is provided, it may be a configured name or an explicit YAML path.

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
      -  Optional strict xTB project settings.
   -  -  ``-f, --filename``
      -  molecular artifact
      -  Required geometry source.
   -  -  ``-i, --index``
      -  1-based selection
      -  Select one or more structures from the source.
   -  -  ``-c, --charge``
      -  integer
      -  Molecular charge.
   -  -  ``-m, --multiplicity``
      -  positive integer
      -  Multiplicity; ChemSmart renders xTB ``--uhf`` as multiplicity minus one.
   -  -  ``-g, --gfn-version``
      -  gfn0/gfn1/gfn2/gfnff
      -  Maintained xTB method family.
   -  -  ``-sm, --solvent-model``
      -  alpb/gbsa
      -  Implicit-solvent model.
   -  -  ``-si, --solvent-id``
      -  validated identifier
      -  Solvent paired with the selected model and GFN version.
   -  -  ``--grad/--no-grad``
      -  boolean
      -  Gradient request, valid only for ``sp``.
   -  -  ``--optimization-level``
      -  xTB level
      -  ``opt``-only convergence level.

Solvent model and identifier must be supplied together.  The allowed solvent
set depends on the model and GFN version; unknown or incomplete pairs fail
before execution.

**********************
 Project YAML Contract
**********************

The only top-level sections are ``sp``, ``opt``, and ``hess``.  Charge and
multiplicity must be declared together in a section.  If both are omitted,
the selected source structure's electronic state is retained.

.. code:: yaml

   sp:
     jobtype: sp
     gfn_version: gfn2
     charge: 0
     multiplicity: 1
     solvent_model: null
     solvent_id: null
     grad: false
   opt:
     jobtype: opt
     gfn_version: gfn2
     optimization_level: vtight
     charge: 0
     multiplicity: 1
     solvent_model: null
     solvent_id: null
     grad: false
   hess:
     jobtype: hess
     gfn_version: gfn2
     charge: 0
     multiplicity: 1
     solvent_model: null
     solvent_id: null
     grad: false

Unknown keys, a contradictory ``jobtype``, optimization settings outside
``opt``, and an electron-count/multiplicity parity mismatch are rejected.

************************
 Preview and Validation
************************

Fake mode renders an isolated preview and records that no chemistry process
ran.  It cannot create a completion receipt.  Real execution requires xTB
6.7.1, CPU resources, an exact argv list, and a matching environment receipt.
Result validation checks the version, termination, method, job kind, charge,
multiplicity, geometry identity, energy, requested optimization or Hessian
artifacts, and source/project provenance.

GPU execution, arbitrary xcontrol text, molecular dynamics, path following,
unsupported constraints, and unregistered workflow families are not part of
this surface.  ChemSmart does not fall back to native xTB input text.
