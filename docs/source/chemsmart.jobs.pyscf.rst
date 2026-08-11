chemsmart.jobs.pyscf package
############################

The PySCF backend materialises standalone Python calculations and uses a
structured HDF5 file as the numerical result authority.

Core modules
============

``chemsmart.jobs.pyscf.settings``
   Stage-specific scientific settings and validation.

``chemsmart.jobs.pyscf.writer``
   Standalone calculation generation for the selected compute interpreter.

``chemsmart.jobs.pyscf.runner``
   Preview and approved CPU or GPU-engine execution.

``chemsmart.jobs.pyscf.environment``
   Interpreter and compute-environment qualification.

``chemsmart.jobs.pyscf.validation``
   Structured result checks for molecular identity, state and quantities.

The maintained CPU execution surface is single-point energy, geometry
optimisation, Hessian and harmonic frequencies. TD/TDA remains preview-only.
See :doc:`pyscf-cli-options` for the user-facing command and YAML contract.
