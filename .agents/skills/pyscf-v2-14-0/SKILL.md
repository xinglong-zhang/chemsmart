---
name: pyscf-v2-14-0
description: Operate, debug, or improve the PySCF 2.14.0 CPU backend used by ChemSmart for SP, geometry optimisation, analytic Hessian and harmonic frequencies, structured HDF5 results, and bounded TD preview. Use for PySCF project settings, compute interpreters, state mapping, density fitting, solvents, dispersion, derivatives, or result validation.
---

# PySCF 2.14.0 for ChemSmart

## Mission

Act as a computational scientist integrating a Python-library backend through
the same ChemSmart CLI and project model as executable programs. Keep the
compute environment explicit and learn from real numerical behavior.

## Maintained surface

- CPU `sp`, `opt`, and `hess` execution;
- fixed-geometry closed-shell gas-phase singlet TDA/TDDFT preview only;
- stage-specific project sections with no implicit scientific inheritance;
- standalone generated Python calculation; and
- structured `LABEL.h5` as the numeric authority.

PySCF's `mol.spin` is `2S`, so ChemSmart converts multiplicity `2S+1` to
`multiplicity - 1`. Preserve atom order, charge, multiplicity and geometry
through every stage.

## Scientific settings

Use either `ab_initio: hf` or a DFT `functional`, never both. Use PySCF basis
spelling such as `def2-svp`. Configure density fitting before GPU conversion,
attach implicit solvent only when the solvent identifier is explicit, and
probe optimisation dependencies at call time rather than trusting imports.

An explicit `hess` analyses the supplied geometry. In a multi-stage workflow,
feed it the validated optimized geometry. Do not run PySCF's independent
thermochemistry calculator; route frequencies into ChemSmart's cross-program
analysis layer.

## CPU environment

Point `PYSCF.EXEFOLDER` at the `bin` directory whose Python owns PySCF 2.14.0,
NumPy, HDF5 and the selected optimiser. The compute script does not import the
ChemSmart package. Record convergence and applied settings in HDF5; process
exit zero is insufficient.

## Self-improvement cycle

When PySCF behavior surprises you, test reference state, basis spelling, grid,
SCF convergence, optimiser availability, derivative support, solvent wrapper,
compute interpreter and HDF5 parsing separately. Improve the smallest general
mapping or diagnostic and rerun a real molecule with a different scientific
role.

Do not silently map unsupported TS, IRC, scan, NEB, post-HF, double-hybrid or
mixed-basis work to a simpler calculation.
