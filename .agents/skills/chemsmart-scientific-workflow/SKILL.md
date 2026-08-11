---
name: chemsmart-scientific-workflow
description: Design, preview, execute with approval, analyse, or improve Gaussian, ORCA, xTB, PySCF, and GPU4PySCF computational-chemistry workflows through ChemSmart project YAML and CLI. Use for identity, state, methods, basis sets, solvents, constraints, multi-stage DAGs, convergence, frequencies, thermochemistry, and physical result interpretation.
---

# ChemSmart Scientific Workflow

## Mission

Act as a challenging computational scientist. Use ChemSmart as the canonical
interface while retaining freedom to choose any scientifically defensible
calculation and post-processing route allowed by the question.

## Define the problem

Before planning, establish:

- molecule or system identity and geometry role;
- coordinate units, atom order, charge, multiplicity and electronic state;
- constraints and relevant environment;
- requested observable, direction, units and physical conditions;
- method or program facts fixed by the question; and
- whether the work is planning, preview, existing-result analysis or real
  execution.

Ask rather than invent a consequential missing fact.

## Build and run the workflow

1. Inspect the live program capability and exact environment.
2. Encode method, basis/ECP, dispersion, solvent, convergence and stage settings
   in project YAML.
3. Represent necessary producer-consumer relationships in a scientific DAG.
4. Let ChemSmart compile and preview canonical CLI operations.
5. Obtain user authority before real execution or scheduler submission.
6. Validate normal termination, convergence, stationary-point character,
   identity, state and requested properties.
7. Perform unit-aware analysis and state signs and conventions explicitly.

Use the validated optimized geometry for a downstream Hessian, SP, IRC or path
job. Never let a filename or model narrative substitute for the geometry and
state handoff.

## Program principles

- **Gaussian and ORCA:** use project YAML and ChemSmart-generated input.
- **xTB:** use the maintained CPU `sp`, `opt`, and `hess` surface; do not add
  native xcontrol or unsupported paths as a fallback.
- **PySCF:** use stage-specific `sp`, `opt`, `hess`, and bounded preview-only
  `td` settings. Structured HDF5 is the numeric result.
- **GPU4PySCF:** treat GPU as a PySCF engine. Never silently fall back to CPU.

## Self-improvement cycle

Treat every unexpected result as a research opportunity. Check the chemistry,
reference wavefunction, derivative availability, convergence, numerical grid,
program version, environment, parser and analysis convention. Design a probe
that separates the explanations. Improve the smallest reusable ChemSmart layer
or this skill, then rerun a materially relevant calculation.

Do not turn one successful method or one paper workflow into a universal rule.
Record the conditions under which a lesson applies and remove guidance when
new evidence disproves it.

## Interpretation

Distinguish electronic, zero-point, enthalpy and free-energy quantities;
adiabatic and vertical geometries; reactant/product subtraction directions;
minimum and transition-state frequency criteria; and planned versus executed
work. A number without units and a physical convention is not a result.
