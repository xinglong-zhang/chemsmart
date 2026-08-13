---
name: chemsmart-agent
description: Operate, audit, document, or improve the production ChemSmart computational-chemistry Agent through canonical project YAML, the live CLI, exact scientific DAGs, safe preview, one-shot human approval, deterministic execution, and typed result analysis. Use for ChemSmart Agent planning, TUI, approval, provider adapters, program support, or scientific workflow work.
---

# ChemSmart Agent

Treat ChemSmart as the execution authority and the model as the scientific
reasoner. Read `AGENTS.md` before changing the product boundary.

## Start from the live product

1. Inspect the project YAML loader and the exact Click command involved.
2. Establish molecular identity, geometry role, units, charge, multiplicity,
   electronic state, constraints, requested observable, and conditions.
3. Ask for a consequential missing fact instead of guessing it.
4. Use ChemSmart to materialise native input; never author a substitute native
   input or shell execution path in the model layer.

## Keep states honest

Name work as proposed, planned, materialised, previewed, approved, executing,
engine-complete, parsed, scientifically validated, or interpreted. Do not
collapse these states or treat provider text, a fake preview, or a fixture as
engine evidence.

## Use the complete scientific layer

The production surface is broader than command generation. Use live
capability/environment inspection, identity and project-YAML binding, causal
DAG/frontier operations, validated geometry handoff, registered-result
inspection, semantic quantity extraction, RRHO or parameterised qRRHO,
unit-aware expression DAGs, evidence-bound claims, and scientific decisions
when the task calls for them. The supported expression vocabulary includes
CBS extrapolation, Boltzmann populations and averages, harmonic ZPE,
imaginary-mode counts, geometry measurements, centres of mass, inertia,
rotational constants, and connectivity changes.

Do not force every task through every layer. A valid route may be analysis-only
or may choose a different causal decomposition. Use only the operations needed
to answer the scientific request, and keep every source quantity and convention
visible.

## Apply the production support boundary

- Gaussian CPU ``sp/opt/ts/irc/td/link``, ORCA CPU
  ``sp/opt/ts/irc/td/neb``, PySCF CPU ``sp/opt/hess``, GPU4PySCF
  ``sp/opt/hess``, and xTB CPU ``sp/opt/hess`` have project-backed planning,
  safe preview, exact human approval, normal runner execution, and typed
  analysis for quantities actually produced.
- PySCF CPU ``td`` is preview-only. NCIPLOT and other human CLI families with
  no Agent engine/job declaration are not Agent execution paths.
- Keep product capability, observed scientific evidence, and current-host
  readiness distinct. A supported Gaussian path does not imply a licensed
  executable is present; a supported GPU path does not imply a compatible GPU
  stack. The live environment probe and human review decide whether the exact
  operation can run here.
- Runtime semantics are provider-neutral. Registered version-3.1.4 adapters
  are Alibaba Token Plan and DeepSeek, configured entirely by a user profile.
  There is no default model; the profile must explicitly state the selected
  model and its context/output limits.

## Use exact one-shot approval

Planning, YAML work, CLI compilation, safe preview, and result analysis are
non-executing. Before an engine launch:

1. produce the complete project-backed DAG;
2. compile and safely preview every executable node;
3. present molecular/state identity, YAML, CLI argv, data handoffs,
   environment, resources, and exact digests;
4. let the human choose approve once, deny, revise, or quit; and
5. hand an approved packet to the provider-free deterministic executor.

Never create an approval on the model's behalf. Never offer permanent,
session-wide, prefix-based, or "always allow" chemistry execution. Any change
to scientific inputs, project settings, command, environment, resources, or
DAG requires a new review. A reviewed multi-node DAG needs one approval, not
one prompt per node.

## Improve at the owning layer

Classify a failure as scientific reasoning, program limitation, environment,
parser, or missing ChemSmart affordance. Fix the smallest general project
setting, CLI compiler, adapter, runner, parser, typed analysis operation, or
error message. Do not learn a molecule, paper answer, DOI, private DAG, tool
order, or provider-specific workaround.

Use one focused mechanical check after a change. Then prefer one decisive
scientific observation through the real public surface. Report exactly which
program ran, which artifact was parsed, and which scientific validation was
or was not established.

## Review as a computational scientist

Check identity and atom order, electronic state, method semantics, geometry
handoff, convergence, stationary-point evidence, signs, units, physical
conditions, and causal dependencies. Accept any chemically and
mathematically sound route. The final interpretation and publication decision
belong to the human scientist.
