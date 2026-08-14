---
name: chemsmart-agent
description: Operate, audit, document, or improve the production ChemSmart computational-chemistry Agent through canonical project YAML, the live CLI, scientific DAGs, safe preview, visible human approval, deterministic execution, and typed result analysis. Use for ChemSmart Agent planning, TUI, approval, provider adapters, program support, or scientific workflow work.
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

- Gaussian CPU ``sp/opt/ts/irc/td/link`` has project-backed planning, native
  input preview, and typed analysis of supplied completed outputs; do not claim
  Agent execution in version 3.1.4.
- ORCA CPU planning covers ``sp/opt/ts/irc/td/neb``. Release-qualified
  execution covers single-points, optimization/frequency, transition-state,
  excited-state, and serial DAG workflows. Treat ``irc`` and ``neb`` as
  preview paths until the selected target is qualified.
- PySCF CPU ``sp/opt/hess`` and xTB CPU ``sp/opt/hess`` have approved real
  execution paths. PySCF CPU ``td`` is preview-only.
- GPU4PySCF ``sp/opt/hess`` is a PySCF-engine configuration and preview
  surface until a compatible GPU target is qualified. NCIPLOT and other human
  CLI families without an Agent declaration are not Agent execution paths.
- Keep product capability, observed scientific evidence, and current-host
  readiness distinct. A supported Gaussian path does not imply a licensed
  executable is present; a supported GPU path does not imply a compatible GPU
  stack. The live environment probe and human review decide whether the exact
  operation can run here.
- Runtime semantics are provider-neutral. Registered version-3.1.4 adapters
  are Alibaba Token Plan and DeepSeek, configured entirely by a user profile.
  There is no default model; the profile must explicitly state the selected
  model and its context/output limits.

## Use visible one-shot approval

Planning, YAML work, CLI compilation, safe preview, and result analysis are
non-executing. Before an engine launch:

1. produce the complete project-backed DAG;
2. compile and safely preview every executable node;
3. present molecular/state identity, effective YAML, ChemSmart CLI operations,
   data handoffs, environment, and resources in the terminal interface;
4. let the human enter ``/approve`` once, or choose ``/deny`` or ``/revise``;
   and
5. hand the displayed workflow to the provider-free deterministic executor.

Never create an approval on the model's behalf. Never offer permanent,
session-wide, prefix-based, or "always allow" chemistry execution. A revised
scientific input, project, environment, resource request, or DAG requires a new
review. Internal receipts remain provenance; never make a human retype their
digests as a second scientific authority. A reviewed multi-node DAG needs one
approval, not one prompt per node.

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
