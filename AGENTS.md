# ChemSmart Agent Operating Contract

## Mission

Develop ChemSmart as a CLI-first, provider-neutral computational-chemistry
automation agent. A model may plan, ask, explain, and propose typed
ScientificTaskSpec and CommandWorkflowSpec objects. It must not author,
patch, or treat as authoritative Gaussian, ORCA, or xTB native input text, or
the standalone Python scripts emitted for PySCF. Those files are generated
artifacts, never a model-editable interface.
Deterministic ChemSmart code owns CLI semantics, command compilation,
permission policy, execution, scientific validation, and evidence recording.

Do not resume GUI, desktop-app, packaging, Studio, or visual-design work unless
the user explicitly reopens it. Do not treat an installed tool, a valid command,
or a plausible model answer as evidence that a calculation or scientific claim
is correct.

## Authority and scope

Use this precedence order:

1. Explicit user request and approval for the current task.
2. This file and the nearest applicable repository instructions.
3. The checked-out CLI parser, configuration, job, runtime, and test contracts.
4. Generated artifacts, deterministic validators, execution receipts, and
   primary sources.
5. Model output and legacy guidance.

`CLAUDE.md`, `.clinerules/`, and `memory-bank/` retain useful historical or
tool-specific material, but they must not weaken this contract. Treat dated
metrics and collection instructions there as historical unless the user
explicitly reactivates that work.

## Before changing anything

- Inspect the branch, ancestry, dirty state, relevant instructions, source
  contracts, and focused tests.
- Preserve unrelated or partial user work. Never reset, clean, stash, or
  overwrite it without explicit authority.
- State assumptions that materially affect chemistry, execution, cost, or
  scientific interpretation.
- Derive CLI behavior from the current Click parser and generated schema; join
  that observation with the immutable program-capability declaration and hash
  both. Neither source alone establishes installation or execution readiness.
- Treat the command compiler as the only authority that turns a model proposal
  into argv. It must resolve live schema options, trusted project and artifact
  references, canonical long flags, and shell-safe rendering. A model never
  supplies executable shell syntax, arbitrary paths, option ordering, aliases,
  quoting, or a native-engine fallback.
- Do not install dependencies, alter environment pins, contact external
  systems, commit, push, or publish without authority for that action.

## Scientific workflow

Make these facts explicit before a calculation is treated as specified:

- molecule identity and stable artifact identifier;
- exact geometry frame and coordinate units;
- charge, multiplicity, electronic-state assumptions, and constraints;
- requested observable, program, job kind, method, basis/ECP, dispersion,
  solvent, temperature/standard-state convention, and resource target;
- required evidence, diagnostics, and limitations.

Ask instead of inventing a scientifically consequential missing fact. Never
infer geometry identity from a filename alone. Compile the typed intent through
the live schema, trusted project/artifact resolver, safe CLI preview, and
independent parser observation before an action can be called previewed.
Generated native inputs and PySCF scripts are downstream evidence produced by
ChemSmart, not model-editable interfaces.

Use these mutually exclusive outcome labels precisely:

- **planned**: intent is recorded, but no executable input exists;
- **previewed**: input or command was rendered but not run;
- **executed**: an engine or scheduler was invoked and a receipt exists;
- **validated**: required deterministic checks passed;
- **reproduced**: an independently rerun, pinned environment produced the
  declared result within the stated tolerance;
- **waiting for approval**, **blocked**, or **failed**: a required condition
  remains unmet.

Do not call a run complete when a required receipt, validator, artifact, or
approval is absent.

## Approvals and execution

Planning, read-only inspection, deterministic validation, and fixture-based
simulation may proceed within the user's scope. Require explicit approval for:

- real local calculations, scheduler/HPC submission, cancellation, or retry;
- writes outside a disposable task workspace or overwrites of user artifacts;
- paid model or compute use, networked external execution, and publication;
- a material change to the agreed molecular model, method, electronic state,
  resource budget, or scientific claim.

Bind approval to the exact command, inputs, project, executable/environment,
and artifact hashes that already exist. For a multi-stage workflow whose
consumer input does not yet exist, approval may instead bind the exact
producer plus a deterministic output-selection and consumer-materialization
rule. Invalidate it when any scientific setting, resource, producer, selection
rule, validator, or command family changes. Keep secrets out of prompts, logs,
commits, and evidence bundles.

Use a paid DeepSeek model call or literature lookup only through an existing
user-owned quota and a short-lived standard Keychain lease. Never print,
persist, transmit in a prompt, or infer a secret. Record only the provider,
endpoint class, key-validation outcome, quota-sufficiency outcome, and
non-secret error class. Do not top up a quota, change a billing plan, or turn a
provider credential into general network authority. Once a user has authorized
the current development phase, lease-bound calls within its recorded quota may
proceed without per-call reapproval; a new provider, target, quota expansion,
or billing change needs new authority.

## Agent architecture

- Keep provider-private reasoning transport-only. It may be replayed
  ephemerally when a provider protocol requires continuation, but it must not
  be persisted, displayed, graded, or used as scientific evidence. Persist
  visible responses, typed actions, tool calls and results, concise public
  decisions, artifacts, usage, validator findings, and terminal outcomes.
- Treat Gaussian, ORCA, xTB, PySCF, and NCIPLOT as programs with distinct
  declared capabilities. Treat CPU and GPU4PySCF as engines of PySCF. A
  registry entry is not proof that a dependency, executable, GPU stack, or
  scientific method is available.
- Let models propose program candidates only. Deterministic host code owns
  environment observation, requested-versus-selected bindings, substitution
  equivalence, approval requirements, preflight, readiness, and terminal
  state. Never silently substitute one program or engine for another, and
  never fall back from GPU to CPU.
- Gaussian-to-PySCF substitution is eligible only for explicitly registered
  single-point, optimization, fixed-geometry Hessian, or optimization-plus-
  frequencies mappings. It always requires hash-bound user approval. TS, IRC,
  TD, scans, QMMM/ONIOM, NEB, post-HF, double-hybrid, mixed basis/ECP, and
  unsupported constraints must not fall back automatically. Preserve those
  steps and scientifically relevant alternatives in the planning draft as
  explicit capability gaps; block only their materialization or execution
  until a supported path is selected.
- Expose the frontier calculation-preparation surface only as typed project
  operations plus synthesize, repair, inspect, and explain command workflow
  operations. Legacy molecule/settings/job/input/execution builders may remain
  in an explicit compatibility profile, but must be absent and fail closed from
  the command-compiled frontier profile.
- Treat raw legacy direct-string synthesis and compact-v8 conversion as
  baseline or migration inputs only. They are not Frontier Runtime V2
  model-surface authorities and may not bypass typed compilation, preview, or
  evidence gates.
- Let a CommandWorkflowSpec bind a workflow ID, task-spec ID, live CLI-schema
  digest, and ordered immutable command nodes. Nodes contain only trusted
  artifact IDs/hashes, project references, declared intent, dependencies,
  constraint IDs, and expected artifact classes. The compiler performs DAG checking, schema
  resolution, canonical argv rendering, safe preview, parser observation, and
  intent round-trip comparison. A structured counterexample may support at
  most two constrained repairs; it must not silently change an explicit
  program, geometry, charge, multiplicity, method, or constraint.
- Keep exploratory and draft planning available before all execution evidence
  exists. Missing hashes, installed-program evidence, or project artifacts are
  explicit unresolved fields at the planning level; they block only the
  affected materialization, preview, approval, or execution transition. Do not
  discard an otherwise useful computational-chemistry plan merely because it
  is not yet executable.
- Bind every repair to the immediately preceding ScientificTaskSpec and
  preflight-receipt digests. In the active command profile, a terminal success
  requires a deterministic `previewed` receipt; a model assertion, command
  string, or proposed repair is never a completion substitute.
- Give each task the smallest relevant tool surface and explicit token, tool,
  wall-time, and compute budgets.
- Use subagents only for bounded, independently verifiable work with declared
  immutable inputs, expected outputs, allowed tools, owner, and merge rule.
  One agent owns each mutable artifact.
- Use a critic as a fresh, read-only cross-examiner. A critic cannot approve,
  execute, or repair its own finding. Deterministic checks or independent
  computation arbitrate disagreements.
- End every run as planned, complete, failed, blocked, or waiting for approval;
  do not loop indefinitely. `planned` is a useful non-executable result and
  must not be mislabeled as blocked merely because action-grade evidence is
  intentionally unresolved.
- When the current CLI or scientific validator cannot express a requested
  task, preserve the requested step in the scientific plan and localize the
  gap. Use `needs_clarification` or `infeasible` only for the affected
  materialization or action path, with a structured reason; independent plan
  nodes may continue. Neither state permits a native-input fallback.

## Evidence and reporting

Record stable IDs, input and output hashes, engine and environment versions,
commands, working directory, timestamps, exit status, parsed values with
units, validator outputs, approval records, and claim-to-evidence links.

Separate observation, computed result, inference, literature statement, and
unresolved uncertainty. A report, notebook, or chat summary is a rendered view
of evidence, not the evidence source. Use QCSchema-compatible records where
practical, retain native engine artifacts, and make each numerical claim
traceable to a receipt.

## Project-local skills

Use the smallest matching skill set:

- `chemsmart-agent-harness` for provider adapters, tool loops, permissions,
  Runtime V2, task graphs, and harness evaluation;
- `chemsmart-scientific-workflow` for Gaussian, ORCA, xTB, and PySCF task
  intake, preflight, approved execution, and physical validation;
- `chemsmart-evidence-audit` for provenance, claims, citations, reports,
  red-teaming, and evaluation.

## Validation and reporting discipline

- During a milestone, prefer source inspection, schema/receipt checks, and
  narrowly scoped deterministic probes. Do not repeatedly run pytest, Ruff, or
  broad checks after each edit. Run one focused suite when a material milestone
  is complete; allow at most one evidence-driven rerun for that milestone.
- Run full agent tests, read-only Ruff, schema/link/citation/secret checks, and
  diff checks only at the preregistered integration/freeze gate. Do not
  autofix, format, or regenerate snapshots unless separately authorized.
- Keep product, runtime, scientific, and release readiness separate. A focused
  green check is not proof of product or scientific readiness.
- Report green checks, blockers, retired metrics, and unverified claims
  separately.
- Preserve backward replay of existing runtime events when evolving the agent;
  extend the current Runtime V2 nucleus instead of introducing a competing
  runtime.
