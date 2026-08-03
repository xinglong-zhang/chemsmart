# ChemSmart 3.1.4 program and agent architecture

Status: implementation specification. This document describes the integrated
design and its required evidence. It does not claim that the final integration
or real-calculation gates have passed.

For implementation order and current evidence state, see
[`chemsmart-v314-integration-implementation-plan.md`](chemsmart-v314-integration-implementation-plan.md)
and
[`../maintenance/integration-status-v314.md`](../maintenance/integration-status-v314.md).

## 1. Objective and scope

This fork keeps the upstream package version at `3.1.4` and identifies the
integrated build by its Git revision and joined capability/schema digest. It
adds three bounded surfaces to the exact upstream v3.1.4 base:

1. PySCF 2.14.0 for `sp`, `opt`, and `hess`, with CPU and GPU4PySCF 1.8.0
   engine declarations.
2. xTB 6.7.1 for CPU `sp`, `opt`, and `hess` only.
3. A capability-driven Runtime V2 agent that converts typed scientific intent
   into project settings, canonical ChemSmart commands, safe previews, and
   deterministic receipts.

The architecture does not turn model text into native input, generated Python,
shell, readiness, approval, or terminal state. Gaussian, ORCA, xTB, and PySCF
artifacts remain products of checked-in deterministic code. Declaring a program
or command is not evidence that its executable or scientific method is usable.

## 2. Current source-of-truth map

| Concern | Current authoritative implementation |
|---|---|
| Program declarations | `chemsmart/settings/capabilities.py` |
| Live CLI observation | `chemsmart/agent/cli_schema.py` |
| Agent capability projection and environment receipts | `chemsmart/agent/capabilities.py` |
| Program substitution | `chemsmart/agent/knowledge.py` |
| Project rendering and loader validation | `chemsmart/agent/projects.py` |
| Semantic command-workflow drafts and materialized DAG contracts | `chemsmart/agent/workflows.py` |
| Typed command compilation and parser inspection | `chemsmart/agent/commands.py` |
| Isolated fake/test preview | `chemsmart/agent/preview.py` |
| Cross-receipt readiness | `chemsmart/agent/preflight.py` |
| Model-visible tool schemas and host implementation | `chemsmart/agent/tool_specs.py`, `chemsmart/agent/tool_runtime.py` |
| Hash-chained event envelope and replay | `chemsmart/agent/runtime/events.py`, `chemsmart/agent/runtime/reducer.py` |
| PySCF settings, generation, execution, and verification | `chemsmart/jobs/pyscf/`, `chemsmart/settings/pyscf.py`, `chemsmart/io/pyscf/output.py` |
| xTB settings, execution, and verification | `chemsmart/jobs/xtb/`, `chemsmart/settings/xtb.py` |
| Advisory program knowledge | `chemsmart/agent/knowledge_packs/` |
| Boundary-campaign public context | `chemsmart/agent/experiments/program_management_context.py` |

The immutable `PROGRAM_CAPABILITIES` registry declares Gaussian, ORCA, xTB,
PySCF, and NCIPLOT. `load_program_capabilities()` projects that registry into
`ProgramCapabilityRegistryV1`; it does not maintain an independent program
inventory. The public `AgentProgramSupportOverlayV1` name is an alias of the
implemented `ProgramSupportOverlayV1`; it may narrow the projection, but cannot
add a program, job, engine, or project capability.

## 3. Lifecycle and authority

The canonical lifecycle is:

`scientific request -> semantic workflow draft with explicit unknowns -> evidence and identity grounding -> capability observation -> optional substitution assessment -> environment observation -> project materialization and loader validation -> command compilation -> parser inspection -> isolated safe preview -> cross-receipt preflight -> approval -> execution -> result verification -> terminal event`

The first transition is intentionally productive rather than hash-gated. A
model may organize stages, compare scientifically plausible routes, identify
required evidence, and expose capability gaps before exact coordinates,
projects, environments, or artifact bytes exist. It records missing facts in
typed `unresolved_fields`; it does not invent them or claim that a draft is
previewable.

Assurance increases progressively:

| Level | Purpose | Required identity |
|---|---|---|
| P0 semantic draft | scientific decomposition, alternatives, dependencies, evidence needs | semantic IDs, source locators when known, explicit unknowns |
| P1 materialization | exact project and artifact grounding, current CLI compilation, safe preview | exact bytes and hashes for every materialized node |
| P2 action | approval, execution, result verification, scientific claims | immutable command, project, input, environment, approval, and result receipts |

Failure to satisfy P1 never justifies a false P2 claim, but it also does not
make a useful P0 research plan disappear. The primary generality protocol is
defined in
[`../evaluation/progressive-planning-generality-protocol.md`](../evaluation/progressive-planning-generality-protocol.md).

The lifecycle has three related but non-interchangeable state views:

- **workflow-node state**: `planned`, `compiled`, or `previewed`;
- **node-preflight plan state**: `blocked`, `ready_for_safe_preview`, or
  `previewed`;
- **execution state**: not requested, waiting for approval, executed, failed,
  or validated.

`compiled` means only that canonical argv exists for fully resolved inputs. A
future producer edge keeps its consumer `planned`, so that consumer has no
invocation and no node-preflight state yet. `ready_for_safe_preview` is a
preflight finding that preview evidence is still required; it is not execution
readiness.

No model message advances a materialization, preflight, or execution state. The
host derives those states from receipts; the model may only propose the P0
draft and explicit unknowns.

Runtime may terminate a planning-only task as `planned` when a typed workflow
draft was recorded. `complete` remains reserved for the active profile's
required green action-grade gates; `planned` never implies preview or execution
readiness.

This is a progressive-assurance lifecycle, not a requirement that every
research idea arrive fully materialized. A model may create a broad
`planned` workflow with explicit unknowns, alternative program candidates,
unresolved project roles, and semantic artifact slots. Missing hashes or
environment receipts block only the transition to compilation, preview,
approval, or execution. They do not block literature interpretation,
scientific decomposition, project drafting, or an honest research plan.

| Decision | Model | Deterministic host | User approval |
|---|---:|---:|---:|
| Interpret a request and propose a typed candidate | propose | validate shape and ancestry | not normally |
| Bind geometry, charge, multiplicity, and artifact identity | no | yes | when scientific identity changes |
| Select paths, flags, aliases, ordering, or quoting | no | compiler only | no |
| Render project YAML candidate | propose fields | canonical renderer | write requires policy decision |
| Accept project settings | no | checked-out loader and validators | material scientific change requires approval |
| Establish installed program or engine | no | exact environment probe | no |
| Substitute one program for another | propose | closed matrix assessment | always for non-identity mapping |
| Establish preview or execution readiness | no | receipt gate | execution still needs approval |
| Run, retry, or fall back | no | approved dispatcher only | exact node approval or an explicitly producer-bound workflow rule |
| Declare terminal success | no | required-green event gate | no |

## 4. Capability and environment planes

`ProgramCapabilityV1` declares project requirements, job types, project-owned
parameters, and engines. The registry currently exposes:

| Program | CLI leaf jobs | Project policy | Declared engines |
|---|---|---|---|
| Gaussian | existing v3.1.4 Gaussian leaves | required | CPU |
| ORCA | existing v3.1.4 ORCA leaves | required | CPU |
| PySCF | `sp`, `opt`, `hess` | required | CPU, GPU |
| xTB | `sp`, `opt`, `hess` | optional, strict when supplied | CPU |
| NCIPLOT | top-level leaf | unsupported | CPU |

The following grades must remain separate:

1. **declared**: present in `PROGRAM_CAPABILITIES`;
2. **schema-bound**: observed in the checked-out Click tree;
3. **supported by the agent overlay**: component evidence exists for the exact
   job and engine;
4. **environment-observed**: the exact executable or compute interpreter and
   required dependencies were probed;
5. **previewed**: exact compiler-owned argv passed fake/test execution and its
   program validator;
6. **executed**: a process receipt exists;
7. **validated**: all required scientific and provenance checks passed.

`ProgramComponentConformanceReceiptV1` binds compiler, preview, preflight,
verifier, and optional execution observations to explicit
`covered_jobtypes`/`covered_engines`.
`build_command_compiled_preview_overlay()` must not promote unobserved jobs or
engines. An absent execution receipt can produce at most preview-only support.

`ProgramEnvironmentQueryV1` refers to a host-owned `EnvironmentTargetV1`.
`EnvironmentCapabilityReceiptV1` reports observations, not model assertions.
For PySCF, `TrustedComputeEnvironmentReceiptV1` binds the exact interpreter,
dependency versions, solver evidence, and GPU evidence. A controller
environment variable is not compute-interpreter evidence.

## 5. Project and CLI compilation contract

### 5.1 Project settings

`ProjectDocumentV1` and `ProjectRenderReceiptV1` create a deterministic YAML
candidate. Only `validate_project_yaml()` and the program's checked-out loader
can produce `ProjectValidationReceiptV1(status="valid")`.

PySCF projects use exactly three complete top-level sections: `sp`, `opt`, and
`hess`. Every section is independently materialized through
`PySCFJobSettings`; missing or unknown sections and fields fail. Method, basis,
solvent, density fitting, grid, convergence, optimizer, engine, and frequency
intent are stage-owned. There is no scientifically complete implicit PySCF
project.

xTB projects use the same three section names. Project omission intentionally
selects ChemSmart's explicit GFN2 defaults. A supplied but unresolved project
does not fall back to those defaults. An explicit electronic state must provide
charge and multiplicity together. Solvent model and solvent identifier must
also be a complete valid pair.

### 5.2 Commands

The maintained integrated leaves are:

- `chemsmart run pyscf sp|opt|hess`
- `chemsmart sub pyscf sp|opt|hess`
- `chemsmart run xtb sp|opt|hess`
- `chemsmart sub xtb sp|opt|hess`

`CommandProposalV1` contains only node identity, `run|sub`, program, job type,
trusted project/input artifact identifiers, scientific-identity digest, charge,
and multiplicity. `compile_command()` resolves the live Click option scope and
produces `CanonicalCommandInvocationV1`. The compiled receipt binds argv,
schema/capability digests, project/input hashes, identity, engine binding, and
repair ancestry. Its `compiled` status is not a preview.

`CommandInspectionReceiptV1` independently observes the Click parse.
`execute_safe_preview()` then invokes only exact compiler-owned argv in an
isolated directory. `run` requires `--fake --no-scratch`; `sub` additionally
requires `--test`. A successful `SafePreviewReceiptV1` binds emitted artifacts,
output, exit status, input/project hashes, and a program-specific validation
receipt. External artifacts are re-hashed before and after preview.

Repair uses a `CommandCounterexampleV1` bound to the prior task,
invocation, and preflight receipts. Repair is field-local and cannot alter
program, geometry, charge, multiplicity, project, or engine authority.

### 5.3 Multi-command planning and future artifacts

`CommandWorkflowDraftV1` is the broad planning contract for an ordered
calculation chain. It records semantic programs and job types, project roles,
dependencies, input slots, expected artifact classes, and unresolved fields
without requiring execution-grade hashes. Unsupported or not-yet-integrated
steps remain visible as findings instead of being deleted or silently
simplified.

`CommandWorkflowSpecV1` is the later materialized command contract.
`CommandNodeV1` records dependencies, immutable scientific and program
bindings, project identity, input bindings, and expected artifact classes.
`ArtifactBindingV1` distinguishes two cases:

- an external artifact is resolved immediately and therefore carries its exact
  SHA-256; and
- a producer output is a typed DAG edge. At draft time only its producer node,
  output ID, and artifact class are known. A host may later allocate a stable
  artifact ID, but its byte hash and producer receipt remain empty until the
  output is materialized.

An unresolved producer edge is valid only while the consumer remains
`planned`. It cannot be compiled, previewed, preflighted, approved, or
executed. After the producer completes, deterministic host code must resolve
the edge with the observed artifact hash and producer receipt, re-materialize
the consumer's scientific-identity binding, and recompute every downstream
digest. A model cannot predict, invent, or copy a future artifact hash.

The frontier tool surface therefore provides a typed workflow-planning
operation. It accepts semantic node intent, direct dependencies, optional
external artifact IDs, explicit unknown fields, and producer-output edges; the
host constructs `CommandWorkflowDraftV1`. It returns the canonical draft, the
currently groundable nodes, unresolved nodes, and non-terminal findings.
Materialization into `CommandWorkflowSpecV1` is a separate deterministic
transition. Existing
`synthesize_command` and `preview_command` remain single-node operations and
may be called only after the selected node has been grounded.

For `OPT -> HESS`, OPT may be compiled and previewed against the initial
geometry. HESS must remain planned until a real optimized artifact and receipt
exist. A second preview of HESS before that execution would be a false-ready
defect, not evidence of successful planning.

The empty digest fields on an unresolved producer edge are not a placeholder
for executable input. They declare that no executable artifact binding exists
yet. The resolved replacement is a new immutable binding with the actual byte
hash and producer-receipt hash; it changes every downstream invocation,
preflight, and approval digest.

## 6. Gaussian-to-PySCF substitution

Program substitution is an explicit scientific change, not an availability
shortcut. `ProgramCandidateProposalV1` is untrusted. The host constructs a
`ProgramSubstitutionRequestV1`, applies
`assess_typed_program_substitution()`, and records
`ProgramSubstitutionReceiptV1(readiness_authority=False)`.

| Requested Gaussian family | PySCF mapping | Additional conditions | Result without exact approval |
|---|---|---|---|
| single point | `pyscf sp` | HF or DFT; uniform basis; source claim | rejected |
| minimum optimization | `pyscf opt` | HF or DFT; uniform basis; no constraints; source claim | rejected |
| fixed-geometry Hessian/frequencies | `pyscf hess` | HF or DFT; uniform basis; source claim | rejected |
| optimization plus frequencies | ordered `pyscf opt` then `pyscf hess` | optimized HDF5 geometry hash must bind the Hessian input | rejected |

All eligible DFT transfers also require a verified functional-semantics
receipt. This is essential for names such as B3LYP whose program conventions
may differ. The following requests block: transition state, IRC, excited state,
scan, QMMM/ONIOM, NEB, post-HF, double hybrid, mixed basis, ECP, mixed
basis/ECP, and any unsupported constraint. Method-name similarity, missing
Gaussian installation, or a green PySCF environment probe does not prove
equivalence.

Identity mapping does not require substitution approval. Every non-identity
mapping requires `ProgramSubstitutionApprovalV1` bound to the exact request
digest. Project, environment, preview, and result gates still apply after
approval.

### 6.1 No-fallback invariant

A red environment, project, compiler, preview, preflight, or result receipt
terminates the current candidate. It does not select another engine, program,
method, basis, solvent, or electronic state. In particular:

- unavailable Gaussian may prompt a typed PySCF proposal, but the proposal is
  inert until the closed matrix and exact approval pass;
- an unavailable or incompatible GPU blocks the GPU candidate and never creates
  a CPU candidate;
- unsupported PySCF or xTB science never falls back to model-written native
  input;
- an execution failure is preserved as failure and does not authorize retry;
- changing any approved binding invalidates the approval rather than repairing
  the action in place.

## 7. PySCF and GPU4PySCF execution contract

ChemSmart owns the standalone Python script generated by
`chemsmart/jobs/pyscf/writer.py`. The script imports only the standard library,
NumPy, h5py, and PySCF-related packages; the model cannot author or edit it.

Before launch, `PySCFJobSettings.validate()` and scientific `preflight()` must
pass. `probe_compute_environment()` executes in the exact selected interpreter
and records interpreter hash, PySCF 2.14.0, NumPy, h5py, method/basis evidence,
and required solver facts. The input receipt binds script, source geometry,
source artifact, project, requested settings, and environment hashes. A child
process failure propagates as failure even when a scheduler wrapper itself ran.

The authoritative result is structured HDF5, not the text log. Schema 2.0
contains real `spec`, `provenance`, `status`, and `results` groups.
`verify_provenance()` compares requested and applied settings, engine,
interpreter/environment, project, input, script, and receipt ancestry.
Required properties that are omitted become findings. A Hessian that consumes
an optimized structure binds the exact source HDF5 bytes and the molecular
state extracted from that artifact.

For a nonlinear minimum, `frequency_validation_receipt()` requires `3N-6`
finite frequencies and the declared imaginary-mode count. The initial H2O
acceptance slice therefore requires three finite modes and zero imaginary
modes, in addition to SCF/optimization convergence and HDF5 provenance.

GPU4PySCF is an engine of the PySCF program. GPU readiness requires the exact
GPU4PySCF 1.8.0 distribution, CuPy, CUDA driver/runtime, cuTENSOR compatibility,
GPU model/UUID, device count, supported method/basis evidence, and an allowed
stage. The current surface permits `sp|opt|hess` planning and fake preview.
An incomplete or incompatible GPU stack blocks execution; it never selects
CPU.

## 8. Bounded xTB execution contract

`XTBJobSettings` accepts only the `gfn0`, `gfn1`, `gfn2`, or `gfnff` literals
exposed by the v3 xTB I/O vocabulary, strict optimization levels, integer charge, positive
multiplicity, optional SP gradient, and validated solvent pairs. The executable
surface remains CPU `sp|opt|hess`.

The runner must bind the source artifact and, when supplied, the exact project
artifact before rendering. It must bind the actual input bytes used by the
child process, including scratch execution. Stale outputs are quarantined or
rejected so a previous calculation cannot satisfy a new run. Result validation
must compare job kind, method, electronic state, solvent, source/project/input
hashes, exit state, and expected artifacts. Native control text, molecular
dynamics, path following, unsupported constraints, arbitrary shell, and silent
job-type reinterpretation remain outside the surface.

## 9. Runtime V2 contracts and events

The principal additive contracts are:

- `ProgramCapabilityRegistryV1`, `AgentProgramSupportOverlayV1` (alias of
  `ProgramSupportOverlayV1`), `ProgramComponentConformanceReceiptV1`,
  `ProgramCapabilityQueryV1` (alias of `CapabilityQueryV1`), and
  `ProgramCapabilityReceiptV1` (alias of `CapabilityQueryReceiptV1`);
- `EnvironmentTargetV1`, `ProgramEnvironmentQueryV1`,
  `ProgramEnvironmentReceiptV1` (alias of
  `EnvironmentCapabilityReceiptV1`), and
  `TrustedComputeEnvironmentReceiptV1`;
- `ProgramCandidateProposalV1`, `ProgramSubstitutionRequestV1`,
  `ProgramSubstitutionApprovalV1`, and `ProgramSubstitutionReceiptV1`;
- `ResolvedProgramBindingV1` and `ResolvedEngineBindingV1`;
- `ProjectDocumentV1`, `ProjectRenderReceiptV1`, and
  `ProjectValidationReceiptV1`;
- `ScientificIdentityBindingV1`, `CommandProposalV1`,
  `CanonicalCommandInvocationV1`, `CommandInspectionReceiptV1`, and
  `CommandCounterexampleV1`;
- `ArtifactInputIntentV1`, `ArtifactOutputIntentV1`,
  `CommandNodeIntentV1`, `CommandWorkflowDraftV1`, `ArtifactOutputV1`,
  `ArtifactBindingV1`, `CommandNodeV1`, and `CommandWorkflowSpecV1`;
- `SafePreviewReceiptV1`, `ProgramValidatorReceiptV1`,
  `ProgramNodePreflightRequestV1`, and `ProgramNodePreflightReceiptV1`.

Runtime V2 records these through hash-chained events:

`program_capability_queried`, `program_environment_queried`,
`program_substitution_assessed`, `program_binding_resolved`,
`engine_binding_resolved`, `project_validated`, `command_workflow_planned`,
`command_compiled`,
`command_inspected`, `safe_preview_observed`,
`program_validator_observed`, `program_node_preflighted`,
`program_result_verified`, `provider_turn_observed`,
`api_attempt_observed`, and `runtime_terminated`.

Legacy streams replay with empty additive lists. Sequence, previous-event hash,
and idempotency identity are verified. Terminal state is absorbing. A complete
terminal event requires a sorted non-empty set of required receipt digests,
the identical green set, and a host-generated completion-gate digest.

## 10. Model-visible tools and provider boundary

`build_command_compiled_tool_surface()` exposes only:

- `inspect_program_capability`;
- `inspect_program_environment`;
- `assess_program_candidate`;
- project render/read/validate operations;
- a typed command-workflow planning operation that cannot resolve future
  producer hashes;
- `synthesize_command`, `repair_command`, and `preview_command`;
- `preflight_program_node`;
- `inspect_calculation_artifact`.

The H0 surface removes the three program-management tools and uses host-bound
program, capability, and environment evidence. Both profiles retain identical
deterministic gates.

For paired experiments, relevant host-bound evidence must also be made visible
through a typed, canonical public context packet. Event-store seeding alone
does not give a model the references needed by action-grade tool calls. The
packet is identical across H0/H1, contains no paths or authority fields, and
distinguishes semantic planning context from materialized evidence. A missing
or unmatched action-grade reference invalidates that action and its benchmark
grade, not the model's ability to submit a broader workflow draft. H1 adds only
the ability to query the same program-management facts independently.

Provider-specific continuation stays inside its adapter. For the current
DeepSeek V4 Flash adapter, thinking-mode continuation may be preserved in
memory across tool turns, while persisted history and receipts contain only
visible response text, tool calls/results, usage, model identity, finish state,
and public decision summaries. Provider continuation is not scientific
evidence.

## 11. Advisory knowledge packs

`AdvisoryProgramKnowledgePackV1` packages exist for Gaussian, ORCA, xTB,
PySCF CPU, and GPU4PySCF. `KnowledgePackActivationReceiptV1` records scoped
activation and exclusion. Packs may explain conventions or propose fields,
but `readiness_authority` is permanently false.

Pack output must be checked against the current project loader, live CLI,
environment probe, method/basis evidence, and deterministic validators. Packs
cannot broaden the substitution matrix or convert reference-only xTB/GPU
knowledge into execution support. Their benefit must be measured against the
single-agent baseline rather than assumed from expert-sounding prose.

## 12. Ownership, sequencing, and acceptance

| Work package | Primary owner | Produced evidence |
|---|---|---|
| Program CLI/settings integration | program maintainer | run/sub schema projection, strict-loader receipts |
| PySCF generation and execution | PySCF maintainer | input, environment, run, HDF5, provenance, frequency receipts |
| xTB bounded execution | xTB maintainer | source/project/input bindings, result receipt |
| Capability and command plane | agent-runtime maintainer | registry join, overlay, invocation, inspection, preview receipts |
| Runtime and provider adapter | runtime maintainer | replay results, public provider receipts, tool-loop outcomes |
| Scientific review | independent read-only reviewer | typed findings; no repairs or approvals |
| Integration release | repository maintainer | migration ledger, full gates, exact published revision |

Implementation order is fixed: establish strict program semantics; join them to
the live CLI; implement capability and environment observation; compile and
preview typed commands; add substitution; add Runtime V2 events and provider
continuation; run focused integration checks; run approved bounded experiments;
then run the single freeze gate.

Acceptance requires:

- exact upstream v3.1.4 behavior plus registered PySCF and xTB leaves;
- run/sub schema parity for the integrated program leaves;
- no native-input, generated-script, shell, readiness, approval, terminal, or
  silent engine-fallback authority exposed to a model;
- fail-closed project, identity, substitution, environment, preview, and result
  validation;
- replay of existing Runtime V2 streams and deterministic rerendering;
- zero successful terminal states while a required receipt is red;
- for the separately approved CPU H2O slice, exact
  `SP(initial) -> OPT(initial) -> HESS(optimized HDF5)` artifact handoff and all
  declared scientific validators green.

One workflow-level pre-execution approval can bind SP and OPT exactly and bind
HESS to the validated optimized-geometry output selected from that exact OPT
node. After OPT, the host materializes and validates the artifact, then
compiles and previews HESS. The approval remains applicable only while the
producer rule, scientific settings, resources, validators, and command family
are unchanged; otherwise the workflow pauses for new approval.

The H2O calculation, any GPU rerun, scheduler submission, broader engine test,
publication, and comparative performance claim remain separate actions. A
green preview or focused suite is not execution or scientific validation.
