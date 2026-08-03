# Capability-driven program management

This is the maintainer playbook for the integrated program-management plane.
The normative architecture and phase ownership are defined in
[`../design/chemsmart-v314-program-agent-architecture.md`](../design/chemsmart-v314-program-agent-architecture.md)
and
[`../design/chemsmart-v314-integration-implementation-plan.md`](../design/chemsmart-v314-integration-implementation-plan.md).
Use
[`../design/chemsmart-progressive-assurance.md`](../design/chemsmart-progressive-assurance.md)
to keep planning requirements separate from preview and execution gates.

## Authority chain

Keep these questions and evidence grades separate:

1. `chemsmart/settings/capabilities.py` declares the maintained program, job,
   project, and engine inventory.
2. `build_live_click_schema()` observes what this checkout exposes.
3. `AgentProgramSupportOverlayV1` (the public alias of
   `ProgramSupportOverlayV1`) narrows the registry to component-observed job
   and engine coverage.
4. `EnvironmentCapabilityReceiptV1` observes the selected executable or exact
   compute interpreter and dependencies.
5. Project loaders and program validators decide whether the scientific request
   is representable.
6. The command compiler and safe preview decide whether canonical argv creates
   valid downstream artifacts.
7. Cross-receipt preflight derives plan state.
8. A separate exact approval decides whether execution may begin.
9. Program-specific result verification determines whether an executed result
   is validated.

A green item never implies the next item.

## Program inventory maintenance

`PROGRAM_CAPABILITIES` is the only immutable declaration. Agent code must call
`load_program_capabilities()` and must not carry a second program enum.

When adding or changing a program:

1. Register the Click group through `chemsmart/cli/subcommands.py` so `run` and
   `sub` share it.
2. Update the capability declaration with sorted job, parameter, engine, and
   loader-bounded parameter-domain tuples.
3. Add a strict settings loader and define whether a project is required,
   optional, or unsupported.
4. Add program-specific preview expectations and result verification.
5. Generate `ProgramComponentConformanceReceiptV1` for only the observed jobs
   and engines.
6. Build a support overlay from those receipts. Do not mark unobserved siblings
   operable.
7. Add additive Runtime V2 events and replay defaults if a new receipt family
   is required.

Use these support meanings:

- `available`: compiler, preview, preflight, verifier, and execution evidence;
- `preview_only`: compiler, preview, preflight, and verifier evidence, but no
  execution evidence;
- `reference_only`: declaration or documentation only;
- `disabled`: intentionally absent from the active profile.

The default overlay starts every program at `reference_only`. Registry or CLI
presence cannot promote it. PySCF CPU and xTB may become `preview_only` only
after exact compiler, preview, preflight, and verifier conformance receipts are
bound to the current schema and fixtures. GPU4PySCF `sp|opt|hess` may become
`preview_only` after current fake-preview conformance. Execution additionally
requires a complete CUDA/GPU4PySCF environment observation; the non-CUDA
integration host remains blocked and never falls back to CPU.

## Integrated CLI/project matrix

| Program | Maintained new jobs | Project rule | Engine rule |
|---|---|---|---|
| PySCF | `sp`, `opt`, `hess` | strict `sp`/`opt`/`hess` YAML required | CPU or explicitly selected GPU |
| xTB | `sp`, `opt`, `hess` | GFN2 defaults if omitted; strict YAML if supplied | CPU only |
| Gaussian | upstream v3.1.4 leaves | required | CPU |
| ORCA | upstream v3.1.4 leaves | required | CPU |
| NCIPLOT | top-level leaf | unsupported | CPU |

The corresponding integrated command paths are `run|sub pyscf sp|opt|hess`
and `run|sub xtb sp|opt|hess`. Model output supplies neither paths nor option
placement; `compile_command()` derives both from the live Click schema.

## Substitution operations

The only non-identity mapping is Gaussian to PySCF for:

- single point;
- minimum optimization;
- fixed-geometry Hessian/frequencies;
- ordered minimum optimization followed by a Hessian over the exact optimized
  HDF5 geometry.

Use `build_program_substitution_request()` followed by
`assess_typed_program_substitution()`. Eligibility requires HF or DFT, a
uniform basis, no unsupported constraint, source evidence, and verified
functional semantics for DFT. A non-identity result remains rejected until an
exact `ProgramSubstitutionApprovalV1` is supplied.

Block materialization of a Gaussian-to-PySCF substitution for TS, IRC, TD,
scan, QMMM/ONIOM, NEB, post-HF, double hybrid, mixed basis, ECP, mixed
basis/ECP, and unsupported constraints. Keep those requested scientific stages
in the workflow draft with an explicit capability gap and possible program
candidates. Missing Gaussian is not evidence that PySCF is equivalent. Never
choose CPU when GPU was selected or choose another program because the
selected environment is red.

## PySCF operational checks

Before preview or execution:

- load one complete `sp`/`opt`/`hess` project section;
- validate method, basis, solvent, state, stage, solver, and engine;
- probe the exact compute interpreter;
- require PySCF 2.14.0, NumPy, and h5py evidence;
- require basis, functional, and optimizer evidence applicable to the request;
- for GPU, require GPU4PySCF 1.8.0, its installed distribution, CuPy, CUDA,
  cuTENSOR, device model/UUID, and allowed-stage evidence.

The generated script, environment receipt, source geometry/artifact, project,
and requested settings are hash-bound before launch. The result HDF5 is the
numeric authority and contains `spec`, `provenance`, `status`, and `results`.
After launch, require exact receipt ancestry, requested/applied settings
agreement, child exit status, required-property presence, and applicable
frequency validation. Text logs are diagnostic views only.

## xTB operational checks

Accept only CPU `sp`, `opt`, and `hess` with a validated
`gfn0`/`gfn1`/`gfn2`/`gfnff` method, electronic state, optimization level,
optional SP gradient, and
complete solvent pair. Preserve source charge/multiplicity unless both fields
are explicitly overridden by an approved layer.

Bind source and optional project before rendering, then bind the exact child
input, including its scratch path. Quarantine stale outputs. Compare the final
result with job kind, method, state, solvent, source, project, input, process,
and expected-artifact receipts. Treat a missing or mutated binding as red.

## Runtime/tool maintenance

The normal command-compiled profile exposes semantic workflow planning,
program capability/environment inspection, candidate assessment, typed project operations, command
synthesis/repair/preview, program-node preflight, and result inspection. The H0
profile removes only the three program-management operations and uses
host-prebound evidence.

Never expose native-input builders, generated-script editing, arbitrary
filesystem paths, shell strings, executable claims, approval, readiness, or
terminal state. A successful `command_compiled` event is not preview evidence.
A useful semantic DAG remains in product state `planned`. A planning-only
Runtime V2 task terminates `planned` when its workflow-draft receipt is
observed; that label is not command, execution, or scientific success. The
product must not be downgraded to `blocked` solely because later action-grade
inputs are unresolved.
A complete terminal event requires every declared receipt gate to be present
and green.

For a multi-command workflow, an unresolved producer-output edge is planning
data only. It keeps the consumer `planned` and cannot enter command synthesis,
preview, preflight, or a node-specific exact-command approval. After the
producer artifact is observed, the host creates a new exact binding and
refreshes downstream command evidence. A workflow-level approval may cover the
consumer when it explicitly binds the producer node, output-selection rule,
scientific settings, resource limits, validators, and invalidation conditions.
A changed semantic field, producer, or command family pauses the workflow for
new approval.

When adding a receipt event, update both `runtime/events.py` validation and
`runtime/reducer.py`. Older event streams must replay with empty additive
state. Preserve sequence, hash-chain, idempotency, and absorbing terminal
semantics.

## Knowledge-pack maintenance

Built-in packs for Gaussian, ORCA, xTB, PySCF CPU, and GPU4PySCF are advisory.
They may propose or explain settings, but must have
`readiness_authority=False`. Activation is scoped by program/engine terms and
exclusions and is recorded in `KnowledgePackActivationReceiptV1`.

Do not add a fact merely because it appeared in an earlier model response.
Bind sources, version the pack, and add a deterministic validator or
environment observation when the fact affects readiness. Evaluate pack-on and
pack-off behavior under the same prompts, tools, sources, and graders.

## Failure routing

| Failure | Required disposition |
|---|---|
| Unknown program/job/engine or Click mismatch | retain the draft and localize a capability gap; block only materialization |
| Environment target missing or mismatched | retain the draft; block preview/execution for that candidate and do not fall back |
| Project loader invalid | retain the semantic node; return field findings and block its command materialization |
| Functional equivalence unverified | retain both program identities; block substitution materialization |
| Input, geometry, project, script, environment, or result hash mismatch | fail current attempt; invalidate approval |
| Safe preview artifact missing or program validator red | block execution |
| Child process nonzero | fail and preserve exact exit status |
| Required property omitted | fail validation with explicit finding |
| GPU evidence incomplete | block GPU; never use CPU implicitly |
| Provider/tool transport fault | record separately; do not change scientific readiness |

## Change acceptance checklist

- Registry and live Click projection agree for the changed program.
- `run` and `sub` reconstruct equivalent program settings.
- Project omission and invalid-project behavior match the declared policy.
- Exact interpreter/executable evidence exists before availability is claimed.
- Planning tools contain semantic IDs and explicit unknowns; materialization
  tools use the exact typed references they consume, never paths or shell.
- Preview uses fake, no-scratch, and scheduler-test modes as applicable.
- Project, source, input, environment, and result mutation are detected.
- No program or engine silently replaces a blocked one.
- Legacy event replay remains valid.
- Success is impossible while any required deterministic receipt is red.

Record focused implementation evidence first. Defer broad integration QA until
the architecture milestone is stable, then run the preregistered focused gate
once and the full freeze gate once.
