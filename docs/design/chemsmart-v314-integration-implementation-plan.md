# ChemSmart 3.1.4 integration implementation plan

Status: active implementation plan. Integration QA and publication gates are
intentionally deferred until the implementation milestones below are complete.

Current milestone state is tracked in
[`../maintenance/integration-status-v314.md`](../maintenance/integration-status-v314.md).
The planning-versus-execution assurance boundary is defined in
[`chemsmart-progressive-assurance.md`](chemsmart-progressive-assurance.md).
The live provider procedure is preregistered in
[`../evaluation/deepseek-v4-flash-h0-h1-protocol.md`](../evaluation/deepseek-v4-flash-h0-h1-protocol.md),
and the release procedure is fixed by
[`../maintenance/release-source-custody-runbook.md`](../maintenance/release-source-custody-runbook.md).

## 1. Delivery target

The release line starts from upstream ChemSmart v3.1.4 revision
`5d486b34501b9c040ad1db30d41a0eaa3850b15f`. The package version remains
`3.1.4`; the fork build is identified by Git revision and joined live-CLI and
capability digests.

The release is delivered as six separately reviewable changes:

1. PySCF 2.14.0 and GPU4PySCF 1.8.0 workflows.
2. Bounded xTB 6.7.1 execution.
3. Capability-driven agent program management.
4. DeepSeek V4 thinking/tool continuation.
5. Program and agent integration tests.
6. Migration and reproducibility documentation.

No benchmark claim or model-quality result is part of the implementation
definition. Those are later empirical outcomes.

## 2. Workstream rules

- Preserve the exact upstream base and the source-custody archive.
- Adapt changes to v3.1.4 conventions; do not copy historical schema digests,
  hard-coded agent program lists, private experiment evidence, or secrets.
- Keep `chemsmart/settings/capabilities.py` as the sole immutable program
  declaration and derive agent projections from it.
- Register integrated programs through the shared `subcommands` list so `run`
  and `sub` observe the same Click groups.
- Require typed settings, trusted artifacts, canonical argv, and deterministic
  receipts at each boundary.
- Do not replace upstream Gaussian, ORCA, or xTB I/O parsers where the v3.1.4
  implementation is already authoritative.
- Do not run broad QA while implementation surfaces are still changing. Use
  source inspection and narrow deterministic probes, then one focused suite at
  the M3 entry gate and one full freeze gate in M4.

## 3. Milestones

### M0 — source custody and neutral integration lineage

Owner: repository maintainer.

Deliverables:

- verified bundles for committed refs;
- binary patches for tracked modifications;
- path-safe archives for non-secret untracked evidence;
- SHA-256 manifest, exclusions, source metadata, and restoration drill;
- path-level migration ledger with `ported`, `adapted`, `archive_only`, or
  `rejected` disposition;
- a neutral integration branch from the exact upstream revision.

Exit evidence: restoration succeeds in a disposable directory and every source
path has one disposition. Worktrees remain retained.

### M1A — PySCF and GPU4PySCF program implementation

Owner: PySCF maintainer.

Implementation:

- preserve `run/sub pyscf sp|opt|hess` and standalone generated scripts;
- enforce strict three-section project YAML and mandatory
  `PySCFJobSettings.validate()`;
- make scientific preflight and post-result `verify_provenance()` mandatory;
- bind exact compute interpreter and PySCF 2.14.0, NumPy, h5py, solver, basis,
  functional, and optional GPU evidence;
- propagate generated child exit status through local/scheduler wrappers;
- bind script, source, geometry, project, requested/applied settings,
  environment, HDF5, and result hashes;
- make requested properties explicit and convert omission to findings;
- materialize the requested title and reject inherited unsupported settings;
- bind optimized-HDF5 geometry when a later Hessian consumes it;
- allow GPU `sp|opt|hess` planning and fake preview, require a complete
  CUDA/GPU4PySCF environment for execution, and never fall back to CPU.

Exit evidence: all program paths can produce internally consistent preview
receipts and all negative boundaries fail closed. Real execution is not part of
this milestone.

### M1B — bounded xTB execution implementation

Owner: xTB maintainer.

Implementation:

- register CPU `sp|opt|hess` under both `run` and `sub`;
- retain the v3 xTB parser and molecule I/O;
- enforce strict method, optimization, job, state, gradient, and solvent
  settings;
- support optional project YAML while rejecting unresolved explicit projects;
- pre-bind source and project bytes, then bind the actual child input bytes;
- quarantine/reject stale outputs and verify fresh result identity;
- reject arbitrary control input, MD, path following, unsupported constraints,
  and conflicting job semantics.

Exit evidence: source/project/input/result bindings are stable for normal and
scratch paths, and the preview validator cannot accept stale or mismatched
artifacts.

### M1C — canonical capability plane

Owner: settings maintainer.

Implementation:

- retain a single immutable registry for Gaussian, ORCA, xTB, PySCF, and
  NCIPLOT;
- derive job/engine/project projections from that registry;
- bind the registry to the observed Click schema and record the joined digest;
- distinguish declaration from agent conformance and environment evidence.

Exit evidence: no independent agent program inventory and no capability status
that infers installation.

### M2A — capability-driven Runtime V2

Owner: agent-runtime maintainer.

Implementation:

- port Runtime V2 events, reducer, permissions, command compilation, safe
  preview, project tools, and artifact/result inspection;
- implement coverage-bound `AgentProgramSupportOverlayV1` (the public alias of
  `ProgramSupportOverlayV1`) and component conformance;
- implement capability and exact-environment query receipts;
- implement requested/selected program and engine bindings;
- implement cross-receipt `ProgramNodePreflightRequestV1` and receipt;
- make legacy event replay default additive program state to empty;
- require receipt-complete host gating for terminal success.

Exit evidence: model output cannot set executable status, readiness, approval,
or terminal state; unknown or unsupported cases block with typed rule IDs.

### M2B — project, command, and substitution tools

Owner: command-plane maintainer.

Implementation:

- expose only host-bound program inspection, typed project operations,
  command-workflow planning, command synthesis/repair/preview, preflight, and
  result inspection;
- add an evidence-aware `ScientificPlanningDraftV1` and typed mutually
  exclusive candidate branches so alternatives, conflicts, and unknowns do
  not have to masquerade as execution-grade bindings;
- keep model arguments free of arbitrary paths, shell, option placement, and
  native input;
- compile through the live Click tree and independently inspect the parse;
- implement the closed Gaussian-to-PySCF matrix for SP, minimum optimization,
  fixed Hessian, and ordered optimization plus Hessian;
- construct multi-command DAGs with stable producer-output IDs while keeping a
  consumer `planned` until the producer's actual artifact hash and receipt are
  available; never preview a Hessian against an invented future optimized
  geometry;
- materialize one selected node at a time and return localized findings without
  erasing the rest of the semantic plan;
- keep semantic workflow drafting available when capability, project,
  geometry, or environment evidence is incomplete; record those gaps as
  unresolved fields and require exact bindings only when materializing an
  executable command node;
- require uniform basis, HF/DFT method family, source evidence, verified DFT
  semantics, no unsupported constraints, and exact user approval;
- block TS, IRC, TD, scan, QMMM/ONIOM, NEB, post-HF, double hybrid, mixed basis,
  and ECP requests.

Exit evidence: missing Gaussian availability never causes automatic PySCF
selection, and no approved substitution alone establishes readiness.

### M2C — advisory program knowledge

Owner: scientific-knowledge maintainer.

Implementation:

- version scoped packs for Gaussian, ORCA, xTB, PySCF CPU, and GPU4PySCF;
- bind activation/exclusion and source IDs in a receipt;
- keep `readiness_authority=False` and xTB/GPU reference status explicit;
- measure pack-on versus pack-off behavior using the same deterministic gates.

Exit evidence: removing a pack changes explanation or proposal quality only;
it never removes a safety, loader, environment, or validator requirement.

### M2D — provider continuation

Owner: provider-adapter maintainer.

Implementation:

- use the official DeepSeek V4 Flash endpoint through a one-use secret lease;
- enable thinking-mode tool calls at provider-supported limits;
- preserve required reasoning continuation only in memory between tool turns;
- persist visible English text, tool calls/results, public decisions, usage,
  latency, errors, and model identity;
- sanitize private continuation from events and artifacts;
- make provider transport failure independent from scientific readiness.

Exit evidence: the real `UnifiedSessionRunner -> ToolLoopRunner ->
CommandCompiledToolHostV1` path can continue across tools without credential or
private-continuation persistence.

### M3A — focused integration gate

Owner: integration maintainer.

Run once after M1 and M2 stop changing materially. Cover PySCF/xTB CLI and
settings, run/sub reconstruction, HDF5 provenance, capability/schema joining,
legacy replay, project tools, substitution, provider continuation, malformed
tool calls, shell/path rejection, and false completion. Permit one
evidence-driven rerun after identified defects are corrected.

Exit evidence: focused results and exact residual defects are recorded. A green
focused suite is implementation evidence, not scientific validation.

### M3B — positive-first real provider experiment

Owner: experiment maintainer; scientific-coverage graders own outcomes.

This is the primary agent-development signal and follows
[`../evaluation/progressive-planning-generality-protocol.md`](../evaluation/progressive-planning-generality-protocol.md).
Run ordinary chemistry tasks through semantic workflow drafting and node-local
materialization before adversarial cases. Cover closed-shell HF/DFT,
optimization/frequencies, open-shell state, solvent/thermochemistry, xTB to
ab-initio staging, explicit program alternatives, and multi-stage handoff.
Measure visible English planning quality, scientific-stage coverage,
false-block rate, unnecessary clarification, and the fraction of eligible
nodes that materialize successfully. Missing noncritical evidence may remain
explicit without terminally blocking the whole draft.

Exit evidence: common tasks produce useful drafts, and fully grounded eligible
nodes traverse project loading, compilation, and safe preview without
case-specific rules.

### M3C — H0/H1 boundary regression

Owner: experiment maintainer; deterministic graders own outcomes.

This regression is run after M3B and is not used as the primary generality or
scientific-planning score.

Run paired H0/H1 English natural-language cases through the real active path.
Before scheduling, materialize a typed public context packet that separates
semantic planning facts from action-grade artifact IDs and receipt references.
Bind the same canonical packet into both arms' prompts and plan identities;
event-only evidence does not satisfy the H0 contract. Missing action-grade
evidence makes the affected command or preview inconclusive, but must not
prevent a scientifically useful workflow draft with explicit unknowns.
The context-packet factory and typed workflow-planning operation are exposed in
the integration source. The campaign cannot start until their active-session
behavior and arm equality are observed at the deferred focused gate.
Cases cover capability queries, PySCF SP/OPT/HESS planning, unavailable-Gaussian
substitution, Gaussian-only blocking, unsupported TD/TS/IRC, functional
semantics, ECP rejection, open-shell state, GPU evidence, no fallback, shell
rejection, false completion, and HDF5 mismatch.

Each request records a unique hypothesis, one changed factor, expected
observation, deterministic oracle, tool/configuration hashes, usage, latency,
and non-secret error class. Attempts are bounded by authorized account quota
and wall time, not by an arbitrary call count. No top-up occurs.

Exit evidence: visible provider responses and typed receipts are graded
separately for wire validity, tool validity, scientific-state preservation,
preview outcome, and terminal honesty.

The OPT-to-HESS case is successful when the host accepts an ordered workflow,
previews the resolved OPT node, and leaves HESS planned with an unresolved
typed producer edge. Requiring two previews before OPT execution is retired as
an invalid oracle because the optimized-coordinate hash does not yet exist.

### M3D — approved H2O CPU slice

Owner: execution maintainer; user owns approval.

Prepare exact approval material for neutral singlet H2O, gas-phase
B3LYP/def2-SVP, CPU, four cores, 4 GB, and no scratch:

`SP(initial XYZ) -> OPT(initial XYZ) -> HESS(exact optimized HDF5 geometry)`

Use one workflow-level approval that exactly binds SP and OPT, the scientific
settings, source, interpreter/environment, resource limits, ordered stages,
and the deterministic rule "HESS consumes the validated optimized geometry
emitted by this OPT node." HESS remains `planned` until that artifact exists.
After OPT, the host resolves the artifact, compiles and previews HESS, and may
continue under the same approval only when every approved semantic field and
the producer-selection rule remain unchanged. Any method, basis, state,
resource, producer, validator, or command-family change pauses for new
approval. No retry or engine substitution is automatic.

After approval, validation requires SCF convergence, optimization convergence,
valid structured HDF5, requested/applied provenance agreement, three finite
nonlinear-water modes, and zero imaginary modes.

### M4 — freeze, commits, and publication of the repository branch

Owner: repository maintainer.

Before the report checkpoint, run the selected focused suite and read-only
lint once, followed by credential, path, diff, and public-tree scans. Defer the
full repository suite and remote release gates until a separately authorized
publication checkpoint. Do not autofix or regenerate snapshots.

Release vetoes are any native-input/shell authority bypass, false-ready state,
artifact substitution, secret leak, silent engine fallback, or success with a
red deterministic gate. If provider cases fail while deterministic code is
sound, keep the new program-management profile disabled by default.

Create three neutral local commits for programs, agent, and tests. Verify exact
authorship and stop for review before any fetch, push, pull request, or remote
replacement. Do not alter upstream or remove archived refs.

## 4. Milestone report schema

Every milestone report must contain:

1. changed files and contracts;
2. observable result and receipt identifiers;
3. failures and causal classification;
4. unknown or untested claims;
5. retain, revise, disable, or reject decision;
6. next owner, prerequisites, and approvals.

Use `planned`, `previewed`, `executed`, `validated`, `blocked`, and `failed`
literally. Never infer a higher state from a lower-state receipt.
