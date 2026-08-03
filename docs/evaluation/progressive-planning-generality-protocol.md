# Progressive planning generality protocol

## Purpose

This is the primary model-development protocol for the integrated ChemSmart
agent. It asks whether the agent can turn varied computational-chemistry goals
into useful, scientifically organized ChemSmart workflow drafts. It does not
require every task to arrive with materialized coordinates, project files, or
artifact hashes.

The separate fourteen-case DeepSeek H0/H1 protocol is a later boundary
regression. It tests typed tools, false completion, fallback, and provenance;
it is not the main measure of scientific generality.

## Progressive-assurance contract

Planning proceeds through three assurance levels.

### P0 — semantic research draft

Required output:

- scientific objective and requested observables;
- molecular systems and electronic states, including explicit unknowns;
- proposed program and job families with alternatives where justified;
- stage-specific method, basis/ECP, solvent, constraint, temperature, and
  standard-state requirements;
- ordered `CommandWorkflowDraftV1` nodes, dependencies, project roles,
  expected artifact classes, analysis steps, and failure branches;
- evidence needed to resolve every consequential unknown;
- classification of each stage as CLI-expressible, capability-gap, or awaiting
  evidence.

P0 may use human-readable source locators and stable semantic IDs. Missing
coordinates, a future optimized geometry, an unselected project, or an
unobserved environment stays explicit in `unresolved_fields`; it does not erase
the workflow or force a terminal failure. No executable command, approval, or
readiness claim is made.

### P1 — deterministic materialization

When sufficient evidence exists, the host resolves exact molecule identity,
geometry bytes, charge/multiplicity, project YAML, program/engine, and current
CLI schema. Only here are artifact, project, schema, environment, and receipt
hashes mandatory. The host promotes eligible draft nodes into
`CommandWorkflowSpecV1`, canonical commands, and safe-preview requests.

Nodes that remain unresolved preserve their scientific role and missing-input
requirements. Materializing one branch must not silently delete or simplify
another.

### P2 — action and scientific validation

Preview, approval, execution, result verification, and claims use exact
immutable bindings. These stages retain the existing permission, provenance,
and deterministic-validator contracts. P0 flexibility does not weaken P2
evidence requirements.

## Generality task set

Pre-register development and held-out tasks before reading model outcomes.
Include at least two tasks in each family:

1. ground-state structure, single-point energy, and vibrational analysis;
2. reaction mechanism with TS/IRC or an explicit unsupported-program branch;
3. thermochemistry, solvation, temperature, and standard-state corrections;
4. excited-state, photochemical, or spectroscopic planning;
5. open-shell, transition-metal, or basis/ECP-sensitive planning;
6. conformer, noncovalent, cluster, or xTB-assisted multistage planning.

At least half of the held-out tasks must differ from implementation fixtures in
both molecule/system and workflow shape. Do not use hidden expected answers to
change prompts, schemas, validators, knowledge packs, or routing rules.

Tasks intentionally vary evidence completeness:

- complete method and coordinate evidence;
- complete methods but missing coordinates;
- partial method evidence requiring a targeted clarification;
- conflicting program terminology;
- one CLI-expressible branch plus one legitimate capability gap.

The desired behavior is not universal completion. It is maximum useful plan
coverage with exact separation of known, inferred, conflicting, and unknown
facts.

## Compared configurations

Retain a single-coordinator reference path.

- `G0`: typed semantic planning plus the live CLI/capability summary.
- `G1`: `G0` plus scoped program knowledge packs and the common typed workflow
  planner.
- `G2`: `G1` plus bounded evidence/specialist decomposition only where tasks
  have independently verifiable source or species partitions.

Use the same model, prompt template, task evidence, token/wall-time envelope,
task order policy, and semantic graders. Exact materialization tools may be
available only for tasks with supplied artifacts, but that availability is
fixed across compared configurations for the same task.

Knowledge packs and specialists are evaluated for added scientific coverage,
not for producing more authoritative language. A coordinator owns final
molecular identity, electronic state, project roles, DAG integration, and
readiness.

## Primary measurements

Score the P0 draft before attempting P1:

- requested scientific-stage coverage;
- correct program and ChemSmart job-family selection;
- preservation of molecule, state, method, basis/ECP, solvent, constraints,
  temperature, standard state, units, and observable intent;
- correct dependency ordering and expected artifact classes;
- explicit and actionable unknowns rather than fabricated settings;
- useful alternatives when multiple methods are scientifically legitimate;
- CLI-expressibility and capability-gap classification;
- clarification precision: necessary questions asked, unnecessary questions
  avoided;
- domain-expert rubric score, visible response quality, model calls, tokens,
  latency, and cost.

Secondary P1 measurements apply only to resolvable nodes:

- project-loader validity;
- draft-to-materialized intent preservation;
- canonical command and safe-preview success;
- false materialization and false-ready rate;
- fraction of unresolved nodes correctly retained rather than dropped.

Do not lower a P0 score merely because an honest draft awaits user-provided
coordinates. Do lower it for omitting the geometry requirement, guessing the
geometry, or claiming preview readiness.

## Retain/revise rule

Retain a component only when it improves held-out P0 scientific coverage or
clarification quality without increasing fabricated consequential facts or
cross-domain regressions. Report per-family results and negative results; an
aggregate average cannot hide a failure concentrated in one chemistry family.

After the candidate planning profile is frozen, run the fourteen-case boundary
regression and deterministic materialization checks. Boundary failures may
disable action-stage tools, but they do not erase the measured P0 planning
capability. Conversely, strong P0 prose cannot authorize P1/P2 actions.

## Evidence and reporting

Persist visible English responses, typed drafts, explicit unknowns, tool
calls/results, configuration and prompt versions, source locators, usage,
latency, graders, and dispositions. Exact hashes are recorded whenever bytes
exist, but absence of bytes at P0 is represented as a typed unknown rather than
a fabricated digest.

Report separately:

1. semantic planning quality;
2. deterministic materialization coverage;
3. boundary-regression outcomes;
4. execution/scientific validation, if separately approved later.

No result from this protocol alone establishes reproduction, execution, broad
autonomy, or state-of-the-art performance.
