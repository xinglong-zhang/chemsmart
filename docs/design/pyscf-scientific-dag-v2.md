# PySCF Scientific DAG V2

## Purpose

ChemSmart controls PySCF through project YAML and canonical CLI commands.  A
model proposes scientific intent; it does not write the standalone Python
driver, choose executable paths, construct shell, approve execution, validate
its own result, or set terminal state.

The DAG joins the model-visible plan to host-owned materialization, preview,
approval, execution, and evidence without introducing a second runtime beside
Runtime V2.

## Scientific graph

For a minimum-search workflow, SP and OPT consume the same initial geometry
but neither calculation is scientific input to the other.  HESS and the
experimental TD preview consume the validated optimized geometry:

```text
initial geometry --+--> SP(initial)
                  |
                  +--> OPT(initial) --+--> HESS(optimized)
                                      +--> TD/TDA preview(optimized)
```

A **control edge** expresses scheduling only.  A **data edge** names the
producer output, consumer input, artifact class, and validated handoff.  A
consumer cannot be materialized from a future hash invented by a model.

Scientific identity binds only an exact geometry artifact, never a project,
result, or filename-derived guess.  Identity must be green before the host can
construct the canonical scientific plan.  A planning response that contains a
useful draft but `scientific_workflow_plan=null` is not a DAG.  If identity or
another prerequisite is repaired later, the coordinator must invoke planning
again and obtain a non-null typed plan before any workflow-level readiness or
edge claim is allowed.

## Authority layers

1. `ScientificWorkflowPlanV2` records task-bound stages, unresolved fields,
   and typed edges.
2. `MaterializedWorkflowV1` joins resolvable nodes to identity, project,
   environment, command, preview, and preflight receipts.
3. `FrozenWorkflowApprovalV1` can be created only from a materialized workflow
   and exact resources.  Each future producer-edge target is bound to one
   exact environment receipt rather than merely sharing a workflow-wide set.
4. `WorkflowRunStateV1` records launch and node outcomes so replay can recover
   the executable frontier.
5. Program-specific validators determine whether outputs are scientifically
   usable.  Model prose cannot promote an execution state.

## Specialist topology

The default remains one coordinator.  A deterministic complexity gate may
request fresh, isolated specialists for scientific method/state, PySCF project
settings, or DAG structure.  Each worker receives an immutable
`ContextManifestV1` and `SpecialistTaskPacketV1`; it returns typed data and
cannot modify the coordinator's project or workflow.  A fresh critic is
read-only and cannot repair, approve, execute, or set readiness.

Specialist proposals use advisory scientific paths only.  Environment
observations may be summarized, but host-owned keys such as readiness,
execution readiness, approval, validation, and terminal state are invalid in
a specialist proposal.  The coordinator imports validated worker records
before it acts; it does not treat specialist prose as a materialized project
or graph.

For controlled harness experiments, three potential pre-coordinator worker
slots and one critic slot are budgeted invariantly.  Each disposable
participant receives the same fixed fraction whether C is off or on, and an
unused slot is not donated to another participant.  Feedback receipts from
workers and a critic are retained before coordinator receipts are appended, so
participant usage and provider-visible projections can be reconciled without
silently reducing the D arm to coordinator-only evidence.

The topology borrows only bounded principles from frontier agent research:
typed execution graphs and admissible transitions from El Agente Grafico,
specialist roles from El Agente Q, and progressive disclosure from El Agente
Forjador.  ChemSmart deliberately excludes model-authored native input,
runtime tool forging, unrestricted peer swarms, shell authority, and automatic
engine retry.

## Current PySCF boundary

- CPU: SP, unconstrained minimum OPT, and fixed-geometry HESS.
- GPU4PySCF: the same stages may be planned and previewed; execution requires
  a complete matching GPU environment and never falls back to CPU. A fake
  preview may render the exact GPU request without a local allocation, but it
  remains explicitly execution-inadmissible.
- Experimental TD: closed-shell singlet RKS, gas-phase TDA or TDDFT, positive
  explicit state count, preview only. The inert artifact carries a typed
  response-materialization manifest naming the RKS response factory, state
  channel, state count, and operation order; it contains no runnable response
  kernel.
- Charge and multiplicity overrides are an indivisible pair in the loader,
  settings validator, generic preflight, CLI builder, and job constructor.
  The canonical density-fitting default is false in both settings and raw-map
  preflight.
- Unsupported: TS, IRC, scans, constraints, post-HF, double hybrids, mixed
  basis/ECP, excited-state optimization, TD gradients/Hessians, spin flip, and
  SOC.

Unsupported requests remain useful typed plans with an explicit reason.  They
never fall through to a hand-written PySCF script.

## Replay and interrupted dispatch

The replayable workflow object and engine launch state are separate.  Before a
provider episode is dispatched, the campaign host writes an atomic intent
bound to the plan and exact implementation freeze.  A complete evidence
envelope can be restored only when it matches that intent.  An intent without
an envelope after a crash is preserved as an inconclusive interrupted
observation and is not redispatched under the same episode identity.  A source
change starts a new freeze; it does not mutate or reinterpret the earlier
workflow evidence.

This campaign recovery rule is stricter than ordinary planning replay because
an external model call may have occurred even when its result was not durably
captured.  It protects exactly-once observation semantics; it does not make an
interrupted scientific workflow complete.

## Evaluation boundary

The Qwen campaign varies decomposition, provider feedback projection, and a
fresh critic while holding the deterministic host fixed.  Full host results
remain authoritative in every arm; compact feedback is only a provider view.
The primary C factor is descriptive unless a separate seeded-oracle study is
run; critic recall and false-rejection claims cannot be inferred merely from
the presence of findings.

Transfer binds H0 and the selected candidate as whole profiles.  Its typed
terminal receipt covers 80 frozen episodes and 40 within-task pairs, but the
chemistry generalization unit is four untouched task cases.  Ten sessions per
case estimate stochastic variability rather than create 40 independent
chemistry tasks.  A multifactor candidate therefore supports only a
whole-profile comparison, not a decomposition-, feedback-, or critic-specific
causal claim.

The v12 gate demonstrated why leaf previews are insufficient: identity was
first bound to the project artifact, planning returned a null scientific plan,
identity was repaired, and SP/OPT previews were later produced without
replanning.  The terminal was complete and safety was green, but the
deterministic DAG checks failed because no canonical SP/OPT sibling and
OPT-to-HESS data-edge graph existed.  This is negative integration evidence,
not a workflow success.  No harness change is retained from one fixture alone;
a fresh exact-source gate and cross-case evidence are required.

The v13 gate then showed that the geometry-only identity and replan
affordances can produce the canonical SP/OPT/HESS graph and safe sibling
previews while an independent specialist wire-format lapse still invalidates
the D factor. Specialist output now records whether it was raw strict JSON or
a bounded normalization of exactly one JSON object after a plain-prose prefix.
The normalized object still passes the unchanged ownership and host-authority
validators, and its receipt prevents normalized acceptance from being counted
as raw schema success.

The subsequent v14 exact-source gate was the first to realize all three
specialists and the coordinator through a safety-green, deterministically
previewed SP/OPT/HESS graph. All three specialist envelopes were raw strict
JSON. This single-case pass establishes active-path integration only. Later
parser hardening rejects duplicate JSON keys, validates envelope digests at
the factorial gate, and exposes strict-versus-normalized counts; therefore the
factorial must use a newer implementation freeze rather than treating v14 as
its source baseline.

The v15 gate repeated the same safety-green preview outcome while exercising
output-envelope receipt admission and a longer bounded repair path. A rejected
scientific decision could not cite an unapproved identity; the coordinator
later rebound it correctly without changing the approved geometry or
electronic state. The 35-minute latency is a harness observation, not evidence
of better chemistry, and constrains how much transfer evaluation fits in the
Goal window.

The scientific DAG grader treats preview-only excited-state work as a future
artifact boundary: the TD node must receive geometry from OPT through an
explicit `geometry_xyz` data edge and cannot be previewed until that producer
artifact exists. The exact response method and state count remain project
semantics, while the host retains identity, gas-phase, reference, and terminal
authority. Mutation and context-resume testing are not aliases for this
boundary; they need separately typed perturbation and transition receipts.
