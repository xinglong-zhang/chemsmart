# ChemSmart progressive-assurance design

Status: active design contract. This document corrects an overly restrictive
interpretation of provenance controls. It does not weaken execution safety;
it prevents execution-grade evidence requirements from blocking ordinary
computational-chemistry planning.

## 1. Design principle

ChemSmart uses precision when an action becomes consequential and flexibility
while scientific work is still exploratory.

> Plan broadly, expose uncertainty, materialize incrementally, and require
> exact evidence only at the boundary that consumes it.

Hashes are useful identities for immutable bytes and receipts. They are not a
measure of intelligence, chemical plausibility, or plan quality. A missing
hash must never be the sole reason that the agent cannot discuss a method,
draft a project configuration, compare program candidates, or construct a
semantic calculation DAG.

## 2. Assurance levels

| Level | Permitted work | Required grounding | Normal outcome |
| --- | --- | --- | --- |
| A0 Explore | Interpret a paper or user goal, retrieve knowledge, compare methods, identify missing information | Source labels and explicit uncertainty; no artifact hash required | exploratory plan |
| A1 Draft | Propose stages, project roles, program candidates, command families, dependencies, expected artifacts, and validation criteria | Typed semantic fields; unknowns and alternatives remain explicit | `planned` workflow draft |
| A2 Materialize | Select one CLI-expressible node, validate project YAML, bind geometry/state, and compile canonical argv | Exact source/project artifacts and receipts consumed by that node | `compiled` node or localized finding |
| A3 Preview | Generate ChemSmart-owned native input in fake/test mode and inspect scientific semantics | A2 plus live CLI, environment, parser, and preview evidence | `previewed` node |
| A4 Execute | Invoke a chemistry engine and validate results | Exact approval and executable/environment/artifact bindings | `executed`, `validated`, or failed |

Failure at A2-A4 rolls back only the affected node to the highest supported
lower level. It does not erase the research plan or force the entire workflow
to `blocked` when useful planning can continue.

## 3. Planning behavior

At A0 and A1 the agent should:

- retain every scientifically meaningful stage, including stages not yet
  expressible by the installed CLI;
- distinguish `known`, `supported inference`, `candidate default`, `unknown`,
  and `conflict` rather than collapsing them into accepted/rejected;
- use program registries, scoped knowledge packs, and Basis Set Exchange
  evidence when that registry is available to find likely valid settings and
  present ranked alternatives;
- ask only questions that change the workflow materially;
- isolate an unsupported node and continue planning independent nodes;
- preserve requested observables and scientific purpose when proposing a
  program substitution;
- provide validation and failure criteria even when execution cannot yet be
  prepared.

The next additive planning contracts are:

- `ScientificPlanningDraftV1`: goal, observables, systems, electronic-state
  requirements, evidence statements, explicit unknowns/conflicts, workflow
  drafts, analysis plan, and targeted questions;
- `CandidateBranchV1`: a mutually exclusive program/method/workflow candidate,
  its applicability assumptions, selection evidence, and expected tradeoffs;
- `NodeMaterializationRequestV1` and `NodeMaterializationReceiptV1`: promote
  one selected draft node using only the exact evidence that node consumes.

These contracts are **specified, not yet implemented**. Today,
`CommandWorkflowDraftV1` can preserve stages and unresolved fields, but it
cannot by itself represent mutually exclusive scientific alternatives or a
complete evidence-localized research narrative. Do not force those concepts
into hash fields or free-form command text.

An unfamiliar literal is not automatically invalid. The agent first performs
case-insensitive canonical lookup, registered aliases, engine/version-aware
registry lookup, and bounded similarity retrieval. If one exact supported
candidate is found, it may be proposed with its source. If several candidates
remain, the plan records them and asks a targeted clarification. Only a
deterministic loader/compiler rejection makes a selected executable setting
invalid.

## 4. Draft-to-command transition

`CommandWorkflowDraftV1` accepts semantic programs, job types, project roles,
input slots, producer edges, expected outputs, and unresolved fields. It does
not require capability, environment, project, identity, or artifact hashes.
The host may annotate undeclared programs, unsupported jobs, or unresolved
inputs without rejecting the draft.

Materialization is node-local:

`draft node -> setting resolution -> project validation -> identity/artifact binding -> capability/environment observation -> CommandWorkflowSpecV1 node -> canonical command`

Only this transition requires the exact references consumed by the compiler.
If materialization fails, the host returns a typed finding with alternatives
and the draft remains available for repair.

Producer outputs are semantic edges at A1. For example, an OPT node may
declare an `optimized_geometry` output and a HESS node may consume it. The HESS
node remains planned until OPT actually produces bytes. No predicted hash is
requested from the model.

## 5. Capability and substitution behavior

Capability grades must support exploration:

- `declared` means ChemSmart has a maintained program concept;
- `not_declared` means command execution is not currently integrated, not that
  the scientific step is meaningless;
- `environment_unknown` means availability has not been observed;
- `preview_supported` and `execution_supported` require progressively stronger
  evidence.

During planning, the agent may compare Gaussian, ORCA, xTB, PySCF, and other
scientifically relevant candidates. During command materialization it may emit
only a registered ChemSmart command. Program substitution remains explicit,
but an approval requirement must not prevent the agent from explaining and
drafting the candidate workflow.

## 6. Evaluation reset

The existing fourteen H0/H1 cases are retained as provider-continuation and
boundary-regression cases. They are not the primary measure of generality.
The primary study is specified in
[`../evaluation/progressive-planning-generality-protocol.md`](../evaluation/progressive-planning-generality-protocol.md).

The primary development stream should be positive-first and chemistry-diverse:

1. closed-shell HF and DFT single points;
2. optimization and frequency workflows;
3. open-shell states;
4. solvent models and thermochemical conventions;
5. xTB pre-optimization followed by ab initio refinement;
6. program-unavailable planning with explicit alternatives;
7. multi-stage artifact handoff;
8. paper-derived plans with incomplete but noncritical evidence.

Measure at least:

- scientific-stage coverage;
- preservation of molecule, state, method, basis, solvent, and observable;
- useful-plan rate before all execution evidence exists;
- project and command materialization rate for eligible nodes;
- unnecessary clarification and false-block rate;
- false-ready and false-success rate;
- cross-domain transfer after a change.

Safety boundary cases remain mandatory before execution is enabled, but they
do not dominate prompt design, tool design, or model selection. A change is
retained only when it improves ordinary-task coverage without introducing an
execution-boundary regression.

## 7. Implementation sequence

1. Keep `plan_command_workflow` semantic and free of execution-grade digest
   arguments.
2. Implement `ScientificPlanningDraftV1` and typed candidate branches without
   requiring materialized artifact hashes.
3. Add a deterministic node-materialization operation that consumes only the
   evidence needed by the selected node.
4. Return localized findings and ranked setting alternatives instead of
   terminally blocking a whole draft.
5. Add positive-first generalization fixtures before expanding adversarial
   cases.
6. Run the existing H0/H1 boundary campaign only after the common planning and
   materialization path works on ordinary tasks.

The release still requires exact provenance and approval at preview/execution
boundaries. Progressive assurance changes when those controls apply, not
whether they exist.
