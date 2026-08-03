---
name: chemsmart-agent-harness
description: Design, audit, test, or document ChemSmart's provider-neutral command-compiled agent harness, including typed CommandWorkflowSpec proposals, CLI-schema grounding, deterministic compilation and safe preview, provider adapters, tool exposure, permissions, Runtime V2, event replay, task graphs, bounded subagents, and harness evaluations. Use when changing or assessing chemsmart/agent runtime, provider, loop, registry, permission, command synthesis, or agent-architecture behavior.
---

# ChemSmart Agent Harness

Use this skill to keep the agent loop auditable, bounded, and independent of a
single model provider. Read `AGENTS.md` first; it supplies repository-wide
authority, safety, and evidence rules.

## Working procedure

1. Inspect the active branch, dirty state, affected runtime contracts, and
   focused tests before proposing a change.
2. Convert a model proposal into typed ScientificTaskSpec and
   CommandWorkflowSpec data. Derive command and option behavior from the live
   Click schema, never a copied command list. Preserve the existing CLI
   contract unless the task explicitly changes it.
3. Keep provider-specific request, tool-call, and continuation conversion in
   adapters. Normalize only observable decisions, tool calls, artifacts,
   approvals, and outcomes.
4. Put authorization, idempotency, budgets, and validation in deterministic
   code rather than a prompt. Require exact, one-shot approval for material
   execution or state changes.
5. Extend the current Runtime V2 contracts and event stream additively. Keep
   old event logs replayable and do not introduce a parallel orchestration
   system.
6. Use subagents only when a task is independently verifiable or genuinely
   parallel. Give each worker immutable inputs, a limited tool set, a budget,
   an expected artifact, and a deterministic merge check.
7. During implementation, use source inspection and deterministic receipts.
   Run one focused runtime, permission, registry, or CLI-schema suite only
   after a material milestone; allow one evidence-driven rerun. Reserve broad
   suites and read-only lint for the declared integration freeze.

## Command-compiled boundary

The model may return a JSON CommandWorkflowSpec or a bounded repair proposal.
It must not return a shell command for execution, native Gaussian/ORCA/xTB
input text, a model-written PySCF script, arbitrary filesystem paths, shell operators, redirections,
environment assignments, option aliases, or quoting decisions.

The deterministic compiler owns this sequence:

1. validate the immutable command DAG and budget;
2. resolve the live Click path and option scope;
3. resolve trusted project and ArtifactBinding identifiers and hashes;
4. render canonical long-flag argv and a display string;
5. run the safe fake/test CLI preview;
6. obtain an independent parser observation and compare semantic intent;
7. persist a CommandPreflightReceipt with schema, project, input, environment,
   preview-artifact, and counterexample references.

Expose only workspace/project operations plus synthesize, repair, inspect, and
explain command-workflow tools to the frontier profile. Keep legacy
molecule/settings/job/input/execution builders in an explicit
`harness_jobs` compatibility profile only; absent tools must fail closed.
Raw direct-string synthesis and compact-v8 conversion are baseline/migration
inputs, not a Frontier Runtime V2 model surface or a way around typed
compilation.
Limit repair to two structured counterexamples and reject any repair that
changes an explicit program, geometry, charge, multiplicity, method, or
constraint.

## Required boundaries

- Do not make a model assertion, a valid command, or a successful tool call a
  scientific pass condition.
- Do not make a direct-string baseline, a compact-v8 compatibility adapter, or
  a legacy job-builder fallback the command authority. They are migration or
  evaluation inputs only.
- Do not persist hidden reasoning. Preserve provider continuation state only as
  opaque protocol state, never as evidence.
- Do not let a planner, worker, or critic approve its own high-risk action.
- Do not enable autonomous execution, dynamic delegation, or a new provider
  protocol without a frozen single-agent baseline and an explicit evaluation.

## Use the references

- Read [runtime-contract.md](references/runtime-contract.md) before adding
  contracts, command workflow payloads, events, task graphs, provider
  capability fields, or replay logic.
- Read [approval-and-evaluation.md](references/approval-and-evaluation.md)
  before changing permissions, dispatch, budgets, or benchmark gates.

## Examples

Use this skill for: “compile a typed ORCA command workflow through the live
schema,” “add an approval-bound task-graph event,” “audit whether a provider
adapter leaks reasoning state,” or “design a bounded subagent experiment.”

Do not use this skill to let a model hand-author a .com, .gjf, or .inp file,
to select a quantum-chemistry method without scientific review, to validate a
frequency calculation, or to publish a result. Combine it with
`chemsmart-scientific-workflow` or `chemsmart-evidence-audit` when appropriate.
