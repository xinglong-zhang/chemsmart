---
name: chemsmart-evidence-audit
description: Audit ChemSmart command workflows, calculations, agent traces, claims, benchmarks, or reports for reproducibility and evidence coverage. Use when verifying citations, auditing CommandWorkflowSpec and preflight receipts, creating claim-evidence manifests, producing canonical research documents, conducting bounded independent critique, designing adversarial tests, or evaluating subagents, documentation, and review components.
---

# ChemSmart Evidence Audit

Use this skill to make a result inspectable and falsifiable. Read `AGENTS.md`
first. Treat a model-generated narrative as an interpretation to audit, not as
evidence.

## Audit workflow

1. Identify the claim type: observation, computed result, inference,
   literature statement, or unresolved uncertainty.
2. Bind every command claim to its typed IR, schema digest, trusted
   project/artifact bindings, canonical argv, safe-preview receipt, parser
   observation, and intent comparison. Bind every numerical statement to a
   stable artifact, parsed value, unit, validator, and calculation environment.
   Mark missing evidence rather than filling it from prose.
3. Verify scholarly metadata programmatically. Keep DOI/arXiv source, title,
   authors, venue, year, correction/retraction check, and retrieval timestamp
   outside the BibTeX file.
4. Generate reports and notebooks from the evidence manifest. Preserve
   ChemSmart-generated native inputs/outputs and structured records so
   rendering remains reproducible. Do not use a model-written native input as
   evidence or silently repair a command/input while auditing it.
5. Give an independent, read-only critic the artifacts and declared
   assumptions. Require it to cite a check, receipt, or source for every
   finding.
6. Resolve findings with deterministic validation, independent recomputation,
   or a human decision. Do not use majority vote as scientific arbitration.
7. End with supported, qualified, rejected, or unresolved claims and separate
   readiness from raw task completion.

## Red-team rules

Seed realistic failures: wrong charge/multiplicity, malformed geometry,
incompatible settings, nonconvergence, imaginary frequencies, truncated
output, bad units, contradictory literature, unavailable executable, and
approval mismatch. Bound review rounds and retain a blocked outcome.

Also seed command-compiler failures: shell injection, invented or
out-of-scope options, wrong option scope/order, stale schema/project/artifact
binding, charge/multiplicity drift, lost constraints, unsafe direct input
fallback, parser disagreement, and reused approval after a semantic change.
The critic may identify these defects but cannot repair, approve, execute, or
mark the command ready.

## Use the references

- Read [evidence-and-reports.md](references/evidence-and-reports.md) for the
  canonical manifest and reproducible-document requirements.
- Read [red-team-and-ablation.md](references/red-team-and-ablation.md) before
  designing a critique or component evaluation.

## Examples

Use this skill for: “audit a typed command preview receipt,” “audit every
claim in this calculation report,” “verify the BibTeX and provenance ledger,”
or “compare a subagent system to the single-agent baseline.”

Do not use it to execute a calculation, silently alter a method or native input
to pass a check, compile a model-written shell command, or grant approval for
an action.
