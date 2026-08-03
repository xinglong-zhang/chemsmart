# Integration maintenance index

This directory records the maintained ChemSmart 3.1.4 integration contract.
Documents describe authority and intended behavior; machine-readable receipts
and deterministic validation remain the evidence for an observed run.

The organizing principle is that ChemSmart is the canonical multi-program CLI
and project-YAML hub. Agent features extend that hub; they do not introduce a
second program configuration or execution authority.

## Start here

1. [`integration-status-v314.md`](integration-status-v314.md) — current
   milestone state, implementation map, owners, and exit gates.
2. [`v314-integration.md`](v314-integration.md) — upstream lineage, custody
   policy, and non-negotiable behavior.
3. [`program-management.md`](program-management.md) — capability,
   environment, substitution, and readiness authority.
4. [`migration-ledger-v314.json`](migration-ledger-v314.json) —
   machine-readable source dispositions.
5. [`deepseek-v4-flash-validation-v1.json`](deepseek-v4-flash-validation-v1.json)
   — preregistered paired model cases.
6. [`../evaluation/deepseek-v4-flash-h0-h1-protocol.md`](../evaluation/deepseek-v4-flash-h0-h1-protocol.md)
   — retained H0/H1 evaluation protocol; it is not evidence that every
   preregistered case has run.
7. [`release-source-custody-runbook.md`](release-source-custody-runbook.md) —
   bounded validation status, source custody, and public-main replacement
   procedure.
8. [`../design/chemsmart-progressive-assurance.md`](../design/chemsmart-progressive-assurance.md)
   — permissive scientific planning with progressively stricter materialization,
   preview, and execution evidence.
9. [`../evaluation/progressive-planning-generality-protocol.md`](../evaluation/progressive-planning-generality-protocol.md)
   — the positive-first, chemistry-diverse primary model-development study.
10. [`progressive-planning-observations-v1.json`](progressive-planning-observations-v1.json)
    — historical narrow source, semantic-draft, and campaign-fixture
    observations; it is not the current release status.

## Architecture and experiments

See `docs/design/` for source-level architecture and `docs/evaluation/` for
provider, deterministic-validation, and conditional engine protocols. User
commands and project schemas are documented in
`docs/source/pyscf-cli-options.rst` and `docs/source/xtb-cli-options.rst`.
The natural-language entrypoints are `chemsmart agent plan` for planning and
safe preview, and `chemsmart agent run` for the same path plus explicitly
approved host execution.

## Skill ownership

Maintained ChemSmart skills live in the user's global agent-skill registry.
Project-local skill packages are intentionally absent from this public source
tree. The checked-out CLI, project loaders, tests, and generated schema remain
authoritative if a globally installed skill becomes stale.

## Evidence discipline

- Source presence is not test evidence.
- A focused test result is not full integration readiness.
- A provider response is not scientific validation.
- A process exit code is not a validated calculation.
- A historical GPU archive is not evidence for the current build until its
  manifest and requested/applied provenance are revalidated.
- Publication to the fork requires the remote-reference and recoverability
  checks in the release runbook; it does not imply that the full suite ran.
