# Integration maintenance index

This directory records the maintained ChemSmart 3.1.4 integration contract.
Documents describe authority and intended behavior; machine-readable receipts
and deterministic validation remain the evidence for an observed run.

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
   — active-path H0/H1 execution and analysis protocol.
7. [`release-source-custody-runbook.md`](release-source-custody-runbook.md) —
   deferred validation, exact commit series, backup, and public main
   replacement procedure.
8. [`../design/chemsmart-progressive-assurance.md`](../design/chemsmart-progressive-assurance.md)
   — permissive scientific planning with progressively stricter materialization,
   preview, and execution evidence.
9. [`../evaluation/progressive-planning-generality-protocol.md`](../evaluation/progressive-planning-generality-protocol.md)
   — the positive-first, chemistry-diverse primary model-development study.
10. [`progressive-planning-observations-v1.json`](progressive-planning-observations-v1.json)
    — narrow source, semantic-draft, and campaign-fixture observations made
    before deferred integration QA.

## Architecture and experiments

See `docs/design/` for source-level architecture and `docs/evaluation/` for
provider, deterministic-validation, and conditional engine protocols. User
commands and project schemas are documented in
`docs/source/pyscf-cli-options.rst` and `docs/source/xtb-cli-options.rst`.

## Evidence discipline

- Source presence is not test evidence.
- A focused test result is not full integration readiness.
- A provider response is not scientific validation.
- A process exit code is not a validated calculation.
- A historical GPU archive is not evidence for the current build until its
  manifest and requested/applied provenance are revalidated.
- Release readiness requires the explicit freeze gates and remote-reference
  verification in the release runbook.
