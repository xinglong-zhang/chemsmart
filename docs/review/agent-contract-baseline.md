# ChemSmart Agent Contract Baseline

This document freezes the behavior that must remain observable while the agent
codebase is decomposed. It is a characterization baseline, not a claim that the
current implementation is correct or production-ready.

## Scope

- Original pinned `fork/main`: `f5c7a3280dd86bf19d915430a87ee24c0b4783e4`
- Baseline stabilization commit: `95faa53c`
- Latest local-main UX integration: `ea2a368e8bf5c812494fc5e465164b621ced4497`
- Contract snapshot: `tests/agent/contracts/agent_contract_baseline.json`
- Contract SHA-256: `9a5e20c313dac01dc4e53d5c19f49602371781a2e9cbcfd357659ffa1a6bc95b`
- HighRisk48 fixture SHA-256: `88df9955a0e8fc587a074969ce356ba0bb1fd2c0bd7859a18d166b93b05a4dd1`

The snapshot covers 30 CLI help paths, slash commands and keybindings, 30 tool
schemas, public imports, prompt sources and rendered prompt examples, workspace
project-YAML resolution, session artifacts, 166 semantic rule IDs, and the 26
canonical command kinds.

## Reproduction

Run the snapshot check in the isolated worktree and ChemSmart environment:

```bash
env -u PYTHONPATH conda run -n chemsmart python \
  scripts/review/snapshot_agent_contracts.py \
  --repo . \
  --check tests/agent/contracts/agent_contract_baseline.json
```

Reproduce the structural baseline with the frozen coverage result:

```bash
env -u PYTHONPATH conda run -n chemsmart python \
  scripts/review/audit_agent_architecture.py \
  --repo . \
  --coverage-json <COVERAGE_JSON> \
  --json docs/review/agent-architecture-baseline.json \
  --markdown docs/review/agent-maintainability-audit.md
```

The coverage input is a local P0 receipt rather than a portable source
artifact. A fresh coverage run must replace it when the final branch is
evaluated.

## Characterized Failures

All 26 compact SPEC fixtures adapt into commands. That does not prove the
commands preserve runtime or chemical intent. The current baseline records 18
semantic passes and 8 semantic rejects:

| Kind | Failed rule |
|---|---|
| `gaussian.irc` | `cmd.runtime.cli_value_error` |
| `gaussian.crest` | `cmd.runtime.cli_value_error` |
| `gaussian.traj` | `cmd.contract.traj_jobtype_required` |
| `gaussian.qrc` | `cmd.runtime.input_not_found` |
| `orca.ts` | `cmd.runtime.input_not_found` |
| `orca.irc` | `cmd.runtime.cli_value_error` |
| `orca.scan` | `cmd.contract.scan_required_parameters` |
| `orca.qrc` | `input.state.electron_multiplicity_parity` |

The intent layer records 21 passes and 5 rejects:

| Kind | Failed rule |
|---|---|
| `gaussian.freq` | `intent.kind` |
| `orca.ts` | `intent.chemistry.tssearch_type` |
| `orca.freq` | `intent.kind` |
| `orca.irc` | `intent.chemistry.inithess` |
| `orca.scan` | coordinate, step-count, and step-size intent rules |

These failures are frozen so refactoring cannot silently hide or reclassify
them. A later behavior fix requires a dedicated test, a semantic commit, and an
explicitly reviewed contract update.

## Baseline Test Receipt

- `tests/agent`: 1,048 passed in two isolated runs.
- Branch-aware coverage: 74.63%.
- Existing order-dependent TUI and ambient-HOME test failures were corrected by
  characterization tests before this snapshot was created.
- After integrating the latest local-main TUI commits, the frozen contract
  remained byte-identical and the complete `tests/agent/tui` slice passed in
  the isolated `chemsmart` Conda environment.

The final regression audit must run the suite in three deterministic order
permutations and replay Full26 and HighRisk48. This P0 receipt does not satisfy
those final gates by itself.
