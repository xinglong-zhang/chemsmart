# ChemSmart 3.1.4 integration record

For the live milestone ledger, implementation ownership, and remaining gates,
see [`integration-status-v314.md`](integration-status-v314.md).

## Scope

This fork is based on upstream ChemSmart revision
`5d486b34501b9c040ad1db30d41a0eaa3850b15f`. The package version remains
`3.1.4`; the exact Git revision and the joined live-CLI-schema/capability digest
identify a fork build.

The integration adds three maintained layers:

1. PySCF 2.14.0 as a first-class `sp|opt|hess` program, with GPU4PySCF 1.8.0
   as a distinct PySCF execution engine.
2. A bounded xTB 6.7.1 execution plane while retaining the upstream v3 xTB
   parser and molecule I/O.
3. A provider-neutral agent runtime in which models propose typed intent and
   deterministic ChemSmart code owns commands, environment evidence,
   approvals, validation, and terminal state.

## Source-custody policy

The integration was prepared from independently preserved source snapshots.
Committed refs, tracked modifications, and non-secret untracked evidence were
archived with SHA-256 manifests and restoration drills before any port began.
Credentials, caches, environments, request headers, and provider-private
reasoning were excluded. The private custody archive is intentionally not part
of this public tree.

Migration decisions are classified as:

- `ported`: copied with semantics retained;
- `adapted`: rewritten for the v3.1.4 contracts;
- `archive_only`: retained as historical evidence but absent from the product;
- `rejected`: deliberately excluded because it violated the new contract.

The machine-readable summary is in
[`migration-ledger-v314.json`](migration-ledger-v314.json).

User-facing project and command contracts are documented in
[`../source/pyscf-cli-options.rst`](../source/pyscf-cli-options.rst) and
[`../source/xtb-cli-options.rst`](../source/xtb-cli-options.rst).

## Non-negotiable behavior

- `run` and `sub` share the same program groups and reconstruct equivalent
  program settings.
- A declared capability never implies that an executable, library, GPU stack,
  or scientific method is available.
- The actual compute interpreter and dependencies are observed before PySCF
  readiness; the controller environment is not a substitute.
- GPU4PySCF never falls back silently to CPU.
- PySCF preflight and post-result provenance verification are mandatory.
- A scheduler wrapper propagates the generated child process exit status.
- xTB rejects unknown settings, incomplete solvent pairs, contradictory job
  types, stale mixed outputs, arbitrary control text, and unsupported workflow
  families.
- Models cannot author shell syntax, native-engine input, generated PySCF
  scripts, readiness, approvals, or terminal outcomes.
- Program substitution is host assessed, preserves requested and selected
  identities, and requires explicit hash-bound approval.

## Evidence grades

Keep these claims distinct. In particular, do not confuse the source enum
`EnvironmentStatus.AVAILABLE` with `SupportLevel.AVAILABLE`:

- `declared`: present in the immutable capability registry;
- `schema-bound`: present in the observed Click command tree;
- `configured`: a server/project configuration was loaded;
- `environment-observed`: the exact executable or interpreter and dependencies
  were observed; this does not establish agent support;
- `preview_only`: current compiler, preview, preflight, and verifier
  conformance receipts exist, but execution evidence does not;
- `previewed`: one exact invocation generated the expected downstream artifacts without
  running chemistry;
- `executed`: the engine ran and emitted an execution receipt;
- `operable`: the support overlay has both current conformance and execution
  evidence (`SupportLevel.AVAILABLE`);
- `validated`: all required deterministic and scientific checks passed.

No lower grade entails a higher grade.

## Maintenance sequence

When adding a program or engine:

1. Add the narrow immutable capability declaration.
2. Register the same Click surface under both `run` and `sub`.
3. Add strict project settings and generated-artifact ownership.
4. Add environment observation, preflight, execution, and result verification.
5. Extend the agent overlay without broadening the core declaration.
6. Prove schema agreement, replay compatibility, and fail-closed negative cases.
7. Add real execution evidence only under an exact artifact- and
   environment-bound approval.
