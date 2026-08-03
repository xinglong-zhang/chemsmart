# ChemSmart 3.1.4 integration and main-replacement runbook

## Purpose and current state

This runbook controls source custody, integration, deferred validation,
reviewable commits, and replacement of the public fork's `main` branch. It
does not authorize publication, a pull request, an upstream push, or deletion
of archived worktrees or remote branches.

Fixed release facts:

- upstream base: `5d486b34501b9c040ad1db30d41a0eaa3850b15f`;
- package version: `3.1.4`;
- integration branch: `integration/chemsmart-v314-pyscf-agent`;
- public remote: `https://github.com/Hongjiseung-ROK/chemsmart.git`;
- release identity: `Hongjiseung-ROK <206387986+Hongjiseung-ROK@users.noreply.github.com>`;
- target ref: `refs/heads/main` on the public fork only.

Current evidence as of 2026-08-03:

| Item | State | Evidence boundary |
| --- | --- | --- |
| Exact upstream base and neutral integration worktree | observed | Current branch resolves to the pinned upstream SHA before integration commits. |
| Private custody archive for target, PySCF source, and two frontier worktrees | observed | Git bundles, tracked patches, non-secret path-safe archives, manifests, exclusions, identities, and restoration records exist outside the public tree. |
| PySCF/GPU4PySCF, bounded xTB, Runtime V2, capability plane, and provider-campaign source | observed | Implemented files are present as an uncommitted integration tree. |
| Initial focused validation | supported but not release-green | It exposed defects; bounded source repairs and narrow probes followed. The single evidence-driven focused rerun is deliberately deferred until implementation is complete. |
| Live DeepSeek campaign | planned | No scored official-endpoint episode has been executed from this tree. |
| CPU H2O SP/OPT/HESS execution | waiting for approval | A coordinate artifact is known. SP/OPT can later receive an exact first approval; HESS cannot receive an exact input-bound approval until OPT produces and validates its geometry artifact. |
| Full QA and public `main` replacement | planned | No release gate, commit series, or push has completed. |

Do not reinterpret source presence or a narrow probe as integration, scientific,
or release readiness.

## M0 — Source custody

### Required contents

For every contributing checkout, preserve before integration:

1. repository identity, branch, HEAD, remotes, and porcelain-v2 status;
2. committed refs in a verified Git bundle;
3. tracked changes in a binary-capable patch;
4. non-secret untracked and required ignored evidence in a path-safe archive;
5. file modes, sizes, SHA-256 digests, explicit exclusions, and deleted paths;
6. a migration ledger disposition of `ported`, `adapted`, `archive_only`, or
   `rejected` for each in-scope frontier path.

Exclude credentials, the ignored API secret file, authorization/request
headers, caches, environments, and provider-private reasoning. Keep the two
frontier worktrees in place after archiving.

### Restoration gate

Restore each bundle plus patch plus path-safe archive into a disposable
directory. Verify HEAD, dirty-state inventory, and archive hashes against the
manifest. Record the drill outcome without embedding private absolute paths in
the public repository. A manifest alone is not custody proof.

M0 is currently **observed complete**, but its private manifest and restoration
records must be reverified at M4 because the old public `main` SHA will be
added to a release backup bundle immediately before replacement.

## M1 — Program integration

M1 remains an implementation milestone until the deferred focused suite is
green.

### PySCF 2.14.0 and GPU4PySCF 1.8.0

- Keep `run/sub pyscf sp|opt|hess` settings reconstruction equivalent.
- Generate standalone scripts deterministically; never expose script authorship
  or editing to a model.
- Require settings validation, scientific preflight, exact compute-interpreter
  evidence, and post-result provenance verification.
- Propagate generated child exit status through local and scheduler wrappers.
- Bind input, project, requested/applied settings, script, interpreter,
  environment, and HDF5 result hashes.
- Record GPU identity and dependency-stack evidence when GPU is requested.
- Never fall back from GPU to CPU.
- Convert missing requested properties into explicit findings.
- Materialize titles and reject inherited no-op settings.

GPU support remains limited to archived SP evidence unless a separate future
execution is approved. It is not part of this release's new execution gate.

### xTB 6.7.1

- Port only `sp|opt|hess`, `gfn0|gfn1|gfn2|gfnff`, CPU execution, strict
  optional project YAML, and supported optimization levels.
- Preserve the v3 molecule and output I/O implementation.
- Reject unknown keys, contradictory job types, incomplete solvent pairs,
  arbitrary control text, stale mixed outputs, MD/path workflows, and
  unsupported constraints.
- Bind source, project, actual execution input, output family, and validation
  receipts before a result can be accepted.

### Capability plane

Maintain one immutable declaration registry for Gaussian, ORCA, xTB, PySCF,
and NCIPLOT. Join it with the live Click schema and record the digest. A
registry row proves neither installation nor execution readiness.

No broad validation is run during M1. Source inspection and narrow deterministic
probes may guide implementation; the milestone receives one focused suite only
after M1 and M2 are substantially complete.

## M2 — Capability-driven agent

Extend Runtime V2 rather than introducing a parallel runtime. Preserve replay
of legacy streams with empty additive program-management state.

Required public contracts cover:

- declared versus observed program support;
- live capability and environment queries;
- model-proposed candidates separated from host decisions;
- explicit requested/selected program and engine bindings;
- host-owned substitution findings and approval requirement;
- project, geometry, electronic state, environment, validator, and preflight
  bindings;
- replayable capability, environment, substitution, preflight, and result
  verification events.

The model-visible frontier surface is limited to project operations, command
synthesis/repair/preview/inspection, calculation inspection, and the narrow
capability/environment/candidate/preflight tools. Models cannot supply paths,
shell syntax, native inputs, generated Python, executable status, approvals,
readiness, or terminal state.

Only registered Gaussian-to-PySCF equivalence classes for SP, optimization,
fixed-geometry Hessian, and optimization plus frequencies may produce a
candidate. Explicit substitution always waits for user approval. TS, IRC, TD,
scan, QMMM/ONIOM, NEB, post-HF, double hybrids, mixed basis/ECP, and unsupported
constraints block.

Knowledge packs remain advisory. The live schema, strict project loaders,
environment probes, registries, deterministic validators, and artifact hashes
are authoritative.

M2 source is **observed present** but its integrated behavior and legacy replay
are **unverified** until M3/M4. Semantic workflow planning, the arm-neutral
shared-context bridge, a distinct `planned` Runtime terminal, and
`schedule_completed` campaign termination exist in source. The broader
evidence-aware planning draft, alternative-branch, and node-local
materialization contracts remain specified work. Current source presence does
not replace the deferred focused observation.

## M3 — Deferred focused validation and real experiments

### M3.1 Focused integration gate

After implementation stabilizes, run one focused suite covering:

- PySCF and xTB CLI/settings and `run/sub` reconstruction;
- generated-script child status and HDF5 provenance;
- capability/live-schema binding and coverage-scoped conformance;
- Runtime V2 event replay, permissions, and project tools;
- substitution, preflight, preview, result-verification, and provider
  continuation contracts.

One evidence-driven rerun is allowed. Do not run formatter, autofix, snapshot
regeneration, broad lint, or the full suite here.

### M3.2 DeepSeek H0/H1 campaign

Execute the preregistered fourteen-case paired protocol in
[`../evaluation/deepseek-v4-flash-h0-h1-protocol.md`](../evaluation/deepseek-v4-flash-h0-h1-protocol.md).
The active path, thinking continuation, adaptive attempts, redaction,
deterministic oracles, and stop rules in that document are release gates. A
direct provider response is not equivalent to an active-path result.

If a safety case fails, keep the capability-driven model profile disabled by
default. Persuasive text cannot override a red host gate.

### M3.3 CPU H2O approval boundary

The currently observed coordinate candidate has SHA-256
`42f720e0b1ae9883ad99e814488bf46093068bb386f007358841562629957045`.
Re-observe it before the execution proposal; a changed hash invalidates the
candidate.

Planning and safe preview may prepare this exact DAG:

1. CPU PySCF SP on the initial neutral-singlet geometry;
2. CPU PySCF OPT on the same initial geometry;
3. CPU PySCF HESS bound to the exact optimized-geometry artifact and hash from
   node 2.

All nodes use gas-phase B3LYP/def2-SVP, four cores, 4 GB, and no scratch. Before any engine
call, stop in `waiting for approval` and present an immutable workflow-level
approval manifest containing:

- the exact SP and OPT canonical argv/display commands and command hashes;
- input/project/schema/capability/interpreter/environment hashes;
- requested and applied settings digests;
- initial coordinate hash and the typed rule for binding the future optimized
  geometry into HESS;
- working directory, wall-time/resource bounds, approval expiry, no-retry
  policy, and expected artifacts;
- validators and terminal conditions.

The user must approve that workflow manifest. HESS remains `planned`; before
OPT execution there is no optimized input hash or exact HESS invocation. The
approval therefore binds HESS semantically to the validated geometry selected
from the exact approved OPT node, together with its settings, resources,
validators, and command family. After OPT, the host validates the HDF5 and
geometry, resolves the producer edge, and compiles and previews HESS. It may
continue under the same approval only if every bound semantic field and
producer rule remains unchanged. Otherwise it pauses for new approval. Do not
retry automatically or substitute another engine.

Execution receipts must then show SCF convergence, optimization convergence,
valid structured HDF5, requested/applied provenance equality, and an exact
optimized-geometry handoff. The stationary point is `validated` only if the
nonlinear water Hessian yields three finite vibrational modes and no imaginary
mode under the recorded convention. Successful process exit alone is not
scientific validation.

One bounded xTB GFN2 water SP was authorized and completed for integration.
No GPU, Gaussian, ORCA, scheduler, or HPC calculation was performed, and no
further engine call is authorized by this runbook.

## M4 — Freeze, commit series, and publication to the fork

### Final gates

At a separately authorized release checkpoint, run exactly once after the
experiment tree is frozen:

1. full test suite;
2. read-only Ruff;
3. live-schema/run-sub parity and legacy event replay checks;
4. credential, private absolute-path, forbidden-token, and generated-cache
   scans;
5. archive restoration and bundle verification;
6. migration-ledger hash regeneration and completeness check;
7. `git diff --check`;
8. public tree and release-range commit-metadata sanitization.

Do not autofix, format, or regenerate snapshots. A failure requires an
evidence-linked correction and an explicit decision about whether the affected
gate may be rerun; do not hide the original failure.

Release requires:

- exact upstream v3.1.4 behavior plus registered PySCF and xTB leaves;
- 100% `run/sub` schema parity;
- successful CPU H2O SP/OPT/HESS receipts and scientific validators;
- zero native-input/shell bypass, false-ready state, artifact substitution,
  secret leakage, silent fallback, or success while a required gate is red;
- a clean public-tree scan for the case-insensitive token supplied by the
  release owner outside the repository;
- all public commits authored and committed by the required release identity.

### Reviewable commit series

Set both author and committer identity to the release identity, then create
only these three neutral Conventional Commit subjects, in order:

1. `feat(programs): integrate PySCF GPU and xTB workflows`
2. `feat(agent): add capability-driven execution workflows`
3. `test: validate integrated program and agent workflows`

Stage explicit path groups for each commit. Exclude credentials, caches,
private receipts, hidden reasoning, raw provider responses, private archive
paths, and bytecode. After every commit, inspect the staged diff summary and
the resulting commit metadata; do not rewrite upstream history.

Use these ownership boundaries when staging:

| Commit | Primary path ownership |
| --- | --- |
| 1 | PySCF and xTB CLI/jobs/I/O/settings/templates, shared program registry, and program-facing documentation |
| 2 | provider-neutral agent contracts, capability/command/project/workflow/runtime/tools, DeepSeek continuation, project-local skills, and architecture documentation |
| 3 | all focused/integration fixtures and tests, including shared registration assertions |

If one file contains inseparable changes for two commits, split its hunks or
assign it to the earliest commit that must compile with it; do not duplicate or
temporarily break the public series merely to satisfy directory ownership.

### Exact SHA-bound `main` replacement

Use a clean release shell and explicit 40-character SHAs. The values below are
placeholders that must be copied from fresh command output; do not derive them
from stale notes.

```sh
git remote get-url fork
gh api user --jq .login
git fetch --no-tags fork refs/heads/main:refs/remotes/fork/main
git rev-parse refs/remotes/fork/main
git rev-parse integration/chemsmart-v314-pyscf-agent
```

Record the first SHA as `OLD_FORK_MAIN` and the second as `RELEASE_SHA` in a
private release receipt. Confirm the remote URL and authenticated account are
exactly the intended fork owner. Do not fetch or push the upstream remote.

Create and verify a recoverable bundle containing the freshly fetched old
public ref:

```sh
git bundle create "$PRIVATE_BACKUP_PATH/fork-main-before-replacement.bundle" refs/remotes/fork/main
git bundle verify "$PRIVATE_BACKUP_PATH/fork-main-before-replacement.bundle"
```

Append the bundle hash, `OLD_FORK_MAIN`, timestamp, and verification outcome to
the private custody manifest. Verify `RELEASE_SHA` is the tip of the
three-commit series and that its ancestry contains the pinned upstream base.

Supply the release owner's forbidden case-insensitive token only through an
environment variable outside the public tree, then require both scans to
return no matches:

```sh
test -n "$PUBLIC_FORBIDDEN_TOKEN"
git grep -I -n -i -e "$PUBLIC_FORBIDDEN_TOKEN" "$RELEASE_SHA" -- .
git log --format=fuller "$RELEASE_SHA" | LC_ALL=C grep -i -- "$PUBLIC_FORBIDDEN_TOKEN"
```

For these two negative scans, exit status 1 means no match; any printed match
or other exit status blocks release. The history scan covers every commit
reachable from the new public main, not only the three-release-commit range.
Also scan the tree for credentials and private absolute paths using the release
manifest's private token list.

Perform the dry-run with the exact old SHA as the lease:

```sh
git push --dry-run fork "$RELEASE_SHA:refs/heads/main" --force-with-lease="refs/heads/main:$OLD_FORK_MAIN"
```

Fetch `fork/main` again and require it still equals `OLD_FORK_MAIN`. Only then
perform the one allowed replacement:

```sh
git fetch --no-tags fork refs/heads/main:refs/remotes/fork/main
test "$(git rev-parse refs/remotes/fork/main)" = "$OLD_FORK_MAIN"
git push fork "$RELEASE_SHA:refs/heads/main" --force-with-lease="refs/heads/main:$OLD_FORK_MAIN"
```

Finally, fetch and require local/remote equality:

```sh
git fetch --no-tags fork refs/heads/main:refs/remotes/fork/main
test "$(git rev-parse refs/remotes/fork/main)" = "$RELEASE_SHA"
```

Record the remote response, local and remote SHAs, bundle hash, and gate
manifest hash without credentials or private paths. Do not open a pull request,
push another branch, push upstream, delete historical remote branches, or
remove archived worktrees.

### Recovery boundary

The verified old-main bundle and `OLD_FORK_MAIN` make the replacement
recoverable. Restoration of the old public ref is a separate consequential
operation and requires a new explicit SHA-bound decision; do not perform an
automatic rollback.

## Milestone report template

At M0 through M4, report separately:

1. what changed;
2. observable evidence obtained;
3. failures and causal classification;
4. remaining unknowns;
5. `retain`, `revise`, or `reject` decision;
6. next approval or decision required.

Never report `executed`, `validated`, `reproduced`, release-ready, or generally
capable beyond the exact receipts and gates that exist.
