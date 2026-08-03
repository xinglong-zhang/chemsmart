# ChemSmart 3.1.4 integration and main-replacement runbook

## Purpose and current state

This runbook controls source custody, the committed integration candidate,
bounded validation, and replacement of the public fork's `main` branch. It
does not authorize publication, a pull request, an upstream push, or deletion
of archived worktrees or remote branches.

Fixed release facts:

- upstream base: `5d486b34501b9c040ad1db30d41a0eaa3850b15f`;
- package version: `3.1.4`;
- integration branch: `integration/chemsmart-v314-pyscf-agent`;
- public remote: `https://github.com/Hongjiseung-ROK/chemsmart.git`;
- release identity: `Hongjiseung-ROK <206387986+Hongjiseung-ROK@users.noreply.github.com>`;
- target ref: `refs/heads/main` on the public fork only.

Current evidence as of 2026-08-04:

| Item | State | Evidence boundary |
| --- | --- | --- |
| Exact upstream base and neutral integration worktree | observed | Current branch resolves to the pinned upstream SHA before integration commits. |
| Private custody archive for target, PySCF source, and two frontier worktrees | observed | Git bundles, tracked patches, non-secret path-safe archives, manifests, exclusions, identities, and restoration records exist outside the public tree. |
| PySCF/GPU4PySCF, bounded xTB, Runtime V2, capability plane, and provider integration | committed | Program, agent, and test layers are recorded as reviewable commits; project-local skill packages were removed separately. |
| Bounded validation | supported | The selected run recorded 277 passing tests before the final fixture correction. The corrected xTB provenance test and read-only Ruff over four corrected files subsequently passed. The full suite was not run. |
| Live DeepSeek integration | observed | Official-endpoint planning, approved CPU PySCF and xTB execution, and GPU-unavailable planning episodes traversed the public agent path. |
| CPU H2O SP/OPT/HESS execution | validated for the bounded case | The workflow completed with exact OPT-to-HESS geometry handoff, finite energy, convergence, three finite nonlinear-water modes, and no imaginary mode. |
| Public `main` replacement | pending remote verification | Publication requires a fresh old-ref observation, recoverable bundle, exact lease, and post-push SHA equality. |

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

M1 is implemented and committed. Its evidence remains bounded by the live
water workflows and focused checks described above.

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

M2 is implemented and committed. The public `agent plan` and `agent run`
entrypoints reached the provider session, typed tools, project promotion,
command preview, approval-bound host execution, validation, and geometry
handoff in bounded live episodes. This does not establish paper-level or
cross-domain generality.

## M3 — Bounded validation and live experiments

### M3.1 Focused validation

The selected focused group recorded 277 passing tests before its final xTB
fixture correction. The corrected provenance test subsequently passed as a
targeted check, and read-only Ruff passed on the four corrected Python files.
The complete focused group and full repository suite were not rerun.

### M3.2 Provider observations

Bounded official-endpoint episodes exercised planning, typed tool
continuation, project promotion, safe preview, host approval, CPU execution,
result validation, and honest GPU unavailability. The preregistered H0/H1
protocol remains an evaluation design; its full case matrix was not executed
and is not a publication gate for this bounded integration.

### M3.3 Approved calculation observations

One CPU PySCF water workflow completed
`SP(initial) -> OPT(initial) -> HESS(optimized)` with the Hessian bound to the
validated optimized-geometry artifact. The receipts record finite energy,
SCF and optimization convergence, requested/applied provenance agreement,
three finite nonlinear-water modes, and no imaginary mode. One CPU xTB GFN2
water single point also completed with normal termination and finite energy.

No GPU, Gaussian, ORCA, scheduler, or HPC calculation was performed. These
bounded observations do not establish broad scientific reproduction or agent
generality, and no further engine call is authorized by this runbook.

## M4 — Freeze, commit series, and publication to the fork

### Publication gates

This bounded publication checkpoint does not claim a full-suite release
qualification. The full repository suite remains explicitly unrun. Before
replacing the fork ref, require:

1. the corrected xTB provenance test and read-only Ruff over the four corrected
   files to remain green;
2. credential, private-path, forbidden-token, and generated-cache scans;
3. `git diff --check`, a clean release worktree, and sanitized public commit
   metadata;
4. successful bounded CPU PySCF and xTB receipts, with GPU unavailability
   producing no CPU fallback;
5. a verified recoverable bundle of the freshly fetched old fork ref;
6. an exact lease-bound dry run, followed by the push and a fresh fetch proving
   remote/local SHA equality.

Do not autofix, format, regenerate snapshots, push upstream, or imply that a
successful remote update proves broad scientific or full-suite readiness.

### Reviewable commit series

The integration is recorded in five neutral Conventional Commit subjects, in
order:

1. `feat(programs): integrate PySCF GPU and xTB workflows`
2. `feat(agent): add capability-driven execution workflows`
3. `test: validate integrated program and agent workflows`
4. `chore(skills): remove project-local skill packages`
5. `docs: align integration and publication guidance`

Stage explicit path groups for each commit. Exclude credentials, caches,
private receipts, hidden reasoning, raw provider responses, private archive
paths, and bytecode. After every commit, inspect the staged diff summary and
the resulting commit metadata; do not rewrite upstream history.

Use these ownership boundaries when staging:

| Commit | Primary path ownership |
| --- | --- |
| 1 | PySCF and xTB CLI/jobs/I/O/settings/templates, shared program registry, and program-facing documentation |
| 2 | provider-neutral agent contracts, capability/command/project/workflow/runtime/tools, provider continuation, and architecture documentation |
| 3 | all focused/integration fixtures and tests, including shared registration assertions |
| 4 | project-local skill-package deletion and global-skill guidance |
| 5 | public fundamentals, bounded evidence status, migration record, and publication procedure |

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
the private custody manifest. Verify `RELEASE_SHA` is the tip of the five-commit
series and that its ancestry contains the pinned upstream base.

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
