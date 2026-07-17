# ChemSmart Agent Provenance Manual Review

This review interprets the automated provenance evidence in
`agent-provenance-audit.json`. It is not a legal opinion and does not claim that
an automated comparison can prove independent authorship.

## Evidence Reviewed

- 151 production Python files under `chemsmart/agent`.
- 44,691 blamed production lines at audit HEAD `3b08f476`.
- 176 tracked ChemSmart agent code and documentation files.
- 14,071 tracked text/code documents from five pinned reference snapshots.
- Per-file `git log --follow` history from boundary commit `7f1cb674`.
- Exact 6-line, exact 50-token, and near 100-token similarity results.

All 151 production files were introduced at or after the agent-history
boundary. Their introduction commits are attributed to the same contributor
under three Git author spellings: `Hong jiseung`, `홍지승`, and `hjs0603`.
The line-level blame totals are 10,784, 33,406, and 501 respectively. This is
lineage evidence, not evidence about how an author formed an implementation.

## Documented Influences

### OpenAI Codex

The design lineage is explicit and precedes the first TUI implementation:

1. Codex UI research was committed as `6a55817f` at 2026-05-08 23:02 +09:00.
2. The ChemSmart TUI design decision followed as `aa5051da` at 23:20.
3. The first TUI shell followed as `8328b73a` at 23:38.

`docs/research/ui_codex.md` records the historical Codex source snapshot
`cce059467a...` and cites individual source locations. The design document then
adapts interaction concepts such as inline scrollback, artifact-specific cells,
approval overlays, and cwd-aware resume to chemistry workflows. The automated
comparison included both that historical snapshot and a later Codex snapshot;
it found no production or documentation match at the configured thresholds.

Classification: documented conceptual and interaction-design influence. No
threshold-level source reuse was identified.

### Claude Code

Claude Code UI research was committed as `331109ab` before the TUI design and
implementation. The note states that Claude Code is useful for interaction
patterns and state-model ideas but is not a reusable TUI codebase. It also
states that the public value for ChemSmart is interaction design rather than
renderer code.

The source snapshot is treated conservatively because its license file states
that rights are reserved and use is subject to Anthropic terms. No exact or
near match reached the audit thresholds.

Classification: documented product-behavior research only. Any future matching
implementation block would require behavior characterization and clean-room
replacement before release.

### Textual

Textual is a declared runtime dependency in `pyproject.toml` and the TUI imports
its public `App`, `Screen`, widget, event, binding, reactive, timer, and worker
APIs. Those imports are normal dependency use and were not treated as copied
implementation. The Textual source snapshot is MIT-licensed evidence, and the
automated comparison found no threshold-level source match.

Classification: direct library dependency through public APIs, with no detected
source incorporation. Dependency packaging must continue to preserve the
library's own license metadata.

### oh-my-openagent

`docs/agent-v13_1-local-adapter.md` explicitly says that harness-engineering
patterns from the oh-my-openagent/oh-my-codex investigation informed compact
context, deterministic adapters, real-parser validation, and artifact evidence.
These are architectural contracts specialized to the ChemSmart CLI and
calculation inputs. The pinned repository uses a Sustainable Use License; no
exact or near source match reached the audit thresholds.

Classification: documented architectural influence. Direct reuse from this
source is not accepted without a separate restrictive-license review.

## Automated Leads and Manual Disposition

The final normalized scan produced:

| Match class | Leads | Manual disposition |
|---|---:|---|
| Exact nonblank block, at least 6 lines and 50 lexical tokens | 0 | None to review |
| Exact lexical-token sequence, at least 50 tokens | 0 | None to review |
| Near token unit, at least 100 tokens and similarity at least 0.85 | 0 | None to review |

An initial detector version counted long punctuation-only comment separators as
token matches. Manual inspection showed that these were false positives. The
detector was corrected to count identifiers, numbers, and meaningful operators,
and synthetic exact/near-match tests were added before the final scan.

## Decision

- Preserve the existing design citations and research notes.
- Do not add `THIRD_PARTY_NOTICES.md` on the basis of this audit because no
  incorporated source block was confirmed.
- Keep Textual's license metadata through normal dependency distribution.
- Treat Claude Code and oh-my-openagent as non-copy sources for implementation.
- Re-run the pinned comparison after the structural refactor because moving code
  must not obscure provenance evidence.

The remaining risk is epistemic: thresholded similarity and Git attribution
cannot rule out lower-similarity adaptation or unrecorded sources. The report
therefore supports an engineering review, not a legal conclusion.
