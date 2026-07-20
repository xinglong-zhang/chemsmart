# ChemSmart Agent Provenance Audit

This is an evidence-based source-similarity review, not a legal opinion. A match is a review lead; absence of a match is not proof of independent authorship.

- Audit HEAD: `3b08f4766259dd577ad8f5414186701c61a46a26`
- Agent-history boundary: `7f1cb674bbaeb524a78a1f3939b33d472ea2e94d`
- Production Python files traced: 151
- Candidate code/document files compared: 176

## Pinned Reference Snapshots

| Source | Commit | License evidence | Review posture |
|---|---|---|---|
| openai-codex | `800715d201651a2a07c2706dca10400109dae3d3` | `LICENSE` / `d17f227e4df5...` | Apache-2.0 evidence; preserve applicable notices for confirmed reuse |
| openai-codex-research | `cce059467af64f05d0a1521344847d2b558b6a80` | `LICENSE` / `d17f227e4df5...` | Historical Apache-2.0 snapshot cited by the TUI research note |
| textual | `06dbeef4bb70fb718236aa418ed658ef4667a126` | `LICENSE` / `94f290a76237...` | MIT evidence; preserve copyright and permission notice for substantial reuse |
| claude-code | `c39cb0f14bfe8bb519bae5bfc55add6867c5e2ab` | `LICENSE.md` / `728158fd1037...` | All-rights-reserved evidence; matching code requires clean-room review |
| oh-my-openagent | `5ef852a32c2c433386eb009bd92ca7c07359d0e6` | `LICENSE.md` / `b61ac928f152...` | Sustainable Use License evidence; matching code requires restrictive-license review |

## Automated Similarity Leads

| Match class | Count | Threshold |
|---|---:|---|
| Exact nonblank lines | 0 | >=6 lines and >=50 tokens |
| Exact token sequence | 0 | >=50 normalized tokens |
| Near token unit | 0 | >=100 tokens and ratio >=0.85 |

No automated lead reached the configured thresholds. This does not replace manual review of history, comments, design citations, or lower-similarity adaptations.

## Lineage Summary

- Files whose first recorded commit predates the boundary: 0
- Introduction authors: 홍지승 (113), Hong jiseung (38)
- Full per-file introduction, latest-change, and blame distributions are retained in `agent-provenance-audit.json`.

## Interpretation Rules

- Conceptual similarity is recorded as design influence, not copied code.
- Compatible reuse still requires the applicable notices and attribution.
- Restrictive, unknown, or unattributed matching blocks require behavior tests followed by clean-room replacement before release.
- No legal conclusion is made by this report.
