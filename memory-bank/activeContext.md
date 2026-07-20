# Active Context (update me whenever focus or decisions change)

_Last updated: 2026-07-10 by Claude (initial handoff to Cline)._

## Current focus

1. **Grow trusted `reasoning_synthesis` rows: 38 → 150-300.**
   Method: `reasoning_accum.py --corpus command-hard`, 4-case batches,
   rotating live reasoning teachers (see 02-accumulation-workflow.md).
2. **Grow `wrong_route_contrast` pairs: 14 → ~100** by pairing every
   observed mis-route with a command-hard corrected run of the same request.
3. Keep diversity green (skeleton ratio 0.394 now; floor 0.35): vary
   phrasing structure, kinds (scan/modred/crest/td are under-represented as
   PASSES), charge/mult, and language (EN/KO/ZH) across batches.

## Active decisions (do not relitigate without new evidence)

- Trusted = `synthesis.reasoning` only; assistant reasoning → review branch.
- command-hard is the default corpus; workflow-hard/agentic-style only for
  trajectory variety after trusted volume is healthy.
- No template-loop data generation, ever (v13 lesson).
- compact_spec family stays excluded from the v17 reasoning mix (persona
  separation); it remains exported for the 3B planner lineage.
- Batch cap 4 cases; stop a teacher on first quota-403.

## Recently landed (context for code you will see)

- Reasoning capture end-to-end: synthesis `_last_reasoning` →
  tools_command payload `reasoning` → episode `synthesis.reasoning` →
  exporter trusted/review split.
- Audit SOTA metrics fixed: deterministic sha1 trajectory fingerprint
  (order-preserving; no empty-bigram collapse), schema-efficiency averaged
  over tool-bearing episodes only, multistep-only duplicate counting.
- pydantic tool-schema fix: legacy exec tools (`run_local`, `submit_hpc`,
  `extract_optimized_geometry`) are deprecated JobArgs stubs; their old
  contract tests are intentionally skipped.
- Full agent pytest suite green: 778 passed, 10 skipped (deprecated tools).

## Blocked / waiting

- v17 training (Part D) waits on trusted volume (~150 minimum).
- Teacher quota resets daily-ish; if all reasoning teachers 403, fall back
  to audit/export/documentation work and note it here.
