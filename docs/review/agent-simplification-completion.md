# ChemSmart Agent Simplification — Completion Report

This report records the final measured state of the `agent-codebase-simplification`
campaign against the frozen pre-refactor baseline in
`agent-maintainability-audit.md` / `agent-architecture-baseline.json`
(baseline HEAD `95faa53c`). The baseline documents stay frozen by design;
this file is the completion evidence.

- Measured HEAD: `622ccbcd` (code-final commit; this report is docs-only on top)
- Measurement tool: `scripts/review/audit_agent_architecture.py`
- Date: 2026-07-17

## Gate results

| Gate | Baseline (95faa53c) | Final | Target | Status |
|---|---:|---:|---:|---|
| Files over 800 lines | 14 | 0 | 0 | met |
| Classes over 500 lines | 4 | 0 | 0 | met |
| Functions over 100 lines | 39 | 0 | 0 | met |
| Import cycles | 1 | 0 | 0 | met |
| Pylint duplicate groups | 8 | 0 | 0 | met |
| Ruff default violations | 0 | 0 | 0 | met |
| Ruff C901 violations | 65 | 27 | see waiver | waived at <15 |
| Physical lines (agent) | 44,691 | 49,782 | no increase | waived |

Structural context: python files 151 → 225, functions 1,694 → 2,157,
classes 202 → 253.

## Waivers (approved 2026-07-17)

1. **C901 residual (27 findings, all complexity ≤ 14).** The completion
   scope was revised to eliminate every finding with complexity ≥ 15
   (12 functions, all resolved). The remaining 27 findings are enumerated
   below and are individually small (11–14); they are left as follow-up
   candidates rather than blockers.
2. **LOC increase (+5,091 over baseline, +11.4%).** The no-LOC-increase
   gate conflicts with decomposing 4,600-line screens and 2,600-line
   sessions into single-responsibility modules; the increase is the
   boilerplate cost of 74 additional modules plus one merged TUI-UX
   feature branch (+798 lines at merge `ea2a368e`). Accepted as a
   deliberate trade for the structural gates above.

## Residual C901 findings (all < 15)

| Function | File | Complexity |
|---|---|---:|
| `resolve` | `permissions.py` | 14 |
| `parse_cclib_output` | `runtime/cclib_parser.py` | 14 |
| `_route_opts` | `postprocess_v8.py` | 13 |
| `stream_decision_event` | `services/cli_presenters.py` | 13 |
| `_public_synthesis_trace` | `tools_command.py` | 13 |
| `_run_detail` | `tui/screens/calculations.py` | 13 |
| `disambiguate` | `kind_disambiguator.py` | 12 |
| `_legacy_plan_to_synthesis_result` | `local/adapter.py` | 12 |
| `load_active_provider_config` | `provider_config.py` | 12 |
| `get_provider` | `providers.py` | 12 |
| `_trim_context_to_budget` | `services/conversation_memory.py` | 12 |
| `summarize_tool_result` | `services/memory_summaries.py` | 12 |
| `_render_command_interpretation` | `tui/widgets/cells/command_interpretation.py` | 12 |
| `parse_pbs_pbsnodes_av` | `wizard/parsers.py` | 12 |
| `_probe_sge` | `wizard/survey.py` | 12 |
| `parse_coordinate_literal` | `harness/command_rules/structures.py` | 11 |
| `_kind_from_request` | `harness/intent.py` | 11 |
| `structured_sequence` | `harness/value_equivalence.py` | 11 |
| `consume` | `services/conversation_memory.py` | 11 |
| `parse_gaussian_scan_definition` | `services/scan_directives.py` | 11 |
| `_format_semantic_result` | `tui/chat_helpers.py` | 11 |
| `_handle_project_write_command` | `tui/mixins/project_write.py` | 11 |
| `_tool_success_note` | `tui/mixins/wizard_commands.py` | 11 |
| `render_tool_result_detail` | `tui/tool_meta.py` | 11 |
| `_completed_detail` | `tui/widgets/cells/plan.py` | 11 |
| `_parse_scan_definition` | `v8_adapter.py` | 11 |
| `parse_duration_seconds` | `wizard/parsers.py` | 11 |

## Behavior-preservation evidence

- **Full test suite**: `tests/agent` 1088/1088 green (also green under
  coverage instrumentation). Whole-repo run: 44 failures, all in
  non-agent core chemistry tests (`test_structures` 18, `test_database`
  17, `test_groupers` 7, `test_converter` 1, `cli/test_main` extras-hint
  1); the identical 44 fail at the pre-refactor branch point in the same
  environment (which additionally failed the 2 TUI tests repaired here) —
  the refactor fixed 2 tests and broke 0.
- **Frozen contract snapshot**: `scripts/review/snapshot_agent_contracts.py`
  regenerated at the pre-refactor branch point (`bf3eea67`) and at final
  HEAD *in the same environment*: the only difference is one traceback
  tail-window offset caused by worktree path length (same exception, same
  rule, same message) — zero semantic contract drift from the refactor.
  Aggregate invariants hold exactly: 30 CLI helps all exit 0, 30 unique
  tools, full26 all adapter-valid, semantic verdicts 18 ok / 8 reject,
  intent verdicts 21 ok / 5 reject.
  (The frozen `agent_contract_baseline.json` additionally differs from any
  current-environment regeneration in click version presentation strings
  only — `[COMMAND]` vs `COMMAND` usage lines, "No such option" wording,
  choice-case display; these predate the refactor.)
- **Install boundaries**: `scripts/review/verify_agent_install.py` passes
  for `core`, `agent`, and `tui` profiles resolved against this tree.
- **Branch coverage**: 76% (`pytest tests/agent --cov=chemsmart/agent
  --cov-branch`), above the 74.63% baseline.
- **CLI smoke**: `python -m chemsmart.cli.main agent --help` and
  `agent tools` render the full 30-tool surface.
- Two TUI tests that were failing at the branch point were repaired
  first (`test(tui): stabilize copy-view click and refresh quit snapshot`):
  one was presentation-only snapshot staleness inherited from the TUI-UX
  merge, one was a test-side scroll/reflow race — neither was a refactor
  regression.
