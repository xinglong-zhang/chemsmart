# chemsmart agent CLI/TUI — Orchestratic Design Decision

Decided by: cs-orchestrator
Date: 2026-05-08
Inputs: PR #26 (Claude Code research), PR #27 (Codex research), `/tmp/ui-research/SYNTHESIS.md`, live read of `chemsmart/agent/{cli,core}.py`

## TL;DR

Build a **Python + Textual + Rich** TUI as `chemsmart agent` (no subcommand). Layout: bottom composer + scrollable transcript + footer with phase indicator and job badges. The existing `AgentSession` and `DecisionLog` are not refactored — the TUI **tails `decision_log.jsonl`** to render typed cells and runs `AgentSession.run(...)` in a Textual worker. This keeps v1 strictly additive: every existing one-shot Click command (`chemsmart agent run/resume/doctor/tools`) keeps working, and the TUI is just a richer surface over the same artifacts.

## Decisions (1 line each)

| # | Decision | Choice | Rationale |
|---|---|---|---|
| 1 | Stack | Python + Textual ≥0.80,<1.0 + Rich | In-process access to `Molecule`/`Plan`; no Rust/Node IPC tax; Rich already used by `chemsmart.utils.logger`. |
| 2 | Render mode | Inline scrollback default; alt-screen only for `/jobs`, `/diff`, `/transcript` overlays | Postmortem on HPC SSH/tmux sessions requires real scrollback. (CC:46, CX:52) |
| 3 | Layout | Bottom composer + transcript above + footer | Convergent CC/CX pattern; matches REPL muscle memory. (CC:74, CX:75) |
| 4 | Streaming | Newline-gated, source-retained | Resize-stable; planner/critic outputs are markdown paragraphs. (CX:124) |
| 5 | Slash commands | Phase-tagged `[I/P/D/R/F/A]` — disabled when phase doesn't match | chemsmart's `/submit` is invalid before dry-run; encode that in the UI. (CX:201) |
| 6 | Approval UX | Per-request overlay: `y / n / s / r` (once/no/session/decline-and-revise). **No bypass mode for HPC.** | Decline-and-revise is the chemistry win: "use def2-TZVP instead, keep the rest." (CX:300) |
| 7 | Cells | Typed cells per chemistry artifact | `MoleculeCell`, `MethodCell`, `DryRunInputCell`, `GeometryHandoffCell`, `RuntimeValidationCell`, `CriticVerdictCell`, `SubmissionPreviewCell`, `JobStatusCell`, `RunResultCell`, `WorkflowCell`, `BatchProgressCell`. |
| 8 | Persistence | Reuse existing `~/.chemsmart/agent/<session-id>/decision_log.jsonl`; add `transcript.jsonl` and `state.json` if not already present | Zero schema migration; `DecisionLog` already exists at `core.py:93`. |
| 9 | Event source | TUI tails `decision_log.jsonl` via watchdog | **No refactor of `AgentSession.run`.** The session writes JSONL; the UI reads JSONL. Crash-safe and identical to the one-shot CLI's view of state. |
| 10 | Long-running jobs | First-class `JobsPanel` (`/jobs`) + footer badges | The single biggest chemistry-specific UX gap neither CC nor CX solves. |
| 11 | Ctrl-C | 1st = soft cancel + arm quit (3s); 2nd = quit | Chemistry runs are expensive; accidental session loss is unforgivable. (CC:121, CX:119) |
| 12 | Keymap | Flat (no Vim mode v1); custom file later | Ship working before customizable. (CX has Vim — defer.) |
| 13 | Autocomplete | `@` → file popup biased to `.xyz/.log/.com/.inp`; `/` → slash popup; phase-aware | Chemists reference structure files constantly. (CC:411, CX:398) |
| 14 | Theming | Dark/light auto-detect via Textual `terminal_dark`; ship 2 default themes (dark/light); custom themes later | Terminal palette probe borrowed from CX:71. |
| 15 | Headless fallback | `chemsmart agent --plain` disables alt-screen + mouse + animations | HPC login nodes have flaky TERM. Ship from day 1. |

## Command surface

### `chemsmart agent` (Click subgroup)

```
chemsmart agent                 # NEW: launch TUI
chemsmart agent ask "…"         # NEW: one-shot, streams to stdout, no TUI
chemsmart agent jobs            # NEW: top-level jobs panel without TUI
chemsmart agent run "…"         # PRESERVED: existing one-shot
chemsmart agent resume <id>     # PRESERVED
chemsmart agent doctor          # PRESERVED
chemsmart agent tools           # PRESERVED
```

Every existing entry point is preserved. The TUI is purely additive.

### REPL slash commands (28 total, phase-gated)

```
[A]  /help            command browser
[A]  /jobs            jobs panel (queued/running/done/failed)
[A]  /queue           remote scheduler queue (squeue)
[A]  /molecule <p>    load .xyz/.log/.com and preview
[I,P] /method         method recommendation panel
[I,P] /server         pick HPC server config
[D]  /dryrun          force dry-run regeneration
[D]  /diff            diff dry-run vs previous iteration
[D]  /submit          submit-approval overlay (HPC)
[D]  /run             run-locally approval overlay
[R]  /cancel <id>     scancel a queued/running job
[F]  /extract <step>  extract optimized geometry
[P,D] /critic         critic verdict + issues
[P,D,R] /plan         current plan as a tree
[A]  /rationale       full planner rationale
[I,P] /scan-setup     guided scan-coordinate form
[I]  /resume [id]     resume a prior session
[I]  /sessions        session picker
[I,F] /clear          new conversation
[A]  /transcript      scrollable transcript pager
[A]  /raw             toggle raw output mode
[A]  /copy            copy last assistant cell
[A]  /export <file>   export transcript
[A]  /theme           theme picker
[A]  /config          settings (model, dry-submit default, allow flags)
[A]  /doctor          inline doctor run
[A]  /tools           list registered tools
[A]  /quit, /exit     exit
```

Phase tags: `[I]` idle, `[P]` planning, `[D]` dry-run-ready, `[R]` running, `[F]` finished, `[A]` always.

### Keybindings (flat)

```
Enter            submit
Ctrl+J / S-Enter newline
Ctrl+G           open composer in $EDITOR
Ctrl+C (1st)     soft cancel + arm quit (3s)
Ctrl+C (2nd)     quit
Ctrl+D           quit if composer empty
Esc              cancel modal/popup; else interrupt
Ctrl+T           transcript pager
Ctrl+R           reverse history search
Ctrl+L           redraw
Ctrl+B           background current local run
Tab              accept autocomplete
@                file autocomplete popup
/                slash-command popup
!                shell mode (one-shot)
?                shortcut hints overlay
o (cell focus)   open underlying file ($EDITOR/PyMOL via xdg-open)
d (DryRunInputCell) diff against previous iteration
e (MethodCell)   edit method/basis
y/n/s/r          approval overlay: yes / no / session / decline-and-revise
```

## Architecture

```
chemsmart/
└── agent/
    ├── __init__.py              # existing
    ├── cli.py                   # MODIFY: add bare `agent` (no subcmd) → launch TUI; keep existing subcommands
    ├── core.py                  # UNCHANGED in v1 (DecisionLog already JSONL)
    ├── providers.py             # unchanged
    ├── registry.py              # unchanged
    ├── tools.py                 # unchanged
    ├── transport.py             # unchanged
    ├── prompts/                 # unchanged
    └── tui/                     # NEW
        ├── __init__.py
        ├── app.py               # ChemsmartTuiApp(textual.App)
        ├── styles.tcss          # Textual CSS
        ├── bindings.py          # keymap
        ├── phase.py             # Phase enum + transition rules
        ├── events.py            # parse decision_log.jsonl entries → typed events
        ├── screens/
        │   ├── chat.py          # main screen
        │   ├── sessions.py      # resume picker
        │   ├── jobs.py          # /jobs panel
        │   ├── diff.py          # /diff overlay
        │   ├── transcript.py    # /transcript pager
        │   └── theme.py         # /theme picker
        ├── widgets/
        │   ├── composer.py      # input editor
        │   ├── transcript.py    # cell stack
        │   ├── footer.py        # phase + job badges + hints
        │   ├── popups/          # @-file, /-slash, ?-hints
        │   └── cells/
        │       ├── base.py
        │       ├── user_message.py
        │       ├── agent_message.py
        │       ├── plan.py
        │       ├── method.py
        │       ├── molecule.py
        │       ├── dry_run_input.py
        │       ├── geometry_handoff.py
        │       ├── runtime_validation.py
        │       ├── critic_verdict.py
        │       ├── submission_preview.py
        │       ├── job_status.py
        │       ├── run_result.py
        │       ├── workflow.py
        │       ├── batch_progress.py
        │       └── error.py
        ├── services/
        │   ├── session_runner.py   # Textual @work wrapping AgentSession.run
        │   ├── log_tailer.py       # watchdog tail of decision_log.jsonl
        │   ├── job_poller.py       # background scheduler poller
        │   └── file_index.py       # @-file autocomplete source
        └── theme/
            ├── dark.tcss
            └── light.tcss
```

### Data flow

```
   User types "optimize h2o.xyz" + Enter
            │
            ▼
   ┌─────────────────┐    Textual @work        ┌────────────────┐
   │ Composer widget │ ───────────────────────▶│ session_runner │
   └─────────────────┘                          │ AgentSession   │
            ▲                                   │   .run(...)    │
            │                                   └────────┬───────┘
            │                                            │ writes
            │                                            ▼
   ┌────────────────┐    parse → typed event     ~/.chemsmart/agent/
   │  Transcript    │◀─────────────────────────  <id>/decision_log.jsonl
   │  (cell stack)  │     (log_tailer worker)               ▲
   └────────────────┘                                       │ writes
            ▲                                               │
            │                                       ┌───────┴────────┐
   ┌────────┴────────┐                              │ scheduler      │
   │ Footer badges   │◀──────────  job_poller  ◀───│ (sacct/squeue) │
   └─────────────────┘                              └────────────────┘
```

**Key invariant:** `AgentSession.run` is unchanged. The TUI is a viewer-and-launcher over the same artifacts the one-shot CLI produces. This means the headless `chemsmart agent run` and the interactive TUI cannot diverge.

### `AgentEvent` parser (`tui/events.py`)

Read each JSONL line from `decision_log.jsonl` and dispatch to a typed `AgentEvent`. Initial mapping (extend as `core.py` adds event types):

| `decision_log` event | TUI event class | Cell rendered |
|---|---|---|
| `request_received` | `RequestReceived` | `UserMessageCell` |
| `plan_proposed` | `PlanReady` | `PlanCell` |
| `step_started` | `StepStarted` | child of `WorkflowCell` |
| `tool_call` | `ToolCall` | typed by tool name |
| `dry_run_ready` | `DryRunReady` | `DryRunInputCell` |
| `critic_verdict` | `CriticVerdict` | `CriticVerdictCell` |
| `geometry_extracted` | `GeometryExtracted` | `GeometryHandoffCell` |
| `submit_preview` | `SubmissionPreview` | `SubmissionPreviewCell` |
| `submit_executed` | `JobSubmitted` | `JobStatusCell` |
| `step_completed` | `StepCompleted` | `RunResultCell` (if results present) |
| `session_summary` | `SessionEnded` | summary footer + phase → finished |
| `tool_error` | `ToolError` | `ErrorCell` |

A small adapter (`events.py`) is the only piece that mutates if `core.py` ever changes its event shape. Versioning lives in the JSONL `event_type` field.

## Phase model

```
            ┌──── /clear ────┐
            │                ▼
[idle] ──submit text──▶ [planning] ──plan_proposed──▶ [dry-run-ready]
                                                            │
                                                  /submit or /run
                                                            ▼
[finished] ◀─session_summary─ [running] ◀── approval(yes) ─┤
            │                                               │
            └─── /clear ────▶ [idle]              approval(decline-and-revise)
                                                            │
                                                            ▼
                                                       [planning]
```

The footer always shows the current phase as a colored chip:
`◌ idle` (grey) · `◐ planning` (yellow) · `◑ dry-run-ready` (cyan) · `▶ running` (green) · `✓ finished` (dim green) · `✗ error` (red).

## Risks & mitigations

| Risk | Mitigation |
|---|---|
| Async streaming through Click context | Bare `chemsmart agent` Click command does only `asyncio.run(ChemsmartTuiApp().run_async())`; nothing more. Click ctx is never passed into Textual. |
| HPC headless terminal compatibility | `--plain` mode (no alt-screen, no mouse, no animation). Test on actual cluster login node before merging Phase 1. |
| Background poller leaking | Use Textual `@work(exclusive=True)`. Install SIGTERM/SIGINT handlers. PID file under `~/.chemsmart/agent/sessions/<id>/.tui.pid`. |
| Concurrent TUI + one-shot on same session | `fcntl.flock` on `state.json`. On stale lock (>10 min), warn but allow read-only attach. |
| Textual API churn | Pin `textual >=0.80,<1.0`. Add snapshot tests via `textual.pilot.run_test()` for every cell type. |
| 50MB Gaussian/ORCA logs in transcript | Never render raw output in cells. Always parse via `chemsmart.io` first; raw `.log` only via `o`-key open in `$EDITOR`/`less`. |
| Light/dark contrast | Probe `terminal_dark` at startup; ship dark + light theme files; theme picker side-preview shows a `.com` snippet (not a Rust diff). |

## Non-goals (v1)

No plugin marketplace · No MCP server hosting · No image input · No voice mode · No web/desktop/IDE companion · No Vim mode · No GitHub/PR integration · No remote-control · No batch scan-setup wizard (defer to Phase 4).

## Phased plan

### Phase 1 — MVP TUI (~1 week, small)
**Goal:** Replace `chemsmart agent run` console output with a proper interactive REPL, while leaving `run` itself untouched.

Deliverables:
- `chemsmart/agent/tui/app.py`, `screens/chat.py`, `widgets/{composer,transcript,footer}.py`.
- Cells: `UserMessage`, `AgentMessage`, `Plan`, `DryRunInput` (read-only), `CriticVerdict`, `Error`.
- `services/session_runner.py` (Textual `@work`-wrapped `AgentSession.run`).
- `services/log_tailer.py` (watchdog tail → events).
- Slash commands: `/help`, `/quit`, `/clear`, `/sessions`, `/resume`, `/tools`, `/doctor`.
- Ctrl+C soft cancel.
- `chemsmart agent --plain` headless mode.
- Snapshot tests (`tests/agent/tui/`) for every cell.
- README section + screenshot in `docs/agent-setup.md`.

Acceptance: I can run `chemsmart agent`, type a request, see the plan + dry-run + critic verdict render as cells, Ctrl+C to interrupt, `/quit` to exit. The same request via `chemsmart agent run` still produces identical output to today's master.

### Phase 2 — Approval & artifact UX (~2 weeks, medium)
- Per-request approval overlay (`y/n/s/r`) for `submit_hpc` and `run_local`.
- Cells: `Method` (with `e` to edit), `GeometryHandoff` (collapsed/expanded; `o` opens), `RuntimeValidation`, `SubmissionPreview`.
- `DryRunInput` gains `d`-key diff against previous iteration (use `difflib.unified_diff`).
- Slash: `/dryrun`, `/submit`, `/run`, `/critic`, `/plan`, `/rationale`.
- `@`-file autocomplete biased to chemistry extensions.
- `Ctrl+G` external editor handoff.
- `[Pasted N chars]` placeholder for >10k pastes.
- Footer status badges (counts only; no jobs panel yet).

Acceptance: A full opt+freq workflow runs end-to-end with explicit approval. Decline-and-revise sends a corrective message back to the planner.

### Phase 3 — Jobs panel & HPC integration (~2 weeks, medium)
- `JobsPanel` overlay (`/jobs`) reading scheduler state files + `state.json`.
- `services/job_poller.py` background worker (configurable poll interval, backoff on failure).
- Slash: `/cancel`, `/extract`, `/server`, `/queue`.
- Cells: `JobStatus`, `RunResult` (parsed energies/frequencies via `chemsmart.io`), `Workflow` (parent for opt+freq+IRC).
- cwd-mismatch prompt on resume.
- `/molecule <path>` slash command + `MoleculeCell`.

Acceptance: A submitted SLURM job is visible in `/jobs` with live status; finishing the job auto-renders a `RunResult` cell with parsed SCF energy and lowest 5 frequencies.

### Phase 4 — Polish & batch (~3–4 weeks, large, optional)
- `BatchProgressCell` for `iterate`/grouper.
- `/scan-setup` guided form.
- Theme picker, custom `~/.chemsmart/keybindings.json`.
- `/export` (Markdown + HTML).
- Performance: virtualize `RichLog` for very long sessions.
- File-locking and crash-recovery hardening.

## Implementation hand-off

Phase 1 is small enough to spawn as **one ao session** per the chemsmart waveform pattern. Suggested branch name: `feat/agent-tui-mvp`. Reference this decision doc and the two research PRs (#26, #27) in the spawn prompt.

Critical files for Phase 1 (touch list):
- `chemsmart/agent/cli.py` (lines 14–17: add bare `agent` command launching TUI)
- `pyproject.toml` (add `textual>=0.80,<1.0`, `watchdog`, `pyperclip` to optional `[agent]` extra)
- `chemsmart/agent/tui/**` (all new)
- `tests/agent/tui/**` (all new — Textual `pilot.run_test`)
- `docs/agent-tui.md` (new)

Do NOT touch in Phase 1: `chemsmart/agent/core.py`, `chemsmart/agent/registry.py`, `chemsmart/agent/tools.py`, `chemsmart/agent/providers.py`. The TUI is strictly a layer over these.
