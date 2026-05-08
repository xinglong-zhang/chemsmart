# OpenAI Codex CLI UI patterns

Research scope: OpenAI Codex CLI / TUI only. No Claude Code comparisons beyond brief framing. Source basis: local clone of `openai/codex` at commit `cce059467af64f05d0a1521344847d2b558b6a80` plus chemsmart grounding from local `AGENTS.md` and merged fork PRs #15, #16, #17, and #23.

## Reading notes
- Codex repo clone used for file/line citations: `/tmp/codex-research` at commit `cce059467af64f05d0a1521344847d2b558b6a80`.
- Chemsmart grounding came from `/Users/hongjiseung/developer/chemsmart/AGENTS.md` and fork PRs `Hongjiseung-ROK/chemsmart#15`, `#16`, `#17`, and `#23`.
- Every non-obvious product/UI claim below is tied to code paths, usually `codex-rs/tui/src/...:line-line`.

---

## 1. Tool overview, repo layout, key crates in `codex-rs/`

- Top-level positioning is simple: `README.md` describes Codex CLI as a local coding agent, `codex-cli/` is the npm launcher, and `codex-rs/` contains the native implementation, including the TUI binary and supporting crates. (`README.md:1-8`, `README.md:13-29`, `codex-cli/package.json:1-21`, `codex-rs/Cargo.toml:1-111`)
- The JS package is intentionally thin: `@openai/codex` exports a single `codex` binary entry and ships platform-specific native payloads under `vendor`; the Node wrapper just selects the correct executable, forwards signals, and preserves exit semantics. (`codex-cli/package.json:5-21`, `codex-cli/bin/codex.js:15-22`, `codex-cli/bin/codex.js:175-229`)
- `codex-rs` is a large Cargo workspace, not a one-crate app. The member list includes app-server, rollout, protocol, state, plugins, model-provider, tooling, and the `tui` crate itself. (`codex-rs/Cargo.toml:1-111`)
- The TUI crate is explicitly called out in Codex’s own contributor guidance as a central, high-touch orchestration area; the repo tells contributors to keep `chatwidget.rs`, `bottom_pane/mod.rs`, and `chat_composer.rs` from growing further, which is a good hint about real ownership boundaries. (`AGENTS.md:42-54`)

### 1.1 Practical repo map
| Path | Role | Why it matters for UI research |
| `README.md` | Product framing and install entrypoint | Explains local-agent positioning and distribution story. (`README.md:1-8`, `README.md:13-29`) |
| `AGENTS.md` | Contributor/convention guide | Names TUI hot files and snapshot-test expectations. (`AGENTS.md:47-54`, `AGENTS.md:109-129`) |
| `codex-cli/` | npm wrapper | Shows native binary handoff, package shape, and signal forwarding. (`codex-cli/package.json:5-21`, `codex-cli/bin/codex.js:120-229`) |
| `codex-rs/tui/` | Native TUI crate | Main chat UI, overlays, input state machine, approvals, theme picker. (`codex-rs/tui/Cargo.toml:24-160`) |
| `codex-rs/rollout/` | Rollout/session record plumbing | Referenced by resume/fork fallback paths. (`codex-rs/Cargo.toml:192-205`, `codex-rs/tui/src/session_resume.rs:1-5`) |
| `codex-rs/protocol/` | Core protocol types | Thread items, tool calls, permissions, plans, etc. feed the UI cells. (`codex-rs/Cargo.toml:187-205`) |
| `codex-rs/state/` | Saved thread metadata/state runtime | Used during resume to recover cwd/model from app-server state. (`codex-rs/Cargo.toml:199-205`, `codex-rs/tui/src/session_resume.rs:67-84`) |
| `codex-rs/app-server-*` | Thread/app-server lifecycle | Transcript loading, resume picker, approvals, and thread switching are app-server-backed. (`codex-rs/Cargo.toml:12-17`, `codex-rs/tui/src/resume_picker/transcript.rs:24-37`) |

### 1.2 Key crates worth knowing before reading the UI
| Crate | Observed role | UI-facing consequence |
| `codex-tui` | Ratatui app surface plus modal/overlay/input logic | Owns rendering, key routing, slash commands, approvals, resume picker, transcript overlay. (`codex-rs/tui/Cargo.toml:24-160`) |
| `codex-app-server-client` | Thread/app-server access from TUI | Resume picker and transcript previews read thread state through it. (`codex-rs/tui/Cargo.toml:29-31`, `codex-rs/tui/src/resume_picker/transcript.rs:24-33`) |
| `codex-app-server-protocol` | Thread item and approval request schemas | History cells and approval overlays mirror these protocol types closely. (`codex-rs/tui/Cargo.toml:29-31`, `codex-rs/tui/src/bottom_pane/approval_overlay.rs:308-338`) |
| `codex-file-search` | Fast file match backend | Drives `@path` autocomplete popup with stale-query protection. (`codex-rs/tui/Cargo.toml:43`, `codex-rs/tui/src/bottom_pane/file_search_popup.rs:16-78`) |
| `codex-rollout` | Rollout/session log substrate | Resume/fork can fall back to rollout JSONL when thread metadata is not yet hydrated. (`codex-rs/Cargo.toml:192-205`, `codex-rs/tui/src/session_resume.rs:144-189`) |
| `codex-state` | Persistent thread metadata store | Resume reads saved cwd/model from it before consulting rollout fallback. (`codex-rs/tui/Cargo.toml:56`, `codex-rs/tui/src/session_resume.rs:67-84`, `codex-rs/tui/src/session_resume.rs:114-137`) |
| `codex-protocol` | Message/tool abstractions | Streaming, markdown cells, and MCP/tool rendering depend on these types. (`codex-rs/tui/Cargo.toml:52`, `codex-rs/tui/src/history_cell.rs:1815-1821`) |
| `codex-terminal-detection` | Terminal/multiplexer detection | Used for alt-screen heuristics and terminal-specific shortcut fallbacks. (`codex-rs/tui/Cargo.toml:57`, `codex-rs/tui/src/lib.rs:1584-1615`, `codex-rs/tui/src/chatwidget.rs:208-240`) |
| `codex-realtime-webrtc` | Realtime voice backend | Explains why some input/interrupt paths special-case live conversations. (`codex-rs/tui/Cargo.toml:53`, `codex-rs/tui/src/chatwidget.rs:10211-10223`) |

### 1.3 Architecture takeaway
- The Codex UI is not a thin terminal wrapper around a backend process; it is a first-class client with its own session persistence, history search, approval UX, transcript overlays, theme picker, and editable composer state machine. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:1-133`, `codex-rs/tui/src/pager_overlay.rs:1-16`, `codex-rs/tui/src/session_log.rs:80-118`)
- For chemsmart, that means the “frontend” should probably own its own chat-session ergonomics rather than just printing agent events from Click. Chemsmart already spans local execution and HPC submission, so a richer local client would pay off. (`/Users/hongjiseung/developer/chemsmart/AGENTS.md:23-35`, `/Users/hongjiseung/developer/chemsmart/AGENTS.md:74-77`)

## 2. Tech stack & rendering library (ratatui? crossterm? alt-screen? mouse?)

- The TUI stack is Rust + `ratatui` + `crossterm`. The crate manifest pulls in `ratatui` with unstable widget-ref and scrolling-region features, and `crossterm` with `bracketed-paste` and `event-stream`. (`codex-rs/tui/Cargo.toml:70-87`)
- Syntax highlighting comes from `syntect` plus `two-face`, not from tree-sitter. (`codex-rs/tui/Cargo.toml:114-115`, `codex-rs/tui/src/render/highlight.rs:1-22`)
- The backend is `CrosstermBackend<Stdout>` wrapped in a custom terminal type; terminal events are normalized into `TuiEvent::{Key,Paste,Resize,Draw}` before the app sees them. (`codex-rs/tui/src/tui.rs:28-39`, `codex-rs/tui/src/tui.rs:402-415`)
- Startup enables bracketed paste, raw mode, keyboard-enhancement flags, and focus-change notifications. That directly supports modified-Enter disambiguation, large paste handling, and focus-sensitive notifications. (`codex-rs/tui/src/tui.rs:156-169`)
- Default mode is **inline viewport**, not fullscreen alt-screen: the init comment explicitly says “history stays in normal scrollback.” (`codex-rs/tui/src/tui.rs:329-369`)
- Alt-screen exists, but it is opt-in/contextual. Codex can disable it entirely with `--no-alt-screen`, and `auto` mode disables it specifically under Zellij because of scrollback behavior. (`codex-rs/tui/src/cli.rs:64-70`, `codex-rs/tui/src/lib.rs:1584-1615`)
- When alt-screen is entered, Codex also enables “alternate scroll,” saves the inline viewport, expands to full terminal size, and restores the previous viewport when leaving. (`codex-rs/tui/src/tui.rs:613-648`)
- The event broker explicitly skips mouse events, and I found no mouse-capture enablement path in `tui.rs`; the current TUI appears intentionally keyboard-first. This is an inference from positive evidence, not a product statement. (`codex-rs/tui/src/tui/event_stream.rs:236-245`, `codex-rs/tui/src/tui.rs:16-25`)
- Color depth is terminal-sensitive: Codex detects truecolor vs ANSI-256 vs ANSI-16/unknown, and on Unix it also probes the terminal’s default foreground/background to choose better themes and diff colors. (`codex-rs/tui/src/terminal_palette.rs:4-18`, `codex-rs/tui/src/terminal_palette.rs:71-151`)

### 2.1 Minimal code snippets worth seeing
```rust
EnableBracketedPaste
```
Bracketed paste is part of initial terminal mode setup. (`codex-rs/tui/src/tui.rs:156-169`)

```rust
Viewport::Inline
```
Conceptually, Codex starts in inline scrollback-preserving mode, then uses alt-screen only for overlay-style views. The exact viewport construction is abstracted through `CustomTerminal`, but the inline intent is explicit in the init comment and alt-screen enter/leave helpers. (`codex-rs/tui/src/tui.rs:329-369`, `codex-rs/tui/src/tui.rs:613-648`)

### 2.2 Why this stack matters for chemsmart
- Chemsmart already has a file-oriented workflow and many long-running local/HPC steps; an inline-first terminal UI is attractive because users can keep full scrollback of dry-run inputs, route lines, and scheduler messages. (`/Users/hongjiseung/developer/chemsmart/AGENTS.md:23-35`, `/Users/hongjiseung/developer/chemsmart/AGENTS.md:117-137`, `codex-rs/tui/src/tui.rs:329-369`)
- At the same time, chemistry-specific overlays (job diff, transcript, queue inspector, scheduler output viewer) can still borrow Codex’s alt-screen pager pattern. (`codex-rs/tui/src/pager_overlay.rs:1-16`, `codex-rs/tui/src/app/event_dispatch.rs:351-366`)

## 3. Layout structure (input box, transcript, status bar, side panels, modal screens)

- The main chat surface is a vertical stack: committed transcript/history above, an in-flight active cell in the middle, optional hook cell below that, and the bottom pane last. (`codex-rs/tui/src/chatwidget.rs:6-18`, `codex-rs/tui/src/chatwidget.rs:10956-11020`)
- The bottom pane is itself a stack: optional status row, optional unified exec footer when no status row exists, pending-thread approvals preview, pending-input preview, then the composer. (`codex-rs/tui/src/bottom_pane/mod.rs:220-230`, `codex-rs/tui/src/bottom_pane/mod.rs:1520-1564`)
- The composer is always retained in memory even when another bottom-pane view is covering it; Codex preserves draft input state while modals/pickers are open. (`codex-rs/tui/src/bottom_pane/mod.rs:200-206`)
- Bottom-pane modals are formalized through a `BottomPaneView` trait with completion/cancel semantics, paste handling hooks, approval-request consumption hooks, and a flag for “Action Required” terminal-title state. (`codex-rs/tui/src/bottom_pane/bottom_pane_view.rs:17-138`)
- Transcript and diff are not bottom-pane views; they are pager overlays rendered in alternate screen. The overlay system has two modes: `Transcript` and generic `Static`. (`codex-rs/tui/src/pager_overlay.rs:1-16`, `codex-rs/tui/src/pager_overlay.rs:50-89`)
- Pager overlays render a dim `/ TITLE` header, scrollable content, `~` filler lines, and a bottom bar with percent scrolled. (`codex-rs/tui/src/pager_overlay.rs:174-180`, `codex-rs/tui/src/pager_overlay.rs:207-247`)
- The transcript overlay is special because it renders committed cells plus a cached live tail derived from the active in-flight cell, so users can inspect current tool execution without waiting for finalization. (`codex-rs/tui/src/chatwidget.rs:6-16`, `codex-rs/tui/src/pager_overlay.rs:409-456`)
- The session resume picker is its own interactive screen with toolbar controls, typed search, expandable transcript context, and pagination. (`codex-rs/tui/src/resume_picker.rs:286-297`)
- Theme selection is also a full selection view and can show side-by-side preview when the terminal is wide enough; this is one of the few places where Codex really uses a “side panel” concept. (`codex-rs/tui/src/theme_picker.rs:14-20`, `codex-rs/tui/src/theme_picker.rs:129-146`, `codex-rs/tui/src/theme_picker.rs:314-390`)

### 3.1 Main viewport anatomy
| Layer | What it holds | Key evidence |
| Committed transcript/history | Finalized `HistoryCell`s | `ChatWidget` owns committed cells plus active cell distinction. (`codex-rs/tui/src/chatwidget.rs:6-16`) |
| Active cell | Mutating in-flight assistant/tool/exec group | Inserted above bottom pane and re-rendered in place while streaming. (`codex-rs/tui/src/chatwidget.rs:6-10`, `codex-rs/tui/src/chatwidget.rs:10956-10973`) |
| Active hook cell | Hook-related inline cell | Conditionally rendered between active cell and bottom pane. (`codex-rs/tui/src/chatwidget.rs:10963-10979`) |
| Status row | Task-running indicator and mirrored exec summary | Only shown when available; hidden during some commentary streaming. (`codex-rs/tui/src/chatwidget.rs:24-27`, `codex-rs/tui/src/bottom_pane/mod.rs:220-226`, `codex-rs/tui/src/bottom_pane/mod.rs:1525-1535`) |
| Unified exec footer | Session summary when no status row exists | Mutually exclusive with status-row surfacing. (`codex-rs/tui/src/bottom_pane/mod.rs:223-226`, `codex-rs/tui/src/bottom_pane/mod.rs:1528-1535`) |
| Pending thread approvals preview | Cross-thread pending approval indicators | Shown inline above composer when present. (`codex-rs/tui/src/bottom_pane/mod.rs:229-230`, `codex-rs/tui/src/bottom_pane/mod.rs:1536-1556`) |
| Pending input preview | Queued drafts / pending steers / rejected steers | Shown above composer with spacing rules. (`codex-rs/tui/src/bottom_pane/mod.rs:227-228`, `codex-rs/tui/src/bottom_pane/mod.rs:1537-1563`) |
| Composer | Editable prompt buffer + popups + attachments | Always retained; final row of bottom pane. (`codex-rs/tui/src/bottom_pane/mod.rs:200-206`, `codex-rs/tui/src/bottom_pane/mod.rs:1560-1563`) |

### 3.2 Overlays and modal surfaces
| Surface | Rendering model | Open path |
| Transcript overlay | Alt-screen pager with live-tail cache | `Ctrl+T` opens `Overlay::Transcript`. (`codex-rs/tui/src/app/input.rs:166-174`, `codex-rs/tui/src/pager_overlay.rs:50-89`, `codex-rs/tui/src/pager_overlay.rs:409-456`) |
| Diff overlay | Alt-screen static pager | `/diff` completes into `Overlay::Static` after entering alt-screen. (`codex-rs/tui/src/chatwidget/slash_dispatch.rs:328-335`, `codex-rs/tui/src/app/event_dispatch.rs:351-366`) |
| Approval overlay | Bottom-pane modal list view | Uses `ApprovalOverlay`, not alt-screen, and can switch threads/fullscreen. (`codex-rs/tui/src/bottom_pane/approval_overlay.rs:155-167`, `codex-rs/tui/src/bottom_pane/approval_overlay.rs:508-547`) |
| Theme picker | Selection view with optional side preview | Opened from `/theme`, previews on movement, restores on cancel. (`codex-rs/tui/src/theme_picker.rs:1-20`, `codex-rs/tui/src/theme_picker.rs:314-390`) |
| Resume picker | Dedicated interactive browser | Search/paginate/expand sessions. (`codex-rs/tui/src/resume_picker.rs:286-297`) |
| Onboarding auth | Full-screenish onboarding step widgets | Rounded-border input and hyperlink-marked URLs. (`codex-rs/tui/src/onboarding/auth.rs:56-104`, `codex-rs/tui/src/onboarding/auth.rs:670-684`) |

### 3.3 Small snippet: bottom pane is a real compositional layout
```rust
flex2.push(/*flex*/ 0, RenderableItem::Borrowed(&self.composer));
```
That line is representative: Codex treats the composer as just one renderable in a richer bottom-pane stack, not as the whole footer. (`codex-rs/tui/src/bottom_pane/mod.rs:1560-1563`)

## 4. Streaming / token rendering & interruption (Ctrl-C, ESC, cancel mid-tool)

- Codex streams markdown deltas into a controller that keeps both raw source and wrapped render lines. Live UI gets incremental `HistoryCell`s; finalized transcript gets source-backed re-renderability. (`codex-rs/tui/src/streaming/controller.rs:1-11`, `codex-rs/tui/src/streaming/controller.rs:29-45`, `codex-rs/tui/src/streaming/controller.rs:237-285`)
- Streaming is **newline-gated**. The controller buffers deltas until a newline boundary, then commits the newly renderable source into queued live lines. (`codex-rs/tui/src/streaming/controller.rs:61-76`, `codex-rs/tui/src/streaming/controller.rs:261-266`)
- On resize, Codex re-renders from source and rebuilds only the not-yet-emitted queue, so streaming output stays stable without duplicating already shown lines. (`codex-rs/tui/src/streaming/controller.rs:119-149`, `codex-rs/tui/src/streaming/controller.rs:222-233`)
- While preamble/commentary is streaming, Codex hides the status row to avoid duplicate progress indicators; when commentary drains, it re-shows the status row. (`codex-rs/tui/src/chatwidget.rs:24-27`)
- Interruptive UI events that arrive during streaming (approvals, permissions, request-user-input, item start/end) are queued FIFO and flushed only after the stream completes, preserving event order. (`codex-rs/tui/src/chatwidget/interrupts.rs:16-106`, `codex-rs/tui/src/chatwidget.rs:4304-4327`)
- Esc has multiple roles: if a task is running and the composer has no popup, the bottom pane can send an interrupt; if pending steers exist, Esc interrupts and then re-submits the steers after interruption. (`codex-rs/tui/src/bottom_pane/mod.rs:608-623`, `codex-rs/tui/src/chatwidget.rs:5225-5236`, `codex-rs/tui/src/chatwidget.rs:3154-3190`)
- Ctrl+C is intentionally not instant-exit. First press arms a timed quit shortcut and, if cancellable work exists, also submits an interrupt; second same-key press before expiry quits. (`codex-rs/tui/src/chatwidget.rs:10205-10320`)
- Ctrl+D participates in the same double-press quit logic, but only when the composer is empty and no modal/popup is active. (`codex-rs/tui/src/chatwidget.rs:10265-10294`)
- Mid-tool interruption is surfaced in cells: interrupted MCP tool calls are marked failed with the string `interrupted`. (`codex-rs/tui/src/history_cell.rs:1832-1835`)

### 4.1 Streaming snippet
```rust
if delta.contains('\n')
```
That one condition is the key design choice: Codex does not repaint every token blindly; it waits for newline-complete source for stable markdown rendering. (`codex-rs/tui/src/streaming/controller.rs:67-72`)

### 4.2 Interruption matrix
| Input | Primary effect | Secondary effect / nuance | Evidence |
| `Esc` during running task | Interrupt current task | Only when no popup/modal should consume it first. (`codex-rs/tui/src/bottom_pane/mod.rs:608-623`) |
| `Esc` with pending steers | Interrupt current task | After interrupt, queued steers are restored/resubmitted. (`codex-rs/tui/src/chatwidget.rs:5225-5236`, `codex-rs/tui/src/chatwidget.rs:3154-3190`) |
| First `Ctrl+C` | Arm timed quit shortcut | Also interrupts active task if cancellable work exists. (`codex-rs/tui/src/chatwidget.rs:10207-10215`, `codex-rs/tui/src/chatwidget.rs:10258-10263`) |
| Second `Ctrl+C` before timeout | Quit | Skips another confirmation layer. (`codex-rs/tui/src/chatwidget.rs:10251-10255`) |
| `Ctrl+C` inside modal | Cancel modal or clear draft/history search first | Bottom pane gets first chance to consume it. (`codex-rs/tui/src/bottom_pane/mod.rs:649-682`, `codex-rs/tui/src/chatwidget.rs:10225-10236`) |
| Live voice + `Ctrl+C` | Stop realtime conversation | Realtime stop has precedence over generic quit/interrupt behavior. (`codex-rs/tui/src/chatwidget.rs:10211-10223`) |

## 5. Slash commands (extracted list) & keybindings (full mapping)

- Slash commands are a first-line UX primitive, not a thin alias layer. Their enum order is intentionally presentation order, with common commands first. (`codex-rs/tui/src/slash_command.rs:7-15`)
- The command popup is dynamically feature-gated, and availability is separately tracked for inline arguments, side conversations, and in-progress tasks. (`codex-rs/tui/src/slash_command.rs:146-241`, `codex-rs/tui/src/bottom_pane/slash_commands.rs:13-45`)
- Keybindings are centralized in `keymap.rs` and remappable via `/keymap`; the defaults include cross-terminal compatibility variants when terminals disagree about modifier preservation. (`codex-rs/tui/src/keymap.rs:542-546`, `codex-rs/tui/src/slash_command.rs:122-123`)

### 5.1 Slash commands
| Command | Description | Inline args? | During task? | Side convo? | Visibility / notes |
| `/model` | model/reasoning picker | No | No | No | Popup-first command. (`codex-rs/tui/src/slash_command.rs:15`, `codex-rs/tui/src/slash_command.rs:106-107`, `codex-rs/tui/src/slash_command.rs:180-202`) |
| `/fast` | fast-mode toggle | Yes | No | No | Can take inline args. (`codex-rs/tui/src/slash_command.rs:16`, `codex-rs/tui/src/slash_command.rs:107-109`, `codex-rs/tui/src/slash_command.rs:152-155`) |
| `/ide` | inject IDE context | Yes | Yes | Yes | One of few side-convo-safe commands. (`codex-rs/tui/src/slash_command.rs:17`, `codex-rs/tui/src/slash_command.rs:110-112`, `codex-rs/tui/src/slash_command.rs:174-175`, `codex-rs/tui/src/slash_command.rs:221-224`) |
| `/permissions` | permissions picker | No | No | No | Opens permissions UI. (`codex-rs/tui/src/slash_command.rs:18`, `codex-rs/tui/src/slash_command.rs:121-123`, `codex-rs/tui/src/slash_command.rs:189-190`) |
| `/keymap` | shortcut remapping UI | Yes | No | No | TUI config surface. (`codex-rs/tui/src/slash_command.rs:19`, `codex-rs/tui/src/slash_command.rs:122-123`, `codex-rs/tui/src/slash_command.rs:156-157`, `codex-rs/tui/src/slash_command.rs:190-190`) |
| `/vim` | composer Vim-mode toggle | No | No | No | Editor-mode toggle. (`codex-rs/tui/src/slash_command.rs:20`, `codex-rs/tui/src/slash_command.rs:123-124`, `codex-rs/tui/src/slash_command.rs:191-191`) |
| `/setup-default-sandbox` | elevated sandbox setup | No | No | No | Windows/elevated sandbox flow. (`codex-rs/tui/src/slash_command.rs:21-24`, `codex-rs/tui/src/slash_command.rs:124-125`, `codex-rs/tui/src/slash_command.rs:192-193`) |
| `/sandbox-add-read-dir` | add read-access directory | Yes | No | No | Windows-visible only. (`codex-rs/tui/src/slash_command.rs:23-24`, `codex-rs/tui/src/slash_command.rs:125-127`, `codex-rs/tui/src/slash_command.rs:160-162`, `codex-rs/tui/src/slash_command.rs:237-240`) |
| `/experimental` | feature-flags popup | No | No | No | Experimental control surface. (`codex-rs/tui/src/slash_command.rs:25`, `codex-rs/tui/src/slash_command.rs:128-128`, `codex-rs/tui/src/slash_command.rs:194-194`) |
| `/approve` | retry blocked auto-review once | No | Yes | No | Enum variant `AutoReview`. (`codex-rs/tui/src/slash_command.rs:26-27`, `codex-rs/tui/src/slash_command.rs:129-129`, `codex-rs/tui/src/slash_command.rs:219-220`) |
| `/memories` | memory configuration | No | No | No | State/config surface. (`codex-rs/tui/src/slash_command.rs:28`, `codex-rs/tui/src/slash_command.rs:130-130`, `codex-rs/tui/src/slash_command.rs:195-195`) |
| `/skills` | skills browser/entrypoint | No | Yes | No | Safe while task runs. (`codex-rs/tui/src/slash_command.rs:29`, `codex-rs/tui/src/slash_command.rs:95-95`, `codex-rs/tui/src/slash_command.rs:207-207`) |
| `/hooks` | hook manager | No | Yes | No | Safe while task runs. (`codex-rs/tui/src/slash_command.rs:30`, `codex-rs/tui/src/slash_command.rs:96-96`, `codex-rs/tui/src/slash_command.rs:208-208`) |
| `/review` | request code review | Yes | No | No | Inline args allowed. (`codex-rs/tui/src/slash_command.rs:31`, `codex-rs/tui/src/slash_command.rs:85-85`, `codex-rs/tui/src/slash_command.rs:150-150`, `codex-rs/tui/src/slash_command.rs:196-196`) |
| `/rename` | rename thread | Yes | Yes | No | Editable thread label. (`codex-rs/tui/src/slash_command.rs:32`, `codex-rs/tui/src/slash_command.rs:86-86`, `codex-rs/tui/src/slash_command.rs:151-151`, `codex-rs/tui/src/slash_command.rs:205-205`) |
| `/new` | start fresh chat | No | No | No | Conversation reset path. (`codex-rs/tui/src/slash_command.rs:33`, `codex-rs/tui/src/slash_command.rs:82-82`, `codex-rs/tui/src/slash_command.rs:181-181`) |
| `/resume` | resume prior session | Yes | No | No | Inline args allowed. (`codex-rs/tui/src/slash_command.rs:34`, `codex-rs/tui/src/slash_command.rs:87-87`, `codex-rs/tui/src/slash_command.rs:160-160`, `codex-rs/tui/src/slash_command.rs:182-182`) |
| `/fork` | fork current chat | No | No | No | Conversation fork. (`codex-rs/tui/src/slash_command.rs:35`, `codex-rs/tui/src/slash_command.rs:89-89`, `codex-rs/tui/src/slash_command.rs:183-183`) |
| `/init` | bootstrap AGENTS file | No | No | No | Repo setup helper. (`codex-rs/tui/src/slash_command.rs:36`, `codex-rs/tui/src/slash_command.rs:83-83`, `codex-rs/tui/src/slash_command.rs:184-184`) |
| `/compact` | context compaction | No | No | No | Conversation summarization path. (`codex-rs/tui/src/slash_command.rs:37`, `codex-rs/tui/src/slash_command.rs:84-84`, `codex-rs/tui/src/slash_command.rs:185-185`) |
| `/plan` | switch plan-mode state | Yes | No | No | Inline args allowed. (`codex-rs/tui/src/slash_command.rs:38`, `codex-rs/tui/src/slash_command.rs:116-116`, `codex-rs/tui/src/slash_command.rs:152-153`, `codex-rs/tui/src/slash_command.rs:197-197`) |
| `/goal` | set/view long-running goal | Yes | Yes | No | Available during task. (`codex-rs/tui/src/slash_command.rs:39`, `codex-rs/tui/src/slash_command.rs:117-117`, `codex-rs/tui/src/slash_command.rs:153-154`, `codex-rs/tui/src/slash_command.rs:213-213`) |
| `/collab` | collaboration-mode switch | No | Yes | No | Experimental mode control. (`codex-rs/tui/src/slash_command.rs:40`, `codex-rs/tui/src/slash_command.rs:118-118`, `codex-rs/tui/src/slash_command.rs:229-229`) |
| `/agent` | switch agent thread | No | Yes | No | Shares semantics with `/subagents`. (`codex-rs/tui/src/slash_command.rs:41`, `codex-rs/tui/src/slash_command.rs:119-120`, `codex-rs/tui/src/slash_command.rs:230-230`) |
| `/side` | spawn ephemeral side conversation | Yes | Yes | No | Inline args allowed. (`codex-rs/tui/src/slash_command.rs:42`, `codex-rs/tui/src/slash_command.rs:120-120`, `codex-rs/tui/src/slash_command.rs:159-160`, `codex-rs/tui/src/slash_command.rs:224-224`) |
| `/copy` | copy last assistant markdown | No | Yes | Yes | Hidden on Android. (`codex-rs/tui/src/slash_command.rs:43`, `codex-rs/tui/src/slash_command.rs:91-91`, `codex-rs/tui/src/slash_command.rs:169-175`, `codex-rs/tui/src/slash_command.rs:203-203`, `codex-rs/tui/src/slash_command.rs:238-239`) |
| `/raw` | toggle raw-scrollback rendering | Yes | Yes | Yes | Copy-friendly selection mode. (`codex-rs/tui/src/slash_command.rs:44`, `codex-rs/tui/src/slash_command.rs:92-92`, `codex-rs/tui/src/slash_command.rs:158-159`, `codex-rs/tui/src/slash_command.rs:170-170`, `codex-rs/tui/src/slash_command.rs:204-204`) |
| `/diff` | open repo diff | No | Yes | Yes | Leads to alt-screen pager. (`codex-rs/tui/src/slash_command.rs:45`, `codex-rs/tui/src/slash_command.rs:93-93`, `codex-rs/tui/src/slash_command.rs:171-171`, `codex-rs/tui/src/slash_command.rs:202-202`, `codex-rs/tui/src/app/event_dispatch.rs:351-366`) |
| `/mention` | insert file mention | No | Yes | Yes | Pairs with file popup. (`codex-rs/tui/src/slash_command.rs:46`, `codex-rs/tui/src/slash_command.rs:94-94`, `codex-rs/tui/src/slash_command.rs:172-172`, `codex-rs/tui/src/slash_command.rs:206-206`) |
| `/status` | session-status view | No | Yes | Yes | Side-convo-safe. (`codex-rs/tui/src/slash_command.rs:47`, `codex-rs/tui/src/slash_command.rs:97-97`, `codex-rs/tui/src/slash_command.rs:173-173`, `codex-rs/tui/src/slash_command.rs:209-209`) |
| `/debug-config` | debug config stack | No | Yes | No | Debug surface. (`codex-rs/tui/src/slash_command.rs:48`, `codex-rs/tui/src/slash_command.rs:98-98`, `codex-rs/tui/src/slash_command.rs:210-210`) |
| `/title` | terminal-title setup | No | Yes | No | TUI config surface. (`codex-rs/tui/src/slash_command.rs:49`, `codex-rs/tui/src/slash_command.rs:99-99`, `codex-rs/tui/src/slash_command.rs:217-217`) |
| `/statusline` | status-line setup | No | Yes | No | TUI config surface. (`codex-rs/tui/src/slash_command.rs:50`, `codex-rs/tui/src/slash_command.rs:100-100`, `codex-rs/tui/src/slash_command.rs:218-218`) |
| `/theme` | theme picker | No | No | No | Live-preview dialog. (`codex-rs/tui/src/slash_command.rs:51`, `codex-rs/tui/src/slash_command.rs:101-101`, `codex-rs/tui/src/slash_command.rs:231-231`) |
| `/mcp` | MCP inventory | Yes | Yes | No | Supports verbose flavor. (`codex-rs/tui/src/slash_command.rs:52`, `codex-rs/tui/src/slash_command.rs:131-131`, `codex-rs/tui/src/slash_command.rs:157-158`, `codex-rs/tui/src/slash_command.rs:214-214`) |
| `/apps` | apps manager | No | Yes | No | Connector-facing UI. (`codex-rs/tui/src/slash_command.rs:53`, `codex-rs/tui/src/slash_command.rs:132-132`, `codex-rs/tui/src/slash_command.rs:215-215`) |
| `/plugins` | plugins browser | No | Yes | No | Plugin-facing UI. (`codex-rs/tui/src/slash_command.rs:54`, `codex-rs/tui/src/slash_command.rs:133-133`, `codex-rs/tui/src/slash_command.rs:216-216`) |
| `/logout` | log-out action | No | No | No | Session control. (`codex-rs/tui/src/slash_command.rs:55`, `codex-rs/tui/src/slash_command.rs:134-134`, `codex-rs/tui/src/slash_command.rs:199-199`) |
| `/quit` | exit action | No | Yes | No | Same core behavior as `/exit`. (`codex-rs/tui/src/slash_command.rs:56`, `codex-rs/tui/src/slash_command.rs:90-90`, `codex-rs/tui/src/slash_command.rs:222-223`) |
| `/exit` | exit action | No | Yes | No | Same core behavior as `/quit`. (`codex-rs/tui/src/slash_command.rs:57`, `codex-rs/tui/src/slash_command.rs:90-90`, `codex-rs/tui/src/slash_command.rs:222-223`) |
| `/feedback` | send logs/support package | No | Yes | No | Support surface. (`codex-rs/tui/src/slash_command.rs:58`, `codex-rs/tui/src/slash_command.rs:81-81`, `codex-rs/tui/src/slash_command.rs:220-220`) |
| `/rollout` | print rollout path | No | Yes | No | Debug-only visibility. (`codex-rs/tui/src/slash_command.rs:59`, `codex-rs/tui/src/slash_command.rs:135-135`, `codex-rs/tui/src/slash_command.rs:225-225`, `codex-rs/tui/src/slash_command.rs:239-239`) |
| `/ps` | list background terminals | No | Yes | No | Background terminal manager. (`codex-rs/tui/src/slash_command.rs:60`, `codex-rs/tui/src/slash_command.rs:102-103`, `codex-rs/tui/src/slash_command.rs:211-212`) |
| `/stop` | stop background terminals | No | Yes | No | Canonical name; `/clean` aliases it. (`codex-rs/tui/src/slash_command.rs:61-62`, `codex-rs/tui/src/slash_command.rs:103-103`, `codex-rs/tui/src/slash_command.rs:212-212`, `codex-rs/tui/src/slash_command.rs:261-268`) |
| `/clean` | alias of `/stop` | No | Yes | No | Parser alias only. (`codex-rs/tui/src/slash_command.rs:61-62`, `codex-rs/tui/src/slash_command.rs:261-268`) |
| `/clear` | clear UI and begin new chat | No | No | No | Conversation reset. (`codex-rs/tui/src/slash_command.rs:63`, `codex-rs/tui/src/slash_command.rs:88-88`, `codex-rs/tui/src/slash_command.rs:198-198`) |
| `/personality` | communication-style picker | No | No | No | Model-style config. (`codex-rs/tui/src/slash_command.rs:64`, `codex-rs/tui/src/slash_command.rs:113-113`, `codex-rs/tui/src/slash_command.rs:188-188`) |
| `/realtime` | realtime-voice toggle | No | Yes | No | Experimental audio path. (`codex-rs/tui/src/slash_command.rs:65`, `codex-rs/tui/src/slash_command.rs:114-114`, `codex-rs/tui/src/slash_command.rs:227-227`) |
| `/settings` | audio-device settings | No | Yes | No | Realtime config. (`codex-rs/tui/src/slash_command.rs:66`, `codex-rs/tui/src/slash_command.rs:115-115`, `codex-rs/tui/src/slash_command.rs:228-228`) |
| `/test-approval` | open approval test path | No | Yes | No | Debug-only visibility. (`codex-rs/tui/src/slash_command.rs:67`, `codex-rs/tui/src/slash_command.rs:136-136`, `codex-rs/tui/src/slash_command.rs:226-226`, `codex-rs/tui/src/slash_command.rs:239-239`) |
| `/subagents` | switch agent thread | No | Yes | No | Serialized command name for multi-agent switching. (`codex-rs/tui/src/slash_command.rs:68-69`, `codex-rs/tui/src/slash_command.rs:119-120`, `codex-rs/tui/src/slash_command.rs:230-230`) |
| `/debug-m-drop` | memory debug helper | No | No | No | Explicitly marked internal-only. (`codex-rs/tui/src/slash_command.rs:71-72`, `codex-rs/tui/src/slash_command.rs:104-105`, `codex-rs/tui/src/slash_command.rs:200-201`) |
| `/debug-m-update` | memory debug helper | No | No | No | Explicitly marked internal-only. (`codex-rs/tui/src/slash_command.rs:73-74`, `codex-rs/tui/src/slash_command.rs:105-105`, `codex-rs/tui/src/slash_command.rs:200-201`) |

### 5.2 Slash-command availability heuristics worth copying
- Codex does not just hide unavailable commands after execution; it proactively models three dimensions: visible at all, safe during a running task, and safe in side conversation. (`codex-rs/tui/src/slash_command.rs:165-241`)
- That is a strong pattern for chemsmart because chemistry tasks have many phase-specific actions: e.g. `/submit`, `/resume-job`, `/inspect-geom`, `/irc-template`, `/scan-setup` should not all be visible or runnable in every state. (chemsmart grounding: `/Users/hongjiseung/developer/chemsmart/AGENTS.md:74-77`, `/Users/hongjiseung/developer/chemsmart/AGENTS.md:193-212`, `/Users/hongjiseung/developer/chemsmart/AGENTS.md:226-237`)

### 5.3 Keybindings: app-wide shortcuts
| Category | Action | Default | Notes | Evidence |
| App | Open transcript | `Ctrl+T` | Opens alt-screen transcript pager. | `codex-rs/tui/src/keymap.rs:549-557`, `codex-rs/tui/src/app/input.rs:166-174` |
| App | Open external editor | `Ctrl+G` | Temporarily restores terminal state before launching editor. | `codex-rs/tui/src/keymap.rs:550-550`, `codex-rs/tui/src/tui.rs:528-565` |
| App | Copy last response | `Ctrl+O` | Copies last assistant markdown. | `codex-rs/tui/src/keymap.rs:551-551`, `codex-rs/tui/src/chatwidget.rs:5144-5151` |
| App | Clear terminal | `Ctrl+L` | TUI-wide terminal clear action. | `codex-rs/tui/src/keymap.rs:552-553` |
| App | Toggle Vim mode | _unbound_ | Available for remapping, no default chord. | `codex-rs/tui/src/keymap.rs:554-554` |
| App | Toggle Fast mode | _unbound_ | Available for remapping, no default chord. | `codex-rs/tui/src/keymap.rs:555-555` |
| App | Toggle raw output | `Alt+R` | Same behavior as `/raw`. | `codex-rs/tui/src/keymap.rs:556-556`, `codex-rs/tui/src/chatwidget/slash_dispatch.rs:324-327` |

### 5.4 Keybindings: chat/composer shortcuts
| Category | Action | Default | Notes | Evidence |
| Chat | Decrease reasoning | `Alt+,` | Adjusts reasoning level. | `codex-rs/tui/src/keymap.rs:559-560` |
| Chat | Increase reasoning | `Alt+.` | Adjusts reasoning level. | `codex-rs/tui/src/keymap.rs:559-560` |
| Chat | Edit queued message | `Alt+Up`, `Shift+Left` | Terminal-specific fallback uses `Shift+Left` in Apple Terminal/Warp/VSCode/tmux. | `codex-rs/tui/src/keymap.rs:561-561`, `codex-rs/tui/src/chatwidget.rs:208-252` |
| Composer | Submit | `Enter` | Immediate submit. | `codex-rs/tui/src/keymap.rs:564-564`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:41-42` |
| Composer | Queue | `Tab` | Queues while task running; otherwise submits. | `codex-rs/tui/src/keymap.rs:565-565`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:41-43` |
| Composer | Toggle shortcuts overlay | `?`, `Shift+?` | Compatibility variants included intentionally. | `codex-rs/tui/src/keymap.rs:566-569`, `codex-rs/tui/src/keymap.rs:544-546` |
| Composer | History search previous | `Ctrl+R` | Reverse incremental search mode. | `codex-rs/tui/src/keymap.rs:570-571`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:31-33` |
| Composer | History search next | `Ctrl+S` | Forward match navigation. | `codex-rs/tui/src/keymap.rs:570-571` |

### 5.5 Keybindings: editor mode
| Action | Default | Notes | Evidence |
| Insert newline | `Ctrl+J`, `Ctrl+M`, `Enter`, `Shift+Enter`, `Alt+Enter` | Modified Enter requires keyboard-enhancement flags. | `codex-rs/tui/src/keymap.rs:574-580`, `codex-rs/tui/src/tui.rs:160-167` |
| Move left | `Left`, `Ctrl+B` | Shell-like. | `codex-rs/tui/src/keymap.rs:581-581` |
| Move right | `Right`, `Ctrl+F` | Shell-like. | `codex-rs/tui/src/keymap.rs:582-582` |
| Move up | `Up`, `Ctrl+P` | History/cursor aware. | `codex-rs/tui/src/keymap.rs:583-583` |
| Move down | `Down`, `Ctrl+N` | History/cursor aware. | `codex-rs/tui/src/keymap.rs:584-584` |
| Move word left | `Alt+B`, `Alt+Left`, `Ctrl+Left` | Multiple terminal encodings supported. | `codex-rs/tui/src/keymap.rs:585-589` |
| Move word right | `Alt+F`, `Alt+Right`, `Ctrl+Right` | Multiple terminal encodings supported. | `codex-rs/tui/src/keymap.rs:590-594` |
| Move line start | `Home`, `Ctrl+A` | Shell-like. | `codex-rs/tui/src/keymap.rs:595-595` |
| Move line end | `End`, `Ctrl+E` | Shell-like. | `codex-rs/tui/src/keymap.rs:596-596` |
| Delete backward | `Backspace`, `Shift+Backspace`, `Ctrl+H` | Terminal compatibility variants. | `codex-rs/tui/src/keymap.rs:597-601` |
| Delete forward | `Delete`, `Shift+Delete`, `Ctrl+D` | Only editor-local; app-level `Ctrl+D` quit logic is gated by emptiness. | `codex-rs/tui/src/keymap.rs:602-606`, `codex-rs/tui/src/chatwidget.rs:10265-10294` |
| Delete backward word | `Alt+Backspace`, `Ctrl+Backspace`, `Ctrl+Shift+Backspace`, `Ctrl+W`, `Ctrl+Alt+H` | Aggressive compatibility set. | `codex-rs/tui/src/keymap.rs:607-619` |
| Delete forward word | `Alt+Delete`, `Ctrl+Delete`, `Ctrl+Shift+Delete`, `Alt+D` | Aggressive compatibility set. | `codex-rs/tui/src/keymap.rs:620-628` |
| Kill line start | `Ctrl+U` | Emacs-like kill buffer behavior preserved across submit. | `codex-rs/tui/src/keymap.rs:629-629`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:52-55` |
| Kill whole line | _unbound_ | Remappable, no default. | `codex-rs/tui/src/keymap.rs:630-630` |
| Kill line end | `Ctrl+K` | Emacs-like kill buffer behavior preserved across submit. | `codex-rs/tui/src/keymap.rs:631-631`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:52-55` |
| Yank | `Ctrl+Y` | Restores killed text. | `codex-rs/tui/src/keymap.rs:632-632`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:52-55` |

### 5.6 Keybindings: Vim normal mode
| Action | Default | Evidence |
| Enter insert | `i`, `Insert` | `codex-rs/tui/src/keymap.rs:634-635` |
| Append after cursor | `a` | `codex-rs/tui/src/keymap.rs:636-636` |
| Append line end | `A`, `Shift+A` | `codex-rs/tui/src/keymap.rs:637-640` |
| Insert line start | `I`, `Shift+I` | `codex-rs/tui/src/keymap.rs:641-644` |
| Open line below | `o` | `codex-rs/tui/src/keymap.rs:645-645` |
| Open line above | `O`, `Shift+O` | `codex-rs/tui/src/keymap.rs:646-649` |
| Move left | `h`, `Left` | `codex-rs/tui/src/keymap.rs:650-650` |
| Move right | `l`, `Right` | `codex-rs/tui/src/keymap.rs:651-651` |
| Move up | `k`, `Up` | `codex-rs/tui/src/keymap.rs:652-652` |
| Move down | `j`, `Down` | `codex-rs/tui/src/keymap.rs:653-653` |
| Word forward | `w` | `codex-rs/tui/src/keymap.rs:654-654` |
| Word backward | `b` | `codex-rs/tui/src/keymap.rs:655-655` |
| Word end | `e` | `codex-rs/tui/src/keymap.rs:656-656` |
| Line start | `0` | `codex-rs/tui/src/keymap.rs:657-657` |
| Line end | `$`, `Shift+$` | `codex-rs/tui/src/keymap.rs:658-661` |
| Delete char | `x` | `codex-rs/tui/src/keymap.rs:662-662` |
| Delete to line end | `D`, `Shift+D` | `codex-rs/tui/src/keymap.rs:663-666` |
| Yank line | `Y`, `Shift+Y` | `codex-rs/tui/src/keymap.rs:667-667` |
| Paste after | `p` | `codex-rs/tui/src/keymap.rs:668-668` |
| Start delete operator | `d` | `codex-rs/tui/src/keymap.rs:669-669` |
| Start yank operator | `y` | `codex-rs/tui/src/keymap.rs:670-670` |
| Cancel operator | `Esc` | `codex-rs/tui/src/keymap.rs:671-671` |

### 5.7 Keybindings: Vim operator-pending mode
| Action | Default | Evidence |
| Delete line | `d` | `codex-rs/tui/src/keymap.rs:673-675` |
| Yank line | `y` | `codex-rs/tui/src/keymap.rs:675-675` |
| Motion left | `h` | `codex-rs/tui/src/keymap.rs:676-676` |
| Motion right | `l` | `codex-rs/tui/src/keymap.rs:677-677` |
| Motion up | `k` | `codex-rs/tui/src/keymap.rs:678-678` |
| Motion down | `j` | `codex-rs/tui/src/keymap.rs:679-679` |
| Motion word forward | `w` | `codex-rs/tui/src/keymap.rs:680-680` |
| Motion word backward | `b` | `codex-rs/tui/src/keymap.rs:681-681` |
| Motion word end | `e` | `codex-rs/tui/src/keymap.rs:682-682` |
| Motion line start | `0` | `codex-rs/tui/src/keymap.rs:683-683` |
| Motion line end | `$`, `Shift+$` | `codex-rs/tui/src/keymap.rs:684-687` |
| Cancel | `Esc` | `codex-rs/tui/src/keymap.rs:688-688` |

### 5.8 Keybindings: pager, list, and approval surfaces
| Surface | Action | Default | Evidence |
| Pager | Scroll up | `Up`, `k` | `codex-rs/tui/src/keymap.rs:690-692` |
| Pager | Scroll down | `Down`, `j` | `codex-rs/tui/src/keymap.rs:691-692` |
| Pager | Page up | `PageUp`, `Shift+Space`, `Ctrl+B` | `codex-rs/tui/src/keymap.rs:693-697` |
| Pager | Page down | `PageDown`, `Space`, `Ctrl+F` | `codex-rs/tui/src/keymap.rs:698-702` |
| Pager | Half-page up | `Ctrl+U` | `codex-rs/tui/src/keymap.rs:703-703` |
| Pager | Half-page down | `Ctrl+D` | `codex-rs/tui/src/keymap.rs:704-704` |
| Pager | Jump top | `Home` | `codex-rs/tui/src/keymap.rs:705-705` |
| Pager | Jump bottom | `End` | `codex-rs/tui/src/keymap.rs:706-706` |
| Pager | Close | `q`, `Ctrl+C` | `codex-rs/tui/src/keymap.rs:707-707` |
| Pager | Close transcript | `Ctrl+T` | `codex-rs/tui/src/keymap.rs:708-708` |
| List | Move up | `Up`, `Ctrl+P`, `k` | `codex-rs/tui/src/keymap.rs:710-715` |
| List | Move down | `Down`, `Ctrl+N`, `j` | `codex-rs/tui/src/keymap.rs:716-720` |
| List | Accept | `Enter` | `codex-rs/tui/src/keymap.rs:721-721` |
| List | Cancel | `Esc` | `codex-rs/tui/src/keymap.rs:722-722` |
| Approval | Open fullscreen request | `Ctrl+A`, `Ctrl+Shift+A` | `codex-rs/tui/src/keymap.rs:724-731`, `codex-rs/tui/src/bottom_pane/approval_overlay.rs:508-520` |
| Approval | Open source thread | `o` | `codex-rs/tui/src/keymap.rs:732-732`, `codex-rs/tui/src/bottom_pane/approval_overlay.rs:522-530` |
| Approval | Approve once | `y` | `codex-rs/tui/src/keymap.rs:733-733` |
| Approval | Approve for session | `a` | `codex-rs/tui/src/keymap.rs:734-734` |
| Approval | Approve for prefix / host | `p` | `codex-rs/tui/src/keymap.rs:735-735` |
| Approval | Deny | `d` | `codex-rs/tui/src/keymap.rs:736-736` |
| Approval | Decline / cancel-and-steer | `Esc`, `n` | `codex-rs/tui/src/keymap.rs:737-737` |
| Approval | Cancel | `c` | `codex-rs/tui/src/keymap.rs:738-738` |

## 6. Tool-call & diff display: how file edits, shell commands, tool outputs render; approval/permission UX

- Codex has at least three distinct execution/result renderers: shell/exec cells, MCP tool-call cells, and diff/file-change renderers. They do not share one generic “tool output block.” (`codex-rs/tui/src/exec_cell/render.rs:195-246`, `codex-rs/tui/src/history_cell.rs:1868-1962`, `codex-rs/tui/src/diff_render.rs:1-32`)
- Shell command transcript rendering uses syntax-highlighted command lines prefixed with `$`, then raw/wrapped output, then a result footer with ✓/✗ and duration. (`codex-rs/tui/src/exec_cell/render.rs:204-246`)
- Compact inline shell output is aggressively truncated: default tool-output previews keep only a head/tail window, insert an ellipsis, and point users to `Ctrl+T` for the full transcript. (`codex-rs/tui/src/exec_cell/render.rs:32-35`, `codex-rs/tui/src/exec_cell/render.rs:103-183`, `codex-rs/tui/src/exec_cell/render.rs:253-260`)
- Read/list/search shell activity can collapse into an “Exploring/Explored” presentation rather than raw shell transcripts. (`codex-rs/tui/src/exec_cell/render.rs:262-359`)
- MCP tool calls use a different visual grammar: bullet or spinner, `Calling` vs `Called`, inline invocation if it fits, otherwise wrapped invocation under a `└` tree branch, then dim result blocks. (`codex-rs/tui/src/history_cell.rs:1868-1962`)
- Diff rendering is sophisticated: per-file or multi-file summaries, rename arrows, line numbers, gutter signs, optional syntax highlighting, theme-aware add/delete backgrounds, hunk-level parser-state preservation, and `⋮` separators between hunks. (`codex-rs/tui/src/diff_render.rs:1-32`, `codex-rs/tui/src/diff_render.rs:403-465`, `codex-rs/tui/src/diff_render.rs:548-734`)
- Display paths are normalized relative to cwd/repo/home to avoid ugly absolute paths in diffs and image tool calls. (`codex-rs/tui/src/diff_render.rs:739-763`)
- Approvals are bottom-pane modal list views with keyboard shortcuts, source-thread switching, optional fullscreen detail, and separate option sets for exec, patch, permissions, and elicitation. (`codex-rs/tui/src/bottom_pane/approval_overlay.rs:155-167`, `codex-rs/tui/src/bottom_pane/approval_overlay.rs:230-298`, `codex-rs/tui/src/bottom_pane/approval_overlay.rs:508-547`)

### 6.1 Shell / exec rendering
| UI element | Observed behavior | Evidence |
| Command line | Syntax-highlighted bash, transcript uses `$ ` initial indent and `    ` continuation indent. | `codex-rs/tui/src/exec_cell/render.rs:210-218` |
| Output preview | Dimmed, wrapped, head/tail truncation when long. | `codex-rs/tui/src/exec_cell/render.rs:103-183` |
| Ellipsis hint | Shows omitted line count plus transcript shortcut hint. | `codex-rs/tui/src/exec_cell/render.rs:253-260` |
| Completion footer | Green `✓` or red `✗ (exit)` plus dim duration. | `codex-rs/tui/src/exec_cell/render.rs:229-242` |
| Read/list/search grouping | Turns multiple low-risk commands into human-semantic “Read/List/Search” blocks. | `codex-rs/tui/src/exec_cell/render.rs:262-359` |

### 6.2 MCP tool-call rendering
| UI element | Observed behavior | Evidence |
| Status marker | Green bullet for success, red for failure, animated indicator while running. | `codex-rs/tui/src/history_cell.rs:1871-1881` |
| Verb | Uses `Calling` while active and `Called` after completion. | `codex-rs/tui/src/history_cell.rs:1882-1886` |
| Invocation layout | Inline if it fits; otherwise wraps under a `└` prefixed child line. | `codex-rs/tui/src/history_cell.rs:1888-1909` |
| Result rendering | Tool result blocks are dimmed and wrapped line by line. | `codex-rs/tui/src/history_cell.rs:1915-1948` |
| Interrupt state | Interrupted calls become failed with error text `interrupted`. | `codex-rs/tui/src/history_cell.rs:1832-1835` |

### 6.3 Diff rendering
| UI element | Observed behavior | Evidence |
| Header summary | Single-file case says Added/Deleted/Edited + path + +/- counts; multi-file case summarizes file count + totals. | `codex-rs/tui/src/diff_render.rs:415-437` |
| Per-file chunk header | Multi-file diffs use `  └ ` headers with rename arrows. | `codex-rs/tui/src/diff_render.rs:439-460` |
| Add/Delete files | Highlights whole file by extension when possible. | `codex-rs/tui/src/diff_render.rs:475-547` |
| Update hunks | Parses unified diff, preserves parser state within hunk, inserts `⋮` between hunks. | `codex-rs/tui/src/diff_render.rs:548-734` |
| Large diffs | Skips syntax highlighting when patch size exceeds guardrails. | `codex-rs/tui/src/diff_render.rs:582-589`, `codex-rs/tui/src/render/highlight.rs:20-22` |
| Path normalization | Prefers cwd-relative or repo-relative display. | `codex-rs/tui/src/diff_render.rs:739-763` |

### 6.4 Approval / permission UX
| Request type | Header/title behavior | Choice design | Evidence |
| Exec approval | Title asks whether to run command or approve host access; header may include thread, reason, permission rule, and highlighted command. | Choices include once/session/prefix-or-host/deny/cancel-and-steer. | `codex-rs/tui/src/bottom_pane/approval_overlay.rs:237-258`, `codex-rs/tui/src/bottom_pane/approval_overlay.rs:626-666`, `codex-rs/tui/src/bottom_pane/approval_overlay.rs:790-875` |
| Permissions approval | Title asks whether to grant requested permissions; header can show permission-rule summary. | Choices include turn-only, turn+strict-auto-review, session, deny. | `codex-rs/tui/src/bottom_pane/approval_overlay.rs:260-263`, `codex-rs/tui/src/bottom_pane/approval_overlay.rs:667-692`, `codex-rs/tui/src/bottom_pane/approval_overlay.rs:993-1025` |
| Patch approval | Title asks whether to make edits; header can show thread + reason. | Choices include once/session/cancel-and-steer. | `codex-rs/tui/src/bottom_pane/approval_overlay.rs:264-267`, `codex-rs/tui/src/bottom_pane/approval_overlay.rs:693-720`, `codex-rs/tui/src/bottom_pane/approval_overlay.rs:973-990` |
| Elicitation / MCP form | Server-specific approval title. | Stable Esc-as-cancel semantics. | `codex-rs/tui/src/bottom_pane/approval_overlay.rs:268-271`, `codex-rs/tui/src/bottom_pane/approval_overlay.rs:1027-1035` |

### 6.5 Short approval snippet
```rust
"Would you like to make the following edits?"
```
That literal is emblematic: Codex frames patch approval in user language, not protocol language. (`codex-rs/tui/src/bottom_pane/approval_overlay.rs:264-267`)

### 6.6 What chemsmart should especially notice
- Chemsmart’s dry-run and submit flows would benefit from **two separate surfaces** the way Codex separates shell output from patch approvals: one surface for “here is the generated input / route / scheduler script,” another for “may I submit or overwrite?”. Chemsmart already distinguishes dry-run inputs and critic validation in PR #16. (`https://github.com/Hongjiseung-ROK/chemsmart/pull/16`)
- PR #15’s optimized-geometry handoff suggests a chemistry-specific “artifact block” renderer: optimized geometry, extracted structure label, source job, and downstream consumer, analogous to Codex’s structured MCP/diff cells. (`https://github.com/Hongjiseung-ROK/chemsmart/pull/15`)

## 7. Color / theme / typography, borders, frames, syntax highlighting

- Codex’s visual system is intentionally terminal-adaptive, not hardcoded to one palette. The syntax highlighter can choose from 32 bundled themes and uses terminal background lightness to pick a default (`catppuccin-latte` on light, `catppuccin-mocha` otherwise). (`codex-rs/tui/src/render/highlight.rs:1-22`, `codex-rs/tui/src/render/highlight.rs:135-171`, `codex-rs/tui/src/render/highlight.rs:184-201`)
- Theme selection is live-previewed as the cursor moves and restored on cancel; persisted theme choice writes `[tui] theme = "..."` into config on confirm. (`codex-rs/tui/src/theme_picker.rs:1-20`, `codex-rs/tui/src/theme_picker.rs:314-390`)
- Diff styling separately adapts to dark/light terminals and color depth, with explicit truecolor and ANSI-256 palettes plus foreground-only ANSI-16 fallback. (`codex-rs/tui/src/diff_render.rs:10-21`, `codex-rs/tui/src/diff_render.rs:60-76`, `codex-rs/tui/src/diff_render.rs:251-268`)
- Status-line coloring is semantic rather than positional: model/path/branch/state/usage/limit/mode/thread/progress each map to syntax-theme scopes with fallback colors. (`codex-rs/tui/src/bottom_pane/status_line_style.rs:16-75`, `codex-rs/tui/src/bottom_pane/status_line_style.rs:77-119`)
- Status-line colors are softened rather than simply dimmed, which helps keep rich themes legible but not overwhelming. (`codex-rs/tui/src/bottom_pane/status_line_style.rs:122-173`)
- Typography is plain terminal text plus `Stylize` helpers; the repo’s own TUI style guide explicitly prefers concise ratatui styling (`.dim()`, `.bold()`, `.cyan()`, etc.) over verbose manual styles. (`AGENTS.md:81-108`)
- Rounded borders appear in onboarding/auth input widgets, showing that Codex is willing to depart from minimal line-only terminal UI when a focused form benefits from framing. (`codex-rs/tui/src/onboarding/auth.rs:33-35`, `codex-rs/tui/src/onboarding/auth.rs:677-683`)
- Animation is pluggable: Codex embeds multiple frame sets (`default`, `codex`, `openai`, `blocks`, `dots`, etc.), each with 36 frames and an 80 ms tick default. (`codex-rs/tui/src/frames.rs:3-71`)

### 7.1 Theme / palette details
| Theme concern | Implementation detail | Evidence |
| Bundled theme count | 32 bundled themes plus custom `.tmTheme` discovery under `{CODEX_HOME}/themes/`. | `codex-rs/tui/src/render/highlight.rs:1-22`, `codex-rs/tui/src/theme_picker.rs:1-20` |
| Highlight guardrails | Rejects inputs >512 KB or >10,000 lines for syntax highlighting. | `codex-rs/tui/src/render/highlight.rs:20-22` |
| Adaptive default | Chooses light vs dark default from terminal background probe. | `codex-rs/tui/src/render/highlight.rs:184-201`, `codex-rs/tui/src/terminal_palette.rs:63-69` |
| Custom themes | Looks for `{CODEX_HOME}/themes/{name}.tmTheme` and warns on invalid theme names. | `codex-rs/tui/src/render/highlight.rs:103-132`, `codex-rs/tui/src/render/highlight.rs:174-181` |
| Theme picker preview | Wide side panel if >=44 columns; narrow stacked fallback otherwise. | `codex-rs/tui/src/theme_picker.rs:14-20`, `codex-rs/tui/src/theme_picker.rs:129-146` |
| Status line separator | Uses ` · ` between segments and underlines PR number. | `codex-rs/tui/src/bottom_pane/status_line_style.rs:12-15`, `codex-rs/tui/src/bottom_pane/status_line_style.rs:98-117` |

### 7.2 Visual snippet
```rust
BorderType::Rounded
```
This shows up in the onboarding auth form, paired with a cyan border to make the API-key field feel like a focused input box. (`codex-rs/tui/src/onboarding/auth.rs:677-683`)

## 8. Multi-line input, paste handling, image input, file autocomplete

- Codex’s composer is one of the richest parts of the TUI. Its own module comment reads like a state-machine spec: input editing, popup routing, slash-command promotion, submit-vs-newline behavior, history merging, large paste placeholders, remote-image rows, and non-bracketed paste burst handling. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:1-133`)
- Enter submits by default, while multiline editing comes from modified Enter/newline bindings and the raw editor keymap. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:41-43`, `codex-rs/tui/src/keymap.rs:574-580`)
- Large pastes do not flood the composer. If a paste exceeds the threshold, Codex inserts a placeholder like `[Pasted Content N chars]`, stores the full text off-buffer, and expands it on submit. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:62-71`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:894-929`)
- Paste handling also normalizes CRLF, detects image-file paths, and can attach the image instead of pasting literal path text. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:912-930`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:932-952`)
- Codex explicitly compensates for terminals that do not emit reliable bracketed-paste events, especially on Windows, by detecting “paste bursts” in rapid key streams and converting them into synthetic paste events. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:92-127`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:955-981`)
- Local and remote image attachments are first-class. Remote image URLs render as non-editable `[Image #N]` rows above the textarea; local image placeholders are numbered in the same sequence. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:73-90`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:3034-3066`)
- Image paste from the system clipboard is bound to `Ctrl+V`/`Alt+V` chord logic at the chat widget level, which writes a temp PNG and attaches it. (`codex-rs/tui/src/chatwidget.rs:5177-5202`, `codex-rs/tui/src/clipboard_paste.rs:49-149`)
- File autocomplete is driven by `@token` detection. The file popup has waiting/idle/no-matches states, uses async search results, and on selection can either insert a quoted path or attach an image if the selected file is an image. (`codex-rs/tui/src/bottom_pane/file_search_popup.rs:16-153`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:1969-2058`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:2407-2451`)
- When the cursor is inside an `@token`, file search takes precedence over slash-command popup so users can say things like `/review @docs/...`. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:3807-3815`)

### 8.1 Input behavior table
| Capability | Observed behavior | Evidence |
| Immediate submit | `Enter` submits; no shell-like “always newline” default. | `codex-rs/tui/src/bottom_pane/chat_composer.rs:41-42`, `codex-rs/tui/src/keymap.rs:564-565` |
| Multiline input | Modified Enter/newline bindings stay in editor and require enhanced keyboard flags. | `codex-rs/tui/src/keymap.rs:574-580`, `codex-rs/tui/src/tui.rs:160-167` |
| History recall | Merges persistent text-only history with local in-session history that can rehydrate elements and image attachments. | `codex-rs/tui/src/bottom_pane/chat_composer.rs:19-37` |
| Large paste | Placeholder token stored in `pending_pastes`, expanded later. | `codex-rs/tui/src/bottom_pane/chat_composer.rs:62-71`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:904-918` |
| Image-path paste | Path-like pastes can attach image and insert trailing space. | `codex-rs/tui/src/bottom_pane/chat_composer.rs:906-908`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:919-924`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:932-952` |
| Paste-burst fallback | Rapid key bursts become explicit paste handling instead of literal typing. | `codex-rs/tui/src/bottom_pane/chat_composer.rs:92-127`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:955-981` |
| Remote image rows | Keyboard-selectable/removable rows above textarea. | `codex-rs/tui/src/bottom_pane/chat_composer.rs:73-90`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:3034-3066` |
| File autocomplete | `@token` popup with async file search, selection, quoting if whitespace. | `codex-rs/tui/src/bottom_pane/file_search_popup.rs:16-153`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:2407-2451` |
| Slash popup precedence | Slash popup is suppressed when cursor is on `@token`. | `codex-rs/tui/src/bottom_pane/chat_composer.rs:3807-3815` |
| Mention popup | `$token` mention path exists separately from `@path` file popup. | `codex-rs/tui/src/bottom_pane/chat_composer.rs:2400-2405` |

### 8.2 Small snippet: large paste is treated as an object
```rust
self.pending_pastes.push((placeholder, pasted));
```
That line captures the philosophy: Codex does not force the textarea to carry every pasted byte immediately. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:915-918`)

### 8.3 Chemsmart-specific relevance
- This placeholder approach is a very strong fit for chemsmart dry-run previews, basis-set blocks, SCAN/IRC coordinate specs, and long scheduler scripts. Users often want to paste or generate long structured text without turning the composer into a wall of unreadable input. (chemsmart grounding: `/Users/hongjiseung/developer/chemsmart/AGENTS.md:23-35`, `/Users/hongjiseung/developer/chemsmart/AGENTS.md:193-212`, `https://github.com/Hongjiseung-ROK/chemsmart/pull/23`)
- Image attachment is less central for chemsmart than file attachment, but the same UI slot could hold structure previews, rendered conformer snapshots, or uploaded scan diagrams. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:73-90`, `/Users/hongjiseung/developer/chemsmart/AGENTS.md:290-299`)

## 9. State persistence: session history, resume, transcript export

- Codex supports persistence at multiple levels: thread metadata in app-server/state, rollout JSONL fallback, opt-in TUI session logging, and transcript overlays/history search inside the active client. (`codex-rs/tui/src/session_resume.rs:1-5`, `codex-rs/tui/src/session_log.rs:80-118`, `codex-rs/tui/src/resume_picker.rs:286-297`)
- Opt-in session logging is controlled by `CODEX_TUI_RECORD_SESSION`; default output is `session-<timestamp>.jsonl`, and records include session start metadata, inbound app events, outbound ops, and session end. (`codex-rs/tui/src/session_log.rs:80-118`, `codex-rs/tui/src/session_log.rs:121-215`)
- Resume/fork can prompt when the saved session cwd differs from the current cwd, which is critical for agent tools that depend on repo-relative paths. (`codex-rs/tui/src/session_resume.rs:86-111`, `codex-rs/tui/src/session_resume.rs:140-142`)
- If app-server metadata is unavailable, Codex parses rollout JSONL, prefers latest `turn_context`, and recovers thread ID, cwd, and model from there. (`codex-rs/tui/src/session_resume.rs:54-64`, `codex-rs/tui/src/session_resume.rs:144-189`)
- The resume picker lazily loads transcript previews rather than eagerly rendering every session body. (`codex-rs/tui/src/resume_picker.rs:286-297`, `codex-rs/tui/src/resume_picker/transcript.rs:24-107`)
- Transcript viewing is a first-class overlay triggered by `Ctrl+T`. (`codex-rs/tui/src/app/input.rs:166-174`, `codex-rs/tui/src/pager_overlay.rs:1-16`)
- I do **not** see a first-class built-in “export transcript” slash command. The practical export-ish surfaces are session JSONL logging, transcript overlay viewing, and `/copy` for the last assistant response. That is an inference from the visible slash-command list and persistence modules. (`codex-rs/tui/src/slash_command.rs:12-75`, `codex-rs/tui/src/session_log.rs:80-118`, `codex-rs/tui/src/pager_overlay.rs:1-16`)

### 9.1 Persistence surfaces
| Surface | What is persisted | Evidence |
| Session JSONL log | TUI app events + outbound ops + metadata | `codex-rs/tui/src/session_log.rs:80-118`, `codex-rs/tui/src/session_log.rs:121-215` |
| Thread metadata | Saved cwd/model via state runtime | `codex-rs/tui/src/session_resume.rs:67-84`, `codex-rs/tui/src/session_resume.rs:114-123` |
| Rollout fallback | Thread id, cwd, model from rollout records | `codex-rs/tui/src/session_resume.rs:54-64`, `codex-rs/tui/src/session_resume.rs:164-179` |
| Composer history | Persistent cross-session text history + local in-session full draft history | `codex-rs/tui/src/bottom_pane/chat_composer.rs:19-37` |
| Transcript overlay cache | Committed cells + live-tail cache key | `codex-rs/tui/src/pager_overlay.rs:6-16`, `codex-rs/tui/src/pager_overlay.rs:423-435` |
| Resume picker transcript preview | Lazy transcript cells on demand | `codex-rs/tui/src/resume_picker.rs:286-297`, `codex-rs/tui/src/resume_picker/transcript.rs:24-107` |

### 9.2 Chemsmart relevance
- Chemsmart should copy the cwd-mismatch prompt almost verbatim in spirit. Job submission, output parsing, and geometry handoff are path-sensitive, and the repo explicitly warns that many workflows depend on `~/.chemsmart` plus local environment assumptions. (`/Users/hongjiseung/developer/chemsmart/AGENTS.md:115-137`, `codex-rs/tui/src/session_resume.rs:86-111`)
- A chemistry TUI also needs persistence of **artifacts**, not just chat: last generated input files, extracted optimized geometries, scheduler target, and whether a draft is dry-run-only. Codex’s session-log/rollout layering is a good starting pattern. (`codex-rs/tui/src/session_log.rs:80-118`, `https://github.com/Hongjiseung-ROK/chemsmart/pull/15`, `https://github.com/Hongjiseung-ROK/chemsmart/pull/16`)

## 10. Notable UX details worth borrowing for chemsmart (chemistry domain: jobs, geometry handoff, dry-run inputs, HPC submission, opt+freq/IRC/scan)

- **Inline-first scrollback is a great default for chemistry.** Codex’s main viewport stays in normal scrollback, which would let chemsmart users keep route lines, generated inputs, and scheduler output in the terminal history while still using overlays for focused views. (`codex-rs/tui/src/tui.rs:329-369`)
- **Model “artifacts” as first-class cells.** Codex has distinct cells for exec, MCP, plans, reasoning, and diffs; chemsmart should likewise have cells for generated `.com/.inp` files, extracted geometries, scheduler scripts, parsed warnings, and dry-run results. (`codex-rs/tui/src/history_cell.rs:1775-1962`, `codex-rs/tui/src/exec_cell/render.rs:195-246`, `https://github.com/Hongjiseung-ROK/chemsmart/pull/15`)
- **Make geometry handoff visible.** PR #15 added `extract_optimized_geometry(job)` precisely because multi-program workflows need explicit handoff. Borrow Codex’s structured tool-call rendering and show: source job → extracted geometry → next job target. (`https://github.com/Hongjiseung-ROK/chemsmart/pull/15`, `codex-rs/tui/src/history_cell.rs:1868-1962`)
- **Separate dry-run preview from submit approval.** PR #16 hardens dry-run validation; chemsmart should mirror Codex by showing generated inputs in transcript-style cells, then a separate approval view for actual local run / HPC submission. (`https://github.com/Hongjiseung-ROK/chemsmart/pull/16`, `codex-rs/tui/src/bottom_pane/approval_overlay.rs:237-298`)
- **Borrow the “cancel and steer” pattern.** Codex approvals include “No, and tell Codex what to do differently,” which is ideal for chemistry cases where the user wants to reject an input file but keep the planning context. (`codex-rs/tui/src/bottom_pane/approval_overlay.rs:869-872`, `codex-rs/tui/src/bottom_pane/approval_overlay.rs:986-989`)
- **Use structured composer placeholders for long generated text.** Large-paste placeholders map naturally to SCAN coordinate specs, IRC keywords, and hand-edited route blocks. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:62-71`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:894-929`, `https://github.com/Hongjiseung-ROK/chemsmart/pull/23`)
- **Gate commands by task phase.** Codex’s slash-command availability model would translate well to chemistry phases: preflight (`/dry-run`, `/basis`, `/charge-spin`), generation (`/inspect-input`), execution (`/submit-hpc`, `/cancel-job`), postprocess (`/extract-geom`, `/thermo`, `/plot-irc`). (`codex-rs/tui/src/slash_command.rs:165-241`, `/Users/hongjiseung/developer/chemsmart/AGENTS.md:193-212`, `/Users/hongjiseung/developer/chemsmart/AGENTS.md:226-237`)
- **Make scan/IRC validation UX explicit.** PR #23 taught the planner to decline scans lacking 1-based atom indices and PR #16 rejects IRC dry-runs missing `irc=`. The TUI should surface these as structured critique/warning cells, not buried stderr. (`https://github.com/Hongjiseung-ROK/chemsmart/pull/16`, `https://github.com/Hongjiseung-ROK/chemsmart/pull/23`)
- **Support focused alt-screen inspectors for transcript, diff, and job artifacts.** For chemsmart that could mean full-screen viewers for route cards, scheduler scripts, optimized structures, vibrational summaries, or IRC trajectories. (`codex-rs/tui/src/pager_overlay.rs:1-16`, `codex-rs/tui/src/app/event_dispatch.rs:351-366`, `/Users/hongjiseung/developer/chemsmart/AGENTS.md:297-299`)
- **Normalize paths and preserve cwd on resume.** Chemical workflows often bounce across worktrees, scratch dirs, and HPC staging directories; copy Codex’s resume cwd prompt and relative-path display rules. (`codex-rs/tui/src/session_resume.rs:86-111`, `codex-rs/tui/src/diff_render.rs:739-763`)
- **Treat planner/critic messages as visible product UX.** PR #17 and PR #23 show that chemistry-specific rationale and reproducibility matter; Codex’s explicit plan/reasoning cells suggest chemsmart should expose “why this route/method/workflow” in dedicated readable cells. (`https://github.com/Hongjiseung-ROK/chemsmart/pull/17`, `https://github.com/Hongjiseung-ROK/chemsmart/pull/23`, `codex-rs/tui/src/resume_picker/transcript.rs:68-93`)
- **HPC submission should look like approval + exec + artifact, not one printout.** Chemsmart’s file-oriented run/sub model and user-local config tree make this especially important. (`/Users/hongjiseung/developer/chemsmart/AGENTS.md:93-113`, `/Users/hongjiseung/developer/chemsmart/AGENTS.md:117-137`, `codex-rs/tui/src/bottom_pane/approval_overlay.rs:237-298`, `codex-rs/tui/src/exec_cell/render.rs:195-246`)

### 10.1 Suggested chemistry-specific command/state ideas inspired by Codex
| Chemsmart idea | Codex pattern to borrow | Why |
| `/dry-run` with artifact preview | `/diff` + exec cell split | Preview generated inputs separately from actual execution. (`codex-rs/tui/src/app/event_dispatch.rs:351-366`, `codex-rs/tui/src/exec_cell/render.rs:195-246`) |
| `/handoff-geom` or auto geometry cell | Structured MCP/tool cell | Make geometry provenance explicit across Gaussian/ORCA. (`codex-rs/tui/src/history_cell.rs:1868-1962`, `https://github.com/Hongjiseung-ROK/chemsmart/pull/15`) |
| `/submit-hpc` approval modal | Exec approval overlay | Show cluster, queue, scratch, walltime, and file targets before submit. (`codex-rs/tui/src/bottom_pane/approval_overlay.rs:626-666`, `/Users/hongjiseung/developer/chemsmart/AGENTS.md:115-137`) |
| `/scan-setup` guided form | BottomPaneView modal | Capture 1-based atom indices and scan_definition cleanly. (`codex-rs/tui/src/bottom_pane/bottom_pane_view.rs:17-138`, `https://github.com/Hongjiseung-ROK/chemsmart/pull/23`) |
| `/optfreq` bundled workflow cell | Plan + exec grouping | Reflect PR #17’s single-job opt+freq route semantics. (`https://github.com/Hongjiseung-ROK/chemsmart/pull/17`, `codex-rs/tui/src/chatwidget.rs:6-18`) |
| `/critic` warning cell | Reasoning/plan cell | Surface chemistry plausibility warnings clearly. (`codex-rs/tui/src/resume_picker/transcript.rs:76-93`, `https://github.com/Hongjiseung-ROK/chemsmart/pull/23`) |

## 11. Concrete code references (paths + line numbers), arranged so the orchestrator can dive in directly

### 11.1 Read these first if you want the main TUI shell
| Path | Why start here | High-value lines |
| `codex-rs/tui/src/chatwidget.rs` | Main chat orchestration: transcript vs active cell, stream handling, interrupt/quit logic, bottom-pane composition. | `1-31`, `208-252`, `3154-3190`, `4298-4327`, `5137-5236`, `10205-10330`, `10956-11020` |
| `codex-rs/tui/src/bottom_pane/mod.rs` | Bottom-pane stack, modal hosting, Esc/Ctrl+C routing, status vs composer layout. | `200-230`, `608-623`, `649-690`, `1520-1564` |
| `codex-rs/tui/src/bottom_pane/chat_composer.rs` | Best single-file description of Codex input UX. | `1-133`, `894-981`, `1969-2058`, `2407-2451`, `3034-3066`, `3807-3885` |
| `codex-rs/tui/src/tui.rs` | Terminal setup, inline-vs-alt-screen behavior, restore logic. | `156-169`, `329-369`, `402-415`, `528-565`, `613-648` |
| `codex-rs/tui/src/pager_overlay.rs` | Transcript/diff overlay system and pager chrome. | `1-16`, `50-89`, `116-248`, `409-456` |

### 11.2 Read these for streaming / interruption / event ordering
| Path | What it owns | High-value lines |
| `codex-rs/tui/src/streaming/controller.rs` | Newline-gated streaming, source retention, resize rewrap. | `1-11`, `61-76`, `119-149`, `222-233`, `237-329` |
| `codex-rs/tui/src/chatwidget/interrupts.rs` | Deferred interrupt queue for approvals and tool/user-input events. | `1-106` |
| `codex-rs/tui/src/tui/event_stream.rs` | Terminal-event normalization and mouse-event skipping. | `236-245` |
| `codex-rs/tui/src/app/input.rs` | Global shortcut entrypoints like transcript overlay opening. | `159-174` |

### 11.3 Read these for commands and keyboard control
| Path | What it owns | High-value lines |
| `codex-rs/tui/src/slash_command.rs` | Canonical slash-command list, descriptions, availability model. | `12-138`, `146-241` |
| `codex-rs/tui/src/bottom_pane/slash_commands.rs` | Feature gating and side-conversation popup filtering. | `13-45`, `155-176` |
| `codex-rs/tui/src/chatwidget/slash_dispatch.rs` | Actual slash-command actions like `/raw`, `/diff`, `/logout`, etc. | `301-327` |
| `codex-rs/tui/src/keymap.rs` | Built-in default bindings across all interaction surfaces. | `542-739` |
| `codex-rs/tui/src/chatwidget.rs` | Terminal-specific fallback for queued-message edit binding. | `208-252` |

### 11.4 Read these for approvals, tools, and diffs
| Path | What it owns | High-value lines |
| `codex-rs/tui/src/bottom_pane/approval_overlay.rs` | Exec/patch/permission/elicitation approval UX. | `155-298`, `508-624`, `626-720`, `790-1025` |
| `codex-rs/tui/src/exec_cell/render.rs` | Shell command transcript rendering, truncation, exploring mode. | `32-35`, `103-183`, `195-260`, `262-359` |
| `codex-rs/tui/src/history_cell.rs` | MCP tool-call cell rendering and interrupted state. | `1775-1962` |
| `codex-rs/tui/src/diff_render.rs` | Diff palette, file summaries, syntax-highlighted hunks. | `1-32`, `60-76`, `195-263`, `403-465`, `475-547`, `548-763` |
| `codex-rs/tui/src/app/event_dispatch.rs` | Turns `/diff` results into alt-screen overlay. | `351-366` |

### 11.5 Read these for theming / visuals
| Path | What it owns | High-value lines |
| `codex-rs/tui/src/render/highlight.rs` | Syntax highlighting engine, theme registry, adaptive default. | `1-22`, `48-51`, `135-171`, `184-224` |
| `codex-rs/tui/src/theme_picker.rs` | Theme picker UX and preview layout. | `1-20`, `129-138`, `166-237`, `314-390` |
| `codex-rs/tui/src/bottom_pane/status_line_style.rs` | Semantic status-line coloring. | `1-15`, `16-75`, `77-173` |
| `codex-rs/tui/src/terminal_palette.rs` | Terminal color-level and fg/bg probing. | `4-18`, `71-151` |
| `codex-rs/tui/src/frames.rs` | Animation frame variants and tick duration. | `3-71` |
| `codex-rs/tui/src/onboarding/auth.rs` | Rounded-border input and OSC-8 hyperlink marking. | `56-104`, `670-684` |

### 11.6 Read these for session persistence / resume
| Path | What it owns | High-value lines |
| `codex-rs/tui/src/session_log.rs` | Opt-in JSONL session recording. | `80-118`, `121-215` |
| `codex-rs/tui/src/session_resume.rs` | Resume/fork cwd prompt and rollout fallback. | `1-5`, `54-64`, `86-111`, `144-189` |
| `codex-rs/tui/src/resume_picker.rs` | Searchable/paginated resume picker. | `286-297` |
| `codex-rs/tui/src/resume_picker/transcript.rs` | Lazy transcript preview loading and fallback transcript cells. | `24-107`, `122-168` |
| `codex-rs/tui/src/app/input.rs` | Transcript overlay hotkey path. | `166-174` |

### 11.7 Read these for packaging / runtime shell behavior
| Path | What it owns | High-value lines |
| `codex-cli/package.json` | npm package entry + files + repo metadata. | `1-21` |
| `codex-cli/bin/codex.js` | Native binary launcher, package selection, signal forwarding. | `15-22`, `78-118`, `120-229` |
| `codex-rs/Cargo.toml` | Workspace member/crate map. | `1-111`, `123-205` |
| `codex-rs/tui/Cargo.toml` | TUI dependencies/features. | `24-160` |
| `README.md` | Product framing / docs map. | `1-8`, `13-29`, `53-60` |
| `AGENTS.md` | Contributor TUI conventions and snapshot-test expectations. | `77-129` |

### 11.8 If I had only one reading order for an implementer
- 1. `chatwidget.rs`
- 2. `bottom_pane/chat_composer.rs`
- 3. `bottom_pane/mod.rs`
- 4. `streaming/controller.rs`
- 5. `pager_overlay.rs`
- 6. `bottom_pane/approval_overlay.rs`
- 7. `keymap.rs` + `slash_command.rs`
- 8. `diff_render.rs` + `exec_cell/render.rs` + `history_cell.rs`
- 9. `theme_picker.rs` + `render/highlight.rs`
- 10. `session_resume.rs` + `session_log.rs` + `resume_picker.rs`

## 12. 5–10 bullet "what to copy for chemsmart" recommendations

- Copy **inline-first main chat + alt-screen inspectors**, so chemistry users keep scrollback for generated inputs and scheduler messages while still having focused viewers for transcript, diff, and artifacts. (`codex-rs/tui/src/tui.rs:329-369`, `codex-rs/tui/src/pager_overlay.rs:1-16`)
- Copy **structured artifact cells** instead of one generic log stream: chemsmart should have separate cell types for plan, dry-run input, geometry handoff, submit script, runtime output, and critic warning. (`codex-rs/tui/src/exec_cell/render.rs:195-246`, `codex-rs/tui/src/history_cell.rs:1868-1962`, `https://github.com/Hongjiseung-ROK/chemsmart/pull/15`)
- Copy **approval overlays with “cancel and steer” choices** for expensive/destructive actions like HPC submit, overwrite, delete scratch, or rerun with changed settings. (`codex-rs/tui/src/bottom_pane/approval_overlay.rs:790-875`, `codex-rs/tui/src/bottom_pane/approval_overlay.rs:973-1025`)
- Copy **newline-gated streaming with source retention**, because chemistry plans and explanations should rewrap cleanly on resize and not become unstable token soup. (`codex-rs/tui/src/streaming/controller.rs:1-11`, `codex-rs/tui/src/streaming/controller.rs:61-76`, `codex-rs/tui/src/streaming/controller.rs:119-149`)
- Copy **phase-aware slash-command gating** so chemsmart only offers the right actions for draft, dry-run, queued, running, failed, and finished job states. (`codex-rs/tui/src/slash_command.rs:165-241`)
- Copy **large-paste-as-placeholder** handling for long route sections, basis/ECP blocks, scan coordinates, and scheduler templates. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:62-71`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:894-929`)
- Copy **cwd-aware resume** and relative-path display, because chemsmart’s workflows are path-sensitive and depend on user-local config plus file artifacts. (`codex-rs/tui/src/session_resume.rs:86-111`, `codex-rs/tui/src/diff_render.rs:739-763`, `/Users/hongjiseung/developer/chemsmart/AGENTS.md:117-137`)
- Copy **rich diff/preview styling** for generated input changes between iterations; chemistry users often care about one changed keyword or coordinate line. (`codex-rs/tui/src/diff_render.rs:403-465`, `codex-rs/tui/src/diff_render.rs:548-734`)
- Copy **dedicated resume/transcript browsing** so users can jump back into prior calculations, not just prior chats. (`codex-rs/tui/src/resume_picker.rs:286-297`, `codex-rs/tui/src/resume_picker/transcript.rs:24-107`)
- Do **not** copy mouse dependence or heavy fullscreen-only assumptions; Codex’s keyboard-first, scrollback-friendly design feels more compatible with research-terminal habits and remote HPC use. (`codex-rs/tui/src/tui/event_stream.rs:236-245`, `codex-rs/tui/src/lib.rs:1584-1615`)

---

## Final synthesis
- Codex CLI’s strongest UI idea is not any single widget; it is the disciplined separation between inline conversation flow, focused overlays, structured tool/artifact cells, and explicit approval states.
- For chemsmart, that separation maps unusually well onto chemistry workflows where “generate”, “inspect”, “validate”, “submit”, “handoff geometry”, and “postprocess” are distinct user intents.
- If I were implementing chemsmart next, I would start by copying: inline scrollback main view, bottom-pane composer with placeholders/autocomplete, transcript overlay, artifact-specific history cells, and approval overlays for submit/overwrite.

## Appendix A. `/keymap` action inventory exposed to users

This appendix is useful because the default-keybinding tables above answer “what keys ship by default,” while the `/keymap` catalog answers “what actions does the product consider stable, user-facing remap targets?” (`codex-rs/tui/src/keymap_setup/actions.rs:1-12`, `codex-rs/tui/src/keymap_setup/actions.rs:87-177`)

| Context | Action id | Description | Evidence |
| Global | `open_transcript` | show transcript view | `codex-rs/tui/src/keymap_setup/actions.rs:89-89` |
| Global | `open_external_editor` | edit current draft externally | `codex-rs/tui/src/keymap_setup/actions.rs:90-90` |
| Global | `copy` | copy latest assistant answer | `codex-rs/tui/src/keymap_setup/actions.rs:91-91` |
| Global | `clear_terminal` | wipe terminal UI area | `codex-rs/tui/src/keymap_setup/actions.rs:92-92` |
| Global | `toggle_vim_mode` | switch composer Vim mode | `codex-rs/tui/src/keymap_setup/actions.rs:93-93` |
| Global | `toggle_fast_mode` | switch Fast mode | `codex-rs/tui/src/keymap_setup/actions.rs:94-95` |
| Global | `toggle_raw_output` | switch raw scrollback mode | `codex-rs/tui/src/keymap_setup/actions.rs:95-95` |
| Chat | `decrease_reasoning_effort` | lower reasoning level | `codex-rs/tui/src/keymap_setup/actions.rs:96-96` |
| Chat | `increase_reasoning_effort` | raise reasoning level | `codex-rs/tui/src/keymap_setup/actions.rs:97-97` |
| Chat | `edit_queued_message` | re-open newest queued draft | `codex-rs/tui/src/keymap_setup/actions.rs:98-98` |
| Composer | `submit` | send current draft | `codex-rs/tui/src/keymap_setup/actions.rs:99-99` |
| Composer | `queue` | queue draft during active task | `codex-rs/tui/src/keymap_setup/actions.rs:100-100` |
| Composer | `toggle_shortcuts` | show/hide shortcut hints | `codex-rs/tui/src/keymap_setup/actions.rs:101-101` |
| Composer | `history_search_previous` | open or move backward in history search | `codex-rs/tui/src/keymap_setup/actions.rs:102-102` |
| Composer | `history_search_next` | move forward in history search | `codex-rs/tui/src/keymap_setup/actions.rs:103-103` |
| Editor | `insert_newline` | insert line break | `codex-rs/tui/src/keymap_setup/actions.rs:104-104` |
| Editor | `move_left` | cursor left | `codex-rs/tui/src/keymap_setup/actions.rs:105-105` |
| Editor | `move_right` | cursor right | `codex-rs/tui/src/keymap_setup/actions.rs:106-106` |
| Editor | `move_up` | cursor up | `codex-rs/tui/src/keymap_setup/actions.rs:107-107` |
| Editor | `move_down` | cursor down | `codex-rs/tui/src/keymap_setup/actions.rs:108-108` |
| Editor | `move_word_left` | jump to previous word start | `codex-rs/tui/src/keymap_setup/actions.rs:109-109` |
| Editor | `move_word_right` | jump to next word end | `codex-rs/tui/src/keymap_setup/actions.rs:110-110` |
| Editor | `move_line_start` | jump to line start | `codex-rs/tui/src/keymap_setup/actions.rs:111-111` |
| Editor | `move_line_end` | jump to line end | `codex-rs/tui/src/keymap_setup/actions.rs:112-112` |
| Editor | `delete_backward` | delete grapheme left | `codex-rs/tui/src/keymap_setup/actions.rs:113-113` |
| Editor | `delete_forward` | delete grapheme right | `codex-rs/tui/src/keymap_setup/actions.rs:114-114` |
| Editor | `delete_backward_word` | delete previous word | `codex-rs/tui/src/keymap_setup/actions.rs:115-115` |
| Editor | `delete_forward_word` | delete next word | `codex-rs/tui/src/keymap_setup/actions.rs:116-116` |
| Editor | `kill_line_start` | delete to line start | `codex-rs/tui/src/keymap_setup/actions.rs:117-117` |
| Editor | `kill_whole_line` | delete whole line | `codex-rs/tui/src/keymap_setup/actions.rs:118-118` |
| Editor | `kill_line_end` | delete to line end | `codex-rs/tui/src/keymap_setup/actions.rs:119-119` |
| Editor | `yank` | paste kill buffer | `codex-rs/tui/src/keymap_setup/actions.rs:120-120` |
| Vim normal | `enter_insert` | enter insert at cursor | `codex-rs/tui/src/keymap_setup/actions.rs:121-121` |
| Vim normal | `append_after_cursor` | insert after cursor | `codex-rs/tui/src/keymap_setup/actions.rs:122-122` |
| Vim normal | `append_line_end` | insert at line end | `codex-rs/tui/src/keymap_setup/actions.rs:123-123` |
| Vim normal | `insert_line_start` | insert at first nonblank | `codex-rs/tui/src/keymap_setup/actions.rs:124-124` |
| Vim normal | `open_line_below` | open line below + insert | `codex-rs/tui/src/keymap_setup/actions.rs:125-125` |
| Vim normal | `open_line_above` | open line above + insert | `codex-rs/tui/src/keymap_setup/actions.rs:126-126` |
| Vim normal | `move_left` | move left | `codex-rs/tui/src/keymap_setup/actions.rs:127-127` |
| Vim normal | `move_right` | move right | `codex-rs/tui/src/keymap_setup/actions.rs:128-128` |
| Vim normal | `move_up` | move up / older history | `codex-rs/tui/src/keymap_setup/actions.rs:129-129` |
| Vim normal | `move_down` | move down / newer history | `codex-rs/tui/src/keymap_setup/actions.rs:130-130` |
| Vim normal | `move_word_forward` | next word start | `codex-rs/tui/src/keymap_setup/actions.rs:131-131` |
| Vim normal | `move_word_backward` | previous word start | `codex-rs/tui/src/keymap_setup/actions.rs:132-132` |
| Vim normal | `move_word_end` | word end motion | `codex-rs/tui/src/keymap_setup/actions.rs:133-133` |
| Vim normal | `move_line_start` | line start motion | `codex-rs/tui/src/keymap_setup/actions.rs:134-134` |
| Vim normal | `move_line_end` | line end motion | `codex-rs/tui/src/keymap_setup/actions.rs:135-135` |
| Vim normal | `delete_char` | delete under cursor | `codex-rs/tui/src/keymap_setup/actions.rs:136-136` |
| Vim normal | `delete_to_line_end` | delete to end of line | `codex-rs/tui/src/keymap_setup/actions.rs:137-137` |
| Vim normal | `yank_line` | yank whole line | `codex-rs/tui/src/keymap_setup/actions.rs:138-138` |
| Vim normal | `paste_after` | paste after cursor | `codex-rs/tui/src/keymap_setup/actions.rs:139-139` |
| Vim normal | `start_delete_operator` | begin delete operator | `codex-rs/tui/src/keymap_setup/actions.rs:140-140` |
| Vim normal | `start_yank_operator` | begin yank operator | `codex-rs/tui/src/keymap_setup/actions.rs:141-141` |
| Vim normal | `cancel_operator` | cancel pending operator | `codex-rs/tui/src/keymap_setup/actions.rs:142-142` |
| Vim operator | `delete_line` | whole-line delete operator | `codex-rs/tui/src/keymap_setup/actions.rs:143-143` |
| Vim operator | `yank_line` | whole-line yank operator | `codex-rs/tui/src/keymap_setup/actions.rs:144-144` |
| Vim operator | `motion_left` | operator motion left | `codex-rs/tui/src/keymap_setup/actions.rs:145-145` |
| Vim operator | `motion_right` | operator motion right | `codex-rs/tui/src/keymap_setup/actions.rs:146-146` |
| Vim operator | `motion_up` | operator motion up | `codex-rs/tui/src/keymap_setup/actions.rs:147-147` |
| Vim operator | `motion_down` | operator motion down | `codex-rs/tui/src/keymap_setup/actions.rs:148-148` |
| Vim operator | `motion_word_forward` | operator next-word motion | `codex-rs/tui/src/keymap_setup/actions.rs:149-149` |
| Vim operator | `motion_word_backward` | operator previous-word motion | `codex-rs/tui/src/keymap_setup/actions.rs:150-150` |
| Vim operator | `motion_word_end` | operator word-end motion | `codex-rs/tui/src/keymap_setup/actions.rs:151-151` |
| Vim operator | `motion_line_start` | operator line-start motion | `codex-rs/tui/src/keymap_setup/actions.rs:152-152` |
| Vim operator | `motion_line_end` | operator line-end motion | `codex-rs/tui/src/keymap_setup/actions.rs:153-153` |
| Vim operator | `cancel` | cancel operator state | `codex-rs/tui/src/keymap_setup/actions.rs:154-154` |
| Pager | `scroll_up` | row-wise up scroll | `codex-rs/tui/src/keymap_setup/actions.rs:155-155` |
| Pager | `scroll_down` | row-wise down scroll | `codex-rs/tui/src/keymap_setup/actions.rs:156-156` |
| Pager | `page_up` | page-wise up scroll | `codex-rs/tui/src/keymap_setup/actions.rs:157-157` |
| Pager | `page_down` | page-wise down scroll | `codex-rs/tui/src/keymap_setup/actions.rs:158-158` |
| Pager | `half_page_up` | half-page up scroll | `codex-rs/tui/src/keymap_setup/actions.rs:159-159` |
| Pager | `half_page_down` | half-page down scroll | `codex-rs/tui/src/keymap_setup/actions.rs:160-160` |
| Pager | `jump_top` | jump to beginning | `codex-rs/tui/src/keymap_setup/actions.rs:161-161` |
| Pager | `jump_bottom` | jump to end | `codex-rs/tui/src/keymap_setup/actions.rs:162-162` |
| Pager | `close` | close pager | `codex-rs/tui/src/keymap_setup/actions.rs:163-163` |
| Pager | `close_transcript` | close transcript pager | `codex-rs/tui/src/keymap_setup/actions.rs:164-164` |
| List | `move_up` | move selection up | `codex-rs/tui/src/keymap_setup/actions.rs:165-165` |
| List | `move_down` | move selection down | `codex-rs/tui/src/keymap_setup/actions.rs:166-166` |
| List | `accept` | confirm current selection | `codex-rs/tui/src/keymap_setup/actions.rs:167-167` |
| List | `cancel` | dismiss selection view | `codex-rs/tui/src/keymap_setup/actions.rs:168-168` |
| Approval | `open_fullscreen` | expand approval details fullscreen | `codex-rs/tui/src/keymap_setup/actions.rs:169-169` |
| Approval | `open_thread` | jump to source thread | `codex-rs/tui/src/keymap_setup/actions.rs:170-170` |
| Approval | `approve` | take primary approve choice | `codex-rs/tui/src/keymap_setup/actions.rs:171-171` |
| Approval | `approve_for_session` | approve for session scope | `codex-rs/tui/src/keymap_setup/actions.rs:172-172` |
| Approval | `approve_for_prefix` | approve with reusable prefix scope | `codex-rs/tui/src/keymap_setup/actions.rs:173-173` |
| Approval | `deny` | take explicit deny choice | `codex-rs/tui/src/keymap_setup/actions.rs:174-174` |
| Approval | `decline` | decline and steer | `codex-rs/tui/src/keymap_setup/actions.rs:175-175` |
| Approval | `cancel` | cancel elicitation flow | `codex-rs/tui/src/keymap_setup/actions.rs:176-176` |

## Appendix B. End-to-end interaction walkthroughs

These are not speculative product docs; they are research summaries of how the code paths fit together.

### B.1 Opening transcript while work is still running
- Global input handler recognizes `Ctrl+T` from the app keymap. (`codex-rs/tui/src/keymap.rs:549-557`, `codex-rs/tui/src/app/input.rs:166-174`)
- TUI enters alternate screen before showing the overlay. (`codex-rs/tui/src/app/input.rs:166-169`, `codex-rs/tui/src/tui.rs:613-648`)
- Overlay type is `Transcript`, not a generic string pager. (`codex-rs/tui/src/app/input.rs:169-173`, `codex-rs/tui/src/pager_overlay.rs:50-58`)
- Transcript overlay renders committed cells plus optional active-cell live tail. (`codex-rs/tui/src/chatwidget.rs:6-16`, `codex-rs/tui/src/pager_overlay.rs:409-456`)
- Live tail is cached by width, revision, stream-continuation flag, and animation tick. (`codex-rs/tui/src/pager_overlay.rs:423-435`)
- Pager header shows a dim `/ T R A N S C R I P T` title. (`codex-rs/tui/src/pager_overlay.rs:174-180`, `codex-rs/tui/src/pager_overlay.rs:443-449`)
- User navigates with pager bindings (`j/k`, arrows, page up/down, etc.). (`codex-rs/tui/src/keymap.rs:690-708`, `codex-rs/tui/src/pager_overlay.rs:250-260`)
- `Ctrl+T` again closes transcript specifically. (`codex-rs/tui/src/keymap.rs:708-708`)
- Leaving overlay restores inline viewport if alt-screen was active. (`codex-rs/tui/src/tui.rs:636-648`)
- The main benefit is that users can inspect the growing history without destabilizing the main streaming viewport. (`codex-rs/tui/src/chatwidget.rs:6-16`, `codex-rs/tui/src/streaming/controller.rs:1-11`)

### B.2 Large paste into composer
- Composer is the sole place where paste text is integrated. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:894-912`)
- Incoming text is normalized from CRLF/CR to `
`. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:912-914`)
- Character count is measured before deciding integration strategy. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:914-915`)
- If count exceeds `LARGE_PASTE_CHAR_THRESHOLD`, Codex creates a placeholder label. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:904-906`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:915-918`)
- The placeholder is inserted into textarea as an element, not plain text. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:916-918`)
- Full pasted text is stored in `pending_pastes` for later expansion. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:915-918`)
- If the paste looks like an image path and image paste is enabled, the composer attaches the image instead. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:919-924`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:932-952`)
- Otherwise the text is inserted directly. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:924-926`)
- Paste-burst state is cleared after explicit paste so later Enter is not misinterpreted. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:927-929`)
- Popup state is re-synced after paste handling. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:928-929`)
- On submit, placeholders are expanded so ranges align with final text. (`codex-rs/tui/src/bottom_pane/chat_composer.rs:45-50`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:57-60`)
- This design would map well to long Gaussian/ORCA route or scan-definition blocks. (chemsmart grounding: `/Users/hongjiseung/developer/chemsmart/AGENTS.md:193-212`, `/Users/hongjiseung/developer/chemsmart/AGENTS.md:226-237`)

### B.3 `/diff` request lifecycle
- User triggers `/diff`; slash dispatch handles the command. (`codex-rs/tui/src/chatwidget/slash_dispatch.rs:328-335`)
- Chat widget adds in-progress state before diff result arrives. (`codex-rs/tui/src/chatwidget/slash_dispatch.rs:329-335`)
- App event dispatch receives `AppEvent::DiffResult(text)`. (`codex-rs/tui/src/app/event_dispatch.rs:351-351`)
- Chat widget clears bottom-pane in-progress state. (`codex-rs/tui/src/app/event_dispatch.rs:352-353`)
- TUI enters alternate screen. (`codex-rs/tui/src/app/event_dispatch.rs:354-355`)
- If diff text is empty, overlay shows `No changes detected.` instead of blank content. (`codex-rs/tui/src/app/event_dispatch.rs:356-360`)
- Otherwise raw lines are ANSI-parsed into pager lines. (`codex-rs/tui/src/app/event_dispatch.rs:356-360`)
- Overlay type is static pager with title `D I F F`. (`codex-rs/tui/src/app/event_dispatch.rs:361-365`)
- Pager handles scrolling and close semantics uniformly. (`codex-rs/tui/src/pager_overlay.rs:116-260`)
- If the diff is itself rendered elsewhere as `FileChange`, Codex can also use the richer syntax-highlighted diff renderer. (`codex-rs/tui/src/diff_render.rs:1-32`, `codex-rs/tui/src/diff_render.rs:403-734`)
- For chemsmart, the same overlay style could inspect generated `.com` / `.inp` changes between plan revisions. (`codex-rs/tui/src/app/event_dispatch.rs:351-366`, `https://github.com/Hongjiseung-ROK/chemsmart/pull/17`)

### B.4 Exec approval lifecycle
- Approval request enters `ApprovalOverlay` either as current request or queued request. (`codex-rs/tui/src/bottom_pane/approval_overlay.rs:169-195`)
- Overlay builds a request-specific header first. (`codex-rs/tui/src/bottom_pane/approval_overlay.rs:214-228`, `codex-rs/tui/src/bottom_pane/approval_overlay.rs:626-720`)
- Exec titles distinguish generic command execution from host/network approval. (`codex-rs/tui/src/bottom_pane/approval_overlay.rs:237-258`)
- Header may show thread label. (`codex-rs/tui/src/bottom_pane/approval_overlay.rs:637-643`)
- Header may show human-readable reason. (`codex-rs/tui/src/bottom_pane/approval_overlay.rs:644-646`)
- Header may summarize permission rule. (`codex-rs/tui/src/bottom_pane/approval_overlay.rs:648-655`, `codex-rs/tui/src/bottom_pane/approval_overlay.rs:878-940`)
- Command itself is de-escaped and syntax-highlighted as bash. (`codex-rs/tui/src/bottom_pane/approval_overlay.rs:657-665`)
- Option list is built from available decisions, not hardcoded blindly. (`codex-rs/tui/src/bottom_pane/approval_overlay.rs:790-875`)
- Keyboard shortcuts can directly choose approve/session/prefix/deny/decline. (`codex-rs/tui/src/bottom_pane/approval_overlay.rs:532-547`, `codex-rs/tui/src/keymap.rs:724-738`)
- `Ctrl+A` can open fullscreen details. (`codex-rs/tui/src/bottom_pane/approval_overlay.rs:512-520`)
- `o` can jump to source thread when relevant. (`codex-rs/tui/src/bottom_pane/approval_overlay.rs:522-530`)
- `Ctrl+C` cancels the modal but does not instantly kill the whole app. (`codex-rs/tui/src/bottom_pane/approval_overlay.rs:561-564`, `codex-rs/tui/src/chatwidget.rs:10225-10236`)

### B.5 Resume picker lifecycle
- Resume picker is launched against app-server thread listings. (`codex-rs/tui/src/resume_picker.rs:286-318`)
- It supports backend filtering by provider/source/eligible cwd before local search. (`codex-rs/tui/src/resume_picker.rs:299-302`)
- Rows are compact multi-line records with stable metadata first and preview last. (`codex-rs/tui/src/resume_picker.rs:289-292`)
- Pagination is cursor-based and deduplicated across pages. (`codex-rs/tui/src/resume_picker.rs:294-297`)
- `Ctrl+E` expands the selected session for recent transcript context. (`codex-rs/tui/src/resume_picker.rs:291-292`)
- Transcript preview loads lazily on demand. (`codex-rs/tui/src/resume_picker/transcript.rs:24-37`)
- User messages in transcript preview can rehydrate text elements and images. (`codex-rs/tui/src/resume_picker/transcript.rs:47-66`)
- Reasoning items can show summary or raw content depending on visibility mode. (`codex-rs/tui/src/resume_picker/transcript.rs:76-93`)
- Fallback transcript cells exist for command execution, file changes, and MCP tool calls. (`codex-rs/tui/src/resume_picker/transcript.rs:110-168`)
- Resume/fork may still prompt for cwd choice if session cwd differs from current cwd. (`codex-rs/tui/src/session_resume.rs:86-111`)
- For chemsmart, the same pattern could browse prior calculations by molecule/job/workdir. (`codex-rs/tui/src/resume_picker.rs:286-297`, `/Users/hongjiseung/developer/chemsmart/AGENTS.md:23-35`)

### B.6 Theme picker lifecycle
- Opening `/theme` builds a `SelectionViewParams` object. (`codex-rs/tui/src/theme_picker.rs:304-390`)
- Current theme snapshot is stored so cancel can restore it. (`codex-rs/tui/src/theme_picker.rs:319-320`, `codex-rs/tui/src/theme_picker.rs:385-390`)
- Entries include bundled themes plus custom themes. (`codex-rs/tui/src/theme_picker.rs:322-323`, `codex-rs/tui/src/render/highlight.rs:174-181`)
- Selection changes live-preview the theme by swapping runtime theme. (`codex-rs/tui/src/theme_picker.rs:373-380`)
- Confirm persists config through `AppEvent::SyntaxThemeSelected`. (`codex-rs/tui/src/theme_picker.rs:353-363`)
- Cancel restores original theme. (`codex-rs/tui/src/theme_picker.rs:385-390`)
- Wide terminals get side-by-side preview panel. (`codex-rs/tui/src/theme_picker.rs:14-20`, `codex-rs/tui/src/theme_picker.rs:129-138`)
- Narrow terminals get compact stacked preview. (`codex-rs/tui/src/theme_picker.rs:18-20`, `codex-rs/tui/src/theme_picker.rs:59-82`)
- Preview content is a syntax-highlighted Rust diff snippet, not arbitrary placeholder color bars. (`codex-rs/tui/src/theme_picker.rs:84-127`, `codex-rs/tui/src/theme_picker.rs:166-237`)
- This is a subtle but smart UX choice: users preview the exact class of code blocks the product often renders. (`codex-rs/tui/src/theme_picker.rs:31-36`, `codex-rs/tui/src/theme_picker.rs:166-237`)

## Appendix C. Chemsmart workflow borrow map

| Chemsmart workflow | Pain point today | Codex pattern to borrow | Source grounding |
| Single-point follow-up after optimization | Need explicit optimized-geometry provenance | Structured artifact cell + transcript overlay | `https://github.com/Hongjiseung-ROK/chemsmart/pull/15`, `codex-rs/tui/src/history_cell.rs:1868-1962`, `codex-rs/tui/src/pager_overlay.rs:1-16` |
| Dry-run input generation | Users need to inspect generated text before run | Large-text placeholder + exec/artifact cell split | `https://github.com/Hongjiseung-ROK/chemsmart/pull/16`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:62-71`, `codex-rs/tui/src/exec_cell/render.rs:195-246` |
| Opt+freq generation | Need one coherent job representation | Grouped cell + plan rationale surface | `https://github.com/Hongjiseung-ROK/chemsmart/pull/17`, `codex-rs/tui/src/resume_picker/transcript.rs:68-93` |
| IRC requests | Required keywords easy to omit | Approval/critic warning modal or warning cell | `https://github.com/Hongjiseung-ROK/chemsmart/pull/16`, `codex-rs/tui/src/bottom_pane/approval_overlay.rs:993-1025` |
| Gaussian scan setup | 1-based atom-index and coordinate spec are error-prone | Guided modal input + file diff preview | `https://github.com/Hongjiseung-ROK/chemsmart/pull/23`, `codex-rs/tui/src/bottom_pane/bottom_pane_view.rs:17-138`, `codex-rs/tui/src/diff_render.rs:403-465` |
| HPC submission | Submission is path/config sensitive and potentially destructive | Explicit submit approval surface | `/Users/hongjiseung/developer/chemsmart/AGENTS.md:115-137`, `codex-rs/tui/src/bottom_pane/approval_overlay.rs:237-298` |
| Planner reproducibility | Need visible rationale and normalized plan artifacts | Plan cell + transcript persistence + session logging | `https://github.com/Hongjiseung-ROK/chemsmart/pull/23`, `codex-rs/tui/src/session_log.rs:80-118`, `codex-rs/tui/src/resume_picker/transcript.rs:68-93` |
| File-heavy worktrees | Need safe resume when cwd changes | Cwd mismatch prompt + relative path display | `codex-rs/tui/src/session_resume.rs:86-111`, `codex-rs/tui/src/diff_render.rs:739-763` |
| Long route/basis text | Composer becomes unreadable | Placeholderized large paste/object model | `codex-rs/tui/src/bottom_pane/chat_composer.rs:62-71`, `codex-rs/tui/src/bottom_pane/chat_composer.rs:894-929` |
| Trajectory / IRC / movie artifacts | Need focused inspection without cluttering main chat | Alt-screen static/transcript pager | `/Users/hongjiseung/developer/chemsmart/AGENTS.md:297-299`, `codex-rs/tui/src/pager_overlay.rs:1-16` |
| Cross-program Gaussian/ORCA workflows | Need explicit phase/state UX | Phase-aware slash commands + approvals | `/Users/hongjiseung/developer/chemsmart/AGENTS.md:193-212`, `/Users/hongjiseung/developer/chemsmart/AGENTS.md:226-237`, `codex-rs/tui/src/slash_command.rs:165-241` |
| Open-shell / basis plausibility warnings | Warnings should be visible, not buried | Dedicated critic/warning history cells | `https://github.com/Hongjiseung-ROK/chemsmart/pull/23`, `codex-rs/tui/src/history_cell.rs:1868-1962` |


## Appendix D. Symbol-level jump list for implementers

This appendix is intentionally terse: open the symbol/region when you need to steal one behavior.

| Need | Jump to | Why |
| Main render composition | `codex-rs/tui/src/chatwidget.rs:10956-11020` | Shows top-level flex layout. |
| Status row hiding during commentary | `codex-rs/tui/src/chatwidget.rs:24-27` | Avoid duplicate progress indicators. |
| Queued-message terminal fallback | `codex-rs/tui/src/chatwidget.rs:208-252` | Alt+Up vs Shift+Left logic. |
| Interrupt aftermath | `codex-rs/tui/src/chatwidget.rs:3154-3190` | How interrupted turns restore queued steers/drafts. |
| Deterministic interrupt ordering | `codex-rs/tui/src/chatwidget.rs:4304-4317` | Why interrupts queue during streaming. |
| Flush interrupts after stream | `codex-rs/tui/src/chatwidget.rs:4320-4327` | Where queued prompts become visible. |
| Copy-last-response hotkey | `codex-rs/tui/src/chatwidget.rs:5144-5151` | Keyboard path to copy action. |
| Clipboard image paste | `codex-rs/tui/src/chatwidget.rs:5177-5202` | Attach image from clipboard. |
| Esc interrupt with pending steers | `codex-rs/tui/src/chatwidget.rs:5225-5236` | Chemistry steering analogue. |
| Ctrl+C quit/interrupt state machine | `codex-rs/tui/src/chatwidget.rs:10205-10320` | One of the most reusable interaction designs. |
| Bottom-pane modal stack | `codex-rs/tui/src/bottom_pane/mod.rs:204-230` | Where composer and views coexist. |
| Esc interrupt from bottom pane | `codex-rs/tui/src/bottom_pane/mod.rs:608-623` | How running-task Esc is routed. |
| Ctrl+C inside bottom pane | `codex-rs/tui/src/bottom_pane/mod.rs:649-682` | Modal-first cancellation behavior. |
| Bottom-pane render stack | `codex-rs/tui/src/bottom_pane/mod.rs:1520-1564` | Status/footer/previews/composer ordering. |
| Composer architecture comment | `codex-rs/tui/src/bottom_pane/chat_composer.rs:1-133` | Best single-file UX spec. |
| Large-paste logic | `codex-rs/tui/src/bottom_pane/chat_composer.rs:894-929` | Object-like paste handling. |
| Image-path paste | `codex-rs/tui/src/bottom_pane/chat_composer.rs:932-952` | Path-to-image attachment behavior. |
| Disable paste-burst safely | `codex-rs/tui/src/bottom_pane/chat_composer.rs:955-981` | State defusing logic. |
| File selection accept | `codex-rs/tui/src/bottom_pane/chat_composer.rs:1969-2058` | How popup inserts path vs attaches image. |
| Quoted path insertion | `codex-rs/tui/src/bottom_pane/chat_composer.rs:2407-2451` | Whitespace-safe path insertion. |
| Remote image rows | `codex-rs/tui/src/bottom_pane/chat_composer.rs:3034-3066` | Attachment row rendering/removal. |
| Slash vs file popup precedence | `codex-rs/tui/src/bottom_pane/chat_composer.rs:3807-3815` | Great UX detail. |
| File popup synchronization | `codex-rs/tui/src/bottom_pane/chat_composer.rs:3853-3885` | Async search wiring. |
| File popup states | `codex-rs/tui/src/bottom_pane/file_search_popup.rs:16-153` | Idle/loading/no matches handling. |
| Bottom-pane view contract | `codex-rs/tui/src/bottom_pane/bottom_pane_view.rs:17-138` | Reusable modal interface. |
| Approval overlay setup | `codex-rs/tui/src/bottom_pane/approval_overlay.rs:155-298` | Request-specific titles/options. |
| Approval shortcuts | `codex-rs/tui/src/bottom_pane/approval_overlay.rs:508-547` | Keyboard-first action routing. |
| Approval footer hints | `codex-rs/tui/src/bottom_pane/approval_overlay.rs:601-624` | Human hint line composition. |
| Exec approval header | `codex-rs/tui/src/bottom_pane/approval_overlay.rs:626-666` | Thread/reason/permission/command display. |
| Permission approval header | `codex-rs/tui/src/bottom_pane/approval_overlay.rs:667-692` | Requested permissions summary. |
| Patch approval header | `codex-rs/tui/src/bottom_pane/approval_overlay.rs:693-720` | Edit approval framing. |
| Exec approval options | `codex-rs/tui/src/bottom_pane/approval_overlay.rs:790-875` | Best place to borrow phrasing. |
| Patch approval options | `codex-rs/tui/src/bottom_pane/approval_overlay.rs:973-990` | Small, clean option set. |
| Permission approval options | `codex-rs/tui/src/bottom_pane/approval_overlay.rs:993-1025` | Turn vs session semantics. |
| Slash command inventory | `codex-rs/tui/src/slash_command.rs:12-138` | Canonical command list and descriptions. |
| Slash availability rules | `codex-rs/tui/src/slash_command.rs:146-241` | Inline-arg / during-task / side-convo model. |
| Slash feature gating | `codex-rs/tui/src/bottom_pane/slash_commands.rs:13-45` | Visibility depends on enabled features. |
| Slash dispatch for raw/diff | `codex-rs/tui/src/chatwidget/slash_dispatch.rs:321-335` | Hook from command to action. |
| Default keymap | `codex-rs/tui/src/keymap.rs:542-739` | All built-in shortcuts. |
| Remappable action catalog | `codex-rs/tui/src/keymap_setup/actions.rs:87-177` | Stable user-facing action ids. |
| Streaming core overview | `codex-rs/tui/src/streaming/controller.rs:1-45` | Source retention model. |
| Newline gate | `codex-rs/tui/src/streaming/controller.rs:61-76` | When live output becomes renderable. |
| Resize rewrap | `codex-rs/tui/src/streaming/controller.rs:119-149` | Pending queue rebuild logic. |
| Assistant stream controller | `codex-rs/tui/src/streaming/controller.rs:237-329` | Active assistant message path. |
| Pager overlay overview | `codex-rs/tui/src/pager_overlay.rs:1-16` | Why transcript/diff overlays exist. |
| Pager header/footer | `codex-rs/tui/src/pager_overlay.rs:174-247` | `/ TITLE` and percent bar. |
| Transcript overlay live tail | `codex-rs/tui/src/pager_overlay.rs:409-456` | In-flight transcript visibility. |
| MCP tool-call cell | `codex-rs/tui/src/history_cell.rs:1775-1962` | Structured tool-call display. |
| Exec transcript renderer | `codex-rs/tui/src/exec_cell/render.rs:195-246` | Command/output/result layout. |
| Exec exploring mode | `codex-rs/tui/src/exec_cell/render.rs:262-359` | Semantic shell grouping. |
| Diff renderer overview | `codex-rs/tui/src/diff_render.rs:1-32` | Add/delete/update rendering philosophy. |
| Diff per-file summary | `codex-rs/tui/src/diff_render.rs:403-465` | Header and rename display. |
| Diff hunk rendering | `codex-rs/tui/src/diff_render.rs:548-734` | Parser-state-preserving hunk flow. |
| Path normalization | `codex-rs/tui/src/diff_render.rs:739-763` | Relative path display helper. |
| Theme picker overview | `codex-rs/tui/src/theme_picker.rs:1-20` | Live preview/cancel restore/persist model. |
| Theme picker preview rendering | `codex-rs/tui/src/theme_picker.rs:166-237` | How preview actually draws diff snippet. |
| Theme picker actions | `codex-rs/tui/src/theme_picker.rs:314-390` | Selection change and cancel behavior. |
| Highlighter guardrails | `codex-rs/tui/src/render/highlight.rs:1-22` | Theme/highlight limits. |
| Adaptive theme choice | `codex-rs/tui/src/render/highlight.rs:184-201` | Light vs dark default logic. |
| Status-line semantic colors | `codex-rs/tui/src/bottom_pane/status_line_style.rs:16-119` | Theme-scope accents. |
| Animation frame variants | `codex-rs/tui/src/frames.rs:47-71` | Spinner/style options. |
| Terminal setup | `codex-rs/tui/src/tui.rs:156-169` | Bracketed paste/raw mode/focus change. |
| Inline init | `codex-rs/tui/src/tui.rs:329-369` | Scrollback-preserving startup. |
| Alt-screen enter/leave | `codex-rs/tui/src/tui.rs:613-648` | Overlay screen management. |
| Alt-screen policy | `codex-rs/tui/src/lib.rs:1584-1615` | Why Zellij disables auto alt-screen. |
| Terminal restore for external editor | `codex-rs/tui/src/tui.rs:528-565` | Pause events + restore modes safely. |
| Color probing | `codex-rs/tui/src/terminal_palette.rs:12-18` | Color-depth detection. |
| Default fg/bg probing | `codex-rs/tui/src/terminal_palette.rs:71-151` | Adaptive palettes and themes. |
| Session logger init | `codex-rs/tui/src/session_log.rs:80-118` | Opt-in transcript/event logging. |
| Session logger events | `codex-rs/tui/src/session_log.rs:121-215` | What gets recorded. |
| Resume cwd prompt | `codex-rs/tui/src/session_resume.rs:86-111` | Mismatch handling. |
| Rollout resume fallback | `codex-rs/tui/src/session_resume.rs:144-189` | Recover thread/cwd/model from JSONL. |
| Resume picker overview | `codex-rs/tui/src/resume_picker.rs:286-297` | Search/paginate/expand session browser. |
| Transcript preview loader | `codex-rs/tui/src/resume_picker/transcript.rs:24-107` | Lazy transcript conversion. |
| Clipboard image ingest | `codex-rs/tui/src/clipboard_paste.rs:49-149` | Temp-PNG pipeline. |
| Auth hyperlinks | `codex-rs/tui/src/onboarding/auth.rs:56-104` | OSC-8 wrapped URLs. |
| Rounded API-key input | `codex-rs/tui/src/onboarding/auth.rs:677-683` | Focused form styling example. |
| npm launcher | `codex-cli/bin/codex.js:120-229` | Signal forwarding and exit mirroring. |
