
# Claude Code UI patterns research

- Research date: 2026-05-08.
- Scope: Claude Code CLI/TUI patterns only; no Codex CLI comparison here.
- Local version observed: `claude --version` returned `2.1.133 (Claude Code)`. [LOCAL-help]
- Local behavioral observation was limited to safe startup/help/config inspection and a scratch-directory interactive launch; I did **not** attempt to reverse-engineer or disassemble the native binary. [LOCAL-help] [LOCAL-trust]
- chemsmart grounding references: repo `AGENTS.md`, `frontier_improvement_plan.md`, and merged fork PRs #15, #16, #17, and #23. [CHEMSMART-AGENTS] [CHEMSMART-plan] [CHEMSMART-PR15] [CHEMSMART-PR16] [CHEMSMART-PR17] [CHEMSMART-PR23]

## How to read this note

- I mark direct official-documentation facts as ordinary bullets with citations.
- I mark best-effort renderer/stack guesses as **Inference** and keep the evidentiary chain explicit.
- When docs describe newer behavior than the locally installed `2.1.133`, I say so.
- When a feature is not documented for the terminal UI, I say **unknown / not documented** rather than guessing.

## 1. Tool overview & where source lives

### 1.1 What is public

- The public `anthropics/claude-code` repository exposes installer/readme material, examples, plugins, scripts, and a large public `CHANGELOG.md`; it does **not** expose the entire interactive terminal renderer source tree in the way an ordinary open-source TUI project would. [REPO-root] [REPO-readme]
- The repo top level currently contains `.claude-plugin/`, `.claude/commands/`, `examples/`, `plugins/`, `scripts/`, `README.md`, `SECURITY.md`, and `CHANGELOG.md`. That is a real extension/examples repository, but not the full application source for the interactive TUI itself. [REPO-root]
- Anthropic’s official docs cover the CLI, sessions, keyboard shortcuts, permissions, hooks, MCP, plugins, fullscreen rendering, worktrees, and desktop/web/IDE surfaces. [DOC-overview] [DOC-commands] [DOC-interactive] [DOC-keybindings] [DOC-settings] [DOC-permissions] [DOC-hooks] [DOC-sessions] [DOC-fullscreen] [DOC-terminal]
- The public repo also ships concrete examples that are especially useful for UI-adjacent extensibility research: settings JSON, hook scripts, and plugin READMEs that describe slash-command workflows. [REPO-settings-examples] [REPO-hook-example] [REPO-plugins-readme]

### 1.2 What is binary-only or otherwise closed

- On this machine, Claude Code 2.1.133 is installed as a single native executable at `/Users/hongjiseung/.local/share/claude/versions/2.1.133`. The file is not accompanied by an open source renderer tree in the install directory. [LOCAL-binary]
- Anthropic’s docs do not publish the full internal Claude Code system prompt or the full terminal renderer implementation. The configuration docs explicitly say the internal system prompt is not published. [DOC-configuration-system-prompt]
- Therefore the **interactive terminal UI implementation details** below must be split into:
  - directly documented behavior,
  - behavior observed locally, and
  - cautious inference from public issues/changelog entries.

### 1.3 Practical implication for chemsmart

- For chemsmart, Claude Code is a good source of **interaction design patterns** and **state-model ideas**, but not a reusable codebase for the TUI shell itself. [REPO-root] [DOC-interactive] [DOC-fullscreen]
- The reusable public material is mainly in docs and extension examples, not in the renderer core. [REPO-root] [REPO-plugins-readme]

## 2. Tech stack & rendering approach

### 2.1 What I can say confidently

- Current recommended installs are **native installs** (`curl ... | bash`, Homebrew cask, WinGet), and the repo README says npm installation is deprecated. [REPO-readme] [DOC-overview]
- The current CLI supports two renderers from the user perspective: `default` and `fullscreen`, exposed through `/tui [default|fullscreen]`. [DOC-commands] [LOCAL-help]
- Fullscreen rendering is explicitly documented as an alternate-screen-buffer renderer that is “flicker-free”, mouse-aware, and virtualizes the render tree to keep memory flat on long conversations. [DOC-fullscreen] [DOC-terminal]
- The changelog for 2.1.132 adds `CLAUDE_CODE_DISABLE_ALTERNATE_SCREEN=1` to opt out of the fullscreen alternate-screen renderer and keep the conversation in native scrollback. [CHANGELOG-2.1.132]
- The docs state that fullscreen mode “draws the interface on the terminal’s alternate screen buffer, like `vim` or `htop`, and only renders messages that are currently visible.” [DOC-fullscreen]
- The docs and changelog repeatedly discuss terminal-protocol behavior such as bracketed paste, mouse tracking, synchronized output, Kitty keyboard protocol, tmux passthrough, and VS Code integrated terminal quirks. That means the renderer is terminal-protocol-heavy, not just line-by-line stdout printing. [DOC-fullscreen] [DOC-terminal] [DOC-keybindings] [CHANGELOG-2.1.132]

### 2.2 Best-effort inference about the underlying stack

- **Inference:** legacy/open npm-distributed Claude Code almost certainly used a JavaScript/Node terminal stack with Ink in at least some 2025 builds, because issue `#404` includes the exact runtime error `Raw mode is not supported on the current process.stdin, which Ink uses as input stream by default`, plus stack traces into `@anthropic-ai/claude-code/cli.js`. [ISSUE-404]
- **Inference:** current native builds likely package a JS application into a standalone runtime rather than being a conventional Curses application, because issue `#25630` reports a native Windows binary crash that prints Bun standalone-runtime metadata (`Bun Canary ... standalone_executable`) and users explicitly contrast it with the older npm/Node path. [ISSUE-25630]
- **Inference:** the present renderer may still descend from the older Ink/JS architecture, but the fullscreen path now looks custom or heavily adapted, because official docs describe alternate-screen virtualization, in-app selection/copy, click-to-expand tool results, transcript mode, and custom rebindable action contexts in much more depth than a stock Ink app would normally expose. This inference is **plausible** but not directly confirmed. [ISSUE-404] [ISSUE-25630] [DOC-fullscreen] [DOC-keybindings]
- **Inference:** the split between `default` and `fullscreen` renderers suggests two drawing strategies:
  - a classic scrollback-preserving terminal renderer, and
  - an alternate-screen virtualized renderer with richer mouse/search/copy semantics. [DOC-commands] [DOC-fullscreen]

### 2.3 What chemsmart should take from this

- Copy the **renderer split**, not the unknown implementation: keep a plain scrollback-safe mode and a richer fullscreen mode. [DOC-commands] [DOC-fullscreen]
- Prefer a TUI architecture that treats the transcript, footer, dialogs, and diff viewers as composable surfaces instead of a single scrolling text stream. [DOC-keybindings] [DOC-interactive] [DOC-fullscreen]

## 3. Layout structure

### 3.1 Startup / gate screens

- On first entry into a directory, Claude Code shows a **workspace trust / safety check** screen before granting itself read-edit-execute access. In my local observation, the screen displayed the path, a short explanation of access, a security-guide link, and two choices (`Yes, I trust this folder` / `No, exit`). [LOCAL-trust] [DOC-security]
- This trust screen is important: it sets the tone that the CLI is not just a chatbot but an agent with filesystem and shell power. [LOCAL-trust] [DOC-security]

### 3.2 Primary chat layout

- The main chat surface is organized around a bottom input box and an upper transcript region. Official docs repeatedly describe the input staying fixed at the bottom in fullscreen rendering. [DOC-fullscreen]
- The prompt bar itself carries state: session name, prompt-bar color, permission-mode border accents, fast-mode indicator, model picker access, and prompt suggestions. [DOC-commands] [DOC-keybindings] [DOC-terminal]
- The transcript region is not just plain assistant text. It can show tool calls, collapsible tool results, PR badges, task summaries, memory markers, and transcript-mode overlays. [DOC-interactive] [DOC-fullscreen] [DOC-terminal]
- There is a **footer/status area** with lightweight indicators for things like tasks, focus mode, PR status, and diff/session state; the keybinding contexts explicitly include `Footer`, which means Anthropic treats the footer as a first-class navigable UI element rather than decorative text. [DOC-keybindings] [DOC-interactive]

### 3.3 Overlays, pickers, tabs, and modal screens

- Claude Code relies heavily on modal or picker-style sub-surfaces rather than permanent side panels. I found no evidence of an always-on left sidebar in the terminal CLI. [DOC-interactive] [DOC-sessions] [DOC-keybindings]
- Documented modal/picker surfaces include:
  - `/resume` session picker, [DOC-sessions]
  - `/config` tabbed settings interface, [DOC-settings]
  - `/permissions` rule-management dialog, [DOC-permissions]
  - `/diff` interactive diff viewer, [DOC-commands]
  - `/theme` theme picker, [DOC-commands] [DOC-terminal]
  - `/doctor` diagnostics screen, [DOC-commands] [DOC-keybindings]
  - `/help` command browser, [DOC-commands]
  - plugin dialog / library screens, implied by `Plugin` keybinding context, [DOC-keybindings]
  - transcript viewer / transcript mode, toggled with `Ctrl+O`, [DOC-interactive] [DOC-fullscreen]
  - rewind/message selector overlays, implied by `MessageSelector` context, [DOC-keybindings] [DOC-checkpointing]
  - model picker / effort slider, [DOC-commands] [DOC-keybindings]
  - select/list components and confirmation dialogs. [DOC-keybindings]
- The settings screen is explicitly described as **tabbed**. [DOC-settings]
- Permission and confirmation dialogs are explicit enough in the keybinding docs to get their own `Confirmation` context. [DOC-keybindings]
- Session and diff navigation are deep enough to get their own contexts (`MessageSelector`, `DiffDialog`, `HistorySearch`, `Task`, `Doctor`, `ThemePicker`, `Plugin`, `ModelPicker`). That is a strong signal that Claude Code’s TUI is closer to a modal application than to a REPL with occasional prompts. [DOC-keybindings]

### 3.4 Transcript viewer versus normal chat view

- `Ctrl+O` toggles a transcript viewer / verbose transcript. [DOC-keybindings] [DOC-interactive]
- In transcript mode, Claude exposes expanded tool details and MCP calls that may otherwise collapse to one-line summaries. [DOC-interactive]
- In fullscreen rendering, transcript mode adds search/navigation controls (`/`, `n`, `N`, `j`, `k`, `g`, `G`, `Ctrl+u`, `Ctrl+d`, `Ctrl+b`, `Ctrl+f`) plus export-to-scrollback (`[`) and export-to-editor (`v`). [DOC-fullscreen] [DOC-interactive]
- There is also a **focus view** (`/focus`) that reduces the conversation to the last prompt, a one-line tool summary with diffstats, and the final response. [DOC-commands] [DOC-fullscreen]

### 3.5 Side panels?

- **Unknown / not documented:** I did not find evidence for a permanent split-pane sidebar in the terminal CLI analogous to some desktop IDE chat tools. [DOC-interactive] [DOC-keybindings]
- What Claude Code does have is a footer-driven micro-dashboard plus modal/detail views. [DOC-keybindings] [DOC-interactive]

## 4. Streaming / token rendering & interruption

### 4.1 Streaming model

- Claude Code streams normal conversation responses and streams tool execution progress into the transcript. The docs frame the experience as an agentic loop where Claude reads files, edits, runs commands, and stays interruptible. [DOC-how] [DOC-interactive]
- Fullscreen mode exists partly because incremental updates can otherwise flicker or cause scrollback jumps while tool output streams in. [DOC-fullscreen] [DOC-terminal]
- Only visible messages are kept in the fullscreen render tree, which is a virtualization optimization specifically for long streamed sessions. [DOC-fullscreen]

### 4.2 Interrupt semantics

- `Ctrl+C` is the canonical documented interrupt for “cancel current input or generation.” [DOC-interactive] [DOC-keybindings]
- Official “How Claude Code works” docs also say you can type a correction and press Enter while Claude is working; Claude stops what it is doing and adjusts its approach. [DOC-how]
- The changelog shows that **Esc also matters operationally** in current/fullscreen sessions:
  - one entry says a remote-control stop/interrupt should cancel the CLI session “the same way local Esc does”, [CHANGELOG-2.1.133]
  - one earlier entry notes “Esc in INSERT” then “press Esc again to interrupt” in Vim mode, [CHANGELOG-2.1.125-ish]
  - another notes a bug where Escape failed to interrupt a running turn when draft text existed. [CHANGELOG-escape-interrupt]
- My read is: **Ctrl+C is the stable documented interrupt; Esc is renderer/mode-sensitive but clearly part of the current interactive interruption model.** [DOC-interactive] [CHANGELOG-2.1.133] [CHANGELOG-escape-interrupt]

### 4.3 Tool and background-task interruption

- Bash commands and agents can be moved into the background with `Ctrl+B`. [DOC-interactive] [DOC-keybindings]
- Background tasks get unique IDs, write output to files, and are cleaned up automatically on exit; overly large outputs (>5GB) are terminated. [DOC-interactive]
- `Ctrl+X Ctrl+K` kills all background agents, with a second-confirm press within 3 seconds. [DOC-interactive] [DOC-keybindings]
- The changelog repeatedly mentions resumed sessions, interrupted turns, and graceful shutdown fixes, which implies Anthropic treats interruption as a core reliability path rather than a rare edge case. [CHANGELOG-2.1.132] [CHANGELOG-resume-cluster]

### 4.4 Token / context handling

- Claude Code stores messages, tool calls, and results in local JSONL transcripts under `~/.claude/projects/`. [DOC-how] [DOC-claude-dir] [DOC-sessions]
- As context fills, Claude Code automatically clears older tool outputs first and then compacts/summarizes the conversation. [DOC-how]
- `/compact` is the explicit manual control for compaction, and `/context` visualizes what is occupying context space. [DOC-commands] [DOC-context-window]
- Side questions (`/btw`) reuse the parent conversation’s prompt cache and intentionally do **not** add to the main history. [DOC-interactive]

## 5. Slash commands (full list) & keybindings

### 5.1 Slash-command notes

- Claude Code commands are recognized only at the **start** of a message. [DOC-commands]
- Typing `/` opens a filtered command browser. [DOC-commands]
- The command catalog mixes built-in commands and “Skill” entries; bundled skills are prompt-based workflows, not hard-coded imperative commands. [DOC-commands] [DOC-skills]
- MCP prompts can also surface as slash commands using the `/mcp__server__prompt` naming pattern. [DOC-slash] [DOC-mcp]
- Skills and legacy custom commands both surface as slash commands, with skills taking precedence when names collide. [DOC-skills]

### 5.2 Full documented slash-command list (current docs)

| Command | Description | Evidence |
|---|---|---|
| `/add-dir <path>` | Add a working directory for file access during the current session. Most .claude/ configuration is not discovered from the added directory. You can later resume the session from the added directory with --continue or --resume | [DOC-commands] |
| `/agents` | Manage agent configurations | [DOC-commands] |
| `/autofix-pr [prompt]` | Spawn a Claude Code on the web session that watches the current branch’s PR and pushes fixes when CI fails or reviewers leave comments. Detects the open PR from your checked-out branch with gh pr view ; to watch a different PR, check out its branch first. By default the remote session is told to fix every CI failure and review comment; pass a prompt to give it different instructions, for example /autofix-pr only fix lint and type errors . Requires the gh CLI and access to Claude Code on the web | [DOC-commands] |
| `/batch <instruction>` | Skill . Orchestrate large-scale changes across a codebase in parallel. Researches the codebase, decomposes the work into 5 to 30 independent units, and presents a plan. Once approved, spawns one background agent per unit in an isolated git worktree . Each agent implements its unit, runs tests, and opens a pull request. Requires a git repository. Example: /batch migrate src/ from Solid to React | [DOC-commands] |
| `/branch [name]` | Create a branch of the current conversation at this point. Switches you into the branch and preserves the original, which you can return to with /resume . Alias: /fork . When CLAUDE_CODE_FORK_SUBAGENT is set, /fork instead spawns a forked subagent and is no longer an alias for this command | [DOC-commands] |
| `/btw <question>` | Ask a quick side question without adding to the conversation | [DOC-commands] |
| `/chrome` | Configure Claude in Chrome settings | [DOC-commands] |
| `/claude-api [migrate|managed-agents-onboard]` | Skill . Load Claude API reference material for your project’s language (Python, TypeScript, Java, Go, Ruby, C#, PHP, or cURL) and Managed Agents reference. Covers tool use, streaming, batches, structured outputs, and common pitfalls. Also activates automatically when your code imports anthropic or @anthropic-ai/sdk . Run /claude-api migrate to upgrade existing Claude API code to a newer model: Claude asks which files to scan and which model to target, then updates model IDs, thinking configuration, and other parameters that changed between versions. Run /claude-api managed-agents-onboard for an interactive walkthrough that creates a new Managed Agent from scratch | [DOC-commands] |
| `/clear` | Start a new conversation with empty context. The previous conversation stays available in /resume . To free up context while continuing the same conversation, use /compact instead. Aliases: /reset , /new | [DOC-commands] |
| `/color [color|default]` | Set the prompt bar color for the current session. Available colors: red , blue , green , yellow , purple , orange , pink , cyan . Use default to reset, or run with no argument to pick a random color. When Remote Control is connected, the color syncs to claude.ai/code | [DOC-commands] |
| `/compact [instructions]` | Free up context by summarizing the conversation so far. Optionally pass focus instructions for the summary. See how compaction handles rules, skills, and memory files | [DOC-commands] |
| `/config` | Open the Settings interface to adjust theme, model, output style , and other preferences. Alias: /settings | [DOC-commands] |
| `/context` | Visualize current context usage as a colored grid. Shows optimization suggestions for context-heavy tools, memory bloat, and capacity warnings | [DOC-commands] |
| `/copy [N]` | Copy the last assistant response to clipboard. Pass a number N to copy the Nth-latest response: /copy 2 copies the second-to-last. When code blocks are present, shows an interactive picker to select individual blocks or the full response. Press w in the picker to write the selection to a file instead of the clipboard, which is useful over SSH | [DOC-commands] |
| `/cost` | Alias for /usage | [DOC-commands] |
| `/debug [description]` | Skill . Enable debug logging for the current session and troubleshoot issues by reading the session debug log. Debug logging is off by default unless you started with claude --debug , so running /debug mid-session starts capturing logs from that point forward. Optionally describe the issue to focus the analysis | [DOC-commands] |
| `/desktop` | Continue the current session in the Claude Code Desktop app. macOS and Windows only. Alias: /app | [DOC-commands] |
| `/diff` | Open an interactive diff viewer showing uncommitted changes and per-turn diffs. Use left/right arrows to switch between the current git diff and individual Claude turns, and up/down to browse files | [DOC-commands] |
| `/doctor` | Diagnose and verify your Claude Code installation and settings. Results show with status icons. Press f to have Claude fix any reported issues | [DOC-commands] |
| `/effort [level|auto]` | Set the model effort level . Accepts low , medium , high , xhigh , or max ; available levels depend on the model and max is session-only. auto resets to the model default. Without an argument, opens an interactive slider; use left and right arrows to pick a level and Enter to apply. Takes effect immediately without waiting for the current response to finish | [DOC-commands] |
| `/exit` | Exit the CLI. Alias: /quit | [DOC-commands] |
| `/export [filename]` | Export the current conversation as plain text. With a filename, writes directly to that file. Without, opens a dialog to copy to clipboard or save to a file | [DOC-commands] |
| `/extra-usage` | Configure extra usage to keep working when rate limits are hit | [DOC-commands] |
| `/fast [on|off]` | Toggle fast mode on or off | [DOC-commands] |
| `/feedback [report]` | Submit feedback about Claude Code. Alias: /bug | [DOC-commands] |
| `/fewer-permission-prompts` | Skill . Scan your transcripts for common read-only Bash and MCP tool calls, then add a prioritized allowlist to project .claude/settings.json to reduce permission prompts | [DOC-commands] |
| `/focus` | Toggle the focus view, which shows only your last prompt, a one-line tool-call summary with edit diffstats, and the final response. The selection persists across sessions; set viewMode in settings to override it. Only available in fullscreen rendering | [DOC-commands] |
| `/heapdump` | Write a JavaScript heap snapshot and a memory breakdown to ~/Desktop , or your home directory on Linux without a Desktop folder, for diagnosing high memory usage. See troubleshooting | [DOC-commands] |
| `/help` | Show help and available commands | [DOC-commands] |
| `/hooks` | View hook configurations for tool events | [DOC-commands] |
| `/ide` | Manage IDE integrations and show status | [DOC-commands] |
| `/init` | Initialize project with a CLAUDE.md guide. Set CLAUDE_CODE_NEW_INIT=1 for an interactive flow that also walks through skills, hooks, and personal memory files | [DOC-commands] |
| `/insights` | Generate a report analyzing your Claude Code sessions, including project areas, interaction patterns, and friction points | [DOC-commands] |
| `/install-github-app` | Set up the Claude GitHub Actions app for a repository. Walks you through selecting a repo and configuring the integration | [DOC-commands] |
| `/install-slack-app` | Install the Claude Slack app. Opens a browser to complete the OAuth flow | [DOC-commands] |
| `/keybindings` | Open or create your keybindings configuration file | [DOC-commands] |
| `/login` | Sign in to your Anthropic account | [DOC-commands] |
| `/logout` | Sign out from your Anthropic account | [DOC-commands] |
| `/loop [interval] [prompt]` | Skill . Run a prompt repeatedly while the session stays open. Omit the interval and Claude self-paces between iterations. Omit the prompt and Claude runs an autonomous maintenance check, or the prompt in .claude/loop.md if present. Example: /loop 5m check if the deploy finished . See Run prompts on a schedule . Alias: /proactive | [DOC-commands] |
| `/mcp` | Manage MCP server connections and OAuth authentication | [DOC-commands] |
| `/memory` | Edit CLAUDE.md memory files, enable or disable auto-memory , and view auto-memory entries | [DOC-commands] |
| `/mobile` | Show QR code to download the Claude mobile app. Aliases: /ios , /android | [DOC-commands] |
| `/model [model]` | Select or change the AI model. For models that support it, use left/right arrows to adjust effort level . With no argument, opens a picker that asks for confirmation when the conversation has prior output, since the next response re-reads the full history without cached context. Once confirmed, the change applies without waiting for the current response to finish | [DOC-commands] |
| `/passes` | Share a free week of Claude Code with friends. Only visible if your account is eligible | [DOC-commands] |
| `/permissions` | Manage allow, ask, and deny rules for tool permissions. Opens an interactive dialog where you can view rules by scope, add or remove rules, manage working directories, and review recent auto mode denials . Alias: /allowed-tools | [DOC-commands] |
| `/plan [description]` | Enter plan mode directly from the prompt. Pass an optional description to enter plan mode and immediately start with that task, for example /plan fix the auth bug | [DOC-commands] |
| `/plugin` | Manage Claude Code plugins | [DOC-commands] |
| `/powerup` | Discover Claude Code features through quick interactive lessons with animated demos | [DOC-commands] |
| `/pr-comments [PR]` | Removed in v2.1.91. Ask Claude directly to view pull request comments instead. On earlier versions, fetches and displays comments from a GitHub pull request; automatically detects the PR for the current branch, or pass a PR URL or number. Requires the gh CLI | [DOC-commands] |
| `/privacy-settings` | View and update your privacy settings. Only available for Pro and Max plan subscribers | [DOC-commands] |
| `/recap` | Generate a one-line summary of the current session on demand. See Session recap for the automatic recap that appears after you’ve been away | [DOC-commands] |
| `/release-notes` | View the changelog in an interactive version picker. Select a specific version to see its release notes, or choose to show all versions | [DOC-commands] |
| `/reload-plugins` | Reload all active plugins to apply pending changes without restarting. Reports counts for each reloaded component and flags any load errors | [DOC-commands] |
| `/remote-control` | Make this session available for remote control from claude.ai. Alias: /rc | [DOC-commands] |
| `/remote-env` | Configure the default remote environment for web sessions started with --remote | [DOC-commands] |
| `/rename [name]` | Rename the current session and show the name on the prompt bar. Without a name, auto-generates one from conversation history | [DOC-commands] |
| `/resume [session]` | Resume a conversation by ID or name, or open the session picker. Alias: /continue | [DOC-commands] |
| `/review [PR]` | Review a pull request locally in your current session. For a deeper cloud-based review, see /ultrareview | [DOC-commands] |
| `/rewind` | Rewind the conversation and/or code to a previous point, or summarize from a selected message. See checkpointing . Aliases: /checkpoint , /undo | [DOC-commands] |
| `/sandbox` | Toggle sandbox mode . Available on supported platforms only | [DOC-commands] |
| `/schedule [description]` | Create, update, list, or run routines . Claude walks you through the setup conversationally. Alias: /routines | [DOC-commands] |
| `/security-review` | Analyze pending changes on the current branch for security vulnerabilities. Reviews the git diff and identifies risks like injection, auth issues, and data exposure | [DOC-commands] |
| `/setup-bedrock` | Configure Amazon Bedrock authentication, region, and model pins through an interactive wizard. Only visible when CLAUDE_CODE_USE_BEDROCK=1 is set. First-time Bedrock users can also access this wizard from the login screen | [DOC-commands] |
| `/setup-vertex` | Configure Google Vertex AI authentication, project, region, and model pins through an interactive wizard. Only visible when CLAUDE_CODE_USE_VERTEX=1 is set. First-time Vertex AI users can also access this wizard from the login screen | [DOC-commands] |
| `/simplify [focus]` | Skill . Review your recently changed files for code reuse, quality, and efficiency issues, then fix them. Spawns three review agents in parallel, aggregates their findings, and applies fixes. Pass text to focus on specific concerns: /simplify focus on memory efficiency | [DOC-commands] |
| `/skills` | List available skills . Press t to sort by token count. Press Space to hide a skill from Claude or the / menu , then Enter to save | [DOC-commands] |
| `/stats` | Alias for /usage . Opens on the Stats tab | [DOC-commands] |
| `/status` | Open the Settings interface (Status tab) showing version, model, account, and connectivity. Works while Claude is responding, without waiting for the current response to finish | [DOC-commands] |
| `/statusline` | Configure Claude Code’s status line . Describe what you want, or run without arguments to auto-configure from your shell prompt | [DOC-commands] |
| `/stickers` | Order Claude Code stickers | [DOC-commands] |
| `/tasks` | List and manage background tasks. Also available as /bashes | [DOC-commands] |
| `/team-onboarding` | Generate a team onboarding guide from your Claude Code usage history. Claude analyzes your sessions, commands, and MCP server usage from the past 30 days and produces a markdown guide a teammate can paste as a first message to get set up quickly. For claude.ai subscribers on Pro, Max, Team, and Enterprise plans, also returns a share link teammates can open directly in Claude Code | [DOC-commands] |
| `/teleport` | Pull a Claude Code on the web session into this terminal: opens a picker, then fetches the branch and conversation. Also available as /tp . Requires a claude.ai subscription | [DOC-commands] |
| `/terminal-setup` | Configure terminal keybindings for Shift+Enter and other shortcuts. Only visible in terminals that need it, like VS Code, Cursor, Windsurf, Alacritty, or Zed | [DOC-commands] |
| `/theme` | Change the color theme. Includes an auto option that matches your terminal’s light or dark background, light and dark variants, colorblind-accessible (daltonized) themes, ANSI themes that use your terminal’s color palette, and any custom themes from ~/.claude/themes/ or plugins. Select New custom theme… to create one | [DOC-commands] |
| `/tui [default|fullscreen]` | Set the terminal UI renderer and relaunch into it with your conversation intact. fullscreen enables the flicker-free alt-screen renderer . With no argument, prints the active renderer | [DOC-commands] |
| `/ultraplan <prompt>` | Draft a plan in an ultraplan session, review it in your browser, then execute remotely or send it back to your terminal | [DOC-commands] |
| `/ultrareview [PR]` | Run a deep, multi-agent code review in a cloud sandbox with ultrareview . Includes 3 free runs on Pro and Max through May 5, 2026, then requires extra usage | [DOC-commands] |
| `/upgrade` | Open the upgrade page to switch to a higher plan tier | [DOC-commands] |
| `/usage` | Show session cost, plan usage limits, and activity stats. See the cost tracking guide for subscription-specific details. /cost and /stats are aliases | [DOC-commands] |
| `/vim` | Removed in v2.1.92. To toggle between Vim and Normal editing modes, use /config → Editor mode | [DOC-commands] |
| `/voice [hold|tap|off]` | Toggle voice dictation , or enable it in a specific mode. Requires a Claude.ai account | [DOC-commands] |
| `/web-setup` | Connect your GitHub account to Claude Code on the web using your local gh CLI credentials. /schedule prompts for this automatically if GitHub isn’t connected | [DOC-commands] |

### 5.3 Keyboard-shortcut architecture

- Claude Code supports a real keybinding system in `~/.claude/keybindings.json`, with live reload and per-context bindings. [DOC-keybindings]
- The docs say custom keybindings require Claude Code `v2.1.18` or later. [DOC-keybindings]
- The file schema is intentionally editor-friendly and can include both a schema URL and a docs URL. [DOC-keybindings]

```json
{
  "$schema": "https://www.schemastore.org/claude-code-keybindings.json",
  "$docs": "https://code.claude.com/docs/en/keybindings",
  "bindings": [
    {
      "context": "Chat",
      "bindings": {
        "ctrl+e": "chat:externalEditor",
        "ctrl+u": null
      }
    }
  ]
}
```
[DOC-keybindings]

### 5.4 Keybinding contexts

| Context | What it governs | Evidence |
|---|---|---|
| `Global` | App-wide bindings such as interrupt, exit, redraw, task-list toggle, and transcript toggle. | [DOC-keybindings] |
| `Chat` | Main chat input area. | [DOC-keybindings] |
| `Autocomplete` | Slash-command and `@` suggestion menus. | [DOC-keybindings] |
| `Settings` | The `/config` settings UI. | [DOC-keybindings] |
| `Confirmation` | Permission and confirmation dialogs. | [DOC-keybindings] |
| `Tabs` | Tab navigation components. | [DOC-keybindings] |
| `Help` | Help menu visibility and dismissal. | [DOC-keybindings] |
| `Transcript` | Verbose transcript viewer behavior. | [DOC-keybindings] |
| `HistorySearch` | Ctrl+R reverse-history search. | [DOC-keybindings] |
| `Task` | Background-task display and task-running state. | [DOC-keybindings] |
| `ThemePicker` | Theme dialog and syntax-highlighting toggle. | [DOC-keybindings] |
| `Attachments` | Image attachment navigation in selection dialogs. | [DOC-keybindings] |
| `Footer` | Footer indicator navigation for tasks/teams/diff/PR-like pills. | [DOC-keybindings] |
| `MessageSelector` | Rewind/summarize message selection. | [DOC-keybindings] |
| `DiffDialog` | Diff viewer navigation. | [DOC-keybindings] |
| `ModelPicker` | Effort-level slider inside the model picker. | [DOC-keybindings] |
| `Select` | Generic list/select controls. | [DOC-keybindings] |
| `Plugin` | Plugin browse/discover/manage dialog. | [DOC-keybindings] |
| `Scroll` | Fullscreen conversation scrolling and selection. | [DOC-keybindings] |
| `Doctor` | The `/doctor` diagnostics screen. | [DOC-keybindings] |

### 5.5 Most important default shortcuts to copy or adapt

| Shortcut | Default action | Why it matters | Evidence |
|---|---|---|---|
| `Ctrl+C` | Interrupt current input or generation. | Fast stop is a first-class control, not an obscure command. | [DOC-interactive] [DOC-keybindings] |
| `Ctrl+D` | Exit Claude Code. | Keeps Unix shell conventions. | [DOC-interactive] [DOC-keybindings] |
| `Ctrl+G` / `Ctrl+X Ctrl+E` | Open the prompt in the external editor. | Critical for multi-line or carefully edited prompts. | [DOC-interactive] [DOC-keybindings] |
| `Ctrl+L` | Redraw; in fullscreen, double-press clears the conversation. | Helps recover from terminal corruption without losing state. | [DOC-interactive] [DOC-keybindings] [DOC-fullscreen] |
| `Ctrl+O` | Toggle transcript viewer / verbose transcript. | Separates concise chat view from detailed execution trace. | [DOC-interactive] [DOC-keybindings] [DOC-fullscreen] |
| `Ctrl+R` | Reverse search command history. | Makes prompt history act like a shell. | [DOC-interactive] [DOC-keybindings] |
| `Ctrl+V` / `Alt+V` on Windows | Paste image. | Treats screenshots as first-class prompt artifacts. | [DOC-interactive] [DOC-keybindings] |
| `Ctrl+B` | Background current bash task or agent. | Turns long jobs into manageable background work. | [DOC-interactive] [DOC-keybindings] |
| `Ctrl+T` | Toggle task list. | Exposes live task decomposition without leaving chat. | [DOC-interactive] [DOC-keybindings] |
| `Shift+Tab` | Cycle permission modes. | Mode changes are quick and inline. | [DOC-interactive] [DOC-keybindings] |
| `Alt/Option+P` | Open model picker. | Useful because model changes do not require session restart. | [DOC-interactive] [DOC-keybindings] |
| `Alt/Option+T` | Toggle extended thinking. | Makes reasoning depth a lightweight runtime toggle. | [DOC-interactive] [DOC-keybindings] |
| `Alt/Option+O` | Toggle fast mode. | Exposes latency/cost tradeoff inline. | [DOC-interactive] [DOC-keybindings] |
| `Ctrl+J` | Insert newline without submitting. | Portable multiline fallback across terminals. | [DOC-interactive] [DOC-keybindings] [DOC-terminal] |
| `Tab` | Accept autocomplete suggestions. | Important for slash commands, shell-mode recall, and prompt suggestions. | [DOC-keybindings] [DOC-interactive] |
| `[` in transcript mode | Dump full conversation to native scrollback. | Smart escape hatch from alt-screen lock-in. | [DOC-interactive] [DOC-fullscreen] |
| `v` in transcript mode | Open transcript in `$VISUAL` / `$EDITOR`. | Excellent bridge from TUI to normal tools. | [DOC-interactive] [DOC-fullscreen] |
| `Esc` + `Esc` | Open rewind/summarize selector. | Makes checkpointing feel local and lightweight. | [DOC-interactive] [DOC-checkpointing] |

### 5.6 Chord support

- Claude Code absolutely supports multi-stroke chords, not just single keys; examples in the official docs include `Ctrl+X Ctrl+K` and `Ctrl+X Ctrl+E`. [DOC-interactive] [DOC-keybindings]
- That matters for chemsmart because it provides room for richer actions without exhausting single-key shortcuts. [DOC-keybindings]

## 6. Tool-call & diff display

### 6.1 Tool calls are visible, but at adjustable detail levels

- `Ctrl+O` toggles a more verbose transcript view that “shows detailed tool usage and execution”. [DOC-interactive]
- In normal view, some repeated or verbose tool activity collapses; the docs specifically note MCP calls can collapse to a line like “Called slack 3 times”. [DOC-interactive]
- In fullscreen mode, collapsed tool results are clickable and expand/collapse in place. [DOC-fullscreen]
- The changelog also mentions improved handling of collapsed long pastes, collapsed search/read summary badges, and collapsed tool content, which confirms the UI leans hard on progressive disclosure instead of flooding the transcript. [CHANGELOG-pastes] [CHANGELOG-focus-transcript]

### 6.2 Diff display

- `/diff` opens an interactive diff viewer. [DOC-commands]
- The docs say left/right arrows switch between the current git diff and individual Claude turns, while up/down browses files. [DOC-commands]
- The keybinding docs allocate a dedicated `DiffDialog` context with `previousSource`, `nextSource`, `previousFile`, `nextFile`, `viewDetails`, and `back`, which means diff review is a core navigation surface rather than an afterthought. [DOC-keybindings]
- Fullscreen theme tokens include separate colors for added lines, removed lines, dimmed context, and word-level add/remove highlights, which implies visual diff fidelity is part of the design language. [DOC-terminal]
- The fullscreen docs also say file paths printed after `Edit` or `Write` are clickable and open in a default application. [DOC-fullscreen]

### 6.3 Permission prompts

- Permission prompts are a central UI concept, not just a config file. The `/permissions` command opens an interactive dialog where you can inspect rules by scope, add/remove rules, manage working directories, and review recent denials. [DOC-permissions] [DOC-commands]
- The permission model distinguishes read-only tools, Bash, and file modification; read-only access is mostly prompt-free, Bash and edits require approval depending on mode. [DOC-permissions]
- The trust prompt I observed locally is a separate higher-level “is this workspace trusted?” gate that appears before session-level tool prompts. [LOCAL-trust]
- The permission system is modeful: `default`, `acceptEdits`, `plan`, `auto`, `dontAsk`, `bypassPermissions`. [DOC-permissions] [LOCAL-help]
- Claude Code treats approval ergonomics as a UI problem: `Shift+Tab` cycles modes, the footer/borders reflect the active mode, and settings/theme tokens explicitly style plan/accept-edits/auto borders differently. [DOC-interactive] [DOC-terminal] [DOC-keybindings]

### 6.4 Shell-mode display

- Prefixing a prompt with `!` enters a shell mode that runs a shell command directly without routing through the full natural-language planning loop. [DOC-interactive]
- Shell mode streams real-time progress/output, adds the command plus output back into the conversation context, supports `Ctrl+B` backgrounding, and has its own visual border token (`bashBorder`). [DOC-interactive] [DOC-terminal]
- That is a strong UI pattern for chemsmart: fast path for deterministic commands, conversational path for planning. [DOC-interactive]

## 7. Color / theme / typography choices, borders, frames

### 7.1 Themes and accessibility

- `/theme` offers `auto`, light/dark variants, colorblind-accessible “daltonized” variants, ANSI-palette themes, custom themes from `~/.claude/themes/`, and plugin-supplied themes. [DOC-commands] [DOC-terminal]
- The theme system supports live reload. [DOC-terminal]
- The docs treat theming as a tokenized design system, not as ad-hoc ANSI colors. [DOC-terminal]

### 7.2 Border and frame vocabulary

- The docs expose explicit theme tokens for the prompt border and for mode accents: `promptBorder`, `planMode`, `autoAccept`, and `bashBorder`. [DOC-terminal]
- Permission dialogs and pickers also have a `permission` color token, which implies framed dialogs are visually distinct from the baseline transcript. [DOC-terminal]
- Warning state also doubles as the “auto mode border” color token according to the theme docs. [DOC-terminal]
- This is notable: Claude Code uses border color as a **mode signifier**, not just as decoration. [DOC-terminal]

### 7.3 Fullscreen transcript backgrounds

- In fullscreen mode, user messages can have background fills (`userMessageBackground`, `userMessageBackgroundHover`). [DOC-terminal]
- Bash transcript entries and memory transcript entries have separate background tokens (`bashMessageBackgroundColor`, `memoryBackgroundColor`). [DOC-terminal]
- This suggests Claude Code visually separates message kinds by background fill in fullscreen, not just by prefix text. [DOC-terminal]

### 7.4 Typography and textual cues

- Claude Code is still terminal-first, so typography is limited to terminal capabilities: color, boldness, spacing, borders, and background fill. [DOC-fullscreen] [DOC-terminal]
- The changelog contains many fixes for emoji/wrapping/grapheme handling, which indicates the renderer cares about correct visual width and Unicode behavior. [CHANGELOG-2.1.132] [CHANGELOG-unicode]
- The theme docs also mention shimmer variants for spinner gradients and dedicated subagent colors, so animated/status color is part of the information design. [DOC-terminal]

### 7.5 Theme-token groups worth copying conceptually

| Token group | What it controls | Why it matters for chemsmart | Evidence |
|---|---|---|---|
| Text/accent | `claude`, `text`, `inverseText`, `inactive`, `subtle`, `suggestion`. | Core readability and suggestion visibility. | [DOC-terminal] |
| Permission/memory | `permission`, `remember`. | Makes system-state markers visually legible. | [DOC-terminal] |
| Status | `success`, `error`, `warning`, `merged`. | Direct mapping to job validation, critic warnings, and submission outcomes. | [DOC-terminal] |
| Input/modes | `promptBorder`, `planMode`, `autoAccept`, `bashBorder`, `ide`, `fastMode`. | Color-encodes mode changes without extra prose. | [DOC-terminal] |
| Diff | `diffAdded`, `diffRemoved`, dimmed variants, word-level variants. | Needed for dry-run input previews and job-spec diffs. | [DOC-terminal] |
| Fullscreen backgrounds | `userMessageBackground`, `bashMessageBackgroundColor`, `memoryBackgroundColor`, `selectionBg`. | Lets message kinds read as different object types. | [DOC-terminal] |
| Subagent colors | `blue_FOR_SUBAGENTS_ONLY` etc. | Useful if chemsmart later exposes multiple specialist agents. | [DOC-terminal] |

## 8. Multi-line input, paste handling, image input, file/dir autocomplete (`@file`, `#tag`)

### 8.1 Multi-line input

- `Enter` submits by default. [DOC-keybindings] [DOC-terminal]
- `Ctrl+J` is the universal newline insertion fallback. [DOC-keybindings] [DOC-terminal]
- `Shift+Enter` works natively in some terminals and can be installed/configured with `/terminal-setup` in others. [DOC-terminal] [DOC-commands]
- Another portable escape hatch is typing a backslash (`\`) and then pressing Enter. [DOC-terminal]
- For large or careful prompts, `Ctrl+G` / `Ctrl+X Ctrl+E` opens the prompt in an external editor; optionally Claude’s last reply can appear as `#`-commented context above the draft. [DOC-interactive] [DOC-keybindings]

### 8.2 Large paste handling

- If you paste more than 10,000 characters, Claude Code collapses the visible input to a `[Pasted text]` placeholder while still sending the full content on submit. [DOC-terminal]
- The docs explicitly recommend moving very large content into files and asking Claude to read the file instead of pasting entire logs or source files. [DOC-terminal]
- Changelog entries show active polish around paste behavior: large-image downscaling, bracketed-paste escape cleanup, multiline paste fixes, and a “Pasting…” footer hint for image paste. [CHANGELOG-2.1.132] [CHANGELOG-image-paste]

### 8.3 Image input

- Claude Code supports drag-and-drop images, clipboard image paste, and image-path references. [DOC-tutorials-images] [DOC-desktop-quickstart] [DOC-keybindings]
- `Ctrl+V` inserts an `[Image #N]` chip at the cursor, so images become first-class prompt references. [DOC-interactive] [DOC-keybindings]
- The changelog notes pasted/attached images are compressed/downscaled to fit token budgets and avoid session breakage. [CHANGELOG-image-paste] [CHANGELOG-image-budget]

### 8.4 File and directory autocomplete with `@`

- In the terminal UI, `@` is explicitly documented as **file path mention** that triggers file-path autocomplete. [DOC-interactive]
- The settings docs expose a `fileSuggestion` hook/command interface so projects can replace the built-in `@` suggestion source with a custom file-indexing command. [DOC-settings]
- The built-in `@` file suggestion returns up to 15 newline-separated path suggestions from the custom provider example. [DOC-settings]
- IDE docs state `@` mentions support fuzzy matching for files/folders and trailing slashes for folder references; I did not find a separate terminal-only doc that contradicts that, so I treat fuzzy folder matching as likely shared behavior but only **explicitly documented** on IDE surfaces. [DOC-vscode] [DOC-desktop]
- Remote Control clients can also query `@` autocomplete suggestions, per changelog. [CHANGELOG-remote-at]

### 8.5 `#tag`?

- **Unknown / not documented:** I did not find official terminal-CLI documentation for `#tag` autocomplete analogous to `@file`. [DOC-interactive] [DOC-commands] [DOC-keybindings] [LOCAL-help]
- `#` does appear in other UI contexts:
  - external-editor commented context uses `#` prefixes, [DOC-interactive]
  - fullscreen theming includes a `memoryBackgroundColor` token for `#` memory entries in the transcript. [DOC-terminal]
- My recommendation is: do **not** assume Claude Code has a user-facing `#tag` mention syntax in the CLI unless Anthropic documents it later. [DOC-interactive] [DOC-terminal]

## 9. State persistence: session history, resume (`--continue`, `--resume`), `~/.claude/projects/`

### 9.1 On-disk model

- Claude Code saves each message, tool use, and result as plaintext JSONL under `~/.claude/projects/<project>/<session-id>.jsonl`. [DOC-how] [DOC-claude-dir] [DOC-sessions]
- Large tool outputs can spill into `tool-results/` subdirectories. [DOC-claude-dir]
- File checkpoints live under `~/.claude/file-history/<session>/`. [DOC-claude-dir] [DOC-checkpointing]
- Prompt history lives in `~/.claude/history.jsonl`. [DOC-claude-dir]
- Debug logs, plans, session-env data, and tasks also have dedicated directories under `~/.claude/`. [DOC-claude-dir] [LOCAL-claude-dir]
- My local machine’s `~/.claude/` tree matches the official docs closely: I observed `projects/`, `file-history/`, `plans/`, `tasks/`, `debug/`, `history.jsonl`, and related state directories. [LOCAL-claude-dir] [DOC-claude-dir]

### 9.2 Resume semantics

- `claude --continue` resumes the most recent session in the current directory. [DOC-sessions] [LOCAL-help]
- `claude --resume` opens the interactive picker. [DOC-sessions] [LOCAL-help]
- `claude --resume <name>` resumes a named session directly, or opens the picker with pre-filled search when ambiguous. [DOC-sessions]
- `claude --from-pr <number>` resumes the session linked to a PR. [DOC-sessions] [LOCAL-help]
- `/resume` inside an active session switches to another conversation. [DOC-sessions] [DOC-commands]
- `--fork-session` branches the transcript into a new session ID while preserving the old session. [DOC-sessions] [LOCAL-help]

### 9.3 Session picker UX

- Each picker row shows: session name or summary/first prompt, time since activity, message count, and git branch. [DOC-sessions]
- `Ctrl+A` widens from current repo to all projects on the machine. [DOC-sessions]
- `Ctrl+W` widens from current worktree to all worktrees of the current repository. [DOC-sessions]
- `Ctrl+B` filters to the current git branch. [DOC-sessions]
- `Space` previews a session. [DOC-sessions]
- `Ctrl+R` renames a session from inside the picker. [DOC-sessions]
- Pasting a GitHub/GitLab/Bitbucket PR/MR URL into the picker search locates the session that created it. [DOC-sessions] [CHANGELOG-PR-resume-search]

### 9.4 Retention and privacy

- Official docs emphasize that transcripts/history are plaintext and not encrypted at rest; OS file permissions are the main protection. [DOC-claude-dir]
- Default cleanup is 30 days unless `cleanupPeriodDays` is changed. [DOC-sessions] [DOC-claude-dir]
- `claude project purge` can delete project-local Claude Code state. [DOC-cli] [CHANGELOG-project-purge]
- `CLAUDE_CODE_SKIP_PROMPT_HISTORY` or `--no-session-persistence` can suppress transcript writes. [DOC-sessions] [DOC-claude-dir]

## 10. Notable UX details worth borrowing for chemsmart

### 10.1 Why this matters specifically for chemsmart

- chemsmart’s current agent design is already typed, inspectable, sessioned, and HPC-aware. The frontier plan explicitly centers a conservative Plan-Approve-Execute loop with planner, critic, dry-run generation, runtime validation, and resumable session artifacts. [CHEMSMART-plan]
- The merged agent PRs reinforce that theme:
  - PR #15 added `extract_optimized_geometry` for cross-program geometry handoff, [CHEMSMART-PR15]
  - PR #16 made the critic validate **all** dry-run inputs and hard-reject malformed IRC flows, [CHEMSMART-PR16]
  - PR #17 routed Gaussian `opt+freq` through a single coherent job artifact, [CHEMSMART-PR17]
  - PR #23 added ORCA ab-initio routing, explicit scan-coordinate handling, stronger chemistry plausibility warnings, and rationale-depth improvements. [CHEMSMART-PR23]
- So the right Claude Code question is not “How do we mimic a generic coding chatbot?” but “Which terminal interaction patterns best expose a typed chemistry workflow with approval gates and resumable state?” [CHEMSMART-plan] [CHEMSMART-PR15] [CHEMSMART-PR16] [CHEMSMART-PR17] [CHEMSMART-PR23]

### 10.2 Borrowable patterns, mapped to chemsmart

| Claude Code pattern | Why it works in Claude Code | chemsmart adaptation | chemsmart grounding |
|---|---|---|---|
| Bottom-fixed input plus rich transcript | Keeps the operator oriented while streamed work continues. | Show planner/critic/tool events in the transcript while keeping the chemistry instruction box fixed at the bottom. | [CHEMSMART-plan] [DOC-fullscreen] |
| Mode-colored prompt border | Makes plan/auto/bash states instantly legible. | Map borders to `plan`, `dry-run ready`, `critic warn`, `runtime partial`, `submit armed`. | [CHEMSMART-plan] [CHEMSMART-PR16] [DOC-terminal] |
| Transcript verbosity toggle | Lets normal chat stay clean while power users inspect tool detail. | Default to concise step summaries, with `Ctrl+O`-style deep view for dry-run inputs, runtime validation payloads, and critic rationales. | [CHEMSMART-PR16] [CHEMSMART-PR23] [DOC-interactive] |
| Interactive diff viewer | Makes code edits reviewable per turn. | Use the same concept for Gaussian/ORCA input previews, scan-coordinate diffs, and before/after geometry handoff artifacts. | [CHEMSMART-PR15] [CHEMSMART-PR17] [CHEMSMART-PR23] [DOC-commands] |
| Session picker with branch/worktree filters | Makes many parallel efforts manageable. | Provide session resume by molecule/job/request, with filters for project, branch, scheduler target, and calculation family. | [CHEMSMART-plan] [DOC-sessions] |
| Background task model | Long-running tools do not block the conversation. | Run local jobs, queue polls, log-tail watchers, and HPC submission monitors as background tasks with IDs. | [CHEMSMART-plan] [DOC-interactive] |
| Workspace trust gate | Warns before granting agentic access. | Ask for trust/confirmation before allowing remote submission or writing into production scratch/project trees. | [DOC-security] [LOCAL-trust] [CHEMSMART-plan] |
| `!` shell-mode fast path | Lets users run deterministic shell actions without full NL planning. | Expose deterministic chemsmart helper commands (`!show latest opt log`, `!tail slurm`, `!open generated com`) with transcript capture. | [DOC-interactive] [CHEMSMART-AGENTS] |
| Prompt external-editor handoff | Supports long structured prompts. | Make method-constraint or scan-definition authoring editable in `$EDITOR`, especially for 1-based atom-index specifications. | [DOC-interactive] [CHEMSMART-AGENTS] [CHEMSMART-PR23] |
| Image chips + `@file` mentions | Handles screenshots and explicit file grounding cleanly. | Use for spectra screenshots, geometry images, route cards, and direct references to `.xyz`, `.log`, `.out`, `.com`, `.inp`. | [DOC-tutorials-images] [DOC-interactive] [CHEMSMART-AGENTS] |
| Task list and recap footer | Keeps long tasks legible after user absence. | Summarize where a job pipeline stands: build molecule -> recommend method -> build settings -> dry-run -> validate runtime -> critic -> run/submit. | [DOC-interactive] [CHEMSMART-plan] |
| Rewind/checkpoint UI | Encourages exploration without fear. | Allow “rewind to before submit,” “summarize noisy log tail,” and “branch this session for an alternate functional/basis guess.” | [DOC-checkpointing] [DOC-sessions] [CHEMSMART-plan] |

### 10.3 Chemistry-specific UX ideas directly inspired by Claude Code

- **Dry-run as diff-first artifact review.** Claude Code’s `/diff` and transcript expansion model suggest chemsmart should show generated Gaussian/ORCA inputs as reviewable artifacts before execution, not as opaque side effects. That lines up directly with the current `dry_run_input` and critic gating design. [DOC-commands] [DOC-interactive] [CHEMSMART-plan] [CHEMSMART-PR16]
- **Geometry handoff should read like a visible state transition.** PR #15 proves geometry extraction/handoff is a meaningful workflow boundary. Claude Code’s clickable tool-result transcript would map nicely to “optimized geometry extracted” cards that can be inspected or reused. [CHEMSMART-PR15] [DOC-fullscreen]
- **Make chemistry plausibility warnings visible but non-panicky.** PR #23 adds critic plausibility warnings (implausible TS/IRC, missing diffuse functions, missing unrestricted treatment). Claude Code’s warning-colored borders/footer badges are a better model than burying such warnings in plain logs. [CHEMSMART-PR23] [DOC-terminal]
- **Opt+freq and multi-step workflows should collapse intelligently.** PR #17 intentionally normalizes opt+freq into one coherent job. Claude Code’s habit of collapsing repeated tool calls into one-line summaries suggests chemsmart should collapse related substeps into one “composite workflow card” until the user asks for details. [CHEMSMART-PR17] [DOC-interactive]
- **Scan-definition entry needs editor-grade affordances.** PR #23’s explicit 1-based scan-coordinate validation is exactly the kind of error-prone structured input that benefits from Claude Code-style multiline editor handoff and rewind. [CHEMSMART-PR23] [CHEMSMART-AGENTS] [DOC-interactive]
- **HPC submission needs background status surfaces.** Claude Code’s background tasks and footer badges map well to queue submission, polling, and result arrival. The current chemsmart plan already treats remote runtime uncertainty as a gate; the UI should keep that state continuously visible. [DOC-interactive] [CHEMSMART-plan]
- **Session resume should restore chemistry context, not just text.** Claude Code’s on-disk session model plus chemsmart’s own session metadata suggest a chemsmart TUI should resume with molecule path, chosen method, generated input artifacts, runtime validation status, critic verdict, and remote job IDs. [DOC-sessions] [DOC-claude-dir] [CHEMSMART-plan]

## 11. Concrete references (URLs, doc paths, changelog entries)

### 11.1 Official docs

- Claude Code overview: `https://code.claude.com/docs/en/overview`
- CLI reference: `https://code.claude.com/docs/en/cli-usage`
- Commands reference: `https://code.claude.com/docs/en/commands`
- Interactive mode: `https://code.claude.com/docs/en/interactive-mode`
- Customize keyboard shortcuts: `https://code.claude.com/docs/en/keybindings`
- Settings: `https://code.claude.com/docs/en/settings`
- Configuration (system prompt note / settings model): `https://code.claude.com/docs/en/configuration`
- Permissions: `https://code.claude.com/docs/en/permissions`
- Sandboxing: `https://code.claude.com/docs/en/sandboxing`
- Hooks: `https://code.claude.com/docs/en/hooks`
- MCP: `https://code.claude.com/docs/en/mcp`
- Skills / slash-command evolution: `https://code.claude.com/docs/en/slash-commands`
- Sessions: `https://code.claude.com/docs/en/sessions`
- Checkpointing: `https://code.claude.com/docs/en/checkpointing`
- Context window: `https://code.claude.com/docs/en/context-window`
- Fullscreen rendering: `https://code.claude.com/docs/en/fullscreen`
- Terminal configuration: `https://code.claude.com/docs/en/terminal-config`
- Tutorials / images: `https://code.claude.com/docs/en/tutorials`
- `.claude` directory / application data: `https://code.claude.com/docs/en/claude-directory`
- VS Code integration docs (useful for `@` mention semantics): `https://code.claude.com/docs/en/vs-code`
- Desktop docs (useful for attachment semantics): `https://code.claude.com/docs/en/desktop`
- How Claude Code works: `https://code.claude.com/docs/en/how-claude-code-works`

### 11.2 Public repo / examples

- Repo root: `https://github.com/anthropics/claude-code`
- README: `https://github.com/anthropics/claude-code/blob/main/README.md`
- CHANGELOG: `https://github.com/anthropics/claude-code/blob/main/CHANGELOG.md`
- Plugins directory: `https://github.com/anthropics/claude-code/tree/main/plugins`
- Examples/settings README: `https://github.com/anthropics/claude-code/tree/main/examples/settings`
- Hook validator example: `https://github.com/anthropics/claude-code/blob/main/examples/hooks/bash_command_validator_example.py`
- Official plugins README: `https://github.com/anthropics/claude-code/blob/main/plugins/README.md`
- Code-review plugin README: `https://github.com/anthropics/claude-code/blob/main/plugins/code-review/README.md`
- Commit-commands plugin README: `https://github.com/anthropics/claude-code/blob/main/plugins/commit-commands/README.md`

### 11.3 Public issues useful for stack inference

- Ink-era error report: `https://github.com/anthropics/claude-code/issues/404`
- Native Bun-bundled binary crash report: `https://github.com/anthropics/claude-code/issues/25630`

### 11.4 chemsmart references for UI adaptation

- `/Users/hongjiseung/developer/chemsmart/AGENTS.md`
- `/Users/hongjiseung/.agent-orchestrator/projects/chemsmart/frontier_improvement_plan.md`
- Fork PR #15: `https://github.com/Hongjiseung-ROK/chemsmart/pull/15`
- Fork PR #16: `https://github.com/Hongjiseung-ROK/chemsmart/pull/16`
- Fork PR #17: `https://github.com/Hongjiseung-ROK/chemsmart/pull/17`
- Fork PR #23: `https://github.com/Hongjiseung-ROK/chemsmart/pull/23`

## 12. 5–10 bullet “what to copy for chemsmart” recommendations

- Copy Claude Code’s **dual-renderer idea**: keep a plain scrollback-safe mode plus a richer fullscreen mode for long chemistry sessions. [DOC-commands] [DOC-fullscreen]
- Copy the **mode-colored input border** and map it to chemistry workflow states (`plan`, `dry-run ready`, `critic warn`, `runtime partial`, `submit armed`). [DOC-terminal] [CHEMSMART-plan] [CHEMSMART-PR16]
- Copy the **transcript/detail split** (`Ctrl+O`-style verbose toggle) so normal users see concise workflow cards while power users can inspect dry-run inputs, critic payloads, and runtime validation details. [DOC-interactive] [CHEMSMART-plan]
- Copy the **interactive diff viewer pattern**, but aim it at generated `.com` / `.inp` / scan definitions / geometry-handoff artifacts rather than only source-code edits. [DOC-commands] [CHEMSMART-PR15] [CHEMSMART-PR17] [CHEMSMART-PR23]
- Copy **background tasks** for long local runs, scheduler submissions, and queue polling, with footer badges and resumable IDs. [DOC-interactive] [CHEMSMART-plan]
- Copy **session resume + branching**, but make the primary resume keys chemistry-native: molecule path, method, job family, queue target, critic verdict, and artifact set. [DOC-sessions] [DOC-claude-dir] [CHEMSMART-plan]
- Copy **external-editor prompt handoff** for structured chemical inputs like scan coordinates, route constraints, solvent/method notes, and long approval comments. [DOC-interactive] [CHEMSMART-PR23]
- Copy **image/file chips** (`Ctrl+V`, `@file`) so users can drop in spectra screenshots, geometry images, and exact structure/log paths without bloating the main prompt text. [DOC-tutorials-images] [DOC-interactive]
- Copy the **workspace trust / permission framing**, but specialize the scary operations to chemsmart realities: writing over previous job artifacts, submitting to HPC, or reusing potentially stale geometries. [LOCAL-trust] [DOC-security] [CHEMSMART-PR15]
- Do **not** copy Claude Code’s GitHub/PR-specific chrome literally; repurpose those footer/status concepts for chemistry entities such as active job ID, geometry handoff status, dry-run count, remote runtime certainty, and critic severity. [DOC-interactive] [CHEMSMART-plan] [CHEMSMART-PR23]

## Appendix A. Small public-repo snippets worth re-reading

### A.1 Strict permissions example

```json
{
  "permissions": {
    "disableBypassPermissionsMode": "disable",
    "ask": ["Bash"],
    "deny": ["WebSearch", "WebFetch"]
  },
  "allowManagedPermissionRulesOnly": true,
  "allowManagedHooksOnly": true
}
```
[REPO-settings-strict]

### A.2 Hook example: block/shape Bash behavior before execution

```python
if issues:
    for message in issues:
        print(f"• {message}", file=sys.stderr)
    # Exit code 2 blocks tool call and shows stderr to Claude
    sys.exit(2)
```
[REPO-hook-example]

### A.3 Plugin structure summary

```text
plugin-name/
├── .claude-plugin/plugin.json
├── commands/
├── agents/
├── skills/
├── hooks/
├── .mcp.json
└── README.md
```
[REPO-plugins-readme]

## Appendix B. UI-relevant changelog timeline I would re-read before implementation

- `2.1.133`: worktree-base setting, focus-mode improvements, help text updates, subagent skill discovery fixes. [CHANGELOG-2.1.133]
- `2.1.132`: disable-alt-screen env var, image-paste footer hint, slash autocomplete sizing fix, fullscreen blank-screen fix, graceful SIGINT shutdown. [CHANGELOG-2.1.132]
- `2.1.129`: package-manager auto-update, plugin-manifest tweaks, history picker scope changes, session title chip fix, `/branch` session-id fix. [CHANGELOG-2.1.129]
- `2.1.118+`: custom themes and live-reload theming in `~/.claude/themes/`. [DOC-terminal] [CHANGELOG-themes]
- `2.1.110` neighborhood: `/tui fullscreen` and no-flicker renderer maturation. [DOC-fullscreen] [CHANGELOG-tui]
- `2.1.89+`: fullscreen rendering introduced as research preview. [DOC-fullscreen]
- Resume-related clusters: multiple entries across the changelog improve large-session loading, corruption recovery, picker behavior, and cross-worktree resume semantics. That is a signal to budget engineering time for resume robustness early. [CHANGELOG-resume-cluster]
- Paste/image clusters: multiple entries harden bracketed paste, multiline paste, image downscaling, image token budgeting, and image attachment persistence. This is a signal that “rich prompt input” is harder than it looks. [CHANGELOG-image-paste] [CHANGELOG-image-budget]
- Fullscreen clusters: many entries address selection, scrolling, link opening, wrapped URLs, transcript dumping, and tmux/VS Code rendering edge cases. This is a signal that a chemistry TUI should start with a plain renderer and only add a no-flicker path once the state model is stable. [DOC-fullscreen] [CHANGELOG-fullscreen-cluster]

## Appendix C. Citation key
- **DOC-overview** — Official overview: https://code.claude.com/docs/en/overview
- **DOC-cli** — CLI reference: https://code.claude.com/docs/en/cli-usage
- **DOC-commands** — Commands reference: https://code.claude.com/docs/en/commands
- **DOC-interactive** — Interactive mode: https://code.claude.com/docs/en/interactive-mode
- **DOC-keybindings** — Keybindings: https://code.claude.com/docs/en/keybindings
- **DOC-settings** — Settings: https://code.claude.com/docs/en/settings
- **DOC-configuration-system-prompt** — Configuration docs (system prompt not published): https://code.claude.com/docs/en/configuration
- **DOC-permissions** — Permissions: https://code.claude.com/docs/en/permissions
- **DOC-security** — Security docs: https://docs.claude.com/en/docs/claude-code/security
- **DOC-sandbox** — Sandboxing: https://code.claude.com/docs/en/sandboxing
- **DOC-hooks** — Hooks: https://code.claude.com/docs/en/hooks
- **DOC-mcp** — MCP: https://code.claude.com/docs/en/mcp
- **DOC-slash** — Skills / slash commands: https://code.claude.com/docs/en/slash-commands
- **DOC-skills** — Skills page (same canonical slash-command/skills documentation): https://code.claude.com/docs/en/slash-commands
- **DOC-sessions** — Sessions: https://code.claude.com/docs/en/sessions
- **DOC-checkpointing** — Checkpointing: https://code.claude.com/docs/en/checkpointing
- **DOC-context-window** — Context window: https://code.claude.com/docs/en/context-window
- **DOC-fullscreen** — Fullscreen rendering: https://code.claude.com/docs/en/fullscreen
- **DOC-terminal** — Terminal configuration and theme tokens: https://code.claude.com/docs/en/terminal-config
- **DOC-tutorials-images** — Tutorials / images: https://code.claude.com/docs/en/tutorials
- **DOC-claude-dir** — `.claude` directory / application data: https://code.claude.com/docs/en/claude-directory
- **DOC-vscode** — VS Code integration docs: https://code.claude.com/docs/en/vs-code
- **DOC-desktop** — Desktop docs: https://code.claude.com/docs/en/desktop
- **DOC-desktop-quickstart** — Desktop quickstart: https://code.claude.com/docs/en/desktop-quickstart
- **DOC-how** — How Claude Code works: https://code.claude.com/docs/en/how-claude-code-works
- **REPO-root** — Public repo root: https://github.com/anthropics/claude-code
- **REPO-readme** — README: https://github.com/anthropics/claude-code/blob/main/README.md
- **REPO-settings-examples** — Settings examples dir: https://github.com/anthropics/claude-code/tree/main/examples/settings
- **REPO-settings-strict** — Strict settings example: https://github.com/anthropics/claude-code/blob/main/examples/settings/settings-strict.json
- **REPO-hook-example** — Hook example: https://github.com/anthropics/claude-code/blob/main/examples/hooks/bash_command_validator_example.py
- **REPO-plugins-readme** — Plugins README: https://github.com/anthropics/claude-code/blob/main/plugins/README.md
- **ISSUE-404** — Ink/raw-mode issue: https://github.com/anthropics/claude-code/issues/404
- **ISSUE-25630** — Native Bun-bundled binary issue: https://github.com/anthropics/claude-code/issues/25630
- **CHANGELOG-2.1.133** — CHANGELOG v2.1.133: https://github.com/anthropics/claude-code/blob/main/CHANGELOG.md#2133
- **CHANGELOG-2.1.132** — CHANGELOG v2.1.132: https://github.com/anthropics/claude-code/blob/main/CHANGELOG.md#2132
- **CHANGELOG-2.1.129** — CHANGELOG v2.1.129: https://github.com/anthropics/claude-code/blob/main/CHANGELOG.md#2129
- **CHANGELOG-themes** — Theme-related changelog entries in main CHANGELOG: https://github.com/anthropics/claude-code/blob/main/CHANGELOG.md
- **CHANGELOG-tui** — TUI/fullscreen changelog entries in main CHANGELOG: https://github.com/anthropics/claude-code/blob/main/CHANGELOG.md
- **CHANGELOG-project-purge** — Project purge entry in CHANGELOG: https://github.com/anthropics/claude-code/blob/main/CHANGELOG.md
- **CHANGELOG-PR-resume-search** — PR URL resume-search entry in CHANGELOG: https://github.com/anthropics/claude-code/blob/main/CHANGELOG.md
- **CHANGELOG-remote-at** — Remote Control `@` autocomplete entry in CHANGELOG: https://github.com/anthropics/claude-code/blob/main/CHANGELOG.md
- **CHANGELOG-image-paste** — Image/paste hardening entries in CHANGELOG: https://github.com/anthropics/claude-code/blob/main/CHANGELOG.md
- **CHANGELOG-image-budget** — Image token-budget entries in CHANGELOG: https://github.com/anthropics/claude-code/blob/main/CHANGELOG.md
- **CHANGELOG-fullscreen-cluster** — Fullscreen rendering cluster in CHANGELOG: https://github.com/anthropics/claude-code/blob/main/CHANGELOG.md
- **CHANGELOG-pastes** — Collapsed paste / paste UX entries in CHANGELOG: https://github.com/anthropics/claude-code/blob/main/CHANGELOG.md
- **CHANGELOG-focus-transcript** — Focus/transcript related entries in CHANGELOG: https://github.com/anthropics/claude-code/blob/main/CHANGELOG.md
- **CHANGELOG-unicode** — Unicode/wrapping rendering entries in CHANGELOG: https://github.com/anthropics/claude-code/blob/main/CHANGELOG.md
- **CHANGELOG-escape-interrupt** — Escape-interrupt entries in CHANGELOG: https://github.com/anthropics/claude-code/blob/main/CHANGELOG.md
- **CHANGELOG-resume-cluster** — Resume/session-loading robustness entries across CHANGELOG: https://github.com/anthropics/claude-code/blob/main/CHANGELOG.md
- **CHANGELOG-2.1.125-ish** — Vim/Esc interrupt behavior entry in CHANGELOG: https://github.com/anthropics/claude-code/blob/main/CHANGELOG.md
- **LOCAL-help** — Local observation on 2026-05-08: `claude --version`, `claude --help`, local install path `/Users/hongjiseung/.local/share/claude/versions/2.1.133`.
- **LOCAL-binary** — Local observation on 2026-05-08: install directory contains a single native executable at `/Users/hongjiseung/.local/share/claude/versions/2.1.133`.
- **LOCAL-trust** — Local observation on 2026-05-08: scratch-directory interactive launch displayed a workspace trust/safety modal before entering the session.
- **LOCAL-claude-dir** — Local observation on 2026-05-08: `~/.claude/` contains `projects/`, `file-history/`, `plans/`, `tasks/`, `debug/`, `history.jsonl`, and related state directories.
- **CHEMSMART-AGENTS** — chemsmart repo conventions and domain map: `/Users/hongjiseung/developer/chemsmart/AGENTS.md`.
- **CHEMSMART-plan** — chemsmart frontier plan: `/Users/hongjiseung/.agent-orchestrator/projects/chemsmart/frontier_improvement_plan.md`.
- **CHEMSMART-PR15** — Merged fork PR #15: https://github.com/Hongjiseung-ROK/chemsmart/pull/15
- **CHEMSMART-PR16** — Merged fork PR #16: https://github.com/Hongjiseung-ROK/chemsmart/pull/16
- **CHEMSMART-PR17** — Merged fork PR #17: https://github.com/Hongjiseung-ROK/chemsmart/pull/17
- **CHEMSMART-PR23** — Merged fork PR #23: https://github.com/Hongjiseung-ROK/chemsmart/pull/23


## Appendix D. Useful action tables to keep beside the implementation plan

### D.1 Permission modes

| Mode | Summary | Why it matters for chemsmart | Evidence |
|---|---|---|---|
| `default` | Standard prompt-on-first-use behavior. | Good baseline for local chemistry planning sessions. | [DOC-permissions] |
| `acceptEdits` | Auto-accepts file edits and common filesystem commands in the working directory. | Could map to low-risk artifact-generation sessions. | [DOC-permissions] |
| `plan` | Read-only exploration; no source-file edits. | Strong precedent for chemsmart planner-only mode before dry-run approval. | [DOC-permissions] |
| `auto` | Auto-approves with background safety checks; research preview. | Probably too aggressive for HPC submission, but interesting for local read-only prep work. | [DOC-permissions] |
| `dontAsk` | Auto-denies unless explicitly pre-approved. | Useful for locked-down review or teaching demos. | [DOC-permissions] |
| `bypassPermissions` | Skips prompts except root/home delete circuit-breakers. | Only suitable for sandbox/container automation, not ordinary chemsmart HPC use. | [DOC-permissions] |

### D.2 Session picker shortcuts

| Shortcut | Session-picker behavior | Evidence |
|---|---|---|
| `↑` / `↓` | Move between sessions. | [DOC-sessions] |
| `→` / `←` | Expand or collapse grouped session branches. | [DOC-sessions] |
| `Enter` | Resume highlighted session. | [DOC-sessions] |
| `Space` | Preview highlighted session. | [DOC-sessions] |
| `Ctrl+V` | Alternate preview on terminals that do not reserve it for paste. | [DOC-sessions] |
| `Ctrl+R` | Rename highlighted session. | [DOC-sessions] |
| `/` | Enter search mode. | [DOC-sessions] |
| printable character | Also enters search mode. | [DOC-sessions] |
| `Ctrl+A` | Widen to all projects on the machine. | [DOC-sessions] |
| `Ctrl+W` | Widen to all worktrees of the current repository. | [DOC-sessions] |
| `Ctrl+B` | Filter to current git branch. | [DOC-sessions] |
| `Esc` | Exit picker or search mode. | [DOC-sessions] |

### D.3 Transcript-mode shortcuts in fullscreen rendering

| Key | Transcript-mode action | Evidence |
|---|---|---|
| `/` | Open search. | [DOC-fullscreen] |
| `n` / `N` | Next / previous match. | [DOC-fullscreen] |
| `j` / `k` or `↑` / `↓` | Scroll one line. | [DOC-fullscreen] |
| `g` / `G` or `Home` / `End` | Jump to top / bottom. | [DOC-fullscreen] |
| `Ctrl+u` / `Ctrl+d` | Half-page scroll. | [DOC-fullscreen] |
| `Ctrl+b` / `Ctrl+f` or `Space` / `b` | Full-page scroll. | [DOC-fullscreen] |
| `[` | Dump full conversation into native terminal scrollback. | [DOC-fullscreen] [DOC-interactive] |
| `v` | Write transcript to a temp file and open in editor. | [DOC-fullscreen] [DOC-interactive] |
| `Ctrl+O`, `Esc`, or `q` | Exit transcript mode. | [DOC-fullscreen] [DOC-interactive] |

### D.4 Fullscreen mouse behavior

| Interaction | Fullscreen behavior | Evidence |
|---|---|---|
| click in prompt | Reposition cursor inside current draft. | [DOC-fullscreen] |
| click collapsed tool result | Expand/collapse tool call plus result. | [DOC-fullscreen] |
| click URL or file path | Open browser or default application. | [DOC-fullscreen] |
| click-drag | Select transcript text. | [DOC-fullscreen] |
| double-click | Select word. | [DOC-fullscreen] |
| triple-click | Select line. | [DOC-fullscreen] |
| mouse release | Copy selection automatically when copy-on-select is enabled. | [DOC-fullscreen] |
| wheel scroll | Scroll conversation inside Claude Code. | [DOC-fullscreen] |

### D.5 Prompt-entry affordances

| Input affordance | Behavior | Evidence |
|---|---|---|
| `/` at start | Command/skill browser. | [DOC-commands] [DOC-interactive] |
| `!` at start | Shell mode. | [DOC-interactive] |
| `@` | File-path mention / autocomplete trigger. | [DOC-interactive] |
| `Enter` | Submit prompt. | [DOC-keybindings] [DOC-terminal] |
| `Ctrl+J` | Insert newline. | [DOC-keybindings] [DOC-terminal] |
| `Shift+Enter` | Insert newline where supported or after `/terminal-setup`. | [DOC-terminal] |
| `Ctrl+G` / `Ctrl+X Ctrl+E` | Open in external editor. | [DOC-interactive] [DOC-keybindings] |
| `Ctrl+V` | Paste image. | [DOC-keybindings] [DOC-interactive] |
| >10k-char paste | Collapse visible draft to `[Pasted text]` placeholder. | [DOC-terminal] |

## Appendix E. Local-observation notes from 2.1.133

- `claude --help` still reads like a mature CLI application rather than a thin launcher: it exposes session flags, permission modes, plugin controls, worktree creation, remote control, print/stream modes, agent selection, and MCP loading. [LOCAL-help]
- The first interactive screen I hit in a scratch repo was not a marketing splash screen; it was a **trust decision** that foregrounded file-read/edit/execute power. That feels very intentional and worth emulating for chemsmart’s remote-submit path. [LOCAL-trust]
- The local install layout on this machine is extremely sparse at the executable level: one large versioned binary plus ordinary user-state directories under `~/.claude/`. This supports the claim that the interactive renderer is shipped as an opaque native artifact while user state remains transparent and file-based. [LOCAL-binary] [LOCAL-claude-dir]
- The state directory is not hidden behind a database. That makes sessions, tasks, prompt history, and file checkpoints auditable and scriptable, which is directly relevant to a scientific workflow tool. [LOCAL-claude-dir] [DOC-claude-dir]
- The docs and the local state tree line up closely enough that I trust the public docs as the main UI source of truth, while using local observation mostly to confirm tone and filesystem shape. [LOCAL-claude-dir] [DOC-claude-dir] [LOCAL-trust]

## Appendix F. One concrete chemsmart UI storyboard using Claude Code patterns

1. User launches `chemsmart chat` in a project directory containing structures and prior outputs. [CHEMSMART-AGENTS]
2. TUI opens in default mode with a bottom-fixed prompt, a small footer, and a trust banner if HPC submission is configured. [LOCAL-trust] [DOC-fullscreen] [DOC-security]
3. User asks: “optimize this structure and then run ORCA single-point.”
4. Transcript shows a concise workflow card: `build_molecule -> recommend_method -> build_gaussian_settings -> build_job -> dry_run_input`. [CHEMSMART-plan]
5. Prompt border stays in plan color until dry-run artifacts are ready. [DOC-terminal] [CHEMSMART-plan]
6. User hits `Ctrl+O` to expand the generated Gaussian route line and inspect the exact dry-run input. [DOC-interactive] [DOC-keybindings]
7. Critic warning appears as a yellow footer/border badge if runtime fields are partial or chemistry plausibility is weak. [CHEMSMART-PR16] [CHEMSMART-PR23] [DOC-terminal]
8. User hits `Ctrl+G` to open a multi-line editor and add a solvent note plus atom-indexed scan definition. [DOC-interactive] [CHEMSMART-PR23]
9. TUI shows a diff-style preview for the modified `.com` or `.inp` artifact before execution. [DOC-commands] [CHEMSMART-PR17] [CHEMSMART-PR23]
10. Local optimization is started as a background task; footer shows task ID and latest status. [DOC-interactive]
11. When optimization finishes, a geometry-handoff card appears, analogous to a clickable tool result. [CHEMSMART-PR15] [DOC-fullscreen]
12. User accepts ORCA single-point generation from the extracted optimized geometry. [CHEMSMART-PR15]
13. HPC submission is shown as a distinct “dangerous” mode transition with stronger confirmation styling than a normal file edit. [DOC-permissions] [DOC-terminal] [CHEMSMART-plan]
14. Queue polling continues in the background; the user can ask side questions or inspect transcript details without blocking. [DOC-interactive]
15. If the user comes back later, `/resume`-style session picking restores not just text but the exact artifact set, queue IDs, and critic state. [DOC-sessions] [CHEMSMART-plan]


## Appendix G. Implementation cautions chemsmart should learn from Claude Code’s changelog

- Resume is not a “small feature”; Claude Code’s changelog has many independent fixes for corruption recovery, picker behavior, cross-worktree lookup, stale names, large-session performance, and interrupted writes. Build resumability early and test it hard. [CHANGELOG-resume-cluster]
- Rich paste is not a “small feature”; bracketed paste, giant pastes, image paste, and pasted leading slash/escape-sequence handling each generated separate fixes. Treat paste/input as a dedicated subsystem. [CHANGELOG-2.1.132] [CHANGELOG-image-paste]
- Fullscreen alt-screen rendering is not a “small feature”; selection, scrollback handoff, URL opening, mouse capture, tmux, VS Code terminals, and Unicode wrapping all generated follow-up fixes. Do not make fullscreen the only rendering path at first. [DOC-fullscreen] [CHANGELOG-fullscreen-cluster]
- Theme tokens and border semantics pay off because they let the UI communicate mode changes without verbose prose. That will matter in chemsmart when the user is switching between plan, dry-run review, local execution, and remote submit. [DOC-terminal] [CHEMSMART-plan]
- Visible, auditable local state under a user directory is a strength for scientific tooling. It supports provenance, postmortem review, and reproducibility better than opaque in-memory session state. [DOC-claude-dir] [LOCAL-claude-dir]
- The split between concise chat, verbose transcript, and focus view is a strong pattern for scientific agents, where users alternate between “just tell me what’s happening” and “show me the exact artifact.” [DOC-interactive] [DOC-fullscreen]
- Background tasks plus resumable session state are especially relevant for HPC and long-running chemistry jobs; Claude Code’s UI patterns are unusually compatible with that domain. [DOC-interactive] [DOC-sessions] [CHEMSMART-plan]
- The `@file` model should be copied, but only if the TUI also has robust artifact typing and discoverability; chemistry users care about exact `.xyz`, `.log`, `.out`, `.chk`, `.com`, and `.inp` provenance. [DOC-interactive] [DOC-settings] [CHEMSMART-AGENTS]
- I would **not** copy Claude Code’s Git-centric assumptions around PR state, commit attribution, and branch naming into the first chemsmart TUI. Replace those with domain-native objects such as molecule, job, geometry, method, runtime target, and scheduler submission. [DOC-commands] [DOC-settings] [CHEMSMART-AGENTS]
- I would copy Claude Code’s idea that every dangerous action should have a visible surface, a resumable transcript, and an obvious way to inspect the exact artifact being acted upon. That principle fits chemsmart almost perfectly. [DOC-permissions] [DOC-fullscreen] [CHEMSMART-plan]


## Appendix H. Bottom-line summary in one screenful

- Claude Code’s public value for chemsmart is **interaction design**, not reusable renderer code. [REPO-root] [LOCAL-binary]
- The most important patterns are: dual renderer, bottom-fixed input, verbose-transcript toggle, diff-first artifact review, background tasks, resumable sessions, and mode-colored borders. [DOC-fullscreen] [DOC-interactive] [DOC-terminal] [DOC-sessions]
- The most important chemistry adaptations are: surface dry-run inputs, geometry handoff, critic warnings, scan-definition editing, and HPC submit state as first-class UI objects. [CHEMSMART-PR15] [CHEMSMART-PR16] [CHEMSMART-PR17] [CHEMSMART-PR23] [CHEMSMART-plan]
