# chemsmart agent TUI

The Textual TUI is the recommended interactive surface for `chemsmart agent`.
It keeps the CLI workflow visible while adding conversation memory, runtime
evidence, slash commands, and resumable sessions.

## Launch

```bash
chemsmart agent
chemsmart agent --plain
```

`--plain` keeps the UI inline and conservative for terminals such as HPC login
nodes.

![chemsmart agent dry-run TUI](../tests/agent/tui/snapshots/dry_run_input.svg)

## Slash Palette

Type `/` to open the command palette. Continue typing to filter commands; for
example `/d` shows matching commands such as `/doctor`, `/dryrun`, and `/deny`.
Use the arrow keys to move, `Tab` to complete, and `Enter` to run. Commands that
are not valid in the current phase remain visible with an unavailable reason.
The same entries remain selectable with the mouse.

![chemsmart slash command palette](../tests/agent/tui/snapshots/slash_help.svg)

Current high-value commands:

| Command | Purpose |
|---|---|
| `/mode` | Show provider routing information; no mode switch is required |
| `/init` | Build and validate a project YAML from a reported method |
| `/write-project [name]` | Approve writing or replacing the latest validated workspace YAML |
| `/tools` | List registered tools |
| `/doctor` | Run provider/runtime diagnostics |
| `/sessions` | Browse saved sessions |
| `/resume <session-id>` | Resume a saved session |
| `/runs` or `/jobs` | Inspect local calculations and scheduler jobs |
| `/run` | Execute the latest validated `chemsmart run` command |
| `/submit` | Execute the latest validated `chemsmart sub` command |

Provider roles are automatic in the same transcript. Local providers use the
CLI synthesis lane; API providers use the unified tool loop. A command becomes
executable only after semantic and intent gates pass, generated-input evidence
exists, and the workspace project YAML remains unchanged. `/run` or `/submit`
then acts as the user's explicit approval, and the execution tool repeats the
semantic gate before starting the real command.

## Keyboard Workflow

The interactive TUI is keyboard-first without removing existing mouse actions.

| Shortcut | Action |
|---|---|
| `F1` | Open the shortcut reference |
| `Ctrl+O` | Toggle compact and expanded transcript details |
| `Ctrl+T` | Inspect the current turn's tool activity |
| `Ctrl+B` | Open the calculation and scheduler monitor |
| `Shift+Tab` | Inspect the active workspace project YAML |
| `Ctrl+R` | Recall the previous request |
| `@` | Open the workspace file picker |
| `Ctrl+G` | Edit the current draft in `$EDITOR` |
| `Tab` while busy | Queue or restore one follow-up request |
| `y`, `s`, `n`, `r` | Approve once, approve for session, deny, or revise |
| `Esc` | Close an overlay and return to the request composer |

A queued request starts only after the active turn completes or fails. It stays
paused while the agent is waiting for clarification or approval.

## Calculation Monitor

Approved local calculations run in a background worker, independently of the
agent conversation worker. The persistent strip between the transcript and
composer shows the program, job kind, input label, elapsed time, and the latest
observed stage such as an SCF cycle, optimization cycle, frequency analysis, or
scan point. It never invents a percentage.

Press `Ctrl+B`, or use `/runs` or `/jobs`, to open the combined calculation
monitor:

| Key | Action |
|---|---|
| `Up` / `Down` | Select a local calculation or scheduler job |
| `Enter` | Show the structured receipt |
| `L` | Follow the bounded raw log |
| `/` | Search the selected log |
| `E` | Extract the latest chemistry result |
| `C` | Request cancellation |
| `Esc` | Return to the request composer |

The main transcript receives one mutable calculation receipt rather than a raw
stdout dump. A completed ORCA single point reports final energy, SCF cycles,
normal termination, method/basis, elapsed time, and output path. Opt, frequency,
scan, NEB, and QMMM jobs add their relevant convergence or region fields.
Process exit status and Gaussian/ORCA normal termination are separate checks;
an exit code of zero without a real output/termination marker is not reported
as a successful chemistry calculation.

Calculation receipts remain in the turn that started them. A later completion
updates that cell in place and cannot move below a newer conversation turn.
Receipts and bounded logs are stored below the agent session calculation store
and restored after a TUI restart. Ask `diagnose the latest calculation` (or the
equivalent Korean request) to inspect the most recent unambiguous result without
re-entering its output path.

Optional key overrides live in the existing workspace-independent agent config:

```yaml
tui:
  tool_detail: compact
  keybindings:
    show_shortcuts: f1
    toggle_transcript: ctrl+o
    show_activity: ctrl+t
    search_history: ctrl+r
    show_project_yaml: shift+tab
```

Safety keys (`Ctrl+C`, `Ctrl+D`, approval keys) cannot be reassigned. Invalid,
duplicate, or reserved bindings fall back to their defaults and are reported in
the transcript.

## Tool Activity

Each provider tool call owns one transcript cell keyed by its provider call ID.
The cell is updated in place from pending through approval to result, so the
transcript does not repeat the same call three times. Successful read-only calls
stay compact; warnings, rejects, and errors expand automatically. Press `Enter`
or `Space` on a tool or public decision-trace cell to inspect arguments, side
effects, semantic evidence, and failed rule IDs.

While a response is being prepared, tool calls, deterministic parser output,
semantic and intent gates, repair attempts, and intermediate commands remain
visible. After the turn finishes, those cells collapse into one **Tool chain**
row. Press `Enter`, `Space`, or click that row to inspect the complete public
evidence again; activate the same Tool chain row once more to hide it. The last
valid command remains the visible final deliverable;
for non-command turns the last substantive answer remains visible. A reject or
blocked summary replaces any stale command as the final deliverable.

Click a rendered answer or command explanation to open a plain-text selection
view. Drag with the mouse to copy a range, press `A` to select all, `C` to copy
the selection, and `Esc` to return to the conversation.

The footer reports the actual runtime phase, active operation, project/YAML and
server state, permission mode, provider/model, jobs, queued prompt, and measured
usage. It displays `YAML OK` or `YAML MISSING` in text as well as color. Token
usage is omitted until the provider returns real measurements; draft size is
shown separately as an estimate.

## Rendering Model Outputs

For command-generation turns, the TUI must show the generated CHEMSMART CLI
command before generated input contents. This keeps the agent grounded to the
same command shape a chemist can run manually:

```text
chemsmart run gaussian -p test -f examples/h2o.xyz -c 0 -m 1 opt
```

The transcript can then render:

- runtime semantic gate verdict and generated input evidence;
- dry-run `.com` / `.inp` contents;
- deterministic command interpretation;
- public decision trace;
- session completion and blocked/warn/reject summaries.

## Command Interpretation

The deterministic command parser expands a generated command into user-facing
facts:

- workspace;
- run vs submit execution;
- program and job type;
- server and dry-run status;
- input file and runtime-derived label;
- charge/multiplicity;
- project and the meaning of `-p`;
- functional, basis, auxiliary basis, and solvent when present;
- route parameters and job-specific options;
- one final English summary sentence.

This parser is deterministic and independent of the provider model.

## Public Decision Trace

API/frontier-provider turns may include a collapsible decision trace. This is
not hidden chain-of-thought. It records observable routing evidence:

- selected action: synthesize, explain, critique, repair, or clarify;
- confidence;
- target command;
- default project;
- short public reasoning/caveats when provided;
- rejected action classes.

The goal is auditability without relying on uninspectable model reasoning.

## Project YAML Build Mode

`/init` switches into project-YAML build mode. A user can paste a reported
computational-method paragraph and ask for a project file such as `co2.yaml`.
The tool path is:

```text
extract_project_protocol
  -> render_project_yaml
  -> validate_project_yaml
  -> critic_project_yaml
  -> write_project_yaml  # after explicit approval
```

Project YAML creation does not require a structure file. Valid output uses
CHEMSMART project-settings shape with top-level `gas:` and/or `solv:` blocks.

The candidate and the write are intentionally separate. After a validated
candidate is rendered, run `/write-project` (or `/write-project co2`) to open
the approval dialog. If the target already exists, choose `Y` to overwrite the
current settings or `N` to preserve them and save a new numbered project. The
write is blocked if the candidate's matching validation receipt is absent,
rejected, or belongs to another provider call.

## Workspace Project Selection

Agent project YAML is resolved from the launch directory, not from a global
project folder:

```text
<workspace>/.chemsmart/gaussian/<project>.yaml
<workspace>/.chemsmart/orca/<project>.yaml
```

One candidate is loaded automatically. With two or more candidates, the footer
shows `YAML MISSING` until one is selected. Press `Shift+Tab`, move with
`Tab`/`Down` or `Shift+Tab`/`Up`, then press `Enter`. The selected path and hash
are bound to a validated command; modifying the file or changing the working
directory invalidates `/run` and `/submit` until the command is regenerated.
The footer always shows `YAML OK <program>:<project>` or `YAML MISSING`.

The provider definition remains global at
`~/.chemsmart/agent/agent.yaml`. Its optional `project:` value can select a
same-named workspace candidate but cannot substitute for a missing local file.
