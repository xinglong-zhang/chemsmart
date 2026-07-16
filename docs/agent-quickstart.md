# CHEMSMART agent quickstart

This guide starts one workspace-local computational-chemistry session. The
agent always produces a real `chemsmart run ...` or `chemsmart sub ...` command
before execution, and deterministic gates remain the source of truth.

## Install And Check The Provider

```bash
conda activate chemsmart
pip install -e ".[agent-tui]"
chemsmart agent doctor
```

The active provider is configured in `~/.chemsmart/agent/agent.yaml`. API and
OpenAI-compatible providers can use the full tool loop. A local provider is a
command-synthesis specialist; project authoring with `/init` requires a
tool-calling API provider.

## Start From A Research Workspace

Project YAML belongs to the directory where CHEMSMART is launched. Do not start
the TUI in the source repository when the calculation belongs elsewhere.

```bash
mkdir -p /private/tmp/chemsmart-water-study/inputs
cp examples/h2o.xyz /private/tmp/chemsmart-water-study/inputs/h2o.xyz
cd /private/tmp/chemsmart-water-study
chemsmart agent
```

The agent searches only these workspace paths for project settings:

```text
.chemsmart/gaussian/<project>.yaml
.chemsmart/orca/<project>.yaml
```

The footer displays `YAML OK <program>:<project>` when one is active and
`YAML MISSING` otherwise.

## Build A Project From A Method

Enter `/init`, then send a concise method description such as:

```text
Create an ORCA project named water_sp. Use gas-phase PBE0/def2-SVP for
single-point calculations and do not request frequencies.
```

The API tool loop performs:

```text
extract_project_protocol
  -> render_project_yaml
  -> validate_project_yaml
  -> critic_project_yaml
```

This creates a validated candidate but does not write a file. Run:

```text
/write-project
```

The TUI asks for explicit write approval. If a file already exists, choose
`Y` to overwrite it or `N` to keep it and create a new numbered project. In
plain mode the equivalent explicit form is:

```text
/write-project water_sp yes overwrite
```

Press `Shift+Tab` to inspect the YAML. If several projects exist, use
`Tab`/`Down` and `Shift+Tab`/`Up` to move, then `Enter` to select one.

## Prepare A Calculation

Send a natural-language request:

```text
Prepare an ORCA single-point calculation for @inputs/h2o.xyz as a neutral
singlet using the active water_sp project.
```

Typing `@` opens the workspace file picker. The completed turn should end with
a visible command similar to:

```text
chemsmart run orca -p water_sp -f inputs/h2o.xyz -c 0 -m 1 sp
```

Before that command becomes executable, CHEMSMART checks:

1. CLI tokenization and real option placement;
2. workspace project resolution;
3. safe fake/dry-run execution;
4. generated ORCA or Gaussian input;
5. program- and job-specific invariants;
6. preservation of the requested program, kind, path, project, charge, and
   multiplicity.

The deterministic parser and gate evidence remain visible while the agent is
working. Once the turn finishes, they collapse into one **Tool chain** row.
Press `Enter`, `Space`, or click it to inspect the evidence again. The final
command remains visible.

## Execute And Monitor

Run the last validated local command with:

```text
/run
```

Real execution always requires this explicit action even in permissive modes.
The status strip appears as soon as the process starts. Press `Ctrl+B` or use
`/runs` to open the calculation monitor:

| Key | Action |
|---|---|
| `Up` / `Down` | Select a calculation or scheduler job |
| `Enter` | Open its structured receipt |
| `L` | Follow the bounded raw log |
| `/` | Search the log |
| `E` | Extract the latest chemistry result |
| `C` | Request cancellation |
| `Esc` | Return to chat |

The conversation remains usable while a local calculation runs. A completed
single-point receipt reports final energy, SCF cycles, normal termination,
method/basis, elapsed time, and output path. Raw output is kept in the monitor,
not dumped into the main transcript.

After completion, ask:

```text
Diagnose the latest calculation result.
```

When there is one unambiguous recent run, the deterministic diagnostics path
uses its stored output path without asking for the project again.

## Keyboard Reference

| Shortcut | Action |
|---|---|
| `F1` | Show all shortcuts |
| `/` | Open and filter the command palette |
| `@` | Select a workspace file |
| `Shift+Tab` | Inspect or select workspace project YAML |
| `Ctrl+B` | Open calculation and scheduler monitor |
| `Ctrl+T` | Inspect current tool activity |
| `Ctrl+O` | Toggle compact/full transcript detail |
| `Ctrl+R` | Recall the previous request |
| `Ctrl+G` | Edit the draft in `$EDITOR` |
| `Tab` while busy | Queue or restore one follow-up request |
| `Esc` | Close an overlay and return to the composer |

Click a rendered answer to open a selection view. Drag to copy a range, press
`A` to select all, and `C` to copy the current selection.

## Sessions And Audit Evidence

Every session is stored below `~/.chemsmart/agent/sessions/<id>/`. The record
may include:

- `decision_log.jsonl` for requests, tool events, approvals, and verdicts;
- `session_metadata.json` for runtime status and harness summary;
- `harness_result.json` for generated-input invariant evidence;
- `runtime_events.jsonl` and calculation receipts for replay/recovery;
- generated `.com` or `.inp` inputs from safe validation.

Use `/sessions` and `/resume <session-id>` in the TUI, or the matching
`chemsmart agent sessions` and `chemsmart agent resume` CLI commands.

## Common Failures

- `YAML MISSING`: build one with `/init` and `/write-project`, or select among
  existing workspace candidates with `Shift+Tab`.
- `/run` blocked after editing YAML: regenerate the command so its stored YAML
  hash matches the current file.
- Project name requested despite an existing YAML: select it with `Shift+Tab`
  and verify the footer reads `YAML OK` before sending the calculation request.
- Missing input file in `@`: confirm the TUI was launched from the intended
  workspace and the file is below that directory.
- Local provider rejects `/init`: switch to a configured tool-calling API
  provider for project authoring, then return to local synthesis if desired.
