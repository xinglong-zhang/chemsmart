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

![chemsmart slash command palette](../tests/agent/tui/snapshots/slash_help.svg)

Current high-value commands:

| Command | Purpose |
|---|---|
| `/mode ask` | Command synthesis, explanation, critique, and repair |
| `/mode run` | Full tool-loop harness mode |
| `/init` | Build and validate a project YAML from a reported method |
| `/tools` | List registered tools |
| `/doctor` | Run provider/runtime diagnostics |
| `/sessions` | Browse saved sessions |
| `/resume <session-id>` | Resume a saved session |
| `/run`, `/submit` | Continue a prepared session toward execution/submission |

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

## Active Project Display

The footer shows provider mode/model information and the active project when
available from `~/.chemsmart/agent/agent.yaml`. The same project should be used
for generated commands unless the user explicitly selects another project.
