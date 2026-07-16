# chemsmart agent

`chemsmart.agent` is the natural-language interface for CHEMSMART. It connects
provider models, deterministic command synthesis, runtime harness checks, and a
Textual TUI. The design goal is not free-form shell automation; it is
auditable, CLI-grounded computational-chemistry workflow setup.

## User Interface

Launch the interactive UI with:

```bash
chemsmart agent
```

Use `chemsmart agent --plain` for terminals where mouse support or animations
are undesirable.

![chemsmart agent dry-run TUI](../../tests/agent/tui/snapshots/dry_run_input.svg)

Typing `/` opens the command palette. Continue typing, such as `/d`, to filter
commands.

![chemsmart slash command palette](../../tests/agent/tui/snapshots/slash_help.svg)

Common TUI commands:

| Command | Purpose |
|---|---|
| `/mode ask` | Command synthesis, command explanation, critique, and repair |
| `/mode run` | Full tool-loop harness mode |
| `/init` | Build a project YAML from a reported method |
| `/tools` | List registered tools |
| `/doctor` | Check provider/runtime health |
| `/sessions` | Browse resumable sessions |
| `/resume <session-id>` | Resume a saved session |
| `/run`, `/submit` | Continue a prepared session toward local run or submission |

## Architecture

```text
user request
  -> provider/router
  -> runtime task envelope + phase-scoped tool catalog
  -> compact SPEC or CLI-grounded command decision
  -> deterministic postprocessor / adapter
  -> real chemsmart run|sub command
  -> real parser + semantic fake/dry-run gate
  -> TUI evidence cells + session artifacts
```

The important separation is:

- model/API providers infer user intent and produce a compact plan or a command
  decision;
- deterministic code owns project defaults, command rendering, parser checks,
  generated-input inspection, and final warnings/rejections;
- risky tools stay behind permission gates and explicit approval paths.

For job-creation requests, the primary user-facing artifact must be a real
CHEMSMART command, not an internal `build_molecule`/`build_job` tool path.
Internal tool calls remain useful for validation and evidence.

### Runtime v2

The provider-independent runtime kernel lives in `runtime/`. It adds a typed
task envelope, phase-scoped tool exposure, lifecycle hooks, artifact receipts,
a bounded repair policy, and an append-only hash-chained event log. Runtime
state is reconstructed from `runtime_events.jsonl`; `runtime_state.json` is an
atomic snapshot, not a second source of truth.

Migration is feature-flagged:

```bash
CHEMSMART_AGENT_RUNTIME_V2=shadow chemsmart agent  # observe, do not restrict
CHEMSMART_AGENT_RUNTIME_V2=active chemsmart agent  # enforce phase tool sets
```

`shadow` keeps the legacy tool surface and records exposure violations.
`active` exposes at most five real tools for the current phase plus the virtual
`ask_user` tool. The flag remains off by default until the frozen HighRisk48
regression matrix passes.

The autonomy boundary is fixed in runtime policy:

| Operation | Runtime policy |
|---|---|
| Read/validate project YAML | Automatic |
| Synthesize/repair and semantic-gate a command | Automatic |
| Execute with `test=True` or preview `submit_hpc(execute=False)` | Automatic safe path |
| Write/update project YAML | Explicit approval every time |
| Real local execution or HPC submission | Explicit approval every time |

`bypass` and legacy `yolo` do not override the final two rows.

New project authoring and project writes are separate runtime phases. Authoring
may extract, render, and validate a candidate automatically; it never writes a
workspace file. A successful authoring response includes a deterministic notice
to use `/write-project`, which enters the approval-gated write path. A rendered
candidate cannot complete unless its latest runtime-loader verdict is `ok` or
`warn`.

## Provider Modes

`~/.chemsmart/agent/agent.yaml` is the source of truth for the active provider.

Supported provider types:

- `openai` / `anthropic`: API-backed frontier providers. These may route between
  command synthesis, command explanation, critique, repair, and clarification.
- `local`: the local fine-tuned CHEMSMART model path. This uses the compact SPEC
  adapter and can run with PyTorch or Apple Silicon MLX, depending on provider
  configuration. Runtime v2 treats it as a synthesis specialist and exposes
  only `synthesize_command` and `repair_command`; orchestration remains
  deterministic.

Provider identity is independent from the wire protocol. OpenAI-compatible
providers such as DeepSeek keep their own identity for routing and receipts but
declare `wire_protocol = "openai"` for tool-call parsing and tool-result messages.

If the provider config sets `project: test`, that project is attached
deterministically. The model should not invent project, functional, basis, or
solvent defaults when runtime configuration owns those fields.

## Grounding Components

| Component | Role |
|---|---|
| `synthesis.py` | One-shot command synthesis, API routing, repair, and command-answer orchestration |
| `command_answerer.py` | Lets a frontier provider phrase an answer from deterministic command facts, then rejects contradictions |
| `model_command_parser.py` | Parses generated `chemsmart run|sub ...` commands into user-friendly facts |
| `postprocess_v8.py` / `v8_adapter.py` | Normalize compact local-model SPEC and render real CLI commands |
| `harness/command_semantics.py` | Safe command execution and generated-input evidence |
| `harness/runner.py` | Runtime harness result aggregation and session artifacts |
| `harness/invariants/gaussian_ts.py` | Gaussian TS route invariants, including duplicate/leaked TS-token rejection |
| `harness/basis_sets/` | BSE-backed basis-set catalog, resolver, and top-k search tool |
| `project_yaml.py` | Extract, render, validate, critique, and write CHEMSMART project YAML |
| `runtime/` | Task contracts, event replay, tool exposure, lifecycle, receipts, and repair policy |
| `tui/` | Textual UI, slash palette, transcript cells, footer/header state, session workers |

## TUI Evidence Cells

Recent TUI updates add:

- generated CHEMSMART command display before dry-run input content;
- deterministic command interpretation with workspace, program, job, server,
  dry-run status, project, method/basis, solvent, route options, and summary;
- collapsible public decision trace for API-routed turns;
- active provider/model/project display in the footer;
- slash-command palette with prefix filtering;
- project-YAML build mode through `/init`.

The public decision trace is intentionally not hidden chain-of-thought. It shows
observable routing evidence, the selected action, confidence, rejected action
classes, and caveats so the user can audit why the agent synthesized,
explained, critiqued, repaired, or asked for clarification.

## Basis-Set Search

`search_basis_sets(query, program, limit=8, role="any")` is a read-only tool
backed by a generated Basis Set Exchange catalog. It avoids prompt/catalog
bloat by returning a compact top-k result instead of all basis names.

Examples it handles:

- `"Karlsruhe triple-zeta diffuse"` -> `def2-TZVPD` among top candidates
- `"RI fit for def2-TZVP"` -> `def2-TZVP-RIFIT`
- `"six thirty one star"` -> `6-31G*`

See `harness/basis_sets/DESIGN.md` for the data contract and ranking policy.

## Project YAML Tools

The project-YAML path is:

```text
extract_project_protocol
  -> render_project_yaml
  -> validate_project_yaml
  -> critic_project_yaml
  -> write_project_yaml  # only after approval
```

Generated YAML must use CHEMSMART project-settings shape with top-level `gas:`
and/or `solv:` blocks. It must not invent wrappers such as `gaussian:`,
`project_name:`, or `method:` as top-level YAML.

## Validation

Focused verification used for this README update:

```bash
env HOME=/private/tmp/chemsmart-doc-verify \
  /opt/anaconda3/envs/chemsmart/bin/python -m pytest \
  tests/agent/harness/test_basis_catalog.py \
  tests/agent/test_command_answerer.py \
  tests/agent/test_model_command_parser.py \
  tests/agent/test_project_yaml_tools.py \
  tests/agent/test_synthesis.py \
  tests/agent/tui/test_synthesis_mode.py \
  tests/agent/tui/test_slash_commands.py \
  tests/agent/tui/test_track_b_ux.py -q
```

Result: `74 passed`.

Also verified:

```bash
git diff --check
```

after the README edits.
