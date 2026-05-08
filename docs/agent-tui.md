# chemsmart agent TUI (Phase 1)

Phase 1 adds an interactive Textual interface to `chemsmart agent` while
keeping the existing one-shot agent subcommands unchanged.

## Install

```bash
pip install -e ".[agent-tui]"
```

## Launch

```bash
chemsmart agent
chemsmart agent --plain
```

`--plain` keeps the TUI inline, disables mouse support, and turns off
animations for conservative terminal environments such as HPC login nodes.

## Phase 1 commands

- `/help`
- `/quit`
- `/clear`
- `/sessions`
- `/resume <session-id>`
- `/tools`
- `/doctor`

## Live transcript cells

Phase 1 renders these decision-log artifacts as typed transcript cells:

- user request
- plan
- dry-run input
- critic verdict
- error

The TUI reads the same `decision_log.jsonl` files written by
`AgentSession.run(...)`, so `chemsmart agent run ...` remains unchanged.
