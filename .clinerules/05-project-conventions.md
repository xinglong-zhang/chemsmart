---
paths:
  - "chemsmart/**"
  - "tests/**"
  - "scripts/**"
---

# Project Conventions (active when touching code)

## Identity & scope

chemsmart: open-source computational chemistry planning & HPC automation
toolkit (Python 3.10+). Gaussian/ORCA job generation, submission, analysis,
plus the agent/TUI layer (`chemsmart/agent/`). Local model path (MLX 4-bit,
Apple Silicon) and cloud providers (OpenAI-compatible/Anthropic) both exist.

## Python style

- Python 3.10 target, type hints on params/returns, PEP 8, line length 79.
- f-strings; context managers for resources; docstrings on public API;
  composition over inheritance.
- Format `black` (79) + `isort` (profile=black); lint `ruff`/`flake8`.
- Never print install commands (`pip/brew/apt install`) in code or docs.

## Testing

- `pytest --strict-markers --disable-warnings`; agent tests in
  `tests/agent/`; markers: `slow`, `agent`, `integration`.
- Every behavior change ships with a test. Run the touched test files
  first, then the full `tests/agent/` suite before declaring done.
- Legacy LLM execution tools (`run_local`, `submit_hpc`,
  `extract_optimized_geometry`) are deprecated JobArgs stubs — their old
  contract tests are intentionally skipped; do not "fix" them back.

## Architecture map (agent layer)

- CLI entry: `chemsmart.cli.main:main` (Click).
- Unified agent loop: `chemsmart/agent/core.py` (AgentSession.run_loop) +
  `loop.py` (ToolLoop) + `registry.py` (ToolRegistry, TOOL_GROUPS).
- Command synthesis: `chemsmart/agent/synthesis.py` (SynthesisSession;
  captures `_last_reasoning`); schema pruning `schema_prune.py`.
- Command tools: `tools_command.py` (synthesize/repair/execute + gate cache).
- Training capture: `training_log.py` writes append-only turn snapshots to
  `var/agent-training/` and model-specific `runs/<model>/` stores. The
  exporter reconstructs positive/review session chains and separate repair
  contrasts; the auditor reports terminal-gate, multi-turn-session,
  diversity, and canonical-kind coverage metrics. Implementations are
  `scripts/training/export_sft.py` and
  `scripts/training/audit_dataset.py`.
- System prompt assembly: `prompts/identity.py:build_system_prompt`.
- Workspace project YAML: `./.chemsmart/<program>/<name>.yaml`
  (`chemsmart/settings/workspace_project.py`).

## Key dependencies (pinned facts)

`numpy~=1.26.4` (numpy-2 ABI is a known trap), `pydantic>=2`, `click`,
`transformers==4.56.2` (local path), `textual` (TUI). Dev: `black`,
`isort`, `ruff`, `pre-commit`.
