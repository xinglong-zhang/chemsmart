# Working with CHEMSMART

Read `AGENTS.md` first. It defines the scientific mission, execution boundary,
behavior-first research method, and repository hygiene rules for every coding
agent.

Use the project-local skills under `.agents/skills` when a task matches them.
The checked-out source, project loaders, and live Click schema remain the
operational authority.

## Common commands

```bash
make env            # create/update the Python 3.10 Conda environment
conda activate chemsmart
make install-dev    # editable install with development and documentation deps
make configure      # create local ~/.chemsmart templates
make test           # lint and complete pytest suite
make docs           # Sphinx HTML build
```

For a focused check, invoke `pytest` directly in the active environment. Do
not run broad suites after every edit, and do not use test success as a
substitute for a real computational-chemistry observation.

## Architecture

- `chemsmart/cli`: Click command registry shared by local `run`, scheduler
  `sub`, and the agent command compiler.
- `chemsmart/settings`: project YAML, server configuration, program
  capabilities, and resource selection.
- `chemsmart/jobs`: program settings, generated input, jobs, runners, scratch
  behavior, and execution semantics.
- `chemsmart/io`: program outputs and molecular data readers.
- `chemsmart/analysis`: cross-program quantities, dimensions,
  thermochemistry, and post-processing.
- `chemsmart/agent`: provider-neutral typed planning, project promotion,
  preview, approved execution, result analysis, and event replay.

`chemsmart run` and `chemsmart sub` attach the same program command families.
Any setting added to one route must survive project loading, Click option
resolution, run/sub reconstruction, writer materialisation, parser
observation, and result analysis where applicable.

## Program boundary

Gaussian, ORCA, xTB, NCIPLOT, and PySCF are program identities. GPU4PySCF is a
PySCF engine. External binaries and licensed programs are never bundled.

PySCF uses a dedicated compute interpreter and writes a structured HDF5 result.
xTB is limited to the maintained CPU `sp`, `opt`, and `hess` surface. ORCA
parallel work requires a compatible ORCA/MPI installation; do not treat a
zero process status as normal ORCA termination.

Models express intent through project YAML and typed tools. Do not add native
input or shell-text fallbacks to repair an agent failure.

## Research discipline

When a model or calculation fails, identify whether the cause is chemistry,
method capability, environment, parser, CLI affordance, or model reasoning.
Change the smallest general layer and rerun a relevant real observation. Keep
benchmark-specific answers and one-off analysis reports out of the repository.
