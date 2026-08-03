---
name: chemsmart-v3-1-4
description: Ground truth for the ChemSmart v3.1.4 CLI, project-YAML contract, immutable program capabilities, structured PySCF backend, and job/runner registry. Use when adding or auditing a program group, writing or debugging project YAML, consuming capabilities.py from a future agent harness, tracing ignored CLI options or rejected YAML keys, or reasoning about run/sub parity, scratch, engines, and server configuration.
---

# ChemSmart v3.1.4

Ground truth for the CLI grammar and settings contract of ChemSmart **v3.1.4**
(`pyproject.toml: version = "3.1.4"`, upstream `main`). Read `AGENTS.md` first.

**This supersedes any v2.0.1-era description of ChemSmart.** This integrated
fork adds a deliberately bounded xTB execution surface to the v3.1.4 base.
See [v3-surface-delta.md](references/v3-surface-delta.md) before asserting that
any program or subcommand exists.

The integrated tree adds first-class `pyscf {sp,opt,hess}` and
`xtb {sp,opt,hess}` commands to v3.1.4. Treat GPU4PySCF as the `gpu` execution
engine of `pyscf`, never as a second program. PySCF emits `label.py`,
`label.out`, `label.h5`, and `label.err`; read only the structured HDF5 sibling
for program state and scientific results.

## The two entry points

Everything runs through `chemsmart run <program> ...` (execute here) or
`chemsmart sub <program> ...` (submit to a scheduler). Both attach the *same*
list from `chemsmart/cli/subcommands.py`, so a program registered there exists in
both, and one absent from it exists in neither.

`run` and `sub` must stay behaviourally identical. `sub` reconstructs the
command line with `CtxObjArguments.reconstruct_command_line()` and strips the
options that exist only on `sub` (`time_hours`, `queue`, `verbose`, `test`,
`print_command`). An option that cannot survive that round trip makes a submitted
job diverge from the same job run locally.

## State handoff: group → subcommand

The program group (e.g. `cli/orca/orca.py`) resolves molecules and settings, then
publishes them on `ctx.obj`: `project_settings`, `job_settings`, `keywords`,
`molecules`, `molecule_indices`, `label`, `filename`. Subcommands consume that,
pick their per-jobtype settings, merge, and return a `Job` (or a list of them).
The group's `result_callback` passes the job up to `run`/`sub`, which executes or
submits it.

## The `keywords` whitelist — the most common bug in this codebase

A group-level CLI override reaches the job **only** if its attribute name was
appended to the `keywords` tuple:

```python
if functional is not None:
    job_settings.functional = functional
    keywords += ("functional",)          # <-- without this line, silently lost
```

Subcommands then call
`project_settings.opt_settings().merge(job_settings, keywords=keywords)`.
Setting the attribute without extending `keywords` is a **silent no-op**: no
error, no warning, and the project-YAML value wins. When a CLI flag "does
nothing", check `keywords` before anything else.

## Settings precedence

Project YAML supplies method, basis, phase, and `freq`; the CLI overrides
per-invocation; the jobrunner owns resources. Concretely:

- **Project YAML** (`-p/--project`) → functional/basis/solvent/freq per jobtype,
  split into a `gas:` section (most jobtypes) and a `solv:` section (`sp`).
- **Input file** → charge/multiplicity/geometry when parsed from `.com/.inp/.log/.out`.
- **CLI flags** → override the above, *if whitelisted in `keywords`*.
- **Jobrunner** (`-n/--num-cores`, `-m/--mem-gb`, `-g/--num-gpus`, `--scratch`) →
  never comes from project YAML.

Read [project-yaml-contract.md](references/project-yaml-contract.md) before
writing or debugging a project YAML. Note the two lookalike mechanisms behave
oppositely: an unknown **project-YAML key raises `ValueError`**, while a
non-whitelisted **`keywords` entry is dropped silently**.

## Adding a program

Read [adding-a-program.md](references/adding-a-program.md). It lists every seam,
and — importantly — marks which seams to copy from ORCA and which to **reject**,
because ORCA's shape assumes an external binary that emits text we must parse.
For a Python-library backend (PySCF), copying that shape wholesale adds thousands
of lines of avoidable code.

Declare executable program facts once in `chemsmart/settings/capabilities.py`.
Before wiring an eventual `chemsmart.agent` merge, read
[agent-harness-forward-compatibility.md](references/agent-harness-forward-compatibility.md)
for exact one-line touch points and the boundaries between Click inventory,
agent validation, project ownership, static typing, and CPU/GPU engines.

Two mechanics govern registration:

- `JobRunner.from_job()` selects a runner purely by `job.TYPE in runner.JOBTYPES`,
  tie-broken on `runner.FAKE == fake`. Registration is via `RegistryMixin`, so
  **the module must be imported** for the subclass to be in the registry.
- Scratch resolves in order: explicit `--scratch/--no-scratch` → the program's
  `SCRATCH` key in the server YAML → the runner class's `SCRATCH` default.

## Boundaries

Use this skill for CLI/settings/registry structure. It does not cover scientific
adequacy of a calculation, approval to execute, or evidence standards for a
claim — use `chemsmart-scientific-workflow` and `chemsmart-evidence-audit` for
those. For PySCF and GPU4PySCF specifics, use `pyscf-v2-14-0` and
`gpu4pyscf-v1-8-0`.
