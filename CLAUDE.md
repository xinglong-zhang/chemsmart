# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**Read `AGENTS.md` first** — it is the authoritative operating contract (hard
rules, approval boundaries, PySCF-specific constraints). This file covers commands
and architecture only, and does not repeat it.

Version-pinned ground truth lives in `.agents/skills/` and is loadable as skills:
`chemsmart-v3-1-4` (CLI, project YAML, registry), `pyscf-v2-14-0`,
`gpu4pyscf-v1-8-0`. Read the relevant skill before changing CLI, settings, job, or
runner code.

## Commands

```bash
make env            # create/update the `chemsmart` conda env from environment.yml
make install-dev    # editable install with [dev,test] extras
make configure      # interactive setup of ~/.chemsmart (server + program folders)
make fmt            # isort, then black -l 79
make lint           # ruff check . --fix
make test           # lint + coverage-clean + pytest with coverage over tests/
make docs           # sphinx HTML
```

The Makefile wraps every command in `conda run -n chemsmart` **unless**
`CONDA_DEFAULT_ENV` is already `chemsmart` — it skips the wrapper then because
some conda versions swallow the wrapped exit code, which once let a failing suite
report as a passing CI step.

Single test:

```bash
pytest tests/test_orca_cli.py::TestORCASolventCLISpCommand::test_solvent_model_and_id_injected_into_sp_settings_group_level
```

Line length is **79** (black). Python is `~=3.10`. `numpy` is pinned to
`~=1.26.4` because rdkit and pymol require numpy 1.x — this constrains what any
new backend can share an environment with.

## Architecture

Four layers. The seams between them are what require reading multiple files.

### 1. `chemsmart/cli/` — Click groups

`chemsmart run <program> ...` executes locally; `chemsmart sub <program> ...`
submits to a scheduler. Both attach the same list from `cli/subcommands.py`, so a
program registered there exists in both and one absent exists in neither.

A program group (`cli/orca/orca.py` is the reference) resolves molecules and
settings, then publishes `project_settings`, `job_settings`, `keywords`,
`molecules`, `molecule_indices`, `label`, `filename` on `ctx.obj`. Subcommands
read that, merge their per-jobtype settings, and return a `Job` or list of jobs;
the group's `result_callback` hands it to `run`/`sub`.

`sub` rebuilds argv with `CtxObjArguments.reconstruct_command_line()` and strips
sub-only options (`time_hours`, `queue`, `verbose`, `test`, `print_command`). An
option that cannot survive that round trip makes a submitted job diverge from the
same job run locally.

**The `keywords` tuple is a merge whitelist**, and forgetting it is the most
common bug here — see `AGENTS.md`.

### 2. `chemsmart/settings/` — project and server configuration

`<Program>ProjectSettings.from_project(name)` resolves
`~/.chemsmart/<program>/<name>.yaml`, then the built-in test projects, then
raises. The YAML's `gas:` section feeds most jobtypes and `solv:` feeds `sp`; if
`gas:` is absent, `solv:` feeds everything. `sp` defaults to `freq: False`, but an
explicit `freq:` in `solv:` overrides that default. Details in the
`chemsmart-v3-1-4` skill.

`settings/server.py` + `submitters.py` handle PBS/SLURM/LSF; `executable.py` maps
a `PROGRAM` name to a block in the server YAML for `EXEFOLDER`, `CONDA_ENV`,
`MODULES`, `ENVARS`, `SCRATCH`.

`settings/capabilities.py` is the immutable, agent-independent declaration of
executable program names, immediate Click jobtype children, project requirements,
project-owned CLI parameters, and CPU/GPU execution engines. Its projections have
deliberately narrower meanings than similarly named agent-harness lists. Read the
`chemsmart-v3-1-4` skill's
`references/agent-harness-forward-compatibility.md` before consuming them.

Note `settings/gaussian.py` and `settings/orca.py` are structurally identical
copies. Do not add a third — see the `adding-a-program` reference.

### 3. `chemsmart/jobs/` — jobs, runners, writers

`JobRunner.from_job()` picks a runner purely by `job.TYPE in runner.JOBTYPES`,
tie-broken on `FAKE`. Registration is `RegistryMixin`, so **the module must be
imported** to be in the registry. `Job._job_is_complete()` reads exactly one
thing: `self._output().normal_termination`.

Scratch resolves: explicit `--scratch/--no-scratch` → the program's `SCRATCH` key
in server YAML → the runner class default.

### 4. `chemsmart/io/` — parsers

Readers for Gaussian, ORCA, xTB, PySCF, PDB, XYZ, and the `.db` store. The
external-program readers are large (`io/orca/output.py` is ~4k lines) because we
do not control those programs' output text. `io/pyscf/output.py` instead reads the
structured sibling `.h5` artifact; it never parses the human-facing `.out` log.

## Program status

First-class `run`/`sub` programs are **Gaussian, ORCA, xTB, NCIPLOT, and
PySCF**. Gaussian, ORCA, xTB, and NCIPLOT are external binaries; PySCF is a
library backend whose executable is a Python interpreter. CREST remains an
input and workflow-preparation source rather than a first-class program.
GPU4PySCF is an execution engine of the single `pyscf` program, not a second
program.

The xTB execution surface is deliberately narrower than the xTB reader:
`sp|opt|hess`, `gfn0|gfn1|gfn2|gfnff`, CPU execution, and strict optional
project YAML. Unknown project keys, contradictory job types, incomplete
solvent pairs, stale mixed outputs, arbitrary xcontrol, dynamics, path
following, and unsupported constraints fail closed.

The working tree implements and Stage A runtime-validated this PySCF contract:

- The CLI surface is `chemsmart run|sub pyscf {sp,opt,hess}` and requires a
  PySCF project YAML; there is no `freq` leaf.
- `--gpu/--no-gpu` wins over the default derived from server YAML `NUM_GPUS`.
  The resolved `cpu`/`gpu` engine is part of the settings digest and provenance.
- ChemSmart emits standalone `label.py`, runs it in a subprocess, writes PySCF's
  human log to `label.out` and the sole machine contract to `label.h5`, and
  captures subprocess stderr in `label.err`.
- HDF5 schema 2.0 exposes real `spec/`, `provenance/`, `status/`, and `results/`
  groups. The reader retains strict compatibility with schema 1.0 artifacts.
- An `opt` whose project settings say `freq: True` composes
  `scf -> opt -> hess` in one process. The explicit `hess` leaf analyzes the
  supplied geometry without optimizing it.
- Settings validation rejects known double-hybrid names and families, DFT grids on HF,
  auxiliary bases when density fitting is disabled, and every solvent model
  without an explicit solvent ID. PCM dielectric lookup occurs only inside the
  final child environment.
- `preflight()` and `verify_provenance()` return typed violations without
  executing chemistry. They are callable APIs, not a current CLI gate: the CLI
  does not yet bind injected environment evidence mandatorily to the finalized
  compute interpreter and execution phase. Full agent orchestration is still a
  future merge because this v3.1.4 checkout intentionally has no
  `chemsmart.agent` package.
- The current `settings_digest` identifies the controller-side requested
  configuration. Target-resolved solvent values are recorded in `spec/`, but
  there is not yet a separately recomputed target-side
  `applied_settings_digest`; the HDF5 spec is therefore a resolved receipt,
  not runtime-object attestation.
- Free energies come from `chemsmart/analysis/thermochemistry.py`, not from
  `pyscf.hessian.thermo.thermo()`.

The constraints behind each of these are in `AGENTS.md` and the pyscf skills.
