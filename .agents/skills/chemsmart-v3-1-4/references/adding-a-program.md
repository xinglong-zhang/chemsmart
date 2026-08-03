# Adding a program to ChemSmart

The seam list, derived by tracing ORCA end to end. **ORCA's shape assumes an
external binary whose text output we must reverse-engineer.** For a Python-library
backend (PySCF), some seams are right to copy and some are actively wrong. Each
row is marked.

## Seams to copy

| Seam | ORCA anchor | Dispatch key |
|---|---|---|
| CLI group | `chemsmart/cli/orca/orca.py` | option decorators + `ctx.obj` |
| CLI subcommands | `chemsmart/cli/orca/{opt,singlepoint}.py` | `@<prog>.group(...)` |
| CLI registry | `chemsmart/cli/subcommands.py` | import + append to `subcommands` |
| Program capabilities | `chemsmart/settings/capabilities.py` | one frozen `ProgramCapability` entry |
| Job settings | `chemsmart/jobs/orca/settings.py` | subclass `MolecularJobSettings` |
| Job classes | `chemsmart/jobs/orca/{job,opt,singlepoint}.py` | `TYPE` string |
| Runner | `chemsmart/jobs/orca/runner.py` | `JOBTYPES`, `PROGRAM`, `SCRATCH`, `FAKE` |
| Server YAML block | `settings/templates/server.yaml` | uppercase program name |
| `config` subcommand | `chemsmart/cli/config.py` | per-program folder setup |
| User settings props | `chemsmart/settings/user.py` | `user_<prog>_settings_dir` |
| Project template | `settings/templates/template_orca_simple.yaml` | — |
| Tests | `tests/test_orca_cli.py`, `conftest.py` fixtures | `CliRunner` + mocked job ctor |

Registration mechanics:

- `JobRunner.from_job()` matches `job.TYPE in runner.JOBTYPES`, then prefers the
  runner whose `FAKE` equals the requested `fake`. Registration is `RegistryMixin`
  — **the module must be imported** or the subclass is not in the registry.
- Every `Job` needs a unique `TYPE`; collisions resolve arbitrarily.
- `Job._job_is_complete()` reads exactly one thing: `self._output().normal_termination`.

## Seams to reject or restructure for a library backend

### Do not clone the output parser

`chemsmart/io/orca/output.py` is **3,966 lines** of regex parsing (verified by
`wc -l`). It exists because ORCA's text output is the only interface we get.

If ChemSmart generates the calculation script itself, it controls the output
format too. Emit a structured results file and read it back with a small typed
reader. Current PySCF uses `label.h5` as its sole machine contract, with real
`spec/`, `provenance/`, `status/`, and `results/` groups; `label.out` remains a
human log. Completion comes from per-stage structured convergence and typed
failure state, not process return or log text. Preserve array dtype/shape and
the documented scientific units. A regex parser for output we ourselves
produced is pure maintenance cost, and it loses precision (printed decimals
instead of float64) for nothing.

### Do not create a third copy of the project-settings machinery

`settings/gaussian.py` (532 lines) and `settings/orca.py` (689 lines) are
structurally identical: the same four classes (`*ProjectSettings`,
`Yaml*ProjectSettings`, `Yaml*ProjectSettingsBuilder`, `*ProjectSettingsManager`),
same method names, differing only in the settings class instantiated and the
jobtype list. A third copy is ~500 duplicated lines.

Prefer a generic `YamlProjectSettings` base parameterised by
`(settings_class, jobtypes)`. Introduce it for the *new* program only; migrating
Gaussian and ORCA onto it is a separate, behaviour-risky change.

Related: the historical default in `read_molecular_job_yaml()` contains
Gaussian/ORCA kinds (`dias`, `resp`, `wbi`, `traj`, `uvvis`, `crest`). Keep that
default for those existing callers, but pass `gas_phase_jobs` and `sp_jobs` for a
narrower backend. PySCF passes `("opt", "hess")` and `("sp",)` and therefore
produces exactly its three settings configurations instead of fabricating
irrelevant ones.

### Keep `EXEFOLDER` optional for a library backend

`Executable.from_servername()` now reads `EXEFOLDER` with `.get()`. Do not
reintroduce a bare subscript: a Python-library backend has no binary folder and a
missing program block may validly mean "use the current interpreter."

`PySCFExecutable.from_servername()` consumes a `PYSCF:` block when present and
falls back to `sys.executable` when it is absent. Within a present block,
`EXEFOLDER` is optional. If configured, it names the environment's `bin/` and
`get_executable()` resolves `<EXEFOLDER>/python`; if omitted, the running
interpreter is used. `CONDA_ENV`, `MODULES`, `ENVARS`, and `SCRATCH` remain useful
independently of `EXEFOLDER`.

## Server YAML

Each configured program gets an uppercase block alongside `SERVER:`. Keys
consumed by `Executable`: optional `EXEFOLDER`, `LOCAL_RUN`, `CONDA_ENV`,
`MODULES`, `SCRIPTS`, `ENVARS`. `SCRATCH` is read separately by
`Executable.program_scratch_from_servername()` for the scratch ladder.
`SERVER.NUM_GPUS` is what makes GPU-capable programs selectable.

The per-program `CONDA_ENV` exists so a program can run in **its own environment**.
Anything the generated job imports must therefore be installable in that env
independently of ChemSmart's own pins (notably `numpy~=1.26.4`, held down by
rdkit/pymol).
