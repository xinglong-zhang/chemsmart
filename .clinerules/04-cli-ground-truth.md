# chemsmart CLI Ground Truth (for judging synthesized commands)

Source-verified facts. When in doubt, trust the semantic gate's verdict and
the real click parser over any model output (including your own reasoning).

## Command shape

```
chemsmart {run|sub} [group opts] {gaussian|orca} [-p PROJECT] -f FILE \
    -c CHARGE -m MULT <jobkind> [job-specific structural opts]
```
- `run` = local execution; `sub` = HPC submission. Server/cores/mem/time
  flags belong on the run/sub GROUP (before the engine), never after it.
- Test-mode injection (only when actually executing): `run` → add
  `--fake --no-scratch`; `sub` → add `--fake --test`.

## Job kinds (canonical, exhaustive)

- gaussian: `com crest dias irc link modred nci opt qrc resp scan sp td
  traj ts userjob wbi` (TDDFT's CLI name is `td`, not `tddft`)
- orca: `inp irc modred neb opt qmmm qrc scan sp ts`
- `qmmm` is a NESTED child under each gaussian jobkind (`… opt qmmm`).
- There is NO `freq` subcommand. Frequency on/off is project-YAML-owned.

## Runtime-owned fields — must NOT appear as CLI flags in trusted rows

functional, basis, ab initio method, aux basis, dispersion (D3/D3BJ),
solvent model/id, freq. These come from the workspace project YAML
(`./.chemsmart/<program>/<project>.yaml`). A command carrying
`--functional`, `--basis`, `-x`, `-b`, or freq smuggled into route params
(`-r '... freq ...'`, `--freq`, `freq=`) is a CANONICAL VIOLATION — the
exporter auto-skips it (`canonical_*` skip reasons); never hand-approve one.

## Structural options the model MAY emit

`-f/--filename`, `-c/--charge`, `-m/--multiplicity`, `-p/--project`,
`-l/--label`, and job-specific structure: scan/modred coordinate lists
(comma/range strings like `"1,2"` — never space-separated), `-s` step size,
`-n` step count, td `-n` nstates, dias fragment-1 indices (flat list),
neb end-file. `additional_route_parameters` must be a STRING when present.

## Project (`-p`) rule

If a workspace project YAML is loaded, the runtime injects the default
project when the model omits `-p`. A command referencing a project that
does not exist in the workspace will be gate-rejected — the fix is the
workspace, not inventing a name.

## Known model failure patterns (watch for these in WRONG rows)

- scan↔modred confusion (both take atom pairs; scan asks drive the bond,
  modred asks freeze it) — the single most common systematic WRONG.
- Router confusion: a plain job request answered by project-YAML authoring
  tools instead of `synthesize_command` (fix: command-hard framing).
- Invented option order: options after the jobkind that belong before it.
