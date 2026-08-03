# Project YAML contract

Ground truth: `read_molecular_job_yaml()` in `chemsmart/jobs/settings.py` plus
the program-specific project-settings builder. A project YAML controls
**functional, basis, phase (solvation), and freq** for its declared jobtypes. It
never controls cores, memory, GPUs, or scratch — those are jobrunner options.

This `gas`/`solv` contract applies to Gaussian, ORCA, and PySCF. xTB uses the
separate strict three-section dialect described below.

## Two sections drive everything

```yaml
gas:
  functional: m062x
  basis: def2-svp
  freq: True
solv:
  functional: m062x
  basis: def2-qzvpd
  freq: False
  solvent_model: smd
  solvent_id: cyclohexane
```

- With the historical Gaussian/ORCA default, **`gas:`** populates `opt, modred,
  ts, irc, scan, nci, crest, dias, resp, set, traj, uvvis, wbi, neb`, and
  **`solv:`** populates `sp`.
- A narrower backend passes explicit jobtype lists. PySCF maps `gas:` to
  `opt,hess` and `solv:` to `sp`.
- **If `gas:` is absent, `solv:` populates every declared jobtype** — the whole
  project becomes implicitly solvated. Verified: with only a `solv:` block,
  `cfg["opt"]["solvent_model"]` is the solvated value. A supported mode, but easy
  to trigger by accident.
- Optional `td:` and `qmmm:` sections override those jobtypes. `qmmm:` is applied
  key-by-key (unrestricted), unlike the others.
- Calls that omit `gas_phase_jobs` and `sp_jobs` retain the historical
  Gaussian/ORCA inventory and produce 15 jobtype configs. A narrower backend
  must pass explicit lists. PySCF passes gas jobs `opt,hess` and solvated job
  `sp`, so it produces exactly three configs.

`sp` gets `freq = False` applied *before* the `solv:` merge, so it is a **default,
not a lock** (verified):

| `solv:` block | resulting `sp` freq |
|---|---|
| omits `freq` | `False` — the default applies |
| says `freq: True` | `True` — the YAML wins |

So "`sp` never runs frequencies" is not guaranteed by the generic loader. If
that matters, assert it on the merged settings. PySCF's `sp` job deliberately
stops after SCF even if an inherited flag is present; use the explicit `hess`
leaf for a fixed-geometry frequency analysis. The shipped templates set
`freq: False` in `solv:` explicitly, which is why the default is rarely observed.

## Unknown keys raise — they are not ignored

Merging uses `update_dict_with_existing_keys(target, source)`. The base layer is a
sibling `defaults.yaml` if present, otherwise
`<Program>JobSettings.default().__dict__`.

**A YAML key that does not already exist on the settings class raises
`ValueError` at load time** (verified), naming the offending key and listing every
valid one:

```
ValueError: Keyword `fucntional` is not in list of keywords
`dict_keys(['ab_initio', 'functional', 'dispersion', 'basis', ...])`
Please double check and rectify!
```

This is a good failure mode — a typo stops the run instead of silently applying a
default. Two consequences:

- The settings class defines the project-YAML schema. A setting cannot be used in
  a YAML before the attribute exists on `<Program>JobSettings.__init__`; adding
  the key first is an immediate hard error, not a quiet no-op.
- The valid-key list in the error is the authoritative surface for that program.
  Read it rather than guessing.

Contrast this with the `keywords` whitelist on `settings.merge()`, which **is**
silent — see the main SKILL.md. The two mechanisms look similar and behave
oppositely.

## xTB's strict optional project dialect

xTB does not require `-p`; omitting it selects ChemSmart's declared GFN2
defaults. If a project is supplied, it contains only `sp`, `opt`, and `hess`
sections. Each section accepts only the fields declared by `XTBJobSettings`;
an unknown section, unknown key, non-mapping section, or a `jobtype` that
disagrees with its section raises before job construction.

```yaml
sp:
  gfn_version: gfn2
  charge: 0
  multiplicity: 1
  solvent_model: null
  solvent_id: null
  grad: false
opt:
  gfn_version: gfn2
  optimization_level: vtight
  charge: 0
  multiplicity: 1
  solvent_model: null
  solvent_id: null
hess:
  gfn_version: gfn2
  charge: 0
  multiplicity: 1
  solvent_model: null
  solvent_id: null
```

The allowed method identifiers are `gfn0`, `gfn1`, `gfn2`, and `gfnff`.
Solvent model and solvent identifier must be supplied together or both omitted.
Arbitrary xcontrol text, dynamics, path-following workflows, and unsupported
constraints are not part of this dialect.

## Resolution order

`<Program>ProjectSettings.from_project(name)` tries, in order:

1. `~/.chemsmart/<program>/<name>.yaml` (user settings; the normal location)
2. `tests/data/<PROGRAM>Tests/project_yaml/<name>.yaml` (built-in test projects)
3. raises `FileNotFoundError` listing the projects it did find

The config directory honours `CHEMSMART_CONFIG_DIR`, defaulting to `~/.chemsmart`
(`CHEMSMARTUserSettings.resolve_config_dir()`).

## Per-program plumbing in user settings

`chemsmart/settings/user.py` has one property set per program:
`user_<prog>_settings_dir`, `<prog>_project_yaml_files`,
`all_available_<prog>_projects`. A new program needs its own set, or
`from_project` cannot find its YAMLs and the error message cannot list them.

## Templates

`chemsmart/settings/templates/template_orca_simple.yaml`,
`template_gaussian_simple.yaml`, `template_pyscf_simple.yaml`, and
`template_xtb_simple.yaml` are the
canonical minimal examples;
`templates/.chemsmart/` holds a full user-config tree used by `chemsmart config`.
A new program should ship an equivalent `template_<prog>_simple.yaml`.
