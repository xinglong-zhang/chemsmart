# What v3.1.4 actually has

Written to correct a v2.0.1-era model of ChemSmart. The package version remains
`3.1.4`; fork-specific capabilities are identified by the Git revision and the
joined live-schema/capability digest.

## xTB and CREST boundaries

The upstream v3.1.4 base contained a rich `chemsmart/io/xtb/` parser but no xTB
execution group. This fork adds only the maintained v3-compatible execution
plane: `run|sub xtb {sp,opt,hess}`, methods `gfn0|gfn1|gfn2|gfnff`, CPU
execution, and strict optional project YAML. The v3 reader remains the output
authority. CREST remains job preparation and an ensemble source rather than a
first-class `run|sub` program.

The integrated tree has five first-class `run`/`sub` program commands:

- **Gaussian**, **ORCA**, **xTB**, and **NCIPLOT** are external-binary backends.
- **PySCF** is a Python-library backend with `sp`, `opt`, and `hess` leaves.
  GPU4PySCF is selected as its `gpu` engine, not registered as another program.

All four have `Executable` subclasses, but `PySCFExecutable` resolves a Python
interpreter and does not require an `EXEFOLDER`. A configured `EXEFOLDER` names
the environment's `bin/`; when omitted, it uses the interpreter running
ChemSmart.

## Top-level groups

From `chemsmart/cli/subcommands.py`: `gaussian`, `orca`, `xtb`, `pyscf`, `mol`,
`grouper`, `nciplot`, `thermochemistry`, `database`, `iterate`, `pka`.

Notable capability present in v3.1.4:

- **`database`** — a ChemSmart `.db` store, with `--ri/--record-index`,
  `--rid/--record-id`, `--sid/--structure-id` selectors usable directly as job
  input on the gaussian/orca groups.
- **`pka`** — including table and CDXML batch modes, which bypass the normal
  filetype validation in the program group.
- **`iterate`**, **`thermochemistry`** (with Boltzmann averaging), **`grouper`**
  (rmsd/tanimoto/isomorphism/formula/energy and more).
- ORCA **`neb`**, and **`qmmm`** as a *nested* subcommand attached under both
  `opt` and `sp` (`chemsmart run orca ... opt qmmm`), not a top-level one —
  `create_orca_qmmm_subcommand()` attaches it, and `cli/orca/__init__.py` keeps
  the direct import commented out.
- ORCA solvation beyond CPCM/SMD: **`cpcmc`** (CPCM with COSMO epsilon) and
  **`cosmors`** (openCOSMO-RS, with a `-sf/--solventfilename` `.cosmorsxyz` path).
- Multi-structure input: `-i/--index` yields one job per selected molecule, with
  original indices tracked in `ctx.obj["molecule_indices"]` for labelling.
- PySCF requires a named project YAML and exposes exactly `sp`, `opt`, and
  `hess`. Its generated `.out` is a human log; schema-2 `.h5` groups are the
  machine contract. `freq: True` on `opt` composes the Hessian after the
  optimization, while `hess` analyzes the supplied fixed geometry.

## Program capability registry

`chemsmart/settings/capabilities.py` is the immutable declaration for the five
executable program commands. It records project requirements, immediate Click
children, project-owned CLI parameter names, and execution engines. Its focused
test compares the declared program/jobtype inventory with the live Click objects.

Use the canonical `PROGRAM_*` names described in
[agent-harness-forward-compatibility.md](agent-harness-forward-compatibility.md).
The old `ENGINE_*`, `PROGRAM_JOBTYPES`, `PROGRAM_ENGINES`, and
`PROJECT_OWNED_PARAMETERS` names are compatibility aliases for the v2 harness,
not additional registries.

## Dependency floor

From `pyproject.toml`: Python `~=3.10`, `click==8.1.8`, **`numpy~=1.26.4`**
(pinned low because rdkit and pymol require numpy 1.x), `ase==3.24.0`,
`rdkit==2025.3.3`, `pymatgen==2025.10.7`, `scipy==1.15.2`, `pandas==2.2.3`.

The numpy pin matters for any new backend: a program whose own stack needs
numpy 2.x cannot share ChemSmart's environment and must rely on the per-program
`CONDA_ENV` block in the server YAML.

## Developer commands

`make env` (conda from `environment.yml`), `make install-dev`, `make configure`,
`make fmt` (black + isort), `make lint` (ruff), `make test` (lint + coverage),
`make docs`. The Makefile wraps commands in `conda run -n chemsmart` unless
`CONDA_DEFAULT_ENV` is already `chemsmart`.

Single test: `pytest tests/test_orca_cli.py::TestClass::test_name`.
