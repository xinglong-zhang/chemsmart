---
name: pyscf-v2-14-0
description: Drive PySCF 2.14.0 through ChemSmart's sp, opt, and hess jobs, including Mole and SCF/DFT construction, project-setting translation, CPU/GPU provenance, structured HDF5 results, typed preflight, and hard capability limits. Use when adding or debugging ChemSmart PySCF jobs, mapping method/basis/solvent/charge/multiplicity settings, reading label.h5, verifying applied settings, or deciding whether a requested calculation is supported.
---

# PySCF 2.14.0 for ChemSmart

A deliberately reduced view: the subset needed to drive PySCF from
`chemsmart run|sub pyscf {sp,opt,hess}`. PySCF's full surface is large and mostly
irrelevant here — see
[capability-boundaries.md](references/capability-boundaries.md) for what is
deferred and what is impossible.

**PySCF is a library, not a binary.** ChemSmart's program conventions assume an
executable that reads an input file. Read
`chemsmart-v3-1-4` → `references/adding-a-program.md` for which of those
conventions to keep and which to drop.

## The whole v1 recipe

Execution model: ChemSmart emits a standalone script whose logic is a fixed
skeleton and whose only variable part is a dict of scalars, then runs it in a
subprocess. The script depends on **pyscf / numpy / h5py only** — never on
chemsmart, so the program's own `CONDA_ENV` can differ from ChemSmart's.

```python
import pyscf
from pyscf import scf, dft, lib

lib.num_threads(NUM_CORES)                    # from the jobrunner, not the YAML
mol = pyscf.M(atom=..., basis=..., charge=..., spin=MULT - 1,
              unit='Angstrom', max_memory=MEM_MB, output='label.out', verbose=4)

if xc is None:
    mf = scf.HF(mol)
else:
    mf = dft.KS(mol, xc=xc)
    mf.grids.atom_grid = (99, 590)                       # or .level = N
if dispersion:  mf.disp = dispersion                    # 'd3bj' | 'd4'
if with_df:  mf = mf.density_fit(auxbasis=auxbasis)     # before .to_gpu()
if engine == 'gpu':  mf = mf.to_gpu()
if solvent_model == 'smd':
    mf = mf.SMD()                                        # after .to_gpu()
    mf.with_solvent.solvent = solvent_id
elif solvent_model:
    mf = mf.PCM()
    mf.with_solvent.method = pcm_method
    mf.with_solvent.eps = target_environment_dielectric
```

Then, by stage:

- **sp** — `e = mf.kernel()`; record `mf.converged` (never assume it).
- **opt** — `from pyscf.geomopt import optimize;
  mol_eq = optimize(mf, maxsteps=..., constraints=...)`.
- **hess** — `hess = mf.Hessian().kernel()` then
  `from pyscf.hessian import thermo; nm = thermo.harmonic_analysis(mf.mol, hess)`.

**Stop at frequencies.** Do *not* call `thermo.thermo()`. ChemSmart already
computes free energies in `chemsmart/analysis/thermochemistry.py`, with its own
`rotational_mode` and quasi-harmonic conventions that Gaussian and ORCA results
already flow through. A second thermochemistry path yields silently disagreeing
free energies for the same molecule across programs.

## Keep the generated artifact contract

Treat these four siblings as one receipted run:

- `label.py` — generated fixed skeleton plus a JSON configuration; never edit it.
- `label.out` — PySCF's human-facing log; never parse it for program state.
- `label.h5` — the sole machine contract.
- `label.err` — subprocess stderr.

Read HDF5 schema 2.0 through real `spec/`, `provenance/`, `status/`, and
`results/` groups. Preserve schema 1.0 only as a strict reader-compatibility
path. Use [chemsmart-settings-map.md](references/chemsmart-settings-map.md) for
the field and unit contract.

## Compose stages in one process

Do not schedule separate `opt` and frequency jobs when one composed run was
requested: that re-converges SCF and re-pays DF factorization and grid
construction. Set `freq: True` on the `opt` project settings to run the ordered
stage list `scf → opt → hess` against one mean-field object. Use the explicit
`hess` leaf only for a supplied geometry that is already stationary at the same
method and basis.

## Validate without executing chemistry

Call `preflight(settings, molecule, environment)` with evidence for the actual
compute interpreter, not facts imported from ChemSmart's possibly different
environment. Collect every returned `PySCFViolation`; do not expect the function
to raise or invoke a solver.

Treat this as a callable API, not a current CLI gate. The CLI does not yet
mandatorily bind injected evidence to the finalized compute interpreter,
configuration identity, and execution phase. Do not wire it into job
construction until the frozen integration contract in
`chemsmart-v3-1-4` →
`references/agent-harness-forward-compatibility.md` is deliberately resolved.

After a run, call `verify_provenance(settings, h5_path)` and require no
violations before claiming that a complete artifact records the requested
settings. It compares the request with `spec/`, the CPU/GPU engine, settings
digest, and available project digest. It accepts both supported HDF5 schema
versions and turns unreadable or mismatched artifacts into typed findings. This
is receipt verification, not target-object introspection: the current digest is
computed from the controller-side request, and no separate target-side applied
digest is frozen yet.

## Reject settings that would overclaim the applied science

- Reject known double-hybrid names and families, including `b2plyp`, `pbe0-2`, and
  `wb97x-2`. This v1 driver constructs only a mean-field object and cannot apply
  the perturbative correlation term. Some hazardous names are nevertheless
  parseable by LibXC, so parser acceptance is not evidence of implementation.
- Reject `defgrid` on an HF calculation, `aux_basis` when `density_fit` is
  false, and a solvent model without an explicit `solvent_id`. Omitting or
  dropping any of those settings would make `spec/` claim a calculation that
  was not actually applied.
- Resolve PCM dielectric data only inside the generated child process. The
  controller and final compute interpreter may use different environments and
  solvent databases.

## Two more things that will silently produce wrong numbers

- **`mol.spin` is not multiplicity.** `spin = 2S = Nα − Nβ`, so
  `spin = multiplicity − 1`. Verified: `pyscf.M(atom='O 0 0 0', spin=2)` gives
  `multiplicity == 3` and `nelec == (5, 3)`. Passing multiplicity straight into
  `spin` silently computes the wrong electronic state.
- **Basis names are hyphenated**: `def2-svp`, `def2-tzvpp` — unlike Gaussian's
  `def2svp`. A wrong name raises, but a *valid but different* name does not.

## Optional dependencies fail at call time, not import time

`from pyscf.geomopt import optimize` succeeds even when neither `geometric` nor
`pyberny` is installed; `optimize(mf)` then raises `ImportError` after the job has
already started. Verified on a host without `geometric`. **An import check is a
false green** — preflight must receive evidence bound to the final compute
interpreter that the selected entry point is callable, and likewise that
`pyscf-dispersion` is available before setting `mf.disp`. Do not add an
unbounded controller-side runtime probe.

## Use the references

- [chemsmart-settings-map.md](references/chemsmart-settings-map.md) — the
  ChemSmart-setting → PySCF-attribute table and the results contract.
- [capability-boundaries.md](references/capability-boundaries.md) — deferred
  modules with import paths, and the hard limits (notably: **no IRC exists**).

For program-name, Click-jobtype, project-ownership, and engine projections, use
`chemsmart-v3-1-4` →
`references/agent-harness-forward-compatibility.md`. Do not reinterpret the
Click inventory as proof that an agent validator or concrete JobKind exists.

For the GPU variant use `gpu4pyscf-v1-8-0`; it is a delta on this skill.
