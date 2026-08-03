# ChemSmart setting → PySCF attribute

## Translation table

| ChemSmart setting | PySCF target | Note |
|---|---|---|
| `functional` | `dft.KS(mol, xc=...)` | libxc aliases in `pyscf/dft/libxc.py`; reject known double-hybrid names and families before construction |
| `ab_initio == 'hf'` | `scf.HF(mol)` | bypasses `dft.KS` entirely |
| `basis` | `pyscf.M(basis=...)` | **hyphenated** (`def2-svp`), unlike Gaussian's `def2svp` |
| `charge` | `pyscf.M(charge=...)` | direct |
| `multiplicity` | `pyscf.M(spin=...)` | **`spin = multiplicity − 1`** (`spin = 2S = Nα−Nβ`) |
| `freq` | Hessian stage | project flag; `True` on `opt` adds `hess`, while the explicit CLI leaf is named `hess` |
| `solvent_model` | `mf.PCM()` / `mf.SMD()` | mapping below |
| `solvent_id` | SMD name / PCM dielectric lookup | required for every solvent model; must be a key of the target environment's `pyscf.solvent.smd.solvent_db` |
| `dispersion` | `mf.disp` | `'d3bj'`, `'d4'`; needs `pyscf-dispersion` |
| `density_fit` | `mf.density_fit(...)` | when false, `aux_basis` must be absent |
| `aux_basis` | `mf.density_fit(auxbasis=...)` | applied only with density fitting enabled |
| `defgrid` | `mf.grids.atom_grid` / `.level` | DFT only; reject on HF; see caveat below |
| `numfreq` | — | no numerical-Hessian analogue; reject at the CLI |
| `additional_route_parameters` | — | no route concept; reject or map explicitly |
| `heavy_elements*`, `gen_genecp_file` | `basis` as a per-element dict | PySCF takes `basis={'default':..., 'Pd':...}`, not GenECP text |
| `num_cores` (jobrunner) | `lib.num_threads(n)` | **not** from project YAML |
| `mem_gb` (jobrunner) | `pyscf.M(max_memory=MB)` | MiB convention: `mem_gb * 1024` |

Restricted vs unrestricted is **derived** from `multiplicity > 1`, not a separate
CLI flag. `scf.HF` / `dft.KS` dispatch on `mol.spin` automatically.

## Solvent mapping

`mf.PCM()` and `mf.SMD()` both exist on SCF/KS objects (verified). PCM's default
`with_solvent.method` is `'C-PCM'` (verified).

| ChemSmart `solvent_model` | PySCF |
|---|---|
| `cpcm` | `mf.PCM()`, `with_solvent.method = 'C-PCM'` |
| `pcm`, `iefpcm` | `mf.PCM()`, `with_solvent.method = 'IEF-PCM'` |
| `cosmo` | `mf.PCM()`, `with_solvent.method = 'COSMO'` |
| `ssvpe` | `mf.PCM()`, `with_solvent.method = 'SS(V)PE'` |
| `smd` | `mf.SMD()`, `with_solvent.solvent = <solvent_id>` |

Valid `with_solvent.method` values are `'C-PCM'`, `'IEF-PCM'`, `'COSMO'`,
`'SS(V)PE'` (`pyscf/solvent/pcm.py`). ChemSmart requires `solvent_id` for every
accepted model and resolves the matching database dielectric explicitly for PCM;
otherwise PySCF silently defaults PCM to water. That lookup occurs only in the
generated child process, against the final interpreter's solvent database. Do
not resolve it by importing PySCF in the controller environment.

For SMD specifically, a fresh `mf.SMD()` has
`with_solvent.solvent == ''` (verified) and fails when built. `solvent_db` has
**179 entries** in 2.14.0 and includes `water`. A ChemSmart `solvent_id` that is
valid for Gaussian or ORCA is not automatically in `solvent_db` — validate at
preflight, do not substitute a near match.

## Applied-science rejection rules

Reject known double-hybrid names and families. The v1 driver builds only an HF/KS
mean-field object and does not apply perturbative correlation; `b2plyp`,
`pbe0-2`, and `wb97x-2` are representative rejects. The last two are especially
hazardous because LibXC can parse their hyphens as subtraction rather than fail.

Also reject `defgrid` with `ab_initio: hf`, `aux_basis` when
`density_fit: false`, and a solvent model without `solvent_id`. These are not
optional cleanups: accepting them while dropping the inapplicable field would
make the HDF5 `spec/` receipt overclaim what the engine applied.

## The `defgrid` caveat

ORCA's `defgrid1/2/3` and PySCF's `atom_grid`/`level` are different schemes with
no official correspondence. Any table mapping them is a **ChemSmart convention,
not a PySCF fact**, and must be documented as such wherever it is defined —
otherwise a "same grid" comparison across programs is not one. `(99, 590)` is the
value gpu4pyscf's own driver uses as its default.

## Results contract

Treat `label.h5` as the sole machine-readable result. `label.out` is PySCF's
human log and is never parsed. Schema 2.0 stores four real top-level HDF5 groups;
the reader also accepts the frozen schema 1.0 layout, where `spec`, `provenance`,
and `status` were JSON datasets.

| HDF5 path | Representative contents | Contract |
|---|---|---|
| `spec/` | program, jobtype, stages, method/xc, basis, charge, spin, multiplicity, target-resolved solvent values, DF/grid/SCF/optimizer settings, engine, settings digest | resolved execution receipt; not target-object introspection |
| `provenance/` | PySCF/GPU4PySCF/libxc/ChemSmart/Python package versions, interpreter, host/platform, thread/CUDA environment, engine, timestamps, wall/core seconds, project and settings digests | reproduction and approval evidence |
| `status/` | `normal_termination`, nested per-stage convergence/iterations, typed `failure` or null | the only completion decision |
| `results/` | energies, final positions, atomic numbers, MO energies/occupancies, optional exact forces, Hessian/frequencies/modes, Mulliken charges, dipole, point group | typed scalar and array results |

Preserve result-array dtype and shape. Store forces as the exact negated final
gradient in Hartree/Bohr when the method provides one; omit the field rather than
inventing a placeholder when gradients are unavailable. Represent explicit nulls
as empty uint8 datasets marked with `chemsmart_is_null=True`, which distinguishes
"known unavailable" from a missing field.

The current `settings_digest` identifies the controller-side requested
configuration. It does not include target-resolved solvent dielectric values or
introspected runtime-object state. A future schema revision must add a separate
target-computed `applied_settings_digest` and recompute it during verification;
until then, do not describe receipt agreement as runtime-object attestation.

Compute `normal_termination` from every requested stage's structured convergence
record. A log file, subprocess return, readable HDF5 file, or converged SCF alone
must not establish completion. Use `verify_provenance()` to compare requested
settings against the stored receipt, provenance engine, and digests.

The public reader exposes total energies in Hartree, positions in Angstrom,
forces in Hartree/Bohr, harmonic frequencies in cm^-1, dipoles in Debye, and
orbital-energy properties in eV to match the Gaussian/ORCA reader boundary.

Rotational symmetry number comes from
`pyscf.hessian.thermo.rotational_symmetry_number(mol)` if a value is needed
directly; ChemSmart's own thermochemistry can also derive it from the geometry.

Store the geometry that produced the frequencies. For an `opt` with
`freq: True`, that is the optimized geometry, not the input one; for the explicit
`hess` leaf, it is the supplied fixed geometry.
