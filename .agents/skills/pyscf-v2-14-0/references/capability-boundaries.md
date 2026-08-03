# Capability boundaries

## Hard limits — cannot be wrapped

- **No IRC.** PySCF 2.14.0 contains **no IRC module**: zero files matching
  `*irc*` and zero mentions of "intrinsic reaction coordinate" in the tree
  (verified). An IRC in ChemSmart-on-PySCF would have to be *implemented* — a
  reaction-path follower on top of gradients — not wrapped. Do not offer `irc` as
  a PySCF jobtype.
- **No native TS optimizer.** Transition-state search exists only by passing
  `transition=True` through `pyscf.geomopt` to geomeTRIC. It carries geomeTRIC's
  requirements (a reasonable Hessian guess) and is only as good as that solver.
- **Optimization does not imply frequencies.** A converged `opt` says nothing
  about whether the structure is a minimum or a saddle point. Run the explicit
  `hess` leaf at the same method and basis, or set `freq: True` on `opt` to chain
  the Hessian stage against the optimized geometry, before any stationary-point
  claim.
- **No double hybrids in the v1 driver.** The driver builds an HF or KS
  mean-field object but never applies a double hybrid's perturbative correlation
  term. All double hybrids remain unsupported. Reject recognized names and
  families, including `b2plyp`, `pbe0-2`, and `wb97x-2`; unfamiliar custom
  expressions require bound functional metadata/classification before they can
  be accepted. LibXC parsing a name does not mean the advertised method ran.

## Applied-setting consistency limits

Reject a DFT `defgrid` on HF, an `aux_basis` when density fitting is disabled,
and every solvent model without an explicit `solvent_id`. Silently discarding
one of these requests would make the applied-settings receipt false. Resolve PCM
dielectric data from `pyscf.solvent.smd.solvent_db` only inside the generated
child process, because its interpreter/environment is authoritative.

## Optional dependencies — require bound evidence, not an import guess

Neither is a core PySCF dependency, and both fail *at call time*:

| Feature | Extra | Failure mode |
|---|---|---|
| `opt` (geomeTRIC) | `pyscf[geomopt]` → `geometric` (or `pyberny`) | `from pyscf.geomopt import optimize` **succeeds**; `optimize(mf)` raises `ImportError` (verified on a host without `geometric`) |
| `mf.disp` | `pyscf[dispersion]` → `pyscf-dispersion` | fails during SCF |

Preflight must receive a callable-entry-point attestation bound to the finalized
compute interpreter/configuration. It inspects that injected evidence but never
invokes a solver or launches a controller-side runtime probe. An import check is
a false green. The current CLI does not supply this mandatory binding, so treat
preflight as a callable surface rather than an active CLI gate.

## Deferred but mapped

Present in 2.14.0, out of the v1 CLI scope. Recorded so a later phase starts from
a known location rather than rediscovering it:

| Capability | Import path |
|---|---|
| TDDFT / TDA excited states | `pyscf.tdscf`, `pyscf.tddft` |
| MP2 | `pyscf.mp` |
| Coupled cluster (CCSD, CCSD(T), EOM) | `pyscf.cc` |
| CASCI / CASSCF / NEVPT2 | `pyscf.mcscf`, `pyscf.mrpt` |
| MC-PDFT family | `pyscf.mcpdft` |
| GW / RPA | `pyscf.gw` |
| ADC | `pyscf.adc` |
| Periodic boundary conditions | `pyscf.pbc` |
| QM/MM | `pyscf.qmmm` |
| Molecular dynamics | `pyscf.md` |
| Localized orbitals | `pyscf.lo` |
| Relativistic (X2C) | `pyscf.x2c` |
| Non-adiabatic couplings | `pyscf.nac` |
| Analytic gradients / Hessians | `pyscf.grad`, `pyscf.hessian` |
| Molden / cube export | `pyscf.tools` |

`td` is the natural next addition: `pyscf.tdscf` has a GPU counterpart, and the
project-YAML `td:` section already exists for Gaussian and ORCA.

## Version note

The pinned integration target is PySCF **2.14.0**, matching the
`pyscf==2.14.0` pin in `gpu4pyscf==1.8.0`'s `requirements.txt`. The CPU
integration and all three ChemSmart leaves have been runtime-exercised with
PySCF 2.14.0; GPU behavior remains source-read until exercised on a CUDA host.

Claims in these skills are marked by how they were checked:

- **Source-verified against 2.14.0** — no-IRC, `solvent_db` size (179),
  PCM method names, the `spin` property definition.
- **Runtime-verified against 2.14.0** — the ChemSmart CPU `sp`, `opt`, and
  `hess` paths, spin/multiplicity relation, solvent attachment, libxc aliasing,
  double-hybrid rejection, call-time optimizer behavior, structured results,
  and local/scheduler command parity. Treat exact numeric references as receipts
  for their recorded environment, not universal expected values.
