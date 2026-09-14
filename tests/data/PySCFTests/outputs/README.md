# Archived PySCF results

Real PySCF 2.14.0 / libxc 7.0.0 runs produced on chemsmart-hpc on
2026-09-12 through the human CLI (`chemsmart run --no-fake --no-scratch
-n 4 -m 8 pyscf -p <project> -f <input> [-c C -m M] -l <label> <stage>`)
at commit b4fce1c7 (result contract v4); the three Hessians were re-run at cb99bbff after the runner stopped producing the `unclassified` word, so their receipts say `validated`. Each directory holds the
artifact (`.h5`), its three receipts (`.receipt.json`, `.input.json`,
`.environment.json`), PySCF's own log (`.out`, never parsed) and
`.reference.json`, PySCF's independent account of the same bytes written
by `reference.py` in the compute env (harmonic analysis and RRHO
thermochemistry from the stored Hessian, point group, isotope-averaged
masses). The generated driver script is deliberately not archived. The
receipts carry this host's paths, as the archived ORCA and xTB logs do;
admission compares digests, never paths.

| directory | what it is | why it is here |
|---|---|---|
| `water_sp` | B3LYP/def2-SVP single point at a distorted water (O–H 1.10 Å, 90°) | a green `sp`; supplied and final structures coincide |
| `water_opt` | the same start optimised (moved 0.083 Å) | supplied ≠ reached; `converged`; `energies` = [E(supplied), E(reached)] |
| `water_hess` | Hessian on the converged geometry (`-f water_opt.h5`) | three real modes; gradient 8e-6 Eh/Bohr; the thermochemistry differential oracle |
| `nh3_planar_opt` | exactly planar D3h ammonia optimised: stays planar | a converged saddle from a symmetric start |
| `nh3_planar_hess` | its Hessian: one imaginary mode at −830 cm⁻¹, two degenerate pairs | `failed_wrong_stationary_point` must be typed; PySCF's `detect_symm` says Cs at 1e-5 Bohr while the structure is D3h to 1e-14 Å in z |
| `hydroxyl_sp` | OH· UKS doublet single point | spin populations [1.026, −0.026]; ⟨S²⟩ 0.7518 |
| `water_opt_maxsteps1` | the distorted water with `opt_maxsteps: 1` | an unconverged optimisation: `results/positions` is the last evaluated geometry (0.066 Å from the input); receipt state `failed` |
| `water_stretched_hess` | Hessian at the converged water with O–H(1) +0.02 Å | three real frequencies at max\|g\| = 0.0185 Eh/Bohr: zero imaginary modes is not stationarity |
| `water_hess_historical_unclassified` | the water Hessian as the pre-cb99bbff runner wrote it: receipt `engine_complete` / `unclassified` | historical receipts keep admitting after the per-plan policy was retired |
| `water_sp_scfmaxiter2` | the distorted water with `scf_maxiter: 2` | a quiet SCF non-convergence: `normal_termination` false, `energies` present |

## Expansion round (2026-09-13, result contract v5)

Produced the same way at commit 9e219c8f (the driver's response,
excited-surface and correlated stages), on the relaxed water of
`inputs/water_relaxed.xyz` (the `water_opt` minimum) unless the row says
otherwise; the electronic state was bound on the command line (`-c 0 -m 1`,
`-c 0 -m 2` for the hydroxyl radical). `reference.py` rebuilds the mean
field from the applied spec and recomputes the TDA/TDDFT roots, oscillator
strengths and transition dipoles, or the MP2/CCSD/CCSD(T) components, with
PySCF's public API; every green fixture agrees with that recomputation to
1e-13 Eh (3e-5 Eh for the UKS hydroxyl case, an SCF re-convergence
difference) and every correlated total to 1e-8 Eh. `orca_differential/`
holds ORCA 6.1.1 runs on the same relaxed water: TDA with `B3LYP/G` (the
VWN3 functional PySCF's `b3lypg` names) agrees with `water_td_singlet` to
1 meV on all three roots and on every oscillator strength; MP2 with ORCA's
default frozen core reproduces `water_mp2_sp_fc1` (frozen_core 1) and
`NoFrozenCore` reproduces `water_mp2_sp` (PySCF's all-electron default),
each to 5e-8 Eh in the correlation energy -- the frozen-core divergence is
a convention, and both programs agree once it is named.

| directory | what it is | why it is here |
|---|---|---|
| `water_td_singlet` | TDA-B3LYP(G)/def2-SVP, three singlet roots | ascending roots, oscillator strengths, transition dipoles (Debye), per-root convergence; the ORCA differential |
| `water_td_triplet` | the triplet manifold | multiplicity 3 records, zero oscillator strengths |
| `water_td_rpa` | full TDDFT (RPA) singlets | differs from TDA on the same reference; both ascending |
| `hydroxyl_td_unrestricted` | OH· UKS-TDA, the one unrestricted manifold | no multiplicities, oscillator strengths served, no per-root ⟨S²⟩ |
| `water_td_cpcm_toluene` | TDA singlets under C-PCM toluene | the artifact records `static_eps_applied` 2.3741 and `response_eps_applied` 1.78: PySCF's non-equilibrium response uses water's optical dielectric for every solvent |
| `water_td_unconverged` | `td_max_cycle: 1` through the CLI | every root unconverged; `stages/td/converged` false; receipt `failed` -- the typed ending, with the per-root flags and energies still inspectable |
| `formaldehyde_s1_opt` | H₂CO from a symmetry-broken start (`inputs/formaldehyde_bent_start.xyz`), root 1 of **one** requested root | `excited_state_root == nstates`, the case where PySCF 2.14's scanner `converged` property raises; converged, final gradient 2.5e-5 Eh/Bohr, gap to the ground state 3.05 eV at the reached geometry |
| `formaldehyde_s1_opt_planar` | H₂CO from the exactly planar start, root 1 of three | a symmetric excited stationary point the host cannot characterise (no Hessian): planar stays planar, emission 3.406 eV, neighbour gap 4.59 eV |
| `water_s1_opt_degenerate` | water, root 1 of three | the followed root ends **degenerate with its neighbour** (gap 3e-6 eV) after a 0.52 Å move: root identity is undecidable from the index, which the gap sensor facts say |
| `formaldehyde_s1_td` | TDA singlets, three roots, on `inputs/formaldehyde_s1_reached.xyz` -- the XYZ `bind_reached_geometry` wrote from `formaldehyde_s1_opt` | the response consumer of an excited-surface producer over real bytes: supplied == reached to 5e-11 Å, root 1 is the producer's end gap (3.052 eV, the emission), the reference energies agree to 3e-12 Eh |
| `water_mp2_sp`, `water_mp2_sp_fc1` | MP2, all-electron and `frozen_core: 1` | `reference_energy`, `correlation_energy`, `total_energy`; the frozen-core convention against ORCA |
| `water_ccsd_sp`, `water_ccsdt_sp` | CCSD; CCSD(T) with `frozen_core: auto` (1 orbital applied) | `ccsd_correlation_energy`, `triples_correction`, `correlation_energy` = their sum |
| `hydroxyl_ump2_sp` | OH· UMP2 on a UHF reference | an open-shell correlated result |
| `water_mp2_opt`, `water_ccsd_opt` | the distorted water optimised on the MP2 and on the CCSD surface | correlated structure producers: amplitude and Λ convergence, final gradient 6e-6 Eh/Bohr, components at the reached geometry |
| `water_ccsd_unconverged` | `cc_max_cycle: 1` through the CLI | amplitudes unconverged; `stages/corr/converged` false; receipt `failed` |

## Surface round (2026-09-14, result contract v6)

Produced the same way in the compute env, on the distorted water of
`inputs/water_distorted.xyz`, with `reference.py` beside each. These
four are the first artifacts to record a `surface`: the electronic
surface a result's geometry and total energy belong to. They exist as a
matched and a mismatched pair, because the question a Hessian answers is
about one surface and the host had no way to ask which.

| directory | what it is | why it is here |
|---|---|---|
| `water_mp2_opt_v6` | MP2/def2-SVP optimisation | surface `mp2:rhf:-:def2-svp:0:1:-:-:-:-:0:-:-:-`; all-electron, so the frozen-core count is 0 rather than absent |
| `water_dft_hess_on_mp2_geometry` | B3LYP Hessian on that MP2 minimum (`-f water_mp2_opt.h5`) | the mismatch: a validated Hessian that characterises a *different* surface from the geometry it consumed |
| `water_dft_opt_v6` | B3LYP/def2-SVP optimisation of the same start | surface `dft:rks:b3lypg:def2-svp:0:1:-:-:-:-:-:-:-:-` |
| `water_dft_hess_on_dft_geometry` | B3LYP Hessian on that B3LYP minimum | the control: same surface, so the Hessian characterises the structure |

One measurement worth keeping beside them. A B3LYP Hessian at the
geometry the S1 (TDA root 1) optimisation of formaldehyde reached is
refused by the validator as `pyscf.result.hessian_invalid`: PySCF's
analytic Hessian there carries a raw antisymmetry of 5.0e-05 Eh/Bohr^2
against a limit calibrated near stationarity, because the ground state
at that geometry is not stationary at all -- its gradient is 0.0998
Eh/Bohr. A cross-surface Hessian is not merely uninformative; at a
geometry far from its own stationary point it is numerically worse too.

