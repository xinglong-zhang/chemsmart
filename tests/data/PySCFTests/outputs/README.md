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
