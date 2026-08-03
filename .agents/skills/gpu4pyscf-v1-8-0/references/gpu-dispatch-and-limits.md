# GPU dispatch and limits

All source-read from the GPU4PySCF 1.8.0 checkout; not runtime-verified (no CUDA
on this host).

## Must produce preflight violations — no silent CPU fallback

From the upstream README's Limitations section:

| Limitation | Consequence for ChemSmart |
|---|---|
| Atomic basis up to **g** orbitals | reject a basis with h/i functions |
| Auxiliary basis up to **i** orbitals | reject the auxbasis, not the basis |
| **Double-hybrid functionals unsupported** | reject; do not drop the MP2 part |
| meta-GGA with density Laplacian | reject |
| Hessian of TDDFT unsupported | `freq` on a TDDFT state is CPU-only |
| DF bounded ≈168 atoms at def2-tzvpd | bounded by **host** RAM, not VRAM |

A request that trips one of these must yield a typed violation and be rejected
once the validation surface is safely wired as a gate. Quietly running it on CPU
changes the engine behind the user's back, and quietly dropping a term (a double
hybrid's MP2 contribution) changes the *method* — both are scientific changes
disguised as fallbacks. B3 is currently callable but is not a CLI gate.

The ≈168-atom bound is worth calling out: it scales with **CPU** memory, so a
larger GPU does not lift it and a GPU-memory check will not predict it.

## Solvent parity is partial

`gpu4pyscf/solvent/` provides `pcm.py` and `smd.py` only. PySCF's `ddcosmo`,
`ddpcm`, and `pol_embed` have **no GPU counterpart**.

A `ddcosmo` request must not silently become PCM — they are different cavity
models and give different solvation energies. Either run that job on CPU
(explicitly) or reject it.

## Production vs experimental

The README separates these. Only the first group should back a ChemSmart jobtype
without an explicit caveat.

**Production:** DF and direct SCF; SCF, analytic gradients, analytic Hessian for
HF and DFT; LDA/GGA/mGGA/hybrid/range-separated via libXC; TDA/TDDFT
(spin-conserved and spin-flip); geometry optimization and TS search via geomeTRIC;
DFTD3/DFTD4 dispersion; VV10 nonlocal gradients and Hessian; GPU ECP; PCM with
gradients and Hessian; SMD; unrestricted HF/DFT; CHELPG/ESP/RESP charges.

**Experimental** — do not wire into a CLI leaf without saying so: MP2/DF-MP2 and
CCSD; polarizability, IR, NMR shielding; Raman; QM/MM with PBC; multi-GPU;
periodic SCF/DFT; TDDFT non-adiabatic coupling; energy decomposition analysis.

Note `gpu4pyscf/properties/` (`ir.py`, `raman.py`, `shielding.py`,
`polarizability.py`, `c6.py`, `eda.py`) and `gpu4pyscf/pop/esp.py` exist as
modules — existence is not production readiness.

## Numerical agreement

GPU and CPU results are not bit-identical: different contraction order, and
GPU4PySCF's defaults differ from PySCF's in places (`lebedev_order = 29` on
solvent, `direct_scf_tol = 1e-14`, `scf_conv_tol = 1e-10` in the reference
driver).

Consequences for ChemSmart:

- Never mix engines within a comparison set. A reaction energy where reactant ran
  on GPU and product on CPU is not a controlled comparison.
- Record the engine and versions per job so a set can be audited for mixing.
- When reproducing an archived number, reproduce the engine too.

## Reference implementations

Worth reading before writing the ChemSmart writer — they are the validated
config→object mapping:

- `gpu4pyscf/tools/method_config.py` — `get_default_config()` and
  `method_from_config()`; the canonical construction order.
- `gpu4pyscf/drivers/dft_driver.py` — energy → gradient → Hessian → thermo in one
  process, with results written to HDF5. The stage-composition pattern.
- `gpu4pyscf/drivers/opt_driver.py` — geomeTRIC via the CPU solver, and the
  `constraints.txt` file convention geomeTRIC expects.
