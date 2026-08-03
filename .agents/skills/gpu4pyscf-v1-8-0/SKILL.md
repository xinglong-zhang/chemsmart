---
name: gpu4pyscf-v1-8-0
description: How GPU4PySCF 1.8.0 differs from PySCF when driven by ChemSmart - the .to_gpu() switch and its ordering, the CUDA/CuPy/cuTENSOR version lock, GPU/CPU dispatch from server settings, and the features that must produce typed violations rather than silently fall back. Use when enabling GPU for a chemsmart pyscf job, diagnosing a GPU/CPU discrepancy, or checking whether a method is GPU-supported.
---

# GPU4PySCF 1.8.0 for ChemSmart

A **delta** on `pyscf-v2-14-0`. The Mole, settings translation, results contract,
and stage composition are identical — use that skill for all of it. Only what
changes on GPU is here.

## The switch

One call converts a built method object:

```python
mf = dft.KS(mol, xc=xc)
if with_df: mf = mf.density_fit(auxbasis=auxbasis)   # BEFORE .to_gpu()
mf.grids.atom_grid = (99, 590)                        # BEFORE .to_gpu()
if use_gpu:  mf = mf.to_gpu()
if solvent:  mf = mf.PCM()                            # AFTER .to_gpu()
             mf.with_solvent.lebedev_order = 29
```

**Ordering is not cosmetic.** Density fitting and grids are configured on the CPU
object; solvent is attached to the already-converted GPU object. This is the
order in `gpu4pyscf/tools/method_config.py::method_from_config`, the reference
implementation. GPU4PySCF sets `lebedev_order = 29` on the solvent object where
the CPU path leaves the default.

`pyscf.lib.misc.to_gpu` rewrites the module path `pyscf.*` → `gpu4pyscf.*` and
prefers a `from_cpu` classmethod when the GPU class defines one. It raises
`ImportError` with install instructions when `gpu4pyscf` is absent.

## Geometry optimization has no GPU driver

`gpu4pyscf/geomopt/` contains only `ase_solver.py`. GPU4PySCF's own
`opt_driver.py` imports **`pyscf.geomopt.geometric_solver.kernel`** — the CPU
driver — and hands it a GPU method object; only the gradient evaluations run on
GPU.

So `opt` emits the *same* geomeTRIC call on both paths and only `mf` differs.
Write one code path, not two.

## Version lock

`gpu4pyscf==1.8.0` pins **`pyscf==2.14.0`** exactly (`requirements.txt`), along
with `cupy-cuda12x==13.4.1`, `cutensor-cu12==2.2.0`,
`gpu4pyscf-libxc-cuda12x==0.8.1`, `pyscf-dispersion==1.5.0`, `geometric==1.1.0`.
A ChemSmart pin must preserve the PySCF equality, and the CUDA-suffixed packages
must match the host toolkit generation (`cuda11x`/`cuda12x`/`cuda13x`).

CuPy and cuTENSOR are interdependent; the upstream-recommended pairs are
CuPy 13.3.0 + cuTENSOR 2.0.2, or CuPy 13.4.1 + cuTENSOR 2.2.0.

**A mismatch degrades silently.** It surfaces as a `UserWarning` — "using cupy as
the tensor contraction engine" — not an exception, so the job still produces
numbers, just slowly. This needs an explicit typed preflight rule fed by bound
environment evidence; a `try/except` around the import will not catch it. B3 is
currently a callable validation surface, not a wired CLI gate.

Because these pins conflict with ChemSmart's own `numpy~=1.26.4`, the GPU stack
belongs in its own environment, selected via the `CONDA_ENV` key of the server
YAML program block. This is why the generated script must not import chemsmart.

## Dispatch

Resolve in the same shape as the scratch ladder:

1. explicit `--gpu` / `--no-gpu`
2. otherwise, server YAML `SERVER.NUM_GPUS > 0`
3. otherwise CPU

Record the resolved engine and the package versions **in the emitted script and
in the results file**, so an archived artifact states which engine produced the
number. The same command on two hosts must not silently mean two different
calculations.

## Limits

Read [gpu-dispatch-and-limits.md](references/gpu-dispatch-and-limits.md) before
enabling GPU for a job. Several methods are unsupported and must produce typed
preflight violations, then be rejected once a safe gate is wired, rather than
fall back to CPU — an unrequested engine change is a scientific change, not an
optimization.

## Verification status

Everything here is **source-read** from the GPU4PySCF 1.8.0 checkout and its
README. None of it is runtime-verified: this host has no CUDA and `import
gpu4pyscf` fails at `import cupy`. Treat GPU claims as documented behaviour to
be confirmed on a GPU host, never as observed behaviour.
