---
name: gpu4pyscf-v1-8-0
description: Configure, qualify, debug, or improve GPU4PySCF 1.8.0 as the GPU execution engine of ChemSmart's PySCF program. Use for CUDA, CuPy, cuTENSOR, GPU conversion order, density fitting, solvents, GPU-supported methods, CPU/GPU discrepancies, or fail-closed GPU availability.
---

# GPU4PySCF 1.8.0 for ChemSmart

## Mission

Act as a computational scientist qualifying a GPU engine, not as a dispatcher
that assumes a declared GPU is usable. GPU4PySCF is a PySCF engine and shares
the same molecule, project settings, generated script and HDF5 result contract.

## Execution order

Build the PySCF method, configure density fitting and grids, convert with
`to_gpu()`, then attach supported solvent behavior. Geometry optimisation uses
the PySCF optimisation driver with GPU energy and gradient evaluations.

The CUDA toolkit, NVIDIA driver, CuPy, cuTENSOR, GPU4PySCF library and PySCF
versions must be mutually compatible. Confirm imports, device visibility,
memory and a real kernel. Never fall back silently to CPU.

## Scientific limits

GPU support is method and property specific. Validate the requested SCF/DFT,
gradient, Hessian, solvent and density-fitting path rather than inferring
support from package import. Compare CPU and GPU energies under identical
settings before using a new pathway for production research.

The x86 CPU server workflow may preview GPU intent but must keep GPU execution
unavailable when no qualified device exists.

## Self-improvement cycle

For a GPU failure or discrepancy, separate driver/runtime, package ABI, device
memory, conversion order, unsupported method and numerical convergence. Run a
minimal CPU/GPU matched calculation, improve the smallest general preflight,
mapping or diagnostic, and challenge it with another molecule or property.

Update this skill when observed version behavior changes. Do not encode one
machine, molecule or benchmark result as a universal GPU rule.
