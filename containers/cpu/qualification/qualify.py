#!/usr/bin/env python3
"""Qualify the ChemSmart CPU image through its public YAML/CLI path.

The calculations in this module are deliberately launched only through
``chemsmart run``.  Python is used after execution to inspect ChemSmart's
structured artifacts and to write one small, public scientific summary.
"""

from __future__ import annotations

import argparse
import importlib
import importlib.metadata
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import re
import shutil
import subprocess
import sys

import numpy as np
from ase.io import read as ase_read

from chemsmart.agent._contracts import TrustedArtifactRefV1, file_sha256
from chemsmart.agent.postprocessing import derive_trusted_thermochemistry
from chemsmart.agent.skills import resolve_skill
from chemsmart.io.pyscf.output import read_pyscf_h5
from chemsmart.io.xtb.output import XTBOutput


PYTHON_PACKAGE_VERSIONS = {
    "chemsmart": "3.1.4",
    "pyscf": "2.14.0",
    "numpy": "1.26.4",
    "geometric": "1.1.1",
    "pyscf-dispersion": "1.5.0",
    "h5py": "3.16.0",
}

CONDA_PACKAGE_VERSIONS = {
    "xtb": "6.7.1",
    "libopenblas": "0.3.34",
    "pymol-open-source": "3.1.0",
    "ffmpeg": "7.1.1",
}

QUALIFICATION_NUM_CORES = 2
OPENBLAS_OPENMP_WARNING = (
    "OpenBLAS Warning: Detect OpenMP Loop and this application may hang."
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def run(command: list[str], *, cwd: Path) -> subprocess.CompletedProcess[str]:
    """Run one public ChemSmart command without persisting its verbose log."""
    completed = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if completed.returncode != 0:
        tail = "\n".join(completed.stdout.splitlines()[-80:])
        raise RuntimeError(
            f"Command failed with exit {completed.returncode}: "
            f"{' '.join(command)}\n{tail}"
        )
    return completed


def run_expected_failure(
    command: list[str], *, cwd: Path, required_text: str
) -> None:
    completed = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    require(
        completed.returncode != 0, f"Command unexpectedly succeeded: {command}"
    )
    require(
        required_text.lower() in completed.stdout.lower(),
        f"Expected refusal containing {required_text!r}: {completed.stdout[-2000:]}",
    )


def package_versions() -> dict[str, str]:
    require(
        sys.version_info[:2] == (3, 10),
        f"Python 3.10 required, observed {sys.version_info.major}.{sys.version_info.minor}",
    )
    observed = {
        name: importlib.metadata.version(name)
        for name in PYTHON_PACKAGE_VERSIONS
    }
    for name, expected in PYTHON_PACKAGE_VERSIONS.items():
        require(
            observed[name] == expected,
            f"{name} version {observed[name]!r} != {expected!r}",
        )
    conda_meta = Path(sys.prefix) / "conda-meta"
    for name, expected in CONDA_PACKAGE_VERSIONS.items():
        records = sorted(conda_meta.glob(f"{name}-*.json"))
        require(
            len(records) == 1,
            f"Expected one conda record for {name}: {records}",
        )
        with records[0].open(encoding="utf-8") as handle:
            version = str(json.load(handle).get("version"))
        require(
            version == expected, f"{name} version {version!r} != {expected!r}"
        )
        observed[name] = version

    version_probe = subprocess.run(
        ["xtb", "--version"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=True,
    ).stdout
    match = re.search(r"\b6\.7\.1\b", version_probe)
    require(
        match is not None,
        f"xTB executable did not report 6.7.1: {version_probe}",
    )
    return observed


def assert_linear_algebra_runtime() -> dict[str, str]:
    """Prove xTB and OpenBLAS resolve one ABI-compatible OpenMP runtime."""

    conda_meta = Path(sys.prefix) / "conda-meta"

    def record(name: str) -> dict:
        matches = sorted(conda_meta.glob(f"{name}-*.json"))
        require(len(matches) == 1, f"Expected one {name} record: {matches}")
        return load_receipt(matches[0])

    openblas = record("libopenblas")
    mutex = record("_openmp_mutex")
    require(
        openblas.get("build") == "openmp_hd680484_0",
        "The CPU image must use the locked OpenMP OpenBLAS build, observed "
        f"{openblas.get('build')!r}.",
    )
    require(
        mutex.get("build") == "7_kmp_llvm",
        "The locked OpenMP mutex must select one LLVM runtime, observed "
        f"{mutex.get('build')!r}.",
    )

    libdir = Path(sys.prefix) / "lib"
    implementations = sorted(libdir.glob("libopenblasp-r*.so"))
    require(
        len(implementations) == 1,
        f"Expected one OpenBLAS implementation: {implementations}",
    )
    implementation = implementations[0].resolve(strict=True)
    blas = (libdir / "libblas.so.3").resolve(strict=True)
    lapack = (libdir / "liblapack.so.3").resolve(strict=True)
    require(
        blas == implementation and lapack == implementation,
        "BLAS and LAPACK do not resolve to the locked OpenBLAS implementation: "
        f"blas={blas}, lapack={lapack}, openblas={implementation}",
    )

    xtb = Path(shutil.which("xtb") or "")
    require(xtb.is_file(), "xTB executable is unavailable for ABI inspection")

    def linked_openmp(path: Path) -> Path:
        completed = subprocess.run(
            ["ldd", str(path)],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        require(
            completed.returncode == 0 and "not found" not in completed.stdout,
            f"Unresolved shared library for {path}: {completed.stdout}",
        )
        match = re.search(
            r"^\s*libgomp\.so\.1\s+=>\s+(\S+)",
            completed.stdout,
            flags=re.MULTILINE,
        )
        require(match is not None, f"No OpenMP runtime resolved for {path}")
        return Path(match.group(1)).resolve(strict=True)

    xtb_openmp = linked_openmp(xtb)
    openblas_openmp = linked_openmp(implementation)
    selected_openmp = (libdir / "libgomp.so.1").resolve(strict=True)
    require(
        xtb_openmp == openblas_openmp == selected_openmp,
        "xTB and OpenBLAS loaded different OpenMP runtimes: "
        f"xtb={xtb_openmp}, openblas={openblas_openmp}, "
        f"selected={selected_openmp}",
    )
    require(
        selected_openmp.name == "libomp.so",
        f"OpenMP mutex did not resolve libgomp ABI to libomp: {selected_openmp}",
    )
    return {
        "openblas_build": str(openblas["build"]),
        "openmp_mutex_build": str(mutex["build"]),
        "blas_lapack_implementation": str(implementation),
        "shared_openmp_runtime": str(selected_openmp),
        "xtb_dynamic_linkage": "resolved",
    }


def import_runtime_dependencies() -> dict[str, str]:
    modules = {
        "ase": "ase",
        "geometric": "geometric",
        "h5py": "h5py",
        "openbabel": "openbabel.openbabel",
        "pymol": "pymol",
        "pyscf": "pyscf",
    }
    imported = {}
    for label, module_name in modules.items():
        module = importlib.import_module(module_name)
        imported[label] = getattr(module, "__version__", "importable")
    import h5py

    imported["hdf5"] = h5py.version.hdf5_version
    conventions = resolve_skill("scientific-conventions")
    require(
        conventions is not None and conventions.origin == "builtin",
        "The packaged scientific-conventions skill is unavailable.",
    )
    imported["scientific-conventions"] = conventions.skill_version
    return imported


def assert_runtime_boundary() -> dict[str, object]:
    require(os.geteuid() != 0, "The research runtime must not run as root.")
    stack_soft, stack_hard = resource.getrlimit(resource.RLIMIT_STACK)
    require(
        stack_soft == resource.RLIM_INFINITY,
        "xTB qualification requires an unlimited process stack; run the "
        "container with --ulimit stack=-1:-1.",
    )
    require(
        importlib.util.find_spec("gpu4pyscf") is None,
        "GPU4PySCF must not be present in the CPU image.",
    )
    unavailable = {
        name: shutil.which(name)
        for name in ("g16", "g09", "gaussian", "orca", "nvidia-smi")
    }
    require(
        all(value is None for value in unavailable.values()),
        f"Excluded engine or GPU executable found: {unavailable}",
    )
    return {
        "uid": os.geteuid(),
        "stack_limit": "unlimited",
        "gpu4pyscf_available": False,
        "gaussian_available": False,
        "orca_available": False,
        "nvidia_runtime_available": False,
    }


def prepare_directory(path: Path) -> Path:
    require(
        not path.exists() or not any(path.iterdir()),
        f"Qualification workspace must be new or empty: {path}",
    )
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_receipt(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    require(isinstance(value, dict), f"Receipt is not a mapping: {path}")
    return value


def one_file(folder: Path, pattern: str) -> Path:
    matches = sorted(folder.glob(pattern))
    require(
        len(matches) == 1, f"Expected one {pattern} in {folder}: {matches}"
    )
    return matches[0]


def one_recursive_file(folder: Path, pattern: str) -> Path:
    matches = sorted(folder.rglob(pattern))
    require(
        len(matches) == 1,
        f"Expected one recursive {pattern} below {folder}: {matches}",
    )
    return matches[0]


def finite_last(values, *, label: str) -> float:
    require(values is not None and len(values), f"Missing {label} values.")
    value = float(values[-1])
    require(math.isfinite(value), f"Non-finite {label}: {value}")
    return value


def pyscf_command(
    fixtures: Path,
    *,
    project: str,
    source: Path,
    label: str,
    jobtype: str,
) -> list[str]:
    return [
        "chemsmart",
        "run",
        "--server",
        "local",
        "--num-cores",
        str(QUALIFICATION_NUM_CORES),
        "--num-gpus",
        "0",
        "--mem-gb",
        "6",
        "--no-scratch",
        "pyscf",
        "--project",
        str(fixtures / project),
        "--filename",
        str(source),
        "--label",
        label,
        jobtype,
    ]


def xtb_command(
    fixtures: Path,
    *,
    source: Path,
    label: str,
    jobtype: str,
) -> list[str]:
    return [
        "chemsmart",
        "run",
        "--server",
        "local",
        "--num-cores",
        str(QUALIFICATION_NUM_CORES),
        "--num-gpus",
        "0",
        "--mem-gb",
        "6",
        "xtb",
        "--project",
        str(fixtures / "xtb-water.yaml"),
        "--filename",
        str(source),
        "--label",
        label,
        jobtype,
    ]


def inspect_pyscf_result(path: Path) -> tuple[dict, dict, dict, dict]:
    spec, provenance, status, results = read_pyscf_h5(path)
    require(
        status.get("normal_termination") is True, f"Incomplete PySCF: {path}"
    )
    require(
        status.get("engine_complete") is True,
        f"Incomplete PySCF stages: {path}",
    )
    require(int(spec.get("charge")) == 0, f"Wrong PySCF charge in {path}")
    require(
        int(spec.get("multiplicity")) == 1,
        f"Wrong PySCF multiplicity in {path}",
    )
    return spec, provenance, status, results


def qualify_pyscf(workspace: Path, fixtures: Path) -> dict[str, object]:
    water = fixtures / "water.xyz"
    cases = {
        "sp": workspace / "pyscf-sp",
        "d3bj": workspace / "pyscf-d3bj",
        "opt": workspace / "pyscf-opt",
        "hess": workspace / "pyscf-hess",
    }
    for folder in cases.values():
        folder.mkdir()

    gpu_refusal = workspace / "pyscf-gpu-refusal"
    gpu_refusal.mkdir()
    gpu_command = pyscf_command(
        fixtures,
        project="pyscf-water.yaml",
        source=water,
        label="pyscf_water_gpu_refusal",
        jobtype="sp",
    )
    gpu_command.insert(gpu_command.index("sp"), "--gpu")
    run_expected_failure(
        gpu_command,
        cwd=gpu_refusal,
        required_text="GPU",
    )
    require(
        not list(gpu_refusal.glob("*.h5")),
        "Rejected GPU request unexpectedly produced a PySCF result.",
    )

    run(
        pyscf_command(
            fixtures,
            project="pyscf-water.yaml",
            source=water,
            label="pyscf_water_sp",
            jobtype="sp",
        ),
        cwd=cases["sp"],
    )
    run(
        pyscf_command(
            fixtures,
            project="pyscf-water-d3bj.yaml",
            source=water,
            label="pyscf_water_d3bj",
            jobtype="sp",
        ),
        cwd=cases["d3bj"],
    )
    run(
        pyscf_command(
            fixtures,
            project="pyscf-water.yaml",
            source=water,
            label="pyscf_water_opt",
            jobtype="opt",
        ),
        cwd=cases["opt"],
    )

    opt_h5 = one_file(cases["opt"], "*.h5")
    run(
        pyscf_command(
            fixtures,
            project="pyscf-water.yaml",
            source=opt_h5,
            label="pyscf_water_hess",
            jobtype="hess",
        ),
        cwd=cases["hess"],
    )

    sp_spec, _, sp_status, sp_results = inspect_pyscf_result(
        one_file(cases["sp"], "*.h5")
    )
    d3_spec, _, _, d3_results = inspect_pyscf_result(
        one_file(cases["d3bj"], "*.h5")
    )
    _, _, opt_status, opt_results = inspect_pyscf_result(opt_h5)
    hess_h5 = one_file(cases["hess"], "*.h5")
    hess_spec, _, _, hess_results = inspect_pyscf_result(hess_h5)

    require(
        d3_spec.get("dispersion") == "d3bj",
        "The PySCF dispersion qualification did not apply D3BJ.",
    )
    require(
        "forces" not in sp_results
        and "forces" not in sp_status.get("properties", {}),
        "Force-free PySCF SP unexpectedly evaluated a nuclear gradient.",
    )

    require(
        opt_status["stages"]["opt"].get("optimizer_converged") is True,
        "PySCF geomeTRIC optimization did not converge.",
    )
    require(
        hess_spec.get("input_artifact_kind") == "pyscf_hdf5",
        "PySCF Hessian did not consume the optimized HDF5 artifact.",
    )
    require(
        np.allclose(
            np.asarray(hess_spec["positions"], dtype=float),
            np.asarray(opt_results["positions"], dtype=float),
            atol=1.0e-12,
            rtol=0.0,
        ),
        "PySCF optimized geometry was not handed to the Hessian job.",
    )
    hessian = np.asarray(hess_results["hessian"], dtype=float)
    require(
        hessian.shape == (3, 3, 3, 3),
        f"Unexpected Hessian shape: {hessian.shape}",
    )
    symmetry_error = float(
        np.max(np.abs(hessian - hessian.transpose(1, 0, 3, 2)))
    )
    require(
        symmetry_error < 1.0e-10, f"Asymmetric PySCF Hessian: {symmetry_error}"
    )
    frequencies = np.asarray(
        hess_results["vibrational_frequencies"], dtype=float
    ).reshape(-1)
    require(
        frequencies.size == 3,
        f"Expected 3 water modes, got {frequencies.size}",
    )
    require(np.isfinite(frequencies).all(), "Non-finite PySCF frequencies.")
    require(
        float(np.min(frequencies)) > -20.0,
        f"Consequential imaginary PySCF water mode: {frequencies.tolist()}",
    )

    hess_artifact = TrustedArtifactRefV1(
        artifact_id="qualification.pyscf.water-hessian",
        kind="pyscf_hdf5",
        sha256=file_sha256(hess_h5),
        size_bytes=hess_h5.stat().st_size,
        path=str(hess_h5),
        cli_value=str(hess_h5),
    )
    qh_gibbs = {}
    for cutoff in (50, 100, 200):
        thermo = derive_trusted_thermochemistry(
            artifact=hess_artifact,
            program="pyscf",
            temperature_k=298.15,
            pressure_atm=1.0,
            entropy_method="grimme",
            entropy_cutoff_cm1=cutoff,
            enthalpy_cutoff_cm1=cutoff,
        )
        require(
            thermo.entropy_cutoff_cm1 == cutoff
            and thermo.enthalpy_cutoff_cm1 == cutoff,
            "Typed quasi-harmonic receipt lost its requested cutoff",
        )
        quantity = next(
            item
            for item in thermo.quantities
            if item.quantity_id == "quasi_harmonic_gibbs_free_energy"
        )
        require(
            math.isfinite(float(quantity.value)),
            "Non-finite typed quasi-harmonic Gibbs energy",
        )
        qh_gibbs[str(cutoff)] = float(quantity.value)

    receipts = {}
    for name, folder in cases.items():
        receipt = load_receipt(one_file(folder, "*.receipt.json"))
        require(
            not receipt.get("findings"),
            f"PySCF {name} findings: {receipt.get('findings')}",
        )
        require(
            receipt.get("engine_complete") is True,
            f"PySCF {name} incomplete receipt",
        )
        receipts[name] = receipt.get("state")

    return {
        "method": f"{sp_spec.get('method')}/{sp_spec.get('basis')}",
        "sp_energy_hartree": finite_last(
            sp_results["energies"], label="PySCF SP energy"
        ),
        "d3bj_energy_hartree": finite_last(
            d3_results["energies"], label="PySCF D3BJ energy"
        ),
        "d3bj_applied": d3_spec.get("dispersion") == "d3bj",
        "optimized_energy_hartree": finite_last(
            opt_results["energies"], label="PySCF optimized energy"
        ),
        "frequencies_cm-1": [float(value) for value in frequencies],
        "hessian_shape": list(hessian.shape),
        "hessian_max_antisymmetry": symmetry_error,
        "optimized_hdf5_handoff": True,
        "typed_grimme_qh_gibbs_hartree": qh_gibbs,
        "typed_qh_cutoffs_cm-1": [50, 100, 200],
        "gpu_request_with_zero_devices_refused": True,
        "force_free_sp": True,
        "receipt_states": receipts,
    }


def inspect_xtb_result(
    folder: Path, *, jobtype: str
) -> tuple[XTBOutput, dict]:
    result = XTBOutput(str(folder))
    require(
        result.normal_termination is True,
        f"xTB {jobtype} did not terminate normally",
    )
    require(result.charge == 0, f"Wrong xTB {jobtype} charge: {result.charge}")
    require(
        result.multiplicity == 1,
        f"Wrong xTB {jobtype} multiplicity: {result.multiplicity}",
    )
    receipt = load_receipt(one_file(folder, "*.xtb-result-receipt.json"))
    environment = load_receipt(
        one_file(folder, "*.xtb-environment-receipt.json")
    )
    require(receipt.get("ready") is True, f"xTB {jobtype} receipt is red")
    require(
        not receipt.get("findings"),
        f"xTB {jobtype} findings: {receipt.get('findings')}",
    )
    require(
        receipt.get("engine_completion", {}).get("normal_termination") is True,
        f"xTB {jobtype} receipt lacks normal termination",
    )
    energy = receipt.get("results", {}).get("final_energy_hartree")
    require(
        energy is not None and math.isfinite(float(energy)),
        f"Missing xTB {jobtype} energy",
    )
    resources = environment.get("resources", {})
    for field in (
        "omp_num_threads",
        "mkl_num_threads",
        "openblas_num_threads",
    ):
        require(
            str(resources.get(field)) == str(QUALIFICATION_NUM_CORES),
            f"xTB {jobtype} {field} did not match requested cores: "
            f"{resources}",
        )
    stderr_paths = sorted(folder.glob("*.err"))
    require(
        len(stderr_paths) == 1,
        f"Expected one xTB stderr artifact in {folder}: {stderr_paths}",
    )
    with stderr_paths[0].open(encoding="utf-8", errors="replace") as handle:
        warning_count = sum(
            line.count(OPENBLAS_OPENMP_WARNING) for line in handle
        )
    require(
        warning_count == 0,
        f"xTB {jobtype} emitted {warning_count} OpenBLAS/OpenMP warnings",
    )
    return result, receipt


def qualify_xtb(workspace: Path, fixtures: Path) -> dict[str, object]:
    water = fixtures / "water.xyz"
    cases = {
        "sp": workspace / "xtb-sp",
        "opt": workspace / "xtb-opt",
        "hess": workspace / "xtb-hess",
    }
    for folder in cases.values():
        folder.mkdir()

    run(
        xtb_command(
            fixtures,
            source=water,
            label="xtb_water_sp",
            jobtype="sp",
        ),
        cwd=cases["sp"],
    )
    run(
        xtb_command(
            fixtures,
            source=water,
            label="xtb_water_opt",
            jobtype="opt",
        ),
        cwd=cases["opt"],
    )
    opt_folder = one_recursive_file(
        cases["opt"], "*.xtb-result-receipt.json"
    ).parent
    opt_xyz = opt_folder / "xtbopt.xyz"
    require(opt_xyz.is_file(), "xTB optimization did not produce xtbopt.xyz")
    run(
        xtb_command(
            fixtures,
            source=opt_xyz,
            label="xtb_water_hess",
            jobtype="hess",
        ),
        cwd=cases["hess"],
    )

    sp_folder = one_recursive_file(
        cases["sp"], "*.xtb-result-receipt.json"
    ).parent
    sp, sp_receipt = inspect_xtb_result(sp_folder, jobtype="sp")
    opt, opt_receipt = inspect_xtb_result(opt_folder, jobtype="opt")
    hess_folder = one_recursive_file(
        cases["hess"], "*.xtb-result-receipt.json"
    ).parent
    hess, hess_receipt = inspect_xtb_result(hess_folder, jobtype="hess")
    require(
        opt.geometry_optimization_converged is True,
        "xTB optimization did not converge",
    )

    staged_hess_xyz = hess_folder / "xtb_water_hess.xyz"
    require(staged_hess_xyz.is_file(), "Missing staged xTB Hessian geometry")
    optimized_atoms = ase_read(opt_xyz)
    hess_input_atoms = ase_read(staged_hess_xyz)
    require(
        optimized_atoms.get_chemical_symbols()
        == hess_input_atoms.get_chemical_symbols(),
        "xTB Hessian atom order differs from optimized geometry.",
    )
    require(
        np.allclose(
            optimized_atoms.get_positions(),
            hess_input_atoms.get_positions(),
            atol=1.0e-10,
            rtol=0.0,
        ),
        "xTB optimized geometry was not handed to the Hessian job.",
    )
    frequencies = np.asarray(
        hess.vibrational_frequencies, dtype=float
    ).reshape(-1)
    require(
        frequencies.size == 3,
        f"Expected 3 xTB water modes, got {frequencies.size}",
    )
    require(np.isfinite(frequencies).all(), "Non-finite xTB frequencies.")
    require(
        float(np.min(frequencies)) > -20.0,
        f"Consequential imaginary xTB water mode: {frequencies.tolist()}",
    )
    require(
        hess.hessian_file is not None, "xTB Hessian artifact was not parsed"
    )

    return {
        "method": "GFN2-xTB",
        "sp_energy_hartree": float(
            sp_receipt["results"]["final_energy_hartree"]
        ),
        "optimized_energy_hartree": float(
            opt_receipt["results"]["final_energy_hartree"]
        ),
        "hessian_energy_hartree": float(
            hess_receipt["results"]["final_energy_hartree"]
        ),
        "frequencies_cm-1": [float(value) for value in frequencies],
        "optimized_geometry_handoff": True,
        "requested_threads": QUALIFICATION_NUM_CORES,
        "openblas_openmp_warning_count": 0,
        "receipt_states": {
            "sp": sp_receipt["validation_state"],
            "opt": opt_receipt["validation_state"],
            "hess": hess_receipt["validation_state"],
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workspace", type=Path, default=Path("/workspace/qualification")
    )
    parser.add_argument(
        "--fixtures", type=Path, default=Path("/opt/chemsmart/qualification")
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("/workspace/qualification-summary.json"),
    )
    args = parser.parse_args()

    workspace = prepare_directory(args.workspace.resolve())
    fixtures = args.fixtures.resolve()
    require(fixtures.is_dir(), f"Fixture directory not found: {fixtures}")
    args.summary.parent.mkdir(parents=True, exist_ok=True)

    summary = {
        "schema_version": "chemsmart.cpu-image-qualification.v1",
        "status": "qualified",
        "platform": "linux/amd64",
        "packages": package_versions(),
        "linear_algebra_runtime": assert_linear_algebra_runtime(),
        "imports": import_runtime_dependencies(),
        "runtime_boundary": assert_runtime_boundary(),
        "pyscf": qualify_pyscf(workspace, fixtures),
        "xtb": qualify_xtb(workspace, fixtures),
        "scientific_scope": (
            "Installed-engine qualification on neutral-singlet water; not a "
            "large benchmark or reproduced literature result."
        ),
    }
    args.summary.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"CPU image qualification failed: {exc}", file=sys.stderr)
        raise
