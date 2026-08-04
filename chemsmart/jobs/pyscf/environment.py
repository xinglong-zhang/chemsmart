"""Deterministic evidence for the exact PySCF compute interpreter.

The ChemSmart controller may run in a different environment from PySCF or
GPU4PySCF.  This module therefore probes the interpreter selected by
``PySCFExecutable`` in a short child process before a calculation starts.
It imports packages and basis metadata only; it never launches an SCF cycle.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path


ENVIRONMENT_RECEIPT_SCHEMA_VERSION = "chemsmart.pyscf-environment.v1"
_PROBE_MARKER = "CHEMSMART_PYSCF_ENVIRONMENT="


class PySCFEnvironmentError(RuntimeError):
    """Raised when exact compute-environment evidence is unavailable."""

    def __init__(self, message, *, receipt=None):
        super().__init__(message)
        self.receipt = receipt or {}


def sha256_file(path):
    """Return the SHA-256 digest of a file without loading it all at once."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_sha256(value):
    """Return a stable digest for JSON-compatible evidence."""
    body = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def write_json_receipt(path, payload):
    """Atomically write a private JSON receipt without following a target link."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_symlink():
        raise ValueError(f"Refusing symlink receipt destination: {destination}")
    descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
        text=True,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary, 0o600)
        if destination.is_symlink():
            raise ValueError(
                f"Refusing symlink receipt destination: {destination}"
            )
        os.replace(temporary, destination)
    finally:
        if temporary.exists():
            temporary.unlink()


def resolve_compute_interpreter(interpreter):
    """Resolve and validate the exact executable used by the job child."""
    raw = os.path.expanduser(str(interpreter or "").strip())
    candidate = raw if os.path.isabs(raw) else shutil.which(raw)
    if not candidate:
        raise PySCFEnvironmentError(
            f"PySCF compute interpreter cannot be resolved: {raw!r}"
        )
    path = Path(candidate).resolve()
    if not path.is_file() or not os.access(path, os.X_OK):
        raise PySCFEnvironmentError(
            f"PySCF compute interpreter is not executable: {path}"
        )
    return path


def probe_compute_environment(
    *, interpreter, job, env, receipt_path, timeout=60
):
    """Probe *interpreter*, persist the evidence, and return its receipt."""
    request = _probe_request(job)
    request_sha256 = canonical_sha256(request)
    try:
        executable = resolve_compute_interpreter(interpreter)
    except PySCFEnvironmentError as exc:
        return _finalize_receipt(
            {
                "schema_version": ENVIRONMENT_RECEIPT_SCHEMA_VERSION,
                "observed_at": datetime.now(timezone.utc).isoformat(),
                "interpreter_requested": str(interpreter),
                "request": request,
                "request_sha256": request_sha256,
                "status": "unavailable",
                "error": {
                    "type": type(exc).__name__,
                    "message": str(exc)[:500],
                },
            },
            receipt_path,
        )
    base = {
        "schema_version": ENVIRONMENT_RECEIPT_SCHEMA_VERSION,
        "observed_at": datetime.now(timezone.utc).isoformat(),
        "interpreter": str(executable),
        "interpreter_sha256": sha256_file(executable),
        "request": request,
        "request_sha256": request_sha256,
        "status": "unavailable",
    }
    try:
        completed = subprocess.run(
            [str(executable), "-c", _PROBE_SCRIPT],
            input=json.dumps(request, sort_keys=True),
            text=True,
            capture_output=True,
            env=dict(env),
            timeout=timeout,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        base["error"] = {
            "type": type(exc).__name__,
            "message": str(exc)[:500],
        }
        return _finalize_receipt(base, receipt_path)

    base["probe_returncode"] = int(completed.returncode)
    if completed.returncode != 0:
        base["error"] = {
            "type": "probe_process_failed",
            "message": completed.stderr[-500:],
        }
        return _finalize_receipt(base, receipt_path)

    payload = _extract_probe_payload(completed.stdout)
    if payload is None:
        base["error"] = {
            "type": "probe_output_invalid",
            "message": "exact interpreter did not emit the receipt marker",
        }
        return _finalize_receipt(base, receipt_path)

    base.update(payload)
    base["status"] = "available"
    base["required_stages"] = request["stages"]
    base["cuda_visible_devices"] = env.get("CUDA_VISIBLE_DEVICES")
    base["omp_num_threads"] = env.get("OMP_NUM_THREADS")
    return _finalize_receipt(base, receipt_path)


def environment_blockers(receipt, *, engine):
    """Return hard blockers not covered by scientific preflight rules."""
    blockers = []
    if receipt.get("status") != "available":
        blockers.append(
            {
                "rule_id": "pyscf.environment.probe_unavailable",
                "field": "environment.status",
                "expected": "available",
                "observed": receipt.get("status"),
                "evidence_ref": "environment:status",
            }
        )
        return blockers
    dependencies = receipt.get("dependencies") or {}
    expected_interpreter = receipt.get("interpreter")
    observed_interpreter = receipt.get("interpreter_observed")
    try:
        interpreter_matches = (
            Path(expected_interpreter).resolve()
            == Path(observed_interpreter).resolve()
        )
    except (OSError, TypeError, ValueError):
        interpreter_matches = False
    if not interpreter_matches:
        blockers.append(
            {
                "rule_id": "pyscf.environment.interpreter_mismatch",
                "field": "environment.interpreter_observed",
                "expected": expected_interpreter,
                "observed": observed_interpreter,
                "evidence_ref": "environment:interpreter_observed",
            }
        )
    for name in ("pyscf", "numpy", "h5py"):
        detail = dependencies.get(name) or {}
        if detail.get("available") is not True:
            blockers.append(
                {
                    "rule_id": "pyscf.environment.dependency_missing",
                    "field": f"environment.dependencies.{name}",
                    "expected": f"{name} importable in compute interpreter",
                    "observed": detail,
                    "evidence_ref": f"environment:dependencies/{name}",
                }
            )
    pyscf_detail = dependencies.get("pyscf") or {}
    if pyscf_detail.get("available") is True and pyscf_detail.get(
        "version"
    ) != "2.14.0":
        blockers.append(
            {
                "rule_id": "pyscf.environment.version_mismatch",
                "field": "environment.dependencies.pyscf.version",
                "expected": "2.14.0",
                "observed": pyscf_detail.get("version"),
                "evidence_ref": "environment:dependencies/pyscf/version",
            }
        )
    request = receipt.get("request") or {}
    basis = request.get("basis")
    basis_available = receipt.get("basis_available") or {}
    if basis and basis_available.get(str(basis)) is not True:
        blockers.append(
            {
                "rule_id": "pyscf.environment.basis_unavailable",
                "field": "settings.basis",
                "expected": "loadable for every requested element",
                "observed": basis,
                "evidence_ref": f"environment:basis_available/{basis}",
            }
        )
    ecp_elements = receipt.get("basis_ecp_required_elements") or []
    if ecp_elements:
        blockers.append(
            {
                "rule_id": "pyscf.environment.ecp_unmaterialized",
                "field": "settings.basis",
                "expected": (
                    "all-electron basis or a future explicit PySCF ECP field"
                ),
                "observed": {
                    "basis": basis,
                    "elements_requiring_ecp": ecp_elements,
                },
                "evidence_ref": "environment:basis_ecp_required_elements",
            }
        )
    aux_basis = request.get("aux_basis")
    aux_basis_available = receipt.get("aux_basis_available") or {}
    if (
        aux_basis
        and aux_basis_available.get(str(aux_basis)) is not True
    ):
        blockers.append(
            {
                "rule_id": "pyscf.environment.aux_basis_unavailable",
                "field": "settings.aux_basis",
                "expected": "loadable for every requested element",
                "observed": aux_basis,
                "evidence_ref": (
                    f"environment:aux_basis_available/{aux_basis}"
                ),
            }
        )
    functional = request.get("functional")
    functional_available = receipt.get("functional_available") or {}
    if (
        functional
        and functional_available.get(str(functional)) is not True
    ):
        blockers.append(
            {
                "rule_id": "pyscf.environment.functional_unavailable",
                "field": "settings.functional",
                "expected": "parseable by target PySCF/libxc",
                "observed": functional,
                "evidence_ref": (
                    f"environment:functional_available/{functional}"
                ),
            }
        )
    if str(engine).lower() == "gpu":
        blockers.extend(_gpu_environment_blockers(receipt))
    return blockers


def _gpu_environment_blockers(receipt):
    blockers = []

    def require(field, expected, predicate):
        observed = receipt.get(field)
        try:
            accepted = bool(predicate(observed))
        except (TypeError, ValueError):
            accepted = False
        if not accepted:
            blockers.append(
                {
                    "rule_id": "pyscf.gpu.runtime_evidence_incomplete",
                    "field": f"environment.{field}",
                    "expected": expected,
                    "observed": observed,
                    "evidence_ref": f"environment:{field}",
                }
            )

    dependency = (receipt.get("dependencies") or {}).get("gpu4pyscf") or {}
    if dependency.get("available") is not True:
        blockers.append(
            {
                "rule_id": "pyscf.gpu.dependency_missing",
                "field": "environment.dependencies.gpu4pyscf",
                "expected": "GPU4PySCF importable in compute interpreter",
                "observed": dependency,
                "evidence_ref": "environment:dependencies/gpu4pyscf",
            }
        )
    elif dependency.get("version") != "1.8.0":
        blockers.append(
            {
                "rule_id": "pyscf.gpu.version_mismatch",
                "field": "environment.dependencies.gpu4pyscf.version",
                "expected": "1.8.0",
                "observed": dependency.get("version"),
                "evidence_ref": "environment:dependencies/gpu4pyscf/version",
            }
        )

    distribution = receipt.get("gpu4pyscf_distribution")
    if not isinstance(distribution, dict):
        blockers.append(
            {
                "rule_id": "pyscf.gpu.distribution_unresolved",
                "field": "environment.gpu4pyscf_distribution",
                "expected": "gpu4pyscf 1.8.0 distribution",
                "observed": distribution,
                "evidence_ref": "environment:packages",
            }
        )
    elif distribution.get("version") != "1.8.0":
        blockers.append(
            {
                "rule_id": "pyscf.gpu.version_mismatch",
                "field": "environment.gpu4pyscf_distribution.version",
                "expected": "1.8.0",
                "observed": distribution.get("version"),
                "evidence_ref": "environment:gpu4pyscf_distribution",
            }
        )

    require(
        "device_count",
        "positive integer",
        lambda value: (
            not isinstance(value, bool)
            and isinstance(value, int)
            and value > 0
        ),
    )
    require("device_name", "non-empty GPU model", lambda value: bool(value))
    require("device_uuid", "non-empty GPU UUID", lambda value: bool(value))
    require(
        "cuda_driver_version",
        "positive CUDA driver version",
        lambda value: value > 0,
    )
    require(
        "cuda_runtime_version",
        "positive CUDA runtime version",
        lambda value: value > 0,
    )
    require(
        "cupy_distribution",
        "CuPy distribution name and version",
        lambda value: (
            isinstance(value, dict)
            and bool(value.get("name"))
            and bool(value.get("version"))
        ),
    )
    require(
        "cutensor_distribution",
        "cuTENSOR distribution name and version",
        lambda value: (
            isinstance(value, dict)
            and bool(value.get("name"))
            and bool(value.get("version"))
        ),
    )
    require(
        "cutensor_runtime",
        "positive cuTENSOR runtime",
        lambda value: value > 0,
    )
    require(
        "cutensor_compatible",
        "approved CuPy/cuTENSOR pair",
        lambda value: value is True,
    )

    runtime = receipt.get("cuda_runtime_version")
    cuda_major = int(runtime) // 1000 if isinstance(runtime, int) else None
    for field, pattern in (
        ("gpu4pyscf_distribution", r"cuda(\d{2})x"),
        ("cupy_distribution", r"cuda(\d{2})x"),
        ("cutensor_distribution", r"cu(\d{2})"),
    ):
        detail = receipt.get(field)
        name = detail.get("name", "") if isinstance(detail, dict) else ""
        match = re.search(pattern, name)
        suffix_major = int(match.group(1)) if match else None
        if cuda_major is None or suffix_major != cuda_major:
            blockers.append(
                {
                    "rule_id": "pyscf.gpu.cuda_suffix_mismatch",
                    "field": f"environment.{field}.name",
                    "expected": (
                        "distribution suffix matching the observed CUDA "
                        f"runtime major {cuda_major!r}"
                    ),
                    "observed": name,
                    "evidence_ref": f"environment:{field}",
                }
            )
    return blockers


def _probe_request(job):
    settings = job.settings
    molecule = job.molecule
    return {
        "basis": settings.basis,
        "aux_basis": settings.aux_basis if settings.density_fit else None,
        "functional": settings.xc,
        "dispersion": settings.dispersion,
        "solvent_id": settings.solvent_id,
        "opt_solver": settings.opt_solver,
        "response_method": getattr(settings, "response_method", None),
        "state_manifold": getattr(settings, "state_manifold", None),
        "nstates": getattr(settings, "nstates", None),
        "engine": settings.engine,
        "symbols": sorted(set(molecule.chemical_symbols)),
        "stages": list(job.stages),
    }


def _extract_probe_payload(stdout):
    for line in reversed(str(stdout or "").splitlines()):
        if not line.startswith(_PROBE_MARKER):
            continue
        try:
            value = json.loads(line[len(_PROBE_MARKER) :])
        except json.JSONDecodeError:
            return None
        return value if isinstance(value, dict) else None
    return None


def _finalize_receipt(receipt, path):
    payload = dict(receipt)
    payload.pop("receipt_sha256", None)
    payload["receipt_sha256"] = canonical_sha256(payload)
    write_json_receipt(path, payload)
    return payload


_PROBE_SCRIPT = r'''
import importlib
import importlib.metadata
import importlib.util
import json
import os
import platform
import sys

request = json.load(sys.stdin)


def distribution_versions():
    prefixes = (
        "pyscf",
        "numpy",
        "h5py",
        "gpu4pyscf",
        "cupy",
        "cutensor",
    )
    found = {}
    for distribution in importlib.metadata.distributions():
        name = str(distribution.metadata.get("Name") or "").lower()
        if name.startswith(prefixes):
            found[name] = distribution.version
    return dict(sorted(found.items()))


def import_detail(name):
    try:
        module = importlib.import_module(name)
        return {
            "available": True,
            "version": getattr(module, "__version__", None),
        }
    except Exception as exc:
        return {
            "available": False,
            "error_type": type(exc).__name__,
        }


def solver_detail(module_name, dependency_name):
    try:
        module = importlib.import_module(module_name)
        kernel = getattr(module, "kernel", None)
        dependency_present = importlib.util.find_spec(dependency_name) is not None
        if not callable(kernel) or not dependency_present:
            return {
                "callable": False,
                "dependency_present": bool(dependency_present),
                "call_probe": "not_reached",
            }
        try:
            # Enter the real adapter entry point with a sentinel that cannot
            # launch chemistry.  Each supported adapter must reject the
            # sentinel only after its optional dependency and callable path
            # have loaded.  This avoids the known false green from importing
            # pyscf.geomopt alone.
            kernel(object())
        except (
            AttributeError,
            NotImplementedError,
            RuntimeError,
            TypeError,
        ) as exc:
            return {
                "callable": True,
                "dependency_present": True,
                "call_probe": "entrypoint_reached",
                "probe_exception_type": type(exc).__name__,
            }
        except Exception as exc:
            return {
                "callable": False,
                "dependency_present": True,
                "call_probe": "unexpected_failure",
                "probe_exception_type": type(exc).__name__,
            }
        return {
            "callable": False,
            "dependency_present": True,
            "call_probe": "sentinel_unexpectedly_accepted",
        }
    except Exception as exc:
        return {
            "callable": False,
            "dependency_present": False,
            "call_probe": "import_failed",
            "probe_exception_type": type(exc).__name__,
        }


def basis_max_l(name, symbols):
    if not name:
        return None
    try:
        from pyscf import gto

        values = []
        for symbol in symbols:
            shells = gto.basis.load(name, symbol)
            values.extend(
                int(shell[0])
                for shell in shells
                if isinstance(shell, (list, tuple)) and shell
            )
        return max(values) if values else None
    except Exception:
        return None


def basis_ecp_required_elements(name, symbols):
    if not name:
        return []
    try:
        from pyscf import gto
    except Exception:
        return []

    required = []
    for symbol in symbols:
        try:
            ecp = gto.basis.load_ecp(name, symbol)
        except Exception:
            ecp = None
        if ecp:
            required.append(str(symbol))
    return sorted(set(required))


def functional_available(name):
    if not name:
        return None
    try:
        from pyscf.dft import libxc

        libxc.parse_xc(name)
        return True
    except Exception:
        return False


def dispersion_detail(method, literal):
    """Probe one exact dispersion literal without constructing a molecule.

    ``parse_disp`` resolves the requested parameterization and ``check_disp``
    verifies target-version support.  The deliberately uninitialized mean-field
    object supplies only the method label inspected by ``check_disp``; this
    performs no molecule construction, integral evaluation, or SCF work.
    """
    requested_method = str(method or "hf").strip().lower()
    base = {
        "schema_version": "chemsmart.pyscf-dispersion-conformance.v1",
        "requested_method": requested_method,
        "requested_literal": literal,
        "supported": False,
        "method_compatible": False,
    }
    try:
        from pyscf import dft
        from pyscf.scf import dispersion as dispersion_api

        if not isinstance(literal, str) or not literal:
            raise ValueError("dispersion literal must be a non-empty string")
        parsed_method, parsed_version, with_3body = dispersion_api.parse_disp(
            requested_method, literal
        )
        canonical_method, _nlc, _implicit_dispersion = (
            dispersion_api.parse_dft(requested_method)
        )
        canonical_method = dispersion_api.XC_MAP.get(
            canonical_method, canonical_method
        )
        if requested_method == "hf":
            method_probe = object()
        else:
            method_probe = object.__new__(dft.rks.RKS)
            method_probe.xc = requested_method
        supported = dispersion_api.check_disp(method_probe, literal) is True
        compatible = (
            str(parsed_method).strip().lower()
            == str(canonical_method).strip().lower()
        )
        base.update(
            {
                "parsed_method": parsed_method,
                "dispersion_version": parsed_version,
                "with_3body": bool(with_3body),
                "supported": supported,
                "method_compatible": compatible,
                "status": (
                    "supported" if supported and compatible else "incompatible"
                ),
            }
        )
    except Exception as exc:
        base.update(
            {
                "status": "invalid",
                "error_type": type(exc).__name__,
                "error_message": str(exc)[:300],
            }
        )
    return base


packages = distribution_versions()
dependencies = {
    name: import_detail(name) for name in ("pyscf", "numpy", "h5py")
}
try:
    import pyscf
    from pyscf.dft import libxc

    pyscf_runtime_version = str(pyscf.__version__)
    libxc_runtime_version = str(libxc.libxc_version())
except Exception:
    pyscf_runtime_version = None
    libxc_runtime_version = None
try:
    dispersion_module_available = (
        importlib.util.find_spec("pyscf.dispersion") is not None
    )
except (ImportError, ModuleNotFoundError, ValueError):
    dispersion_module_available = False
dependencies["pyscf_dispersion"] = {
    "available": bool(
        dispersion_module_available
        or "pyscf-dispersion" in packages
    )
}
dependencies["gpu4pyscf"] = import_detail("gpu4pyscf")
dependencies["cupy"] = import_detail("cupy")

orbital_basis_l = basis_max_l(
    request.get("basis"), request.get("symbols") or []
)
auxiliary_basis_l = basis_max_l(
    request.get("aux_basis"), request.get("symbols") or []
)
result = {
    "python_version": platform.python_version(),
    "interpreter_observed": sys.executable,
    "pyscf_version": pyscf_runtime_version,
    "libxc_version": libxc_runtime_version,
    "packages": packages,
    "dependencies": dependencies,
    "solver_callables": {
        "geometric": solver_detail(
            "pyscf.geomopt.geometric_solver", "geometric"
        ),
        "berny": solver_detail("pyscf.geomopt.berny_solver", "berny"),
        "ase": solver_detail("pyscf.geomopt.ase_solver", "ase"),
    },
    "basis_max_l": {
        str(request.get("basis")): orbital_basis_l
    },
    "basis_available": {
        str(request.get("basis")): orbital_basis_l is not None
    },
    "basis_ecp_required_elements": basis_ecp_required_elements(
        request.get("basis"), request.get("symbols") or []
    ),
    "aux_basis_max_l": {
        str(request.get("aux_basis")): auxiliary_basis_l
    } if request.get("aux_basis") else {},
    "aux_basis_available": {
        str(request.get("aux_basis")): auxiliary_basis_l is not None
    } if request.get("aux_basis") else {},
    "functional_available": {
        str(request.get("functional")): functional_available(
            request.get("functional")
        )
    } if request.get("functional") else {},
}

dispersion = request.get("dispersion")
if dispersion is not None:
    result["dispersion_conformance"] = dispersion_detail(
        request.get("functional") or "hf", dispersion
    )

try:
    from pyscf.solvent.smd import solvent_db

    solvent_id = str(request.get("solvent_id") or "").strip().lower()
    result["solvent_ids"] = [solvent_id] if solvent_id in solvent_db else []
except Exception:
    result["solvent_ids"] = []

functional = request.get("functional")
if functional:
    normal = str(functional).lower().replace("_", "-")
    markers = (
        "b2plyp", "b2gp-plyp", "b2gpplyp", "pwpb95", "pbe0-2",
        "wb97x-2", "qidh", "-dh", "revdsd", "dsd", "dhdft",
        "xyg3", "xygjos", "double-hybrid",
    )
    laplacian = None
    try:
        from pyscf.dft import libxc

        checker = getattr(libxc, "needs_laplacian", None)
        if callable(checker):
            laplacian = bool(checker(functional))
    except Exception:
        pass
    functional_evidence = {
        "double_hybrid": any(marker in normal for marker in markers),
    }
    if laplacian is not None:
        functional_evidence["laplacian_meta_gga"] = laplacian
    result["functional_metadata"] = {
        str(functional): functional_evidence
    }
else:
    result["functional_metadata"] = {}

gpu_distribution = None
for name in sorted(packages):
    if name == "gpu4pyscf" or name.startswith("gpu4pyscf-cuda"):
        gpu_distribution = {"name": name, "version": packages[name]}
        break
result["gpu4pyscf_distribution"] = gpu_distribution

cupy_distribution = None
cutensor_distribution = None
for name in sorted(packages):
    if cupy_distribution is None and (
        name == "cupy" or name.startswith("cupy-cuda")
    ):
        cupy_distribution = {"name": name, "version": packages[name]}
    if cutensor_distribution is None and name.startswith("cutensor"):
        cutensor_distribution = {"name": name, "version": packages[name]}
result["cupy_distribution"] = cupy_distribution
result["cutensor_distribution"] = cutensor_distribution
recommended_pairs = {
    ("13.3.0", "2.0.2"),
    ("13.4.1", "2.2.0"),
}
observed_pair = (
    cupy_distribution.get("version") if cupy_distribution else None,
    cutensor_distribution.get("version") if cutensor_distribution else None,
)
result["cupy_cutensor_pair"] = list(observed_pair)
result["cutensor_compatible"] = observed_pair in recommended_pairs

result["cuda_available"] = 0
result["device_count"] = 0
try:
    import cupy

    count = int(cupy.cuda.runtime.getDeviceCount())
    result["cuda_available"] = count
    result["device_count"] = count
    result["cuda_driver_version"] = int(
        cupy.cuda.runtime.driverGetVersion()
    )
    result["cuda_runtime_version"] = int(
        cupy.cuda.runtime.runtimeGetVersion()
    )
    if count:
        properties = cupy.cuda.runtime.getDeviceProperties(0)
        name = properties.get("name") if isinstance(properties, dict) else None
        if isinstance(name, bytes):
            name = name.decode("utf-8", errors="replace")
        result["device_name"] = name
        uuid = properties.get("uuid") if isinstance(properties, dict) else None
        if isinstance(uuid, bytes):
            uuid = uuid.hex()
        elif isinstance(uuid, (list, tuple)):
            uuid = "".join("%02x" % int(item) for item in uuid)
        result["device_uuid"] = uuid
    try:
        from cupy_backends.cuda.libs import cutensor

        version = int(cutensor.get_version())
        result["cutensor_runtime"] = version
        result["cutensor_runtime_available"] = version > 0
    except Exception:
        result["cutensor_runtime_available"] = False
except Exception as exc:
    result["gpu_probe_error_type"] = type(exc).__name__

print("CHEMSMART_PYSCF_ENVIRONMENT=" + json.dumps(result, sort_keys=True))
'''


__all__ = [
    "ENVIRONMENT_RECEIPT_SCHEMA_VERSION",
    "PySCFEnvironmentError",
    "canonical_sha256",
    "environment_blockers",
    "probe_compute_environment",
    "resolve_compute_interpreter",
    "sha256_file",
    "write_json_receipt",
]
