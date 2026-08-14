"""Deterministic environment, preview, and result receipts for xTB 6.7.1."""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
from collections.abc import Mapping
from pathlib import Path

import numpy as np


XTB_REQUIRED_VERSION = "6.7.1"
XTB_ENVIRONMENT_SCHEMA_VERSION = "chemsmart.xtb-environment.v1"
XTB_PREVIEW_SCHEMA_VERSION = "chemsmart.xtb-preview.v1"
XTB_RESULT_SCHEMA_VERSION = "chemsmart.xtb-result-validation.v1"
_MISSING = object()
_OPENBLAS_OPENMP_WARNING = (
    "OpenBLAS Warning: Detect OpenMP Loop and this application may hang."
)


class XTBEnvironmentError(RuntimeError):
    """Raised when the exact xTB executable/resource contract is red."""

    def __init__(self, findings, *, receipt_path=None):
        self.findings = tuple(findings)
        self.receipt_path = receipt_path
        super().__init__(
            "xTB environment validation failed: "
            + ", ".join(item["rule_id"] for item in self.findings)
        )


class XTBResultValidationError(RuntimeError):
    """Raised when an executed xTB result cannot be accepted as complete."""

    def __init__(self, findings, *, receipt_path=None, returncode=None):
        self.findings = tuple(findings)
        self.receipt_path = receipt_path
        self.returncode = returncode
        super().__init__(
            "xTB result validation failed: "
            + ", ".join(item["rule_id"] for item in self.findings)
        )


def canonical_sha256(value):
    """Return the SHA-256 of a canonical JSON representation."""

    body = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def sha256_file(path):
    """Hash *path* without loading the whole artifact into memory."""

    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json_receipt(path, payload):
    """Atomically persist a private receipt without following target links."""

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


def finalize_receipt(path, payload):
    """Bind a receipt to its canonical content and write it."""

    body = dict(payload)
    body.pop("receipt_sha256", None)
    body["receipt_sha256"] = canonical_sha256(body)
    write_json_receipt(path, body)
    return body


def _load_bound_receipt(path):
    try:
        with open(path, encoding="utf-8") as handle:
            payload = json.load(handle)
        observed = payload.get("receipt_sha256")
        body = dict(payload)
        body.pop("receipt_sha256", None)
        if observed and observed == canonical_sha256(body):
            return payload
    except (OSError, TypeError, ValueError, json.JSONDecodeError):
        pass
    return None


def _finding(rule_id, field, expected, observed, evidence_ref):
    return {
        "rule_id": rule_id,
        "field": field,
        "expected": expected,
        "observed": observed,
        "evidence_ref": evidence_ref,
    }


def _artifact_record(path):
    path = Path(path)
    if not path.is_file():
        return None
    return {
        "path": path.name,
        "size": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def _openblas_openmp_warning_count(path):
    """Count the xTB/OpenBLAS incompatibility warning without loading stderr."""

    try:
        with open(path, encoding="utf-8", errors="replace") as handle:
            return sum(line.count(_OPENBLAS_OPENMP_WARNING) for line in handle)
    except OSError:
        return 0


def _bound_artifact_record(path, *, role):
    """Bind an explicitly declared artifact to its exact path and bytes."""

    if path is None:
        return None
    declared = Path(os.path.abspath(os.fspath(path)))
    if not declared.is_file():
        raise FileNotFoundError(
            f"Cannot bind missing xTB {role} artifact: {declared}"
        )
    resolved = declared.resolve(strict=True)
    stat = resolved.stat()
    return {
        "role": role,
        "declared_path": str(declared),
        "resolved_path": str(resolved),
        "size": stat.st_size,
        "sha256": sha256_file(resolved),
    }


def bind_xtb_declared_artifacts(job):
    """Capture source and project provenance before rendering or execution."""

    return {
        "source_artifact": _bound_artifact_record(
            getattr(job, "source_filename", None), role="source_geometry"
        ),
        "project_artifact": _bound_artifact_record(
            getattr(job, "project_source_file", None), role="project_yaml"
        ),
    }


def bind_xtb_execution_input(path):
    """Bind the exact rendered XYZ path passed to the xTB child process."""

    return _bound_artifact_record(path, role="execution_input_geometry")


def _verify_bound_artifact(record, *, findings, evidence_ref):
    """Compare a live artifact with a previously captured binding."""

    if record is None:
        return
    declared = Path(record["declared_path"])
    observed = None
    if declared.is_file():
        try:
            resolved = declared.resolve(strict=True)
            observed = {
                "resolved_path": str(resolved),
                "size": resolved.stat().st_size,
                "sha256": sha256_file(resolved),
            }
        except OSError:
            observed = None
    expected = {
        "resolved_path": record["resolved_path"],
        "size": record["size"],
        "sha256": record["sha256"],
    }
    if observed != expected:
        findings.append(
            _finding(
                "xtb.provenance.artifact_mismatch",
                record["role"],
                expected,
                observed,
                evidence_ref,
            )
        )


def verify_xtb_provenance_binding(binding):
    """Return deterministic findings for mutation or loss of bound inputs."""

    findings = []
    for key in (
        "source_artifact",
        "project_artifact",
        "execution_input_artifact",
    ):
        _verify_bound_artifact(
            (binding or {}).get(key),
            findings=findings,
            evidence_ref=f"binding:{key}",
        )
    return findings


def _audit_xtb_execution_input_binding(record, *, receipt_parent, artifacts):
    """Audit a rendered input after an optional validated scratch cleanup.

    The child argv and provenance receipt remain bound to the original input
    path.  While that path exists, it is therefore the only acceptable audit
    source.  If cleanup removed it, the runner's durable copy may stand in for
    the bytes only when both the artifact manifest and the live copy retain
    the exact size and digest captured by the original binding.
    """

    findings = []
    if record is None:
        return findings
    required_fields = {
        "role",
        "declared_path",
        "resolved_path",
        "size",
        "sha256",
    }
    if not required_fields.issubset(record):
        findings.append(
            _finding(
                "xtb.provenance.artifact_mismatch",
                "execution_input_geometry",
                sorted(required_fields),
                sorted(record),
                "binding:execution_input_artifact",
            )
        )
        return findings

    declared = Path(record["declared_path"])
    # A present file, directory, or symlink is never replaced by a fallback.
    # This keeps mutation/substitution detection strict even when an exact
    # durable copy also happens to be available beside the receipt.
    if declared.exists() or declared.is_symlink():
        _verify_bound_artifact(
            record,
            findings=findings,
            evidence_ref="binding:execution_input_artifact",
        )
        return findings

    expected = {
        "size": record["size"],
        "sha256": record["sha256"],
    }
    basename = declared.name
    manifest_record = (
        artifacts.get(basename) if isinstance(artifacts, Mapping) else None
    )
    observed_manifest = (
        {
            "size": manifest_record.get("size"),
            "sha256": manifest_record.get("sha256"),
        }
        if isinstance(manifest_record, Mapping)
        else None
    )
    durable_path = Path(receipt_parent) / basename
    durable_record = (
        None if durable_path.is_symlink() else _artifact_record(durable_path)
    )
    observed_durable = (
        {
            "size": durable_record.get("size"),
            "sha256": durable_record.get("sha256"),
        }
        if isinstance(durable_record, Mapping)
        else None
    )
    if observed_manifest == expected and observed_durable == expected:
        return findings

    findings.append(
        _finding(
            "xtb.provenance.artifact_mismatch",
            record["role"],
            {
                "original_binding": expected,
                "artifact_manifest": expected,
                "durable_copy": expected,
            },
            {
                "original_path": None,
                "artifact_manifest": observed_manifest,
                "durable_copy": observed_durable,
            },
            f"artifact:{basename}",
        )
    )
    return findings


def _settings_payload(settings):
    return {
        name: getattr(settings, name)
        for name in sorted(settings.FIELDS)
    }


def _molecule_payload(molecule):
    return {
        "symbols": list(molecule.symbols),
        "positions_angstrom": np.asarray(
            molecule.positions, dtype=float
        ).tolist(),
        "charge": molecule.charge,
        "multiplicity": molecule.multiplicity,
    }


def probe_xtb_environment(
    *,
    executable,
    num_cores,
    num_gpus,
    mem_gb,
    env,
    receipt_path,
    settings=None,
    timeout=15,
):
    """Probe the exact binary and localize scientific preflight findings."""

    environment_findings = []
    scientific_findings = []
    requested = os.path.expanduser(str(executable or "").strip())
    candidate = requested if os.path.isabs(requested) else shutil.which(requested)
    resolved = None
    if candidate:
        try:
            resolved = Path(candidate).resolve(strict=True)
        except OSError:
            resolved = None
    if resolved is None or not resolved.is_file() or not os.access(resolved, os.X_OK):
        environment_findings.append(
            _finding(
                "xtb.environment.executable_unavailable",
                "environment.executable",
                "an executable file",
                requested or None,
                "environment:executable",
            )
        )

    if isinstance(num_cores, bool) or not isinstance(num_cores, int) or num_cores < 1:
        environment_findings.append(
            _finding(
                "xtb.resources.cores_invalid",
                "resources.num_cores",
                "positive integer",
                num_cores,
                "environment:resources/num_cores",
            )
        )
    if num_gpus != 0:
        environment_findings.append(
            _finding(
                "xtb.resources.gpu_unsupported",
                "resources.num_gpus",
                0,
                num_gpus,
                "environment:resources/num_gpus",
            )
        )

    if (
        isinstance(mem_gb, bool)
        or not isinstance(mem_gb, (int, float))
        or not math.isfinite(mem_gb)
        or mem_gb <= 0
    ):
        environment_findings.append(
            _finding(
                "xtb.resources.memory_invalid",
                "resources.mem_gb",
                "positive finite number",
                mem_gb,
                "environment:resources/mem_gb",
            )
        )

    omp_threads = env.get("OMP_NUM_THREADS")
    if str(omp_threads) != str(num_cores):
        environment_findings.append(
            _finding(
                "xtb.resources.omp_mismatch",
                "environment.OMP_NUM_THREADS",
                str(num_cores),
                omp_threads,
                "environment:variables/OMP_NUM_THREADS",
            )
        )

    payload = {
        "schema_version": XTB_ENVIRONMENT_SCHEMA_VERSION,
        "required_version": XTB_REQUIRED_VERSION,
        "status": "unavailable",
        "execution_ready": False,
        "resources": {
            "engine": "cpu",
            "num_cores": num_cores,
            "num_gpus": num_gpus,
            "mem_gb": mem_gb,
            "omp_num_threads": omp_threads,
            "mkl_num_threads": env.get("MKL_NUM_THREADS"),
            "openblas_num_threads": env.get("OPENBLAS_NUM_THREADS"),
        },
        "environment_findings": environment_findings,
        "scientific_findings": scientific_findings,
        "findings": [],
    }
    if resolved is not None:
        payload["executable"] = str(resolved)
        payload["executable_sha256"] = sha256_file(resolved)

    if resolved is not None and not environment_findings:
        try:
            completed = subprocess.run(
                [str(resolved), "--version"],
                capture_output=True,
                text=True,
                env=dict(env),
                timeout=timeout,
                check=False,
            )
            payload["probe_returncode"] = int(completed.returncode)
            version_text = f"{completed.stdout}\n{completed.stderr}"
            match = re.search(
                r"(?:xtb\s+version|version)\s+v?(\d+\.\d+\.\d+)",
                version_text,
                flags=re.IGNORECASE,
            )
            observed_version = match.group(1) if match else None
            payload["observed_version"] = observed_version
            if completed.returncode != 0:
                environment_findings.append(
                    _finding(
                        "xtb.environment.version_probe_failed",
                        "environment.probe_returncode",
                        0,
                        completed.returncode,
                        "environment:version_probe",
                    )
                )
            if observed_version != XTB_REQUIRED_VERSION:
                environment_findings.append(
                    _finding(
                        "xtb.environment.version_mismatch",
                        "environment.version",
                        XTB_REQUIRED_VERSION,
                        observed_version,
                        "environment:version",
                    )
                )
        except (OSError, subprocess.TimeoutExpired) as exc:
            environment_findings.append(
                _finding(
                    "xtb.environment.version_probe_failed",
                    "environment.version_probe",
                    "successful bounded argv probe",
                    type(exc).__name__,
                    "environment:version_probe",
                )
            )

    solvent_capability = None
    if settings is not None:
        solvent_capability = settings.solvent_capability().as_dict()
        if solvent_capability["status"] == "unknown":
            scientific_findings.append(
                _finding(
                    "xtb.solvent.capability_unknown",
                    "settings.solvent",
                    "confirmed support for the pinned xTB environment",
                    {
                        "method": solvent_capability["method"],
                        "model": solvent_capability["model"],
                        "identifier": solvent_capability["identifier"],
                    },
                    "settings:solvent_capability",
                )
            )
        elif solvent_capability["status"] == "confirmed_unsupported":
            scientific_findings.append(
                _finding(
                    "xtb.solvent.capability_unsupported",
                    "settings.solvent",
                    "a supported solvent combination",
                    solvent_capability,
                    "settings:solvent_capability",
                )
            )

    findings = [*environment_findings, *scientific_findings]
    payload["solvent_capability"] = solvent_capability
    payload["environment_findings"] = environment_findings
    payload["scientific_findings"] = scientific_findings
    payload["findings"] = findings
    payload["status"] = (
        "available" if not environment_findings else "unavailable"
    )
    payload["preflight_state"] = (
        "ready"
        if not findings
        else (
            "needs_clarification"
            if not environment_findings and scientific_findings
            else "blocked"
        )
    )
    payload["execution_ready"] = not findings
    return finalize_receipt(receipt_path, payload)


def build_preview_receipt(*, job, command, preview_dir, receipt_path):
    """Create an explicit non-executed preview receipt."""

    input_record = _artifact_record(Path(preview_dir) / f"{job.label}.xyz")
    solvent_capability = job.settings.solvent_capability().as_dict()
    payload = {
        "schema_version": XTB_PREVIEW_SCHEMA_VERSION,
        "state": "previewed",
        "executed": False,
        "execution_ready": False,
        "program": "xtb",
        "required_version": XTB_REQUIRED_VERSION,
        "jobtype": job._TYPE_TO_JOBTYPE[job.TYPE],
        "canonical_argv": list(command),
        "command_sha256": canonical_sha256(list(command)),
        "settings": _settings_payload(job.settings),
        "settings_sha256": canonical_sha256(_settings_payload(job.settings)),
        "solvent_capability": solvent_capability,
        "clarification_required": solvent_capability["status"] == "unknown",
        "input_artifact": input_record,
        "preview_directory": os.path.relpath(preview_dir, job.folder),
        "finding": {
            "rule_id": "xtb.preview.not_executed",
            "message": "No xTB process or chemistry engine was invoked.",
        },
    }
    return finalize_receipt(receipt_path, payload)


def _route_solvent_id(route):
    model = route.solvent_model
    if model is None:
        return None
    flags = (f"--{model}", "-g") if model == "gbsa" else (f"--{model}",)
    for flag in flags:
        if flag in route.route_inputs:
            index = route.route_inputs.index(flag) + 1
            if index < len(route.route_inputs):
                return route.route_inputs[index]
    return None


def _route_settings(
    route, *, method=_MISSING, charge=_MISSING, multiplicity=_MISSING
):
    return {
        "jobtype": route.jobtype,
        "gfn_version": (
            route.gfn_version if method is _MISSING else method
        ),
        "charge": route.charge if charge is _MISSING else charge,
        "multiplicity": (
            multiplicity
            if multiplicity is not _MISSING
            else (None if route.uhf is None else route.uhf + 1)
        ),
        "optimization_level": route.optimization_level,
        "solvent_model": route.solvent_model,
        "solvent_id": _route_solvent_id(route),
        "grad": route.grad,
    }


def _expected_mode_count(molecule):
    count = molecule.num_atoms
    if count == 1:
        return 0
    coordinates = np.asarray(molecule.positions, dtype=float)
    centered = coordinates - coordinates.mean(axis=0)
    singular_values = np.linalg.svd(centered, compute_uv=False)
    scale = max(float(singular_values[0]), 1.0)
    linear = bool(
        len(singular_values) < 2
        or float(singular_values[1]) <= 1.0e-8 * scale
    )
    return 3 * count - (5 if linear else 6)


def _compare(findings, *, rule_id, field, expected, observed, evidence_ref):
    if observed != expected:
        findings.append(
            _finding(rule_id, field, expected, observed, evidence_ref)
        )


def _safe_value(obj, attribute, findings, *, evidence_ref):
    try:
        return getattr(obj, attribute)
    except Exception as exc:
        findings.append(
            _finding(
                "xtb.result.semantic_observation_failed",
                f"result.{attribute}",
                "deterministically parseable value",
                type(exc).__name__,
                evidence_ref,
            )
        )
        return None


def validate_xtb_result(
    *,
    job,
    command,
    environment_receipt,
    provenance_binding,
    returncode,
    receipt_path,
):
    """Validate exact output semantics and persist the completion receipt."""

    from chemsmart.io.xtb.output import XTBOutput
    from chemsmart.io.xtb.route import XTBRoute

    findings = verify_xtb_provenance_binding(provenance_binding)
    warnings = []
    requested = _settings_payload(job.settings)
    command_route = XTBRoute(" ".join(command))
    command_settings = _route_settings(command_route)
    _compare(
        findings,
        rule_id="xtb.command.input_mismatch",
        field="command.input_geometry",
        expected=(
            None
            if not provenance_binding.get("execution_input_artifact")
            else provenance_binding["execution_input_artifact"][
                "declared_path"
            ]
        ),
        observed=(
            None if len(command) < 2 else os.path.abspath(command[1])
        ),
        evidence_ref="command:argv/1",
    )
    for field in ("jobtype", "gfn_version", "charge", "multiplicity"):
        _compare(
            findings,
            rule_id="xtb.command.settings_mismatch",
            field=f"command.{field}",
            expected=requested[field],
            observed=command_settings[field],
            evidence_ref=f"command:settings/{field}",
        )
    if requested["jobtype"] == "opt":
        _compare(
            findings,
            rule_id="xtb.command.settings_mismatch",
            field="command.optimization_level",
            expected=requested["optimization_level"],
            observed=command_settings["optimization_level"],
            evidence_ref="command:settings/optimization_level",
        )
    for field in ("solvent_model", "solvent_id", "grad"):
        _compare(
            findings,
            rule_id="xtb.command.settings_mismatch",
            field=f"command.{field}",
            expected=requested[field],
            observed=command_settings[field],
            evidence_ref=f"command:settings/{field}",
        )
    execution_input = provenance_binding.get("execution_input_artifact")
    input_record = (
        None
        if execution_input is None
        else _artifact_record(execution_input["declared_path"])
    )
    if input_record is None or input_record.get("size", 0) == 0:
        findings.append(
            _finding(
                "xtb.result.input_missing",
                "result.input_geometry",
                "non-empty rendered XYZ artifact",
                None,
                "artifact:input_geometry",
            )
        )
    stderr_path = Path(job.folder) / f"{job.label}.err"
    openblas_warning_count = _openblas_openmp_warning_count(stderr_path)
    if openblas_warning_count:
        warnings.append(
            _finding(
                "xtb.environment.openblas_openmp_warning",
                "stderr.openblas_openmp_warning_count",
                0,
                openblas_warning_count,
                f"artifact:{stderr_path.name}",
            )
        )

    output = None
    main = None
    observed_settings = None
    normal_termination = None
    final_energy = None
    try:
        output = XTBOutput(folder=job.folder)
        main = output.main_out
    except Exception as exc:
        findings.append(
            _finding(
                "xtb.result.output_unreadable",
                "result.output",
                "one readable xTB main output",
                type(exc).__name__,
                "artifact:main_output",
            )
        )

    _compare(
        findings,
        rule_id="xtb.result.process_failed",
        field="process.returncode",
        expected=0,
        observed=returncode,
        evidence_ref="process:returncode",
    )
    if not environment_receipt or not environment_receipt.get("execution_ready"):
        findings.append(
            _finding(
                "xtb.result.environment_unbound",
                "environment.execution_ready",
                True,
                bool(
                    environment_receipt
                    and environment_receipt.get("execution_ready")
                ),
                "receipt:environment",
            )
        )

    if main is None:
        unreadable = any(
            item["rule_id"] == "xtb.result.output_unreadable"
            for item in findings
        )
        if not unreadable:
            findings.append(
                _finding(
                    "xtb.result.output_missing",
                    "result.output",
                    "xTB main output",
                    None,
                    "artifact:main_output",
                )
            )
    else:
        route = _safe_value(
            main,
            "route_object",
            findings,
            evidence_ref="output:route",
        )
        input_geometry = _safe_value(
            output,
            "input_geometry",
            findings,
            evidence_ref="artifact:input_geometry",
        )
        expected_positions = np.asarray(job.molecule.positions, dtype=float)
        observed_positions = (
            None
            if input_geometry is None
            else np.asarray(input_geometry.positions, dtype=float)
        )
        positions_match = bool(
            observed_positions is not None
            and observed_positions.shape == expected_positions.shape
            and np.allclose(
                observed_positions,
                expected_positions,
                rtol=0.0,
                atol=1.0e-10,
            )
        )
        if not positions_match:
            findings.append(
                _finding(
                    "xtb.result.input_geometry_mismatch",
                    "result.input_positions_angstrom",
                    canonical_sha256(expected_positions.tolist()),
                    (
                        None
                        if observed_positions is None
                        else canonical_sha256(observed_positions.tolist())
                    ),
                    "artifact:input_geometry",
                )
            )
        _compare(
            findings,
            rule_id="xtb.result.input_identity_mismatch",
            field="result.input_symbols",
            expected=list(job.molecule.symbols),
            observed=(
                None
                if input_geometry is None
                else list(input_geometry.symbols)
            ),
            evidence_ref="artifact:input_geometry",
        )
        _compare(
            findings,
            rule_id="xtb.result.version_mismatch",
            field="result.version",
            expected=XTB_REQUIRED_VERSION,
            observed=_safe_value(
                main,
                "version",
                findings,
                evidence_ref="output:version",
            ),
            evidence_ref="output:version",
        )
        normal_termination_observation = _safe_value(
            main,
            "normal_termination",
            findings,
            evidence_ref="output:normal_termination",
        )
        normal_termination = (
            None
            if normal_termination_observation is None
            else bool(normal_termination_observation)
        )
        _compare(
            findings,
            rule_id="xtb.result.termination_missing",
            field="result.normal_termination",
            expected=True,
            observed=normal_termination is True,
            evidence_ref="output:normal_termination",
        )
        unpaired = _safe_value(
            main,
            "unpaired_electrons",
            findings,
            evidence_ref="output:settings/multiplicity",
        )
        if route is None:
            observed_settings = {
                field: None
                for field in (
                    "jobtype",
                    "gfn_version",
                    "charge",
                    "multiplicity",
                    "optimization_level",
                    "solvent_model",
                    "solvent_id",
                    "grad",
                )
            }
        else:
            observed_settings = _route_settings(
                route,
                method=_safe_value(
                    main,
                    "method",
                    findings,
                    evidence_ref="output:settings/gfn_version",
                ),
                charge=_safe_value(
                    main,
                    "net_charge",
                    findings,
                    evidence_ref="output:settings/charge",
                ),
                multiplicity=None if unpaired is None else unpaired + 1,
            )
        for field in ("jobtype", "gfn_version", "charge", "multiplicity"):
            _compare(
                findings,
                rule_id="xtb.result.settings_mismatch",
                field=f"settings.{field}",
                expected=requested[field],
                observed=observed_settings[field],
                evidence_ref=f"output:settings/{field}",
            )
        if requested["jobtype"] == "opt":
            _compare(
                findings,
                rule_id="xtb.result.settings_mismatch",
                field="settings.optimization_level",
                expected=requested["optimization_level"],
                observed=observed_settings["optimization_level"],
                evidence_ref="output:settings/optimization_level",
            )
        for field in ("solvent_model", "solvent_id", "grad"):
            _compare(
                findings,
                rule_id="xtb.result.settings_mismatch",
                field=f"settings.{field}",
                expected=requested[field],
                observed=observed_settings[field],
                evidence_ref=f"output:settings/{field}",
            )
        final_energy = _safe_value(
            output,
            "final_energy",
            findings,
            evidence_ref="output:final_energy",
        )
        if final_energy is None or not math.isfinite(final_energy):
            findings.append(
                _finding(
                    "xtb.result.energy_missing",
                    "result.final_energy_hartree",
                    "finite value",
                    final_energy,
                    "output:final_energy",
                )
            )

        if requested["jobtype"] == "opt":
            optimization_converged = _safe_value(
                output,
                "geometry_optimization_converged",
                findings,
                evidence_ref="output:optimization",
            )
            if not optimization_converged:
                findings.append(
                    _finding(
                        "xtb.result.optimization_not_converged",
                        "result.geometry_optimization_converged",
                        True,
                        False,
                        "output:optimization",
                    )
                )
            optimized = _safe_value(
                output,
                "xtbopt_geometry",
                findings,
                evidence_ref="artifact:optimized_geometry",
            )
            observed_symbols = None if optimized is None else list(optimized.symbols)
            _compare(
                findings,
                rule_id="xtb.result.optimized_geometry_missing",
                field="result.optimized_symbols",
                expected=list(job.molecule.symbols),
                observed=observed_symbols,
                evidence_ref="artifact:optimized_geometry",
            )

        if requested["jobtype"] == "hess":
            hessian = _artifact_record(Path(job.folder) / "hessian")
            if hessian is None or hessian.get("size", 0) == 0:
                findings.append(
                    _finding(
                        "xtb.result.hessian_missing",
                        "result.hessian",
                        "non-empty hessian artifact",
                        None,
                        "artifact:hessian",
                    )
                )
            frequencies = _safe_value(
                output,
                "vibrational_frequencies",
                findings,
                evidence_ref="output:vibrational_frequencies",
            )
            expected_modes = _expected_mode_count(job.molecule)
            if frequencies is None or len(frequencies) != expected_modes:
                findings.append(
                    _finding(
                        "xtb.result.frequency_mode_count",
                        "result.vibrational_frequencies_cm-1",
                        expected_modes,
                        None if frequencies is None else len(frequencies),
                        "output:vibrational_frequencies",
                    )
                )
            elif not all(math.isfinite(float(value)) for value in frequencies):
                findings.append(
                    _finding(
                        "xtb.result.frequency_nonfinite",
                        "result.vibrational_frequencies_cm-1",
                        "all finite",
                        list(frequencies),
                        "output:vibrational_frequencies",
                    )
                )

        if requested["grad"]:
            gradient_paths = (
                Path(job.folder) / "gradient",
                Path(job.folder) / f"{job.label}.engrad",
            )
            gradient_available = any(
                path.is_file() and path.stat().st_size
                for path in gradient_paths
            )
            if not gradient_available:
                findings.append(
                    _finding(
                        "xtb.result.gradient_missing",
                        "result.gradient",
                        "non-empty gradient or engrad artifact",
                        None,
                        "artifact:gradient",
                    )
                )

    artifact_names = {
        f"{job.label}.xyz",
        f"{job.label}.out",
        f"{job.label}.err",
        f"{job.label}.engrad",
        "gradient",
        "hessian",
        "vibspectrum",
        "g98.out",
        "xtbopt.xyz",
    }
    artifacts = {}
    for name in sorted(artifact_names):
        record = _artifact_record(Path(job.folder) / name)
        if record is not None:
            artifacts[name] = record

    source_record = provenance_binding.get("source_artifact")
    project_record = provenance_binding.get("project_artifact")
    if returncode is None:
        engine_state = "indeterminate"
    elif returncode != 0:
        engine_state = "failed"
    elif main is None:
        engine_state = "indeterminate"
    elif normal_termination is True:
        engine_state = "completed"
    else:
        engine_state = "incomplete"
    validation_state = "validated" if not findings else "failed"
    payload = {
        "schema_version": XTB_RESULT_SCHEMA_VERSION,
        # ``state`` remains as a compatibility view.  New consumers must use
        # the two orthogonal fields below rather than treating a zero process
        # exit as scientific validation.
        "state": validation_state,
        "engine_state": engine_state,
        "validation_state": validation_state,
        "engine_completion": {
            "state": engine_state,
            "returncode": returncode,
            "normal_termination": normal_termination,
        },
        "program": "xtb",
        "jobtype": requested["jobtype"],
        "ready": not findings,
        "returncode": returncode,
        "canonical_argv": list(command),
        "command_sha256": canonical_sha256(list(command)),
        "requested_settings": requested,
        "requested_settings_sha256": canonical_sha256(requested),
        "requested_molecule": _molecule_payload(job.molecule),
        "requested_molecule_sha256": canonical_sha256(
            _molecule_payload(job.molecule)
        ),
        "command_settings": command_settings,
        "command_settings_sha256": canonical_sha256(command_settings),
        "applied_settings": observed_settings,
        "applied_settings_sha256": (
            None
            if observed_settings is None
            else canonical_sha256(observed_settings)
        ),
        "results": {"final_energy_hartree": final_energy},
        "source_artifact": source_record,
        "source_index": job.source_index,
        "project_reference": job.project_reference,
        "project_artifact": project_record,
        "execution_input_artifact": execution_input,
        "provenance_binding_sha256": canonical_sha256(provenance_binding),
        "environment_receipt_sha256": (
            None
            if not environment_receipt
            else environment_receipt.get("receipt_sha256")
        ),
        "artifacts": artifacts,
        "warnings": warnings,
        "findings": findings,
    }
    return finalize_receipt(receipt_path, payload)


def load_validated_result_receipt(job):
    """Return a still-bound green receipt, otherwise ``None``.

    Completion is invalidated by a settings change, receipt substitution, or
    mutation/removal of any recorded output artifact.
    """

    result_path = Path(job.folder) / f"{job.label}.xtb-result-receipt.json"
    environment_path = (
        Path(job.folder) / f"{job.label}.xtb-environment-receipt.json"
    )
    result = _load_bound_receipt(result_path)
    environment = _load_bound_receipt(environment_path)
    if not result or result.get("schema_version") != XTB_RESULT_SCHEMA_VERSION:
        return None
    if result.get("ready") is not True or result.get("state") != "validated":
        return None
    if result.get("validation_state", "validated") != "validated":
        return None
    if result.get("engine_state", "completed") != "completed":
        return None
    requested = _settings_payload(job.settings)
    if result.get("requested_settings_sha256") != canonical_sha256(requested):
        return None
    if result.get("requested_molecule_sha256") != canonical_sha256(
        _molecule_payload(job.molecule)
    ):
        return None
    if not environment or environment.get("execution_ready") is not True:
        return None
    if (
        result.get("environment_receipt_sha256")
        != environment.get("receipt_sha256")
    ):
        return None
    for name, record in (result.get("artifacts") or {}).items():
        if Path(name).name != name:
            return None
        path = Path(job.folder) / name
        if not path.is_file():
            return None
        if path.stat().st_size != record.get("size"):
            return None
        if sha256_file(path) != record.get("sha256"):
            return None
    for attribute, key in (
        ("source_filename", "source_artifact"),
        ("project_source_file", "project_artifact"),
    ):
        record = result.get(key)
        path = getattr(job, attribute, None)
        if path and record is None:
            return None
        if record is None:
            continue
        if not path or not os.path.isfile(path):
            return None
        if os.path.abspath(path) != record.get("declared_path"):
            return None
        resolved = Path(path).resolve(strict=True)
        if str(resolved) != record.get("resolved_path"):
            return None
        if resolved.stat().st_size != record.get("size"):
            return None
        if sha256_file(resolved) != record.get("sha256"):
            return None
    return result


def audit_xtb_result_receipt(
    receipt_path,
    *,
    expected_jobtype,
    expected_charge,
    expected_multiplicity,
    expected_settings=None,
    expected_source_sha256="",
    expected_project_sha256="",
    audit_mode="strict",
):
    """Independently audit a completed xTB receipt for the agent host.

    The xTB runner creates the primary scientific validation receipt.  The
    approved-execution host must nevertheless verify its digest, ancestry,
    settings, environment, and artifact bindings rather than trusting three
    status fields.  This second check never reconstructs or reruns the job.

    ``strict`` is the default and remains required for runtime reuse.  The
    explicit ``archive`` mode is analysis-only: it may accept a missing
    original source/project path or executable after relocation, but still
    requires the bound receipt, requested settings and molecule, normal
    termination, durable execution input, and every local artifact byte to
    validate.  Any original path that still exists is checked strictly.
    """

    if audit_mode not in {"strict", "archive"}:
        raise ValueError("audit_mode must be 'strict' or 'archive'")

    path = Path(receipt_path)
    rule_ids = []
    provenance_limitations = []
    observation = {
        "receipt_path": path.name,
        "expected_jobtype": expected_jobtype,
        "audit_mode": audit_mode,
    }
    try:
        with path.open(encoding="utf-8") as handle:
            result = json.load(handle)
    except (OSError, TypeError, ValueError, json.JSONDecodeError):
        return observation, ("xtb.result.receipt_unreadable",)
    if not isinstance(result, Mapping):
        return observation, ("xtb.result.receipt_not_object",)

    receipt_body = dict(result)
    observed_receipt_sha256 = str(
        receipt_body.pop("receipt_sha256", "")
    )
    if observed_receipt_sha256 != canonical_sha256(receipt_body):
        rule_ids.append("xtb.result.receipt_digest_mismatch")

    observation.update(
        {
            "schema_version": result.get("schema_version"),
            "program": result.get("program"),
            "jobtype": result.get("jobtype"),
            "engine_state": result.get("engine_state"),
            "validation_state": result.get("validation_state"),
            "ready": result.get("ready"),
            "requested_settings_sha256": result.get(
                "requested_settings_sha256"
            ),
            "requested_molecule_sha256": result.get(
                "requested_molecule_sha256"
            ),
            "environment_receipt_sha256": result.get(
                "environment_receipt_sha256"
            ),
        }
    )
    expected_fields = {
        "schema_version": XTB_RESULT_SCHEMA_VERSION,
        "program": "xtb",
        "jobtype": expected_jobtype,
        "engine_state": "completed",
        "validation_state": "validated",
        "ready": True,
        "returncode": 0,
    }
    for field, expected in expected_fields.items():
        if result.get(field) != expected:
            rule_ids.append(f"xtb.result.{field}_mismatch")
    completion = result.get("engine_completion")
    if not isinstance(completion, Mapping) or (
        completion.get("state"),
        completion.get("returncode"),
        completion.get("normal_termination"),
    ) != ("completed", 0, True):
        rule_ids.append("xtb.result.engine_completion_mismatch")
    if result.get("findings") not in ([], ()):
        rule_ids.append("xtb.result.findings_not_empty")

    result_values = result.get("results")
    final_energy = (
        result_values.get("final_energy_hartree")
        if isinstance(result_values, Mapping)
        else None
    )
    energy_evidence = "receipt"
    if final_energy is None:
        # Older v1 receipts validated this value but did not persist it. Read
        # the still-bound main output without rerunning the engine.
        try:
            from chemsmart.io.xtb.output import XTBOutput

            final_energy = XTBOutput(folder=path.parent).final_energy
            energy_evidence = "main_output"
        except Exception:
            final_energy = None
    try:
        energy_is_finite = math.isfinite(float(final_energy))
    except (TypeError, ValueError):
        energy_is_finite = False
    observation["final_energy_hartree"] = final_energy
    observation["final_energy_evidence"] = energy_evidence
    if not energy_is_finite:
        rule_ids.append("xtb.result.final_energy_nonfinite_or_missing")

    requested = result.get("requested_settings")
    if not isinstance(requested, Mapping):
        rule_ids.append("xtb.result.requested_settings_missing")
        requested = {}
    elif result.get("requested_settings_sha256") != canonical_sha256(
        requested
    ):
        rule_ids.append("xtb.result.requested_settings_digest_mismatch")
    required_settings = {
        "jobtype": expected_jobtype,
        "charge": expected_charge,
        "multiplicity": expected_multiplicity,
    }
    required_settings.update(dict(expected_settings or {}))
    for field, expected in sorted(required_settings.items()):
        if requested.get(field, _MISSING) != expected:
            rule_ids.append(f"xtb.result.requested_settings.{field}_mismatch")

    command_settings = result.get("command_settings")
    if not isinstance(command_settings, Mapping):
        rule_ids.append("xtb.result.command_settings_missing")
        command_settings = {}
    elif result.get("command_settings_sha256") != canonical_sha256(
        command_settings
    ):
        rule_ids.append("xtb.result.command_settings_digest_mismatch")
    applied = result.get("applied_settings")
    if not isinstance(applied, Mapping):
        rule_ids.append("xtb.result.applied_settings_missing")
        applied = {}
    elif result.get("applied_settings_sha256") != canonical_sha256(applied):
        rule_ids.append("xtb.result.applied_settings_digest_mismatch")
    semantic_fields = (
        "jobtype",
        "gfn_version",
        "charge",
        "multiplicity",
        "solvent_model",
        "solvent_id",
        "grad",
    )
    if expected_jobtype == "opt":
        semantic_fields += ("optimization_level",)
    for field in semantic_fields:
        expected = requested.get(field, _MISSING)
        if command_settings.get(field, _MISSING) != expected:
            rule_ids.append(f"xtb.result.command_settings.{field}_mismatch")
        if applied.get(field, _MISSING) != expected:
            rule_ids.append(f"xtb.result.applied_settings.{field}_mismatch")

    molecule = result.get("requested_molecule")
    if not isinstance(molecule, Mapping):
        rule_ids.append("xtb.result.requested_molecule_missing")
        molecule = {}
    elif result.get("requested_molecule_sha256") != canonical_sha256(
        molecule
    ):
        rule_ids.append("xtb.result.requested_molecule_digest_mismatch")
    if (molecule.get("charge"), molecule.get("multiplicity")) != (
        expected_charge,
        expected_multiplicity,
    ):
        rule_ids.append("xtb.result.requested_molecule_state_mismatch")
    symbols = molecule.get("symbols")
    positions = molecule.get("positions_angstrom")
    coordinates_valid = bool(
        isinstance(symbols, list)
        and symbols
        and all(isinstance(symbol, str) and symbol for symbol in symbols)
        and isinstance(positions, list)
        and len(positions) == len(symbols)
    )
    if coordinates_valid:
        try:
            coordinates_valid = all(
                len(row) == 3
                and all(math.isfinite(float(value)) for value in row)
                for row in positions
            )
        except (TypeError, ValueError):
            coordinates_valid = False
    if not coordinates_valid:
        rule_ids.append("xtb.result.requested_molecule_geometry_invalid")

    canonical_argv = result.get("canonical_argv")
    if not (
        isinstance(canonical_argv, list)
        and len(canonical_argv) >= 2
        and all(isinstance(item, str) and item for item in canonical_argv)
    ):
        rule_ids.append("xtb.result.canonical_argv_invalid")
        canonical_argv = []
    elif result.get("command_sha256") != canonical_sha256(canonical_argv):
        rule_ids.append("xtb.result.command_digest_mismatch")

    binding = {
        "source_artifact": result.get("source_artifact"),
        "project_artifact": result.get("project_artifact"),
        "execution_input_artifact": result.get("execution_input_artifact"),
    }
    if result.get("provenance_binding_sha256") != canonical_sha256(binding):
        rule_ids.append("xtb.result.provenance_binding_digest_mismatch")
    verifiable_binding = {}
    for role, record in binding.items():
        required_fields = {
            "role",
            "declared_path",
            "resolved_path",
            "size",
            "sha256",
        }
        if record is not None and (
            not isinstance(record, Mapping)
            or not required_fields.issubset(record)
            or not isinstance(record.get("declared_path"), str)
            or not record.get("declared_path")
            or not isinstance(record.get("resolved_path"), str)
            or not record.get("resolved_path")
            or isinstance(record.get("size"), bool)
            or not isinstance(record.get("size"), int)
            or record.get("size") < 0
            or not isinstance(record.get("sha256"), str)
            or re.fullmatch(r"[0-9a-f]{64}", record.get("sha256")) is None
        ):
            rule_ids.append(f"xtb.result.{role}_record_invalid")
            verifiable_binding[role] = None
        else:
            verifiable_binding[role] = record
    strict_binding = {
        "source_artifact": verifiable_binding["source_artifact"],
        "project_artifact": verifiable_binding["project_artifact"],
    }
    for role, record in strict_binding.items():
        if record is None:
            if audit_mode == "archive":
                provenance_limitations.append(f"{role}_not_declared")
            continue
        declared = Path(record["declared_path"])
        original_path_present = declared.exists() or declared.is_symlink()
        if audit_mode == "archive" and not original_path_present:
            provenance_limitations.append(
                f"{role}_original_path_unavailable"
            )
            continue
        for finding in verify_xtb_provenance_binding({role: record}):
            rule_ids.append(str(finding.get("rule_id")))
    execution_input = verifiable_binding["execution_input_artifact"]
    if audit_mode == "archive" and execution_input is None:
        rule_ids.append("xtb.result.execution_input_artifact_record_invalid")
    for finding in _audit_xtb_execution_input_binding(
        execution_input,
        receipt_parent=path.parent,
        artifacts=result.get("artifacts"),
    ):
        rule_ids.append(str(finding.get("rule_id")))
    if audit_mode == "archive" and execution_input is not None:
        declared = Path(execution_input["declared_path"])
        if not (declared.exists() or declared.is_symlink()):
            provenance_limitations.append(
                "execution_input_original_path_unavailable_durable_copy_verified"
            )
    source = binding["source_artifact"]
    project = binding["project_artifact"]
    if expected_source_sha256 and (
        not isinstance(source, Mapping)
        or source.get("sha256") != expected_source_sha256
    ):
        rule_ids.append("xtb.result.source_artifact_mismatch")
    if expected_project_sha256 and (
        not isinstance(project, Mapping)
        or project.get("sha256") != expected_project_sha256
    ):
        rule_ids.append("xtb.result.project_artifact_mismatch")

    environment_path = None
    suffix = ".xtb-result-receipt.json"
    if path.name.endswith(suffix):
        environment_path = path.with_name(
            path.name[: -len(suffix)] + ".xtb-environment-receipt.json"
        )
    environment = (
        _load_bound_receipt(environment_path)
        if environment_path is not None
        else None
    )
    if not isinstance(environment, Mapping):
        rule_ids.append("xtb.result.environment_receipt_unavailable")
    else:
        if result.get("environment_receipt_sha256") != environment.get(
            "receipt_sha256"
        ):
            rule_ids.append("xtb.result.environment_receipt_mismatch")
        if (
            environment.get("schema_version")
            != XTB_ENVIRONMENT_SCHEMA_VERSION
            or environment.get("status") != "available"
            or environment.get("preflight_state") != "ready"
            or environment.get("execution_ready") is not True
            or environment.get("required_version") != XTB_REQUIRED_VERSION
            or environment.get("observed_version") != XTB_REQUIRED_VERSION
            or environment.get("findings") not in ([], ())
        ):
            rule_ids.append("xtb.result.environment_receipt_red")
        executable = environment.get("executable")
        executable_sha256 = environment.get("executable_sha256")
        executable_declared = (
            Path(executable)
            if isinstance(executable, str) and executable
            else None
        )
        executable_present = bool(
            executable_declared is not None
            and (
                executable_declared.exists()
                or executable_declared.is_symlink()
            )
        )
        if audit_mode == "archive" and not executable_present:
            if (
                executable_declared is None
                or not isinstance(executable_sha256, str)
                or re.fullmatch(r"[0-9a-f]{64}", executable_sha256) is None
            ):
                rule_ids.append("xtb.result.environment_executable_mismatch")
            else:
                provenance_limitations.append(
                    "execution_environment_executable_unavailable"
                )
                if canonical_argv and canonical_argv[0] != executable:
                    rule_ids.append("xtb.result.command_executable_mismatch")
        else:
            try:
                executable_path = executable_declared.resolve(strict=True)
            except (AttributeError, OSError, RuntimeError):
                executable_path = None
            if (
                executable_path is None
                or not executable_path.is_file()
                or sha256_file(executable_path) != executable_sha256
            ):
                rule_ids.append("xtb.result.environment_executable_mismatch")
            elif canonical_argv:
                try:
                    command_executable = Path(canonical_argv[0]).resolve(
                        strict=True
                    )
                except (OSError, RuntimeError):
                    command_executable = None
                if command_executable != executable_path:
                    rule_ids.append("xtb.result.command_executable_mismatch")

    artifacts = result.get("artifacts")
    if not isinstance(artifacts, Mapping):
        rule_ids.append("xtb.result.artifact_manifest_missing")
        artifacts = {}
    for name, record in sorted(artifacts.items()):
        if Path(str(name)).name != name or not isinstance(record, Mapping):
            rule_ids.append("xtb.result.artifact_record_invalid")
            continue
        artifact_path = path.parent / name
        if (
            not artifact_path.is_file()
            or artifact_path.stat().st_size != record.get("size")
            or sha256_file(artifact_path) != record.get("sha256")
        ):
            rule_ids.append(f"xtb.result.artifact.{name}_mismatch")
    if not any(str(name).endswith(".out") for name in artifacts):
        rule_ids.append("xtb.result.main_output_missing")
    if not any(str(name).endswith(".xyz") for name in artifacts):
        rule_ids.append("xtb.result.input_geometry_missing")

    provenance_limitations = tuple(sorted(set(provenance_limitations)))
    observation["provenance_limitations"] = provenance_limitations
    observation["provenance_status"] = (
        "archived_validated_native_result_degraded_provenance"
        if provenance_limitations
        else "workspace_exact_validated_native_result"
    )

    return observation, tuple(sorted(set(rule_ids)))


__all__ = [
    "XTB_ENVIRONMENT_SCHEMA_VERSION",
    "XTB_PREVIEW_SCHEMA_VERSION",
    "XTB_REQUIRED_VERSION",
    "XTB_RESULT_SCHEMA_VERSION",
    "XTBEnvironmentError",
    "XTBResultValidationError",
    "audit_xtb_result_receipt",
    "build_preview_receipt",
    "bind_xtb_declared_artifacts",
    "bind_xtb_execution_input",
    "canonical_sha256",
    "finalize_receipt",
    "load_validated_result_receipt",
    "probe_xtb_environment",
    "sha256_file",
    "validate_xtb_result",
    "verify_xtb_provenance_binding",
]
