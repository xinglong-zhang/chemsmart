"""Pure validation for requested and completed PySCF calculations.

The functions in this module never launch a solver, import an optional
dependency, or mutate an artifact.  Optional-package and GPU facts are
supplied as an ``environment`` mapping, which makes preflight deterministic
and lets an agent collect every known blocker in one pass.

Recognised environment evidence is intentionally small and permissive.  The
following keys may be flat or nested under ``dependencies`` / ``gpu``:

``solver_callables``
    A mapping from ``geometric``, ``berny`` or ``ase`` to the probed Python
    entry point (or ``{"callable": True}`` for an external-env attestation).
``solvent_ids`` / ``solvent_db``
    Valid PySCF solvent names as a collection or mapping.
``pyscf_dispersion``
    Whether the dispersion extra is available.
``gpu4pyscf``, ``cupy``, ``cuda_available``, ``cutensor_compatible``
    Required GPU runtime evidence.  ``basis_max_l`` and
    ``aux_basis_max_l`` optionally carry pre-probed angular momenta.

No availability conclusion is inferred by importing from the current
ChemSmart interpreter: the generated calculation may run in a different
``CONDA_ENV``.
"""

from __future__ import annotations

import json
import math
import os
from collections.abc import Collection, Mapping
from dataclasses import dataclass
from numbers import Integral, Real
from typing import Any

import numpy as np
from ase.data import atomic_numbers as ASE_ATOMIC_NUMBERS

from chemsmart.jobs.pyscf.settings import (
    FUNCTIONAL_DIVERGENCES,
    PYSCF_DEFGRIDS,
    PYSCF_ENGINES,
    PYSCF_OPT_SOLVERS,
    PYSCF_SOLVENT_MODELS,
    is_double_hybrid_functional,
)


@dataclass(frozen=True, slots=True)
class PySCFViolation:
    """One stable, machine-repairable PySCF validation finding."""

    rule_id: str
    field: str
    expected: Any
    observed: Any
    evidence_ref: str


RULE_UNSUPPORTED_SETTING = "pyscf.settings.unsupported"
RULE_INVALID_SETTING = "pyscf.settings.invalid_value"
RULE_SOLVENT_MODEL = "pyscf.solvent.model_unsupported"
RULE_SOLVENT_ID_REQUIRED = "pyscf.solvent.id_required"
RULE_SOLVENT_DATABASE = "pyscf.solvent.database_unavailable"
RULE_SOLVENT_ID_UNKNOWN = "pyscf.solvent.id_unknown"
RULE_SOLVER_UNCALLABLE = "pyscf.solver.uncallable"
RULE_DISPERSION_MISSING = "pyscf.dispersion.dependency_missing"
RULE_GPU_DEPENDENCY = "pyscf.gpu.dependency_missing"
RULE_GPU_CUDA = "pyscf.gpu.cuda_unavailable"
RULE_GPU_TENSOR = "pyscf.gpu.cutensor_incompatible"
RULE_GPU_BASIS = "pyscf.gpu.basis_angular_momentum"
RULE_GPU_AUX_BASIS = "pyscf.gpu.aux_basis_angular_momentum"
RULE_GPU_DOUBLE_HYBRID = "pyscf.gpu.double_hybrid"
RULE_GPU_LAPLACIAN = "pyscf.gpu.laplacian_meta_gga"
RULE_GPU_FUNCTIONAL_EVIDENCE = "pyscf.gpu.functional_unverified"
RULE_GPU_SOLVENT = "pyscf.gpu.solvent_unsupported"
RULE_GPU_TD_HESSIAN = "pyscf.gpu.tddft_hessian"
RULE_CHARGE_INVALID = "pyscf.electrons.charge_invalid"
RULE_MULTIPLICITY_INVALID = "pyscf.electrons.multiplicity_invalid"
RULE_ATOMIC_NUMBERS = "pyscf.electrons.atomic_numbers_unavailable"
RULE_ELECTRON_COUNT = "pyscf.electrons.count_inconsistent"
RULE_SPIN_PARITY = "pyscf.electrons.spin_parity"
RULE_SPIN_EXCEEDS = "pyscf.electrons.spin_exceeds_count"
RULE_SPIN_MISMATCH = "pyscf.electrons.spin_mismatch"
RULE_PREFLIGHT_INPUT = "pyscf.preflight.invalid_input"
RULE_PROVENANCE_READ = "pyscf.provenance.artifact_unreadable"
RULE_PROVENANCE_INCOMPLETE = "pyscf.provenance.incomplete_calculation"
RULE_PROVENANCE_SPEC = "pyscf.provenance.spec_mismatch"
RULE_PROVENANCE_ENGINE = "pyscf.provenance.engine_mismatch"
RULE_PROVENANCE_DIGEST = "pyscf.provenance.digest_mismatch"
RULE_PROVENANCE_PROJECT_DIGEST = (
    "pyscf.provenance.project_yaml_digest_mismatch"
)
RULE_PROVENANCE_RECEIPT = "pyscf.provenance.receipt_mismatch"
RULE_PROVENANCE_APPLIED_DIGEST = (
    "pyscf.provenance.applied_settings_digest_mismatch"
)
RULE_FREQUENCY_GEOMETRY = "pyscf.frequency.geometry_invalid"
RULE_FREQUENCY_MODE_COUNT = "pyscf.frequency.mode_count"
RULE_FREQUENCY_NONFINITE = "pyscf.frequency.nonfinite"
RULE_FREQUENCY_IMAGINARY = "pyscf.frequency.imaginary_modes"

FREQUENCY_VALIDATION_SCHEMA_VERSION = (
    "chemsmart.pyscf-frequency-validation.v1"
)
_LINEARITY_RELATIVE_TOLERANCE = 1.0e-8


_MISSING = object()

_SUPPORTED_FIELDS = frozenset(
    {
        "ab_initio",
        "functional",
        "dispersion",
        "basis",
        "aux_basis",
        "defgrid",
        "scf_tol",
        "scf_maxiter",
        "charge",
        "multiplicity",
        "freq",
        "density_fit",
        "opt_solver",
        "opt_maxsteps",
        "engine",
        "jobtype",
        "title",
        "solvent_model",
        "solvent_id",
        "project_yaml_digest",
    }
)

# Defaults distinguish inherited-but-unrequested fields from requests that the
# generated PySCF driver would otherwise ignore.
_UNSUPPORTED_DEFAULTS = {
    "semiempirical": None,
    "numfreq": False,
    "additional_route_parameters": None,
    "route_to_be_written": None,
    "modred": None,
    "gen_genecp_file": None,
    "heavy_elements": None,
    "heavy_elements_basis": None,
    "light_elements_basis": None,
    "custom_solvent": None,
    "forces": False,
    "input_string": None,
}

_SOLVER_ENTRYPOINTS = {
    "geometric": "pyscf.geomopt.geometric_solver.kernel",
    "berny": "pyscf.geomopt.berny_solver.kernel",
    "ase": "pyscf.geomopt.ase_solver.kernel",
}

_ANGULAR_MOMENTA = {
    "s": 0,
    "p": 1,
    "d": 2,
    "f": 3,
    "g": 4,
    "h": 5,
    "i": 6,
    "j": 7,
}

_GPU_SOLVENT_MODELS = frozenset(
    {"pcm", "iefpcm", "cpcm", "cosmo", "ssvpe", "smd"}
)


def preflight(settings, molecule, environment) -> list[PySCFViolation]:
    """Return every known pre-compute violation without running chemistry.

    ``environment`` is evidence gathered for the interpreter that will run
    the generated script.  Solver objects are inspected with ``callable`` but
    never invoked.  Missing evidence fails closed only for a requested
    optional feature (optimization, dispersion, solvent lookup, or GPU).
    """

    violations: list[PySCFViolation] = []
    checks = (
        _check_unsupported_settings,
        _check_setting_values,
        _check_solvent,
        _check_solver,
        _check_dispersion,
        _check_gpu,
        _check_electrons,
    )
    for check in checks:
        try:
            violations.extend(check(settings, molecule, environment))
        except Exception as exc:  # validation must enumerate, never raise
            violations.append(
                PySCFViolation(
                    rule_id=RULE_PREFLIGHT_INPUT,
                    field=check.__name__,
                    expected="valid preflight input",
                    observed=type(exc).__name__,
                    evidence_ref="preflight:input",
                )
            )
    return violations


def verify_provenance(
    settings, h5_path, *, expected_receipt=None
) -> list[PySCFViolation]:
    """Compare requested settings with the receipt stored in ``label.h5``.

    Both the historical schema-1 JSON datasets and schema-2 HDF5 groups are
    accepted.  Artifact/read errors are returned as typed violations rather
    than raised, preserving the same bounded-repair interface as preflight.
    This verifies receipt agreement and completion, not target-object
    introspection; the current digest is the controller-side requested-settings
    identity.
    """

    try:
        spec, provenance, status = _read_artifact_contract(h5_path)
    except Exception as exc:  # corrupt/missing artifacts are findings
        return [
            PySCFViolation(
                rule_id=RULE_PROVENANCE_READ,
                field="h5_path",
                expected="readable PySCF HDF5 provenance",
                observed=type(exc).__name__,
                evidence_ref="h5:/",
            )
        ]

    if not isinstance(spec, Mapping):
        return [
            PySCFViolation(
                rule_id=RULE_PROVENANCE_READ,
                field="spec",
                expected="mapping of applied settings",
                observed=type(spec).__name__,
                evidence_ref="h5:/spec",
            )
        ]

    violations: list[PySCFViolation] = []
    normal_termination = (
        status.get("normal_termination", False)
        if isinstance(status, Mapping)
        else False
    )
    if normal_termination is not True:
        failure = (
            status.get("failure") if isinstance(status, Mapping) else None
        )
        violations.append(
            PySCFViolation(
                rule_id=RULE_PROVENANCE_INCOMPLETE,
                field="status.normal_termination",
                expected=True,
                observed=bool(normal_termination),
                evidence_ref=(
                    "h5:/status/failure"
                    if failure is not None
                    else "h5:/status/normal_termination"
                ),
            )
        )
    requested = _requested_spec(settings)
    for field, expected in requested.items():
        observed = spec.get(field, _MISSING)
        if observed is _MISSING and expected is None:
            observed = None
        if not _equivalent(field, expected, observed):
            violations.append(
                PySCFViolation(
                    rule_id=RULE_PROVENANCE_SPEC,
                    field=field,
                    expected=expected,
                    observed=(
                        "<missing>" if observed is _MISSING else observed
                    ),
                    evidence_ref=f"h5:/spec/{field}",
                )
            )

    if isinstance(provenance, Mapping):
        expected_engine = requested.get("engine", _MISSING)
        observed_engine = provenance.get("engine", _MISSING)
        if expected_engine is not _MISSING and not _equivalent(
            "engine", expected_engine, observed_engine
        ):
            violations.append(
                PySCFViolation(
                    rule_id=RULE_PROVENANCE_ENGINE,
                    field="provenance.engine",
                    expected=expected_engine,
                    observed=(
                        "<missing>"
                        if observed_engine is _MISSING
                        else observed_engine
                    ),
                    evidence_ref="h5:/provenance/engine",
                )
            )

        spec_digest = spec.get("settings_digest", _MISSING)
        provenance_digest = provenance.get("settings_digest", _MISSING)
        if (
            spec_digest is _MISSING
            or provenance_digest is _MISSING
            or spec_digest != provenance_digest
        ):
            violations.append(
                PySCFViolation(
                    rule_id=RULE_PROVENANCE_DIGEST,
                    field="settings_digest",
                    expected=(
                        "present in both spec and provenance"
                        if spec_digest is _MISSING
                        else spec_digest
                    ),
                    observed=(
                        "<missing>"
                        if provenance_digest is _MISSING
                        else provenance_digest
                    ),
                    evidence_ref="h5:/provenance/settings_digest",
                )
            )

        expected_project_digest = _member(
            settings, "project_yaml_digest", _MISSING
        )
        observed_project_digest = provenance.get(
            "project_yaml_digest", _MISSING
        )
        if (
            expected_project_digest not in (_MISSING, None)
            and expected_project_digest != observed_project_digest
        ):
            violations.append(
                PySCFViolation(
                    rule_id=RULE_PROVENANCE_PROJECT_DIGEST,
                    field="project_yaml_digest",
                    expected=expected_project_digest,
                    observed=(
                        "<missing>"
                        if observed_project_digest is _MISSING
                        else observed_project_digest
                    ),
                    evidence_ref="h5:/provenance/project_yaml_digest",
                )
            )
        violations.extend(
            _verify_bound_receipt(
                spec,
                provenance,
                expected_receipt,
            )
        )
    elif expected_receipt:
        violations.append(
            PySCFViolation(
                rule_id=RULE_PROVENANCE_RECEIPT,
                field="provenance",
                expected="mapping with bound receipt hashes",
                observed=type(provenance).__name__,
                evidence_ref="h5:/provenance",
            )
        )
    return violations


def frequency_validation_receipt(
    *, symbols, positions, frequencies, expected_imaginary_modes=0
):
    """Validate a harmonic-frequency result without running chemistry.

    The mode-count oracle is ``3N-5`` for a linear molecule and ``3N-6`` for
    a nonlinear molecule. A single atom has zero vibrational modes. Linearity
    is classified from the result geometry using a recorded relative singular-
    value tolerance. The current PySCF v1 Hessian contract represents a local
    minimum, so callers use ``expected_imaginary_modes=0`` unless a future
    typed stationary-point contract explicitly says otherwise.
    """
    symbol_list = list(symbols or [])
    findings = []
    geometry_class = "unknown"
    expected_mode_count = None

    try:
        coordinates = np.asarray(positions, dtype=float)
    except (TypeError, ValueError):
        coordinates = np.empty((0, 3), dtype=float)
    geometry_valid = (
        coordinates.shape == (len(symbol_list), 3)
        and len(symbol_list) > 0
        and bool(np.isfinite(coordinates).all())
    )
    if not geometry_valid:
        findings.append(
            PySCFViolation(
                rule_id=RULE_FREQUENCY_GEOMETRY,
                field="results.positions",
                expected=(
                    f"finite ({len(symbol_list)}, 3) Cartesian geometry"
                ),
                observed={
                    "shape": list(coordinates.shape),
                    "all_finite": bool(np.isfinite(coordinates).all()),
                },
                evidence_ref="h5:/results/positions",
            )
        )
    elif len(symbol_list) == 1:
        geometry_class = "atomic"
        expected_mode_count = 0
    elif len(symbol_list) == 2:
        geometry_class = "linear"
        expected_mode_count = 1
    else:
        centered = coordinates - coordinates.mean(axis=0)
        singular_values = np.linalg.svd(centered, compute_uv=False)
        scale = float(singular_values[0]) if singular_values.size else 0.0
        if not math.isfinite(scale) or scale <= 0.0:
            findings.append(
                PySCFViolation(
                    rule_id=RULE_FREQUENCY_GEOMETRY,
                    field="results.positions",
                    expected="non-degenerate molecular geometry",
                    observed={"largest_singular_value": scale},
                    evidence_ref="h5:/results/positions",
                )
            )
        else:
            second = (
                float(singular_values[1])
                if singular_values.size > 1
                else 0.0
            )
            linear = second <= _LINEARITY_RELATIVE_TOLERANCE * scale
            geometry_class = "linear" if linear else "nonlinear"
            expected_mode_count = 3 * len(symbol_list) - (5 if linear else 6)

    complex_frequencies = False
    frequency_shape = None
    try:
        raw_frequencies = np.asarray(frequencies)
        frequency_shape = list(raw_frequencies.shape)
        complex_frequencies = bool(np.iscomplexobj(raw_frequencies))
        if complex_frequencies:
            values = np.asarray(raw_frequencies.real, dtype=float).reshape(-1)
        else:
            values = np.asarray(raw_frequencies, dtype=float).reshape(-1)
    except (TypeError, ValueError):
        values = np.asarray([], dtype=float)

    if frequency_shape is None or len(frequency_shape) != 1:
        findings.append(
            PySCFViolation(
                rule_id=RULE_FREQUENCY_MODE_COUNT,
                field="results.vibrational_frequencies",
                expected="one-dimensional frequency vector",
                observed=frequency_shape,
                evidence_ref="h5:/results/vibrational_frequencies",
            )
        )
    if expected_mode_count is not None and len(values) != expected_mode_count:
        findings.append(
            PySCFViolation(
                rule_id=RULE_FREQUENCY_MODE_COUNT,
                field="results.vibrational_frequencies",
                expected=expected_mode_count,
                observed=len(values),
                evidence_ref="h5:/results/vibrational_frequencies",
            )
        )

    finite_mask = np.isfinite(values)
    nonfinite_count = int((~finite_mask).sum())
    if complex_frequencies or nonfinite_count:
        findings.append(
            PySCFViolation(
                rule_id=RULE_FREQUENCY_NONFINITE,
                field="results.vibrational_frequencies",
                expected="finite real wavenumbers",
                observed={
                    "complex": complex_frequencies,
                    "nonfinite_count": nonfinite_count,
                },
                evidence_ref="h5:/results/vibrational_frequencies",
            )
        )

    imaginary_count = int((values[finite_mask] < 0.0).sum())
    if imaginary_count != expected_imaginary_modes:
        findings.append(
            PySCFViolation(
                rule_id=RULE_FREQUENCY_IMAGINARY,
                field="results.vibrational_frequencies",
                expected={
                    "imaginary_mode_count": expected_imaginary_modes,
                    "convention": "negative real wavenumber",
                },
                observed={"imaginary_mode_count": imaginary_count},
                evidence_ref="h5:/results/vibrational_frequencies",
            )
        )

    return {
        "schema_version": FREQUENCY_VALIDATION_SCHEMA_VERSION,
        "state": "validated" if not findings else "failed",
        "atom_count": len(symbol_list),
        "geometry_class": geometry_class,
        "linearity_relative_tolerance": _LINEARITY_RELATIVE_TOLERANCE,
        "expected_mode_count": expected_mode_count,
        "observed_mode_count": len(values),
        "finite_mode_count": int(finite_mask.sum()),
        "expected_imaginary_mode_count": expected_imaginary_modes,
        "observed_imaginary_mode_count": imaginary_count,
        "frequency_unit": "cm^-1",
        "findings": findings,
    }


def _verify_bound_receipt(spec, provenance, expected_receipt):
    """Verify exact controller-to-child hashes when a runner supplies them.

    ``expected_receipt`` is optional to keep historical V1 artifacts readable.
    A hardened run always supplies it and therefore fails closed on absent or
    substituted script, input, environment, geometry, project, or settings
    evidence.
    """
    if not expected_receipt:
        return []

    violations = []
    expected_fields = (
        "run_id",
        "run_nonce",
        "script_sha256",
        "input_receipt_sha256",
        "environment_receipt_sha256",
        "input_geometry_sha256",
        "input_artifact_kind",
        "input_artifact_sha256",
        "requested_settings_sha256",
        "project_yaml_digest",
    )
    spec_mirrors = frozenset(
        {
            "run_id",
            "run_nonce",
            "input_geometry_sha256",
            "input_artifact_kind",
            "input_artifact_sha256",
            "requested_settings_sha256",
        }
    )
    for field in expected_fields:
        expected = expected_receipt.get(field, _MISSING)
        if expected is _MISSING:
            continue
        observed = provenance.get(field, _MISSING)
        if observed is _MISSING or observed != expected:
            violations.append(
                PySCFViolation(
                    rule_id=RULE_PROVENANCE_RECEIPT,
                    field=f"provenance.{field}",
                    expected=expected,
                    observed=(
                        "<missing>" if observed is _MISSING else observed
                    ),
                    evidence_ref=f"h5:/provenance/{field}",
                )
            )
        if field in spec_mirrors:
            spec_observed = spec.get(field, _MISSING)
            if spec_observed is _MISSING or spec_observed != expected:
                violations.append(
                    PySCFViolation(
                        rule_id=RULE_PROVENANCE_RECEIPT,
                        field=f"spec.{field}",
                        expected=expected,
                        observed=(
                            "<missing>"
                            if spec_observed is _MISSING
                            else spec_observed
                        ),
                        evidence_ref=f"h5:/spec/{field}",
                    )
                )

    spec_applied = spec.get("applied_settings_sha256", _MISSING)
    provenance_applied = provenance.get(
        "applied_settings_sha256", _MISSING
    )
    require_applied = bool(
        expected_receipt.get("require_applied_settings_sha256", False)
    )
    if require_applied and spec_applied in (_MISSING, None, ""):
        violations.append(
            PySCFViolation(
                rule_id=RULE_PROVENANCE_APPLIED_DIGEST,
                field="spec.applied_settings_sha256",
                expected="non-empty digest of applied settings",
                observed=(
                    "<missing>" if spec_applied is _MISSING else spec_applied
                ),
                evidence_ref="h5:/spec/applied_settings_sha256",
            )
        )
    if spec_applied != provenance_applied:
        violations.append(
            PySCFViolation(
                rule_id=RULE_PROVENANCE_APPLIED_DIGEST,
                field="provenance.applied_settings_sha256",
                expected=(
                    "same value as spec.applied_settings_sha256"
                    if spec_applied is _MISSING
                    else spec_applied
                ),
                observed=(
                    "<missing>"
                    if provenance_applied is _MISSING
                    else provenance_applied
                ),
                evidence_ref="h5:/provenance/applied_settings_sha256",
            )
        )
    if spec_applied not in (_MISSING, None, ""):
        from chemsmart.jobs.pyscf.writer import (
            APPLIED_SPEC_FIELDS,
            PySCFScriptWriter,
        )

        applied_payload = {
            field: spec.get(field) for field in APPLIED_SPEC_FIELDS
        }
        observed_digest = PySCFScriptWriter.settings_digest(applied_payload)
        if observed_digest != spec_applied:
            violations.append(
                PySCFViolation(
                    rule_id=RULE_PROVENANCE_APPLIED_DIGEST,
                    field="spec.applied_settings_sha256",
                    expected=observed_digest,
                    observed=spec_applied,
                    evidence_ref="h5:/spec/applied_settings_sha256",
                )
            )
    return violations


def _check_unsupported_settings(settings, _molecule, _environment):
    defaults = dict(_UNSUPPORTED_DEFAULTS)
    declared = _member(settings, "UNSUPPORTED", {})
    if isinstance(declared, Mapping):
        defaults.update(declared)

    violations = []
    for field, value in _settings_items(settings):
        if str(field).startswith("_") or field in _SUPPORTED_FIELDS:
            continue
        if field in defaults and _is_unset(value, defaults[field]):
            continue
        violations.append(
            PySCFViolation(
                rule_id=RULE_UNSUPPORTED_SETTING,
                field=str(field),
                expected="supported PySCF setting or its unset default",
                observed=value,
                evidence_ref=f"settings:{field}",
            )
        )
    return violations


def _check_setting_values(settings, _molecule, _environment):
    violations = []
    jobtype = _member(settings, "jobtype", None)
    if jobtype not in {"sp", "opt", "hess"}:
        violations.append(
            PySCFViolation(
                rule_id=RULE_INVALID_SETTING,
                field="jobtype",
                expected=("sp", "opt", "hess"),
                observed=jobtype,
                evidence_ref="settings:jobtype",
            )
        )
    freq = _member(settings, "freq", False)
    if type(freq) is not bool:
        violations.append(
            PySCFViolation(
                rule_id=RULE_INVALID_SETTING,
                field="freq",
                expected="strict boolean",
                observed=freq,
                evidence_ref="settings:freq",
            )
        )
    elif jobtype == "sp" and freq:
        violations.append(
            PySCFViolation(
                rule_id=RULE_INVALID_SETTING,
                field="freq",
                expected=False,
                observed=True,
                evidence_ref="settings:freq",
            )
        )
    elif jobtype == "hess" and not freq:
        violations.append(
            PySCFViolation(
                rule_id=RULE_INVALID_SETTING,
                field="freq",
                expected=True,
                observed=False,
                evidence_ref="settings:freq",
            )
        )
    density_fit = _member(settings, "density_fit", True)
    if type(density_fit) is not bool:
        violations.append(
            PySCFViolation(
                rule_id=RULE_INVALID_SETTING,
                field="density_fit",
                expected="strict boolean",
                observed=density_fit,
                evidence_ref="settings:density_fit",
            )
        )
    scf_tol = _member(settings, "scf_tol", None)
    if scf_tol is not None and (
        isinstance(scf_tol, bool)
        or not isinstance(scf_tol, Real)
        or not math.isfinite(float(scf_tol))
        or float(scf_tol) <= 0
    ):
        violations.append(
            PySCFViolation(
                rule_id=RULE_INVALID_SETTING,
                field="scf_tol",
                expected="finite positive scalar",
                observed=scf_tol,
                evidence_ref="settings:scf_tol",
            )
        )
    for field in ("scf_maxiter", "opt_maxsteps"):
        value = _member(settings, field, 100 if field == "opt_maxsteps" else None)
        if value is not None and (
            isinstance(value, bool)
            or not isinstance(value, Integral)
            or int(value) <= 0
        ):
            violations.append(
                PySCFViolation(
                    rule_id=RULE_INVALID_SETTING,
                    field=field,
                    expected="positive integer",
                    observed=value,
                    evidence_ref=f"settings:{field}",
                )
            )
    engine = _member(settings, "engine", "cpu")
    if engine is not None and str(engine).lower() not in PYSCF_ENGINES:
        violations.append(
            PySCFViolation(
                rule_id=RULE_INVALID_SETTING,
                field="engine",
                expected=tuple(PYSCF_ENGINES),
                observed=engine,
                evidence_ref="settings:engine",
            )
        )

    solver = _member(settings, "opt_solver", "geometric")
    if solver is not None and str(solver).lower() not in PYSCF_OPT_SOLVERS:
        violations.append(
            PySCFViolation(
                rule_id=RULE_INVALID_SETTING,
                field="opt_solver",
                expected=tuple(PYSCF_OPT_SOLVERS),
                observed=solver,
                evidence_ref="settings:opt_solver",
            )
        )

    ab_initio = _member(settings, "ab_initio", None)
    functional = _member(settings, "functional", None)
    if ab_initio is not None and str(ab_initio).lower() != "hf":
        violations.append(
            PySCFViolation(
                rule_id=RULE_INVALID_SETTING,
                field="ab_initio",
                expected="hf or a DFT functional",
                observed=ab_initio,
                evidence_ref="settings:ab_initio",
            )
        )
    if ab_initio is not None and functional is not None:
        violations.append(
            PySCFViolation(
                rule_id=RULE_INVALID_SETTING,
                field="method",
                expected="exactly one of ab_initio or functional",
                observed={
                    "ab_initio": ab_initio,
                    "functional": functional,
                },
                evidence_ref="settings:method",
            )
        )

    if is_double_hybrid_functional(functional):
        violations.append(
            PySCFViolation(
                rule_id=RULE_INVALID_SETTING,
                field="functional",
                expected="non-double-hybrid DFT functional",
                observed=functional,
                evidence_ref="settings:functional",
            )
        )

    if ab_initio is None and not functional:
        violations.append(
            PySCFViolation(
                rule_id=RULE_INVALID_SETTING,
                field="method",
                expected="ab_initio='hf' or a DFT functional",
                observed=None,
                evidence_ref="settings:method",
            )
        )

    basis = _member(settings, "basis", None)
    if not basis:
        violations.append(
            PySCFViolation(
                rule_id=RULE_INVALID_SETTING,
                field="basis",
                expected="non-empty PySCF basis name",
                observed=basis,
                evidence_ref="settings:basis",
            )
        )

    defgrid = _member(settings, "defgrid", None)
    if defgrid is not None and str(defgrid).lower() not in PYSCF_DEFGRIDS:
        violations.append(
            PySCFViolation(
                rule_id=RULE_INVALID_SETTING,
                field="defgrid",
                expected=tuple(sorted(PYSCF_DEFGRIDS)),
                observed=defgrid,
                evidence_ref="settings:defgrid",
            )
        )
    elif defgrid is not None and ab_initio is not None:
        violations.append(
            PySCFViolation(
                rule_id=RULE_INVALID_SETTING,
                field="defgrid",
                expected="unset for ab_initio HF",
                observed=defgrid,
                evidence_ref="settings:defgrid",
            )
        )

    density_fit = _member(settings, "density_fit", True) is True
    aux_basis = _member(settings, "aux_basis", None)
    if not density_fit and aux_basis is not None:
        violations.append(
            PySCFViolation(
                rule_id=RULE_INVALID_SETTING,
                field="aux_basis",
                expected="unset when density_fit is disabled",
                observed=aux_basis,
                evidence_ref="settings:aux_basis",
            )
        )
    return violations


def _check_solvent(settings, _molecule, environment):
    model = _member(settings, "solvent_model", None)
    solvent_id = _member(settings, "solvent_id", None)
    violations = []

    if model is None:
        if solvent_id is not None:
            violations.append(
                PySCFViolation(
                    rule_id=RULE_SOLVENT_MODEL,
                    field="solvent_model",
                    expected="model for the requested solvent_id",
                    observed=None,
                    evidence_ref="settings:solvent_model",
                )
            )
        return violations

    normal_model = str(model).strip().lower()
    if normal_model not in PYSCF_SOLVENT_MODELS:
        violations.append(
            PySCFViolation(
                rule_id=RULE_SOLVENT_MODEL,
                field="solvent_model",
                expected=tuple(sorted(PYSCF_SOLVENT_MODELS)),
                observed=model,
                evidence_ref="settings:solvent_model",
            )
        )

    if normal_model in PYSCF_SOLVENT_MODELS and not solvent_id:
        violations.append(
            PySCFViolation(
                rule_id=RULE_SOLVENT_ID_REQUIRED,
                field="solvent_id",
                expected="non-empty PySCF solvent name",
                observed=solvent_id,
                evidence_ref="settings:solvent_id",
            )
        )
        return violations

    if solvent_id:
        solvent_keys = _solvent_keys(environment)
        if solvent_keys is _MISSING:
            violations.append(
                PySCFViolation(
                    rule_id=RULE_SOLVENT_DATABASE,
                    field="environment.solvent_db",
                    expected="probed PySCF solvent names",
                    observed="<missing>",
                    evidence_ref="environment:solvent_db",
                )
            )
        elif str(solvent_id).strip().lower() not in solvent_keys:
            violations.append(
                PySCFViolation(
                    rule_id=RULE_SOLVENT_ID_UNKNOWN,
                    field="solvent_id",
                    expected="member of PySCF solvent_db",
                    observed=solvent_id,
                    evidence_ref="environment:solvent_db",
                )
            )
    return violations


def _check_solver(settings, _molecule, environment):
    if not _optimization_requested(settings, environment):
        return []

    solver = str(_member(settings, "opt_solver", "geometric")).lower()
    if solver not in _SOLVER_ENTRYPOINTS:
        return []
    probe = _first_environment(
        environment,
        ("solver_callables", solver),
        ("solver_entrypoints", solver),
        ("solvers", solver),
        ("entrypoints", _SOLVER_ENTRYPOINTS[solver]),
        (_SOLVER_ENTRYPOINTS[solver],),
    )
    if _callable_probe(probe):
        return []
    return [
        PySCFViolation(
            rule_id=RULE_SOLVER_UNCALLABLE,
            field="opt_solver",
            expected=_SOLVER_ENTRYPOINTS[solver],
            observed=_describe_probe(probe),
            evidence_ref=f"environment:solver_callables/{solver}",
        )
    ]


def _check_dispersion(settings, _molecule, environment):
    dispersion = _member(settings, "dispersion", None)
    if not dispersion:
        return []
    probe = _dependency_probe(
        environment,
        "pyscf_dispersion",
        "pyscf-dispersion",
        "dispersion_available",
    )
    if _available_probe(probe):
        return []
    return [
        PySCFViolation(
            rule_id=RULE_DISPERSION_MISSING,
            field="dispersion",
            expected="pyscf-dispersion available",
            observed=_describe_probe(probe),
            evidence_ref="environment:dependencies/pyscf-dispersion",
        )
    ]


def _check_gpu(settings, _molecule, environment):
    if str(_member(settings, "engine", "cpu")).lower() != "gpu":
        return []

    violations = []
    for dependency, aliases in (
        ("gpu4pyscf", ("gpu4pyscf", "gpu4pyscf_available")),
        ("cupy", ("cupy", "cupy_available")),
    ):
        probe = _dependency_probe(environment, *aliases)
        if not _available_probe(probe):
            violations.append(
                PySCFViolation(
                    rule_id=RULE_GPU_DEPENDENCY,
                    field=f"environment.{dependency}",
                    expected=f"{dependency} available",
                    observed=_describe_probe(probe),
                    evidence_ref=f"environment:dependencies/{dependency}",
                )
            )

    cuda_probe = _first_environment(
        environment,
        ("cuda_available",),
        ("gpu", "cuda_available"),
        ("gpu", "device_count"),
        ("num_cuda_devices",),
    )
    if not _positive_or_available(cuda_probe):
        violations.append(
            PySCFViolation(
                rule_id=RULE_GPU_CUDA,
                field="environment.cuda_available",
                expected="at least one visible CUDA device",
                observed=_describe_probe(cuda_probe),
                evidence_ref="environment:gpu/cuda_available",
            )
        )

    tensor_probe = _cutensor_compatibility(environment)
    if tensor_probe is not True:
        violations.append(
            PySCFViolation(
                rule_id=RULE_GPU_TENSOR,
                field="environment.cutensor_compatible",
                expected="compatible CuPy/cuTENSOR pair using cuTENSOR",
                observed=_describe_probe(tensor_probe),
                evidence_ref="environment:gpu/cutensor_compatible",
            )
        )

    basis = _member(settings, "basis", None)
    basis_l = _basis_max_l(environment, "basis", basis)
    if basis_l is _MISSING:
        violations.append(
            PySCFViolation(
                rule_id=RULE_GPU_BASIS,
                field="basis",
                expected="probed maximum angular momentum <= g (l=4)",
                observed="<missing>",
                evidence_ref="environment:gpu/basis_max_l",
            )
        )
    elif basis_l > _ANGULAR_MOMENTA["g"]:
        violations.append(
            PySCFViolation(
                rule_id=RULE_GPU_BASIS,
                field="basis",
                expected="maximum angular momentum <= g (l=4)",
                observed=_angular_label(basis_l),
                evidence_ref="environment:gpu/basis_max_l",
            )
        )

    aux_basis = _member(settings, "aux_basis", None)
    aux_l = _basis_max_l(environment, "aux_basis", aux_basis)
    if aux_basis and aux_l is _MISSING:
        violations.append(
            PySCFViolation(
                rule_id=RULE_GPU_AUX_BASIS,
                field="aux_basis",
                expected="probed maximum angular momentum <= i (l=6)",
                observed="<missing>",
                evidence_ref="environment:gpu/aux_basis_max_l",
            )
        )
    elif aux_l is not _MISSING and aux_l > _ANGULAR_MOMENTA["i"]:
        violations.append(
            PySCFViolation(
                rule_id=RULE_GPU_AUX_BASIS,
                field="aux_basis",
                expected="maximum angular momentum <= i (l=6)",
                observed=_angular_label(aux_l),
                evidence_ref="environment:gpu/aux_basis_max_l",
            )
        )

    functional = _member(settings, "functional", None)
    double_hybrid = _functional_evidence(
        environment, functional, "double_hybrid"
    )
    laplacian_meta_gga = _functional_evidence(
        environment, functional, "laplacian_meta_gga"
    )
    name_is_double_hybrid = is_double_hybrid_functional(functional)
    functional_rejected = double_hybrid is True or name_is_double_hybrid
    if functional_rejected:
        violations.append(
            PySCFViolation(
                rule_id=RULE_GPU_DOUBLE_HYBRID,
                field="functional",
                expected="GPU-supported non-double-hybrid functional",
                observed=functional,
                evidence_ref="environment:functional_metadata/double_hybrid",
            )
        )

    if laplacian_meta_gga is True:
        violations.append(
            PySCFViolation(
                rule_id=RULE_GPU_LAPLACIAN,
                field="functional",
                expected="meta-GGA not requiring density Laplacian",
                observed=functional,
                evidence_ref=(
                    "environment:functional_metadata/laplacian_meta_gga"
                ),
            )
        )

    if (
        functional is not None
        and not functional_rejected
        and (double_hybrid is _MISSING or laplacian_meta_gga is _MISSING)
    ):
        missing_flags = []
        if double_hybrid is _MISSING:
            missing_flags.append("double_hybrid")
        if laplacian_meta_gga is _MISSING:
            missing_flags.append("laplacian_meta_gga")
        violations.append(
            PySCFViolation(
                rule_id=RULE_GPU_FUNCTIONAL_EVIDENCE,
                field="functional",
                expected="explicit GPU capability classification",
                observed={"missing": missing_flags},
                evidence_ref="environment:functional_metadata",
            )
        )

    model = _member(settings, "solvent_model", None)
    if model and str(model).lower() not in _GPU_SOLVENT_MODELS:
        violations.append(
            PySCFViolation(
                rule_id=RULE_GPU_SOLVENT,
                field="solvent_model",
                expected=tuple(sorted(_GPU_SOLVENT_MODELS)),
                observed=model,
                evidence_ref="settings:solvent_model",
            )
        )

    jobtype = str(_member(settings, "jobtype", "")).lower()
    stages = _requested_stages(settings)
    if "td" in jobtype and "hess" in stages:
        violations.append(
            PySCFViolation(
                rule_id=RULE_GPU_TD_HESSIAN,
                field="jobtype",
                expected="no GPU TDDFT Hessian",
                observed=jobtype,
                evidence_ref="settings:jobtype",
            )
        )

    return violations


def _check_electrons(settings, molecule, _environment):
    violations = []
    charge_raw = _resolved_value(settings, molecule, "charge")
    multiplicity_raw = _resolved_value(settings, molecule, "multiplicity")
    charge = _integral_value(charge_raw)
    multiplicity = _integral_value(multiplicity_raw)

    if charge is _MISSING:
        violations.append(
            PySCFViolation(
                rule_id=RULE_CHARGE_INVALID,
                field="charge",
                expected="integer resolved charge",
                observed=charge_raw,
                evidence_ref="resolved:charge",
            )
        )
    if multiplicity is _MISSING or multiplicity < 1:
        violations.append(
            PySCFViolation(
                rule_id=RULE_MULTIPLICITY_INVALID,
                field="multiplicity",
                expected="integer >= 1",
                observed=multiplicity_raw,
                evidence_ref="resolved:multiplicity",
            )
        )

    numbers = _atomic_numbers(molecule)
    if numbers is _MISSING:
        violations.append(
            PySCFViolation(
                rule_id=RULE_ATOMIC_NUMBERS,
                field="molecule.atomic_numbers",
                expected="valid atomic numbers or chemical symbols",
                observed="<unavailable>",
                evidence_ref="molecule:atomic_numbers",
            )
        )
        return violations
    if charge is _MISSING:
        return violations

    electron_count = sum(numbers) - charge
    if electron_count < 0:
        violations.append(
            PySCFViolation(
                rule_id=RULE_ELECTRON_COUNT,
                field="electron_count",
                expected=">= 0 from sum(Z) - charge",
                observed=electron_count,
                evidence_ref="derived:electron_count",
            )
        )

    declared_count = _declared_electron_count(molecule)
    if declared_count is not _MISSING:
        converted_count = _integral_value(declared_count)
        if converted_count is _MISSING or converted_count != electron_count:
            violations.append(
                PySCFViolation(
                    rule_id=RULE_ELECTRON_COUNT,
                    field="electron_count",
                    expected=electron_count,
                    observed=declared_count,
                    evidence_ref="molecule:electron_count",
                )
            )

    if multiplicity is _MISSING or multiplicity < 1:
        return violations
    spin = multiplicity - 1
    declared_spin = _stored_member(settings, "spin", _MISSING)
    if declared_spin is not _MISSING:
        converted_spin = _integral_value(declared_spin)
        if converted_spin is _MISSING or converted_spin != spin:
            violations.append(
                PySCFViolation(
                    rule_id=RULE_SPIN_MISMATCH,
                    field="spin",
                    expected=spin,
                    observed=declared_spin,
                    evidence_ref="settings:spin",
                )
            )

    if spin > electron_count:
        violations.append(
            PySCFViolation(
                rule_id=RULE_SPIN_EXCEEDS,
                field="multiplicity",
                expected=f"multiplicity <= {electron_count + 1}",
                observed=multiplicity,
                evidence_ref="derived:electron_count",
            )
        )
    elif (electron_count - spin) % 2:
        violations.append(
            PySCFViolation(
                rule_id=RULE_SPIN_PARITY,
                field="multiplicity",
                expected=(
                    "electron count and multiplicity - 1 with equal parity"
                ),
                observed={
                    "electron_count": electron_count,
                    "spin": spin,
                },
                evidence_ref="derived:electron_count",
            )
        )
    return violations


def _requested_spec(settings):
    requested = {}
    direct_fields = (
        "title",
        "engine",
        "basis",
        "ab_initio",
        "charge",
        "multiplicity",
        "dispersion",
        "density_fit",
        "aux_basis",
        "defgrid",
        "scf_tol",
        "scf_maxiter",
        "solvent_model",
        "solvent_id",
        "jobtype",
    )
    for field in direct_fields:
        value = _member(settings, field, _MISSING)
        if value is _MISSING:
            continue
        if value is None and field in {
            "charge",
            "multiplicity",
            "jobtype",
        }:
            continue
        requested[field] = value

    functional = _member(settings, "functional", _MISSING)
    ab_initio = _member(settings, "ab_initio", _MISSING)
    # Do not evaluate ``PySCFJobSettings.xc`` here: that property may emit a
    # scientific warning for divergent functional names. Verification is a
    # pure comparison, so only a literally stored override is inspected.
    direct_xc = _stored_member(settings, "xc", _MISSING)
    if direct_xc is not _MISSING:
        requested["xc"] = direct_xc
    elif ab_initio is not _MISSING and str(ab_initio).lower() == "hf":
        requested["xc"] = None
    elif functional is not _MISSING and functional is not None:
        key = str(functional).strip().lower()
        requested["xc"] = FUNCTIONAL_DIVERGENCES.get(key, (functional, None))[
            0
        ]

    if ab_initio not in (_MISSING, None):
        requested["method"] = str(ab_initio).lower()
    elif functional not in (_MISSING, None):
        requested["method"] = str(functional).lower()

    multiplicity = _integral_value(_member(settings, "multiplicity", _MISSING))
    if multiplicity is not _MISSING:
        requested["spin"] = multiplicity - 1

    stages = _requested_stages(settings)
    if stages:
        requested["stages"] = stages
    if "opt" in stages:
        requested["opt_solver"] = _member(settings, "opt_solver", "geometric")
        requested["opt_maxsteps"] = _member(settings, "opt_maxsteps", 100)
    return requested


def _requested_stages(settings):
    jobtype = _member(settings, "jobtype", None)
    if not jobtype:
        return []
    normal = str(jobtype).lower()
    if normal.endswith("sp") or normal == "sp":
        return ["scf"]
    if "opt" in normal:
        stages = ["scf", "opt"]
        if bool(_member(settings, "freq", False)):
            stages.append("hess")
        return stages
    if "hess" in normal or "freq" in normal:
        return ["scf", "hess"]
    return []


def _optimization_requested(settings, environment):
    jobtype = str(_member(settings, "jobtype", "")).lower()
    if "opt" in jobtype:
        return True
    stages = _first_environment(environment, ("required_stages",), ("stages",))
    if isinstance(stages, Collection) and not isinstance(stages, str):
        return "opt" in {str(stage).lower() for stage in stages}
    return False


def _solvent_keys(environment):
    source = _first_environment(
        environment,
        ("solvent_ids",),
        ("solvent_db",),
        ("pyscf", "solvent_ids"),
        ("pyscf", "solvent_db"),
        ("capabilities", "solvent_ids"),
    )
    if source is _MISSING or source is None:
        return _MISSING
    if isinstance(source, Mapping):
        source = source.keys()
    elif isinstance(source, str):
        source = (source,)
    if not isinstance(source, Collection):
        return _MISSING
    return {str(key).strip().lower() for key in source}


def _dependency_probe(environment, *names):
    paths = []
    for name in names:
        paths.extend(
            [
                (name,),
                ("dependencies", name),
                ("packages", name),
                ("modules", name),
                ("gpu", name),
            ]
        )
    return _first_environment(environment, *paths)


def _available_probe(probe):
    if probe is _MISSING or probe is None or probe is False:
        return False
    if isinstance(probe, Mapping):
        for key in ("available", "installed", "importable"):
            if key in probe:
                return bool(probe[key])
        return False
    return bool(probe)


def _callable_probe(probe):
    if callable(probe):
        return True
    if probe is _MISSING or probe is None or probe is False:
        return False
    if isinstance(probe, Mapping):
        value = probe.get("callable", _MISSING)
        return callable(value) or value is True
    # A bool is accepted only as an explicit external-environment
    # attestation in the solver_callables slot.  An importable module object
    # is deliberately not accepted: that is the known false-green mode.
    return probe is True


def _positive_or_available(probe):
    if isinstance(probe, bool):
        return probe
    value = _integral_value(probe)
    if value is not _MISSING:
        return value > 0
    return _available_probe(probe)


def _cutensor_compatibility(environment):
    explicit = _first_environment(
        environment,
        ("cutensor_compatible",),
        ("gpu", "cutensor_compatible"),
    )
    if explicit is not _MISSING:
        return explicit is True

    backend = _first_environment(
        environment, ("tensor_backend",), ("gpu", "tensor_backend")
    )
    if backend is not _MISSING:
        return str(backend).strip().lower() == "cutensor"

    cupy_version = _first_environment(
        environment, ("cupy_version",), ("gpu", "cupy_version")
    )
    cutensor_version = _first_environment(
        environment,
        ("cutensor_version",),
        ("gpu", "cutensor_version"),
    )
    known_pairs = {("13.3.0", "2.0.2"), ("13.4.1", "2.2.0")}
    if cupy_version is not _MISSING and cutensor_version is not _MISSING:
        return (str(cupy_version), str(cutensor_version)) in known_pairs
    return _MISSING


def _basis_max_l(environment, role, basis_name):
    names = {
        "basis": ("basis_max_l", "basis_max_angular_momentum"),
        "aux_basis": (
            "aux_basis_max_l",
            "auxbasis_max_l",
            "aux_basis_max_angular_momentum",
        ),
    }[role]
    paths = []
    for name in names:
        paths.extend(((name,), ("gpu", name)))
    value = _first_environment(environment, *paths)
    if isinstance(value, Mapping):
        value = value.get(
            basis_name,
            value.get(str(basis_name).lower(), value.get("default", _MISSING)),
        )
    if value is _MISSING or value is None:
        return _MISSING
    if isinstance(value, str):
        normal = value.strip().lower()
        if normal in _ANGULAR_MOMENTA:
            return _ANGULAR_MOMENTA[normal]
    converted = _integral_value(value)
    return converted


def _functional_evidence(environment, functional, flag):
    direct = _first_environment(environment, (flag,), ("gpu", flag))
    if direct is not _MISSING:
        return direct is True
    metadata = _first_environment(
        environment,
        ("functional_metadata",),
        ("gpu", "functional_metadata"),
    )
    if not isinstance(metadata, Mapping):
        return _MISSING
    item = metadata.get(
        functional,
        metadata.get(str(functional).lower(), metadata.get("default", {})),
    )
    if not isinstance(item, Mapping) or flag not in item:
        return _MISSING
    return item[flag] is True


def _atomic_numbers(molecule):
    raw = _first_member(molecule, "atomic_numbers", "numbers")
    if raw is not _MISSING:
        try:
            values = list(raw)
            converted = [_integral_value(value) for value in values]
            if converted and all(value is not _MISSING for value in converted):
                return converted
        except (TypeError, ValueError):
            pass

    symbols = _first_member(molecule, "chemical_symbols", "symbols")
    if symbols is _MISSING:
        return _MISSING
    try:
        values = [ASE_ATOMIC_NUMBERS[str(symbol)] for symbol in symbols]
    except (KeyError, TypeError):
        return _MISSING
    return values or _MISSING


def _declared_electron_count(molecule):
    direct = _first_member(
        molecule, "electron_count", "num_electrons", "nelectron"
    )
    if direct is not _MISSING:
        return direct
    info = _member(molecule, "info", {})
    if isinstance(info, Mapping):
        for key in ("electron_count", "num_electrons", "nelectron"):
            if key in info:
                return info[key]
    return _MISSING


def _resolved_value(settings, molecule, field):
    value = _member(settings, field, _MISSING)
    if value is not _MISSING and value is not None:
        return value
    return _member(molecule, field, _MISSING)


def _integral_value(value):
    if value is _MISSING or value is None or isinstance(value, bool):
        return _MISSING
    if isinstance(value, Integral):
        return int(value)
    if isinstance(value, Real):
        return int(value) if float(value).is_integer() else _MISSING
    if isinstance(value, str):
        try:
            converted = int(value.strip())
        except ValueError:
            return _MISSING
        return converted if str(converted) == value.strip() else _MISSING
    return _MISSING


def _angular_label(value):
    for label, number in _ANGULAR_MOMENTA.items():
        if value == number:
            return f"{label} (l={number})"
    return f"l={value}"


def _settings_items(settings):
    if isinstance(settings, Mapping):
        return list(settings.items())
    try:
        return list(vars(settings).items())
    except TypeError:
        return []


def _is_unset(value, unset):
    try:
        comparison = value == unset
        return comparison if isinstance(comparison, bool) else False
    except Exception:
        return False


def _first_member(source, *names):
    for name in names:
        value = _member(source, name, _MISSING)
        if value is not _MISSING:
            return value
    return _MISSING


def _member(source, name, default=_MISSING):
    try:
        if isinstance(source, Mapping):
            return source.get(name, default)
        return getattr(source, name, default)
    except Exception:
        return default


def _stored_member(source, name, default=_MISSING):
    """Read mappings/instance storage without evaluating a property."""
    if isinstance(source, Mapping):
        return source.get(name, default)
    try:
        return vars(source).get(name, default)
    except TypeError:
        return default


def _first_environment(environment, *paths):
    for path in paths:
        value = _environment_path(environment, path)
        if value is not _MISSING:
            return value
    return _MISSING


def _environment_path(environment, path):
    current = environment
    for name in path:
        current = _member(current, name, _MISSING)
        if current is _MISSING:
            return _MISSING
    return current


def _describe_probe(probe):
    if probe is _MISSING:
        return "<missing>"
    if callable(probe):
        return "callable"
    if isinstance(probe, Mapping):
        return {
            str(key): ("callable" if callable(value) else bool(value))
            for key, value in probe.items()
            if key in {"available", "installed", "importable", "callable"}
        }
    if isinstance(probe, (str, int, float, bool)) or probe is None:
        return probe
    return type(probe).__name__


def _equivalent(field, expected, observed):
    if observed is _MISSING:
        return False
    expected = _plain_value(expected)
    observed = _plain_value(observed)
    if field in {
        "engine",
        "basis",
        "ab_initio",
        "xc",
        "method",
        "dispersion",
        "defgrid",
        "solvent_model",
        "solvent_id",
        "opt_solver",
        "jobtype",
    }:
        if isinstance(expected, str) and isinstance(observed, str):
            return expected.strip().lower() == observed.strip().lower()
    return expected == observed


def _plain_value(value):
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, tuple):
        return [_plain_value(item) for item in value]
    if isinstance(value, list):
        return [_plain_value(item) for item in value]
    if isinstance(value, Mapping):
        return {str(key): _plain_value(item) for key, item in value.items()}
    return value


def _read_artifact_contract(h5_path):
    # Prefer B2's public schema reader when available.  The fallback keeps
    # this module usable against an in-progress Stage-A tree and recognizes
    # the exact schema-2 null marker.
    try:
        from chemsmart.io.pyscf.output import read_pyscf_h5
    except ImportError:
        read_pyscf_h5 = None

    if read_pyscf_h5 is not None:
        payload = read_pyscf_h5(h5_path)
        if isinstance(payload, Mapping):
            return (
                payload["spec"],
                payload.get("provenance", {}),
                payload.get("status", {}),
            )
        if isinstance(payload, tuple) and len(payload) >= 3:
            return payload[0], payload[1], payload[2]

    import h5py

    with h5py.File(os.fspath(h5_path), "r") as handle:
        return (
            _decode_h5_node(handle["spec"]),
            _decode_h5_node(handle["provenance"]),
            _decode_h5_node(handle["status"]),
        )


def _decode_h5_node(node):
    import h5py

    if isinstance(node, h5py.Group):
        return {key: _decode_h5_node(node[key]) for key in node}
    if bool(node.attrs.get("chemsmart_is_null", False)):
        return None
    return _decode_h5_value(node[()])


def _decode_h5_value(value):
    if hasattr(value, "item") and not hasattr(value, "tolist"):
        value = value.item()
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    if isinstance(value, str):
        stripped = value.strip()
        if stripped == "null" or stripped.startswith(("{", "[", '"')):
            try:
                return json.loads(stripped)
            except json.JSONDecodeError:
                pass
        return value
    if hasattr(value, "tolist"):
        return _plain_value(value.tolist())
    return value


__all__ = [
    "FREQUENCY_VALIDATION_SCHEMA_VERSION",
    "PySCFViolation",
    "frequency_validation_receipt",
    "preflight",
    "verify_provenance",
]
