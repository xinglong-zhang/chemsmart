"""Focused, compute-free tests for PySCF validation."""

import copy
import json
from dataclasses import FrozenInstanceError
from types import SimpleNamespace

import h5py
import numpy as np
import pytest

from chemsmart.jobs.pyscf.validation import (
    PySCFViolation,
    frequency_validation_receipt,
    preflight,
    verify_provenance,
)


class CallableProbe:
    """A callable whose invocation would prove that preflight ran work."""

    def __init__(self):
        self.calls = 0

    def __call__(self, *args, **kwargs):
        del args, kwargs
        self.calls += 1
        raise AssertionError("preflight must not invoke solver entry points")


@pytest.fixture
def water():
    return SimpleNamespace(
        chemical_symbols=["O", "H", "H"],
        charge=0,
        multiplicity=1,
    )


def test_violation_is_immutable():
    finding = PySCFViolation("rule", "field", "expected", "seen", "ref")
    with pytest.raises(FrozenInstanceError):
        finding.field = "changed"


def test_preflight_accepts_valid_evidence_without_calling_solver(water):
    solver = CallableProbe()
    settings = {
        "functional": "pbe",
        "basis": "def2-svp",
        "charge": 0,
        "multiplicity": 1,
        "jobtype": "opt",
        "opt_solver": "geometric",
        "engine": "cpu",
        "dispersion": "d4",
        "solvent_model": "smd",
        "solvent_id": "water",
    }
    environment = {
        "solver_callables": {"geometric": solver},
        "dependencies": {"pyscf-dispersion": True},
        "solvent_ids": {"water", "methanol"},
    }

    assert preflight(settings, water, environment) == []
    assert solver.calls == 0


def test_preflight_reports_import_only_false_green_and_all_blockers(water):
    settings = {
        "functional": "pbe",
        "basis": "def2-svp",
        "charge": 0,
        # Ten electrons cannot have spin=1; this is a parity mismatch.
        "multiplicity": 2,
        "jobtype": "opt",
        "opt_solver": "geometric",
        "engine": "cpu",
        "dispersion": "d4",
        "solvent_model": "smd",
        "solvent_id": "not-a-pyscf-solvent",
        "gen_genecp_file": "basis.genecp",
    }
    environment = {
        "solver_callables": {
            "geometric": {"importable": True},
        },
        "dependencies": {"pyscf-dispersion": False},
        "solvent_db": {"water": object()},
    }

    findings = preflight(settings, water, environment)
    rule_ids = {finding.rule_id for finding in findings}

    assert "pyscf.settings.unsupported" in rule_ids
    assert "pyscf.solver.uncallable" in rule_ids
    assert "pyscf.dispersion.dependency_missing" in rule_ids
    assert "pyscf.solvent.id_unknown" in rule_ids
    assert "pyscf.electrons.spin_parity" in rule_ids


def test_preflight_reports_gpu_limits_and_environment(water):
    settings = {
        "functional": "b2plyp",
        "basis": "mock-orbital-basis",
        "aux_basis": "mock-aux-basis",
        "charge": 0,
        "multiplicity": 1,
        "jobtype": "hess",
        "engine": "gpu",
        "solvent_model": "smd",
        "solvent_id": "water",
    }
    environment = {
        "dependencies": {"gpu4pyscf": True, "cupy": True},
        "solvent_ids": {"water"},
        "gpu": {
            "cuda_available": True,
            "cutensor_compatible": False,
            "basis_max_l": "h",
            "aux_basis_max_l": "j",
            "laplacian_meta_gga": True,
        },
    }

    rule_ids = {
        finding.rule_id for finding in preflight(settings, water, environment)
    }
    assert "pyscf.gpu.cutensor_incompatible" in rule_ids
    assert "pyscf.gpu.basis_angular_momentum" in rule_ids
    assert "pyscf.gpu.aux_basis_angular_momentum" in rule_ids
    assert "pyscf.gpu.double_hybrid" in rule_ids
    assert "pyscf.gpu.laplacian_meta_gga" in rule_ids


def test_gpu_preflight_fails_closed_without_capability_evidence(water):
    settings = {
        "functional": "pbe",
        "basis": "unclassified-basis",
        "charge": 0,
        "multiplicity": 1,
        "jobtype": "hess",
        "engine": "gpu",
    }
    environment = {
        "dependencies": {"gpu4pyscf": True, "cupy": True},
        "gpu": {
            "cuda_available": True,
            "cutensor_compatible": True,
        },
    }

    findings = preflight(settings, water, environment)
    rule_ids = {finding.rule_id for finding in findings}

    assert "pyscf.gpu.basis_angular_momentum" in rule_ids
    assert "pyscf.gpu.functional_unverified" in rule_ids


def test_preflight_returns_typed_method_selection_violations(water):
    settings = {
        "ab_initio": "mp2",
        "functional": "pbe",
        "basis": "def2-svp",
        "charge": 0,
        "multiplicity": 1,
        "jobtype": "sp",
        "engine": "cpu",
    }

    findings = preflight(settings, water, {})
    invalid_fields = {
        finding.field
        for finding in findings
        if finding.rule_id == "pyscf.settings.invalid_value"
    }

    assert invalid_fields == {"ab_initio", "method"}


def test_preflight_returns_typed_missing_method_basis_and_grid(water):
    settings = {
        "charge": 0,
        "multiplicity": 1,
        "jobtype": "sp",
        "engine": "cpu",
        "defgrid": "defgrid99",
    }

    findings = preflight(settings, water, {})
    invalid_fields = {
        finding.field
        for finding in findings
        if finding.rule_id == "pyscf.settings.invalid_value"
    }

    assert invalid_fields == {"method", "basis", "defgrid"}


@pytest.mark.parametrize(
    ("settings", "field"),
    [
        (
            {
                "functional": "b2plyp",
                "basis": "def2-svp",
                "charge": 0,
                "multiplicity": 1,
                "engine": "cpu",
            },
            "functional",
        ),
        (
            {
                "ab_initio": "hf",
                "basis": "def2-svp",
                "defgrid": "defgrid2",
                "charge": 0,
                "multiplicity": 1,
                "engine": "cpu",
            },
            "defgrid",
        ),
        (
            {
                "functional": "pbe0",
                "basis": "def2-svp",
                "density_fit": False,
                "aux_basis": "def2-universal-jkfit",
                "charge": 0,
                "multiplicity": 1,
                "engine": "cpu",
            },
            "aux_basis",
        ),
    ],
)
def test_preflight_types_requests_that_cannot_be_fully_applied(
    water, settings, field
):
    findings = preflight(settings, water, {})

    assert any(
        finding.rule_id == "pyscf.settings.invalid_value"
        and finding.field == field
        for finding in findings
    )


def test_preflight_checks_declared_electron_count_against_charge():
    molecule = SimpleNamespace(
        chemical_symbols=["O", "H", "H"],
        charge=0,
        multiplicity=1,
        electron_count=9,
    )

    findings = preflight(
        {"charge": 0, "multiplicity": 1, "engine": "cpu"},
        molecule,
        {},
    )

    mismatch = [
        finding
        for finding in findings
        if finding.rule_id == "pyscf.electrons.count_inconsistent"
    ]
    assert len(mismatch) == 1
    assert mismatch[0].expected == 10
    assert mismatch[0].observed == 9


def test_preflight_converts_malformed_evidence_to_a_violation(water):
    class BrokenEnvironment:
        @property
        def solvent_ids(self):
            raise RuntimeError("broken probe")

    settings = {
        "charge": 0,
        "multiplicity": 1,
        "engine": "cpu",
        "solvent_model": "smd",
        "solvent_id": "water",
    }

    findings = preflight(settings, water, BrokenEnvironment())
    assert findings
    assert all(isinstance(finding, PySCFViolation) for finding in findings)


def test_frequency_validation_accepts_three_finite_water_modes():
    receipt = frequency_validation_receipt(
        symbols=["O", "H", "H"],
        positions=[
            [0.0, 0.0, 0.0],
            [0.7586, 0.0, 0.5043],
            [-0.7586, 0.0, 0.5043],
        ],
        frequencies=[1595.0, 3657.0, 3756.0],
    )

    assert receipt["state"] == "validated"
    assert receipt["geometry_class"] == "nonlinear"
    assert receipt["expected_mode_count"] == 3
    assert receipt["observed_mode_count"] == 3
    assert receipt["finite_mode_count"] == 3
    assert receipt["observed_imaginary_mode_count"] == 0
    assert receipt["findings"] == []


def test_frequency_validation_rejects_wrong_nonfinite_and_imaginary_modes():
    receipt = frequency_validation_receipt(
        symbols=["O", "H", "H"],
        positions=[
            [0.0, 0.0, 0.0],
            [0.7586, 0.0, 0.5043],
            [-0.7586, 0.0, 0.5043],
        ],
        frequencies=[-101.0, float("nan")],
    )

    rule_ids = {finding.rule_id for finding in receipt["findings"]}

    assert receipt["state"] == "failed"
    assert "pyscf.frequency.mode_count" in rule_ids
    assert "pyscf.frequency.nonfinite" in rule_ids
    assert "pyscf.frequency.imaginary_modes" in rule_ids


REQUESTED_SETTINGS = {
    "functional": "b3lyp",
    "basis": "def2-svp",
    "charge": 0,
    "multiplicity": 1,
    "engine": "cpu",
    "solvent_model": "smd",
    "solvent_id": "water",
}

APPLIED_SPEC = {
    "xc": "b3lyp",
    "method": "b3lyp",
    "basis": "def2-svp",
    "charge": 0,
    "multiplicity": 1,
    "spin": 0,
    "engine": "cpu",
    "solvent_model": "smd",
    "solvent_id": "water",
    "settings_digest": "mock-digest",
}

PROVENANCE = {
    "engine": "cpu",
    "settings_digest": "mock-digest",
}


def _write_artifact(
    path, schema_version, spec, provenance=PROVENANCE, status=None
):
    status = status or {"normal_termination": True, "failure": None}
    with h5py.File(path, "w") as handle:
        handle.attrs["schema_version"] = schema_version
        if schema_version == "1.0":
            handle.create_dataset("spec", data=json.dumps(spec))
            handle.create_dataset("provenance", data=json.dumps(provenance))
            handle.create_dataset("status", data=json.dumps(status))
            handle.create_group("results")
            return

        _write_mapping(handle.create_group("spec"), spec)
        _write_mapping(handle.create_group("provenance"), provenance)
        _write_mapping(
            handle.create_group("status"),
            status,
        )
        handle.create_group("results")


def _write_mapping(group, mapping):
    for key, value in mapping.items():
        if isinstance(value, dict):
            _write_mapping(group.create_group(key), value)
        elif value is None:
            dataset = group.create_dataset(
                key, data=np.empty((0,), dtype=np.uint8)
            )
            dataset.attrs["chemsmart_is_null"] = True
        else:
            group.create_dataset(key, data=value)


@pytest.mark.parametrize("schema_version", ["1.0", "2.0"])
def test_verify_provenance_accepts_matching_v1_and_v2_artifacts(
    tmp_path, schema_version
):
    path = tmp_path / f"matching-{schema_version}.h5"
    _write_artifact(path, schema_version, APPLIED_SPEC)

    assert verify_provenance(REQUESTED_SETTINGS, path) == []


@pytest.mark.parametrize("schema_version", ["1.0", "2.0"])
@pytest.mark.parametrize(
    ("field", "wrong_value"),
    [("xc", "pbe"), ("spin", 2), ("solvent_model", "pcm")],
)
def test_verify_provenance_catches_injected_spec_mismatches(
    tmp_path, schema_version, field, wrong_value
):
    spec = copy.deepcopy(APPLIED_SPEC)
    spec[field] = wrong_value
    path = tmp_path / f"wrong-{field}-{schema_version}.h5"
    _write_artifact(path, schema_version, spec)

    findings = verify_provenance(REQUESTED_SETTINGS, path)
    mismatch = [
        finding
        for finding in findings
        if finding.rule_id == "pyscf.provenance.spec_mismatch"
        and finding.field == field
    ]

    assert len(mismatch) == 1
    assert mismatch[0].observed == wrong_value
    assert mismatch[0].evidence_ref == f"h5:/spec/{field}"


def test_verify_provenance_returns_unreadable_artifact_violation(tmp_path):
    findings = verify_provenance(
        REQUESTED_SETTINGS, tmp_path / "does-not-exist.h5"
    )

    assert len(findings) == 1
    assert findings[0].rule_id == "pyscf.provenance.artifact_unreadable"


def test_verify_provenance_requires_settings_digest_in_both_groups(tmp_path):
    path = tmp_path / "missing-provenance-digest.h5"
    _write_artifact(path, "2.0", APPLIED_SPEC, {"engine": "cpu"})

    findings = verify_provenance(REQUESTED_SETTINGS, path)

    assert len(findings) == 1
    assert findings[0].rule_id == "pyscf.provenance.digest_mismatch"


def test_verify_provenance_checks_project_yaml_digest(tmp_path):
    path = tmp_path / "wrong-project-digest.h5"
    provenance = dict(PROVENANCE, project_yaml_digest="observed-digest")
    _write_artifact(path, "2.0", APPLIED_SPEC, provenance)
    settings = dict(
        REQUESTED_SETTINGS,
        project_yaml_digest="expected-digest",
    )

    findings = verify_provenance(settings, path)

    assert len(findings) == 1
    assert findings[0].rule_id == (
        "pyscf.provenance.project_yaml_digest_mismatch"
    )
    assert findings[0].evidence_ref == ("h5:/provenance/project_yaml_digest")


def test_verify_provenance_never_greenlights_failed_calculation(tmp_path):
    path = tmp_path / "failed-with-matching-spec.h5"
    _write_artifact(
        path,
        "2.0",
        APPLIED_SPEC,
        status={
            "normal_termination": False,
            "failure": {
                "type": "RuntimeError",
                "message": "SCF did not converge",
            },
        },
    )

    findings = verify_provenance(REQUESTED_SETTINGS, path)

    assert any(
        finding.rule_id == "pyscf.provenance.incomplete_calculation"
        for finding in findings
    )


def test_verify_provenance_binds_controller_receipt_hashes(tmp_path):
    path = tmp_path / "substituted-script.h5"
    spec = dict(
        APPLIED_SPEC,
        run_id="run-a",
        run_nonce="nonce-a",
        input_geometry_sha256="geometry-sha",
        input_artifact_kind="pyscf_hdf5",
        input_artifact_sha256="wrong-source-sha",
        requested_settings_sha256="requested-sha",
        applied_settings_sha256=None,
    )
    provenance = dict(
        PROVENANCE,
        run_id="run-a",
        run_nonce="nonce-a",
        script_sha256="wrong-script-sha",
        input_receipt_sha256="input-sha",
        environment_receipt_sha256="environment-sha",
        input_geometry_sha256="geometry-sha",
        input_artifact_kind="pyscf_hdf5",
        input_artifact_sha256="wrong-source-sha",
        requested_settings_sha256="requested-sha",
        applied_settings_sha256=None,
        project_yaml_digest=None,
    )
    _write_artifact(path, "2.0", spec, provenance)

    findings = verify_provenance(
        REQUESTED_SETTINGS,
        path,
        expected_receipt={
            "run_id": "run-a",
            "run_nonce": "nonce-a",
            "script_sha256": "expected-script-sha",
            "input_receipt_sha256": "input-sha",
            "environment_receipt_sha256": "environment-sha",
            "input_geometry_sha256": "geometry-sha",
            "input_artifact_kind": "pyscf_hdf5",
            "input_artifact_sha256": "expected-source-sha",
            "requested_settings_sha256": "requested-sha",
            "project_yaml_digest": None,
            "require_applied_settings_sha256": False,
        },
    )

    assert any(
        finding.rule_id == "pyscf.provenance.receipt_mismatch"
        and finding.field == "provenance.script_sha256"
        for finding in findings
    )
    assert any(
        finding.rule_id == "pyscf.provenance.receipt_mismatch"
        and finding.field == "spec.input_artifact_sha256"
        for finding in findings
    )


def test_verify_provenance_requires_applied_digest_for_real_run(tmp_path):
    path = tmp_path / "missing-applied-digest.h5"
    spec = dict(
        APPLIED_SPEC,
        input_geometry_sha256="geometry-sha",
        requested_settings_sha256="requested-sha",
        applied_settings_sha256=None,
    )
    provenance = dict(
        PROVENANCE,
        script_sha256="script-sha",
        input_receipt_sha256="input-sha",
        environment_receipt_sha256="environment-sha",
        input_geometry_sha256="geometry-sha",
        requested_settings_sha256="requested-sha",
        applied_settings_sha256=None,
        project_yaml_digest=None,
    )
    _write_artifact(path, "2.0", spec, provenance)

    findings = verify_provenance(
        REQUESTED_SETTINGS,
        path,
        expected_receipt={
            "script_sha256": "script-sha",
            "input_receipt_sha256": "input-sha",
            "environment_receipt_sha256": "environment-sha",
            "input_geometry_sha256": "geometry-sha",
            "requested_settings_sha256": "requested-sha",
            "project_yaml_digest": None,
            "require_applied_settings_sha256": True,
        },
    )

    assert any(
        finding.rule_id
        == "pyscf.provenance.applied_settings_digest_mismatch"
        for finding in findings
    )
