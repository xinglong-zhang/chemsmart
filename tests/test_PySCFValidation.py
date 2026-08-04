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
    _independent_mass_weighted_frequencies,
    frequency_validation_receipt,
    preflight,
    validate_pyscf_result,
    verify_provenance,
)
from chemsmart.jobs.pyscf.writer import (
    APPLIED_SPEC_FIELDS,
    RESULT_CONTRACT_VERSION,
    PySCFScriptWriter,
    pyscf_reference_family,
    write_pyscf_h5,
)


class CallableProbe:
    """A callable whose invocation would prove that preflight ran work."""

    def __init__(self):
        self.calls = 0

    def __call__(self, *args, **kwargs):
        del args, kwargs
        self.calls += 1
        raise AssertionError("preflight must not invoke solver entry points")


def test_host_reference_derivation_preserves_one_electron_hf_semantics():
    assert pyscf_reference_family(
        symbols=("H",),
        charge=0,
        multiplicity=2,
        xc=None,
    ) == "rohf"
    assert pyscf_reference_family(
        symbols=("O", "H", "H"),
        charge=0,
        multiplicity=1,
        xc="b3lypg",
    ) == "rks"


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
        "dispersion_conformance": {
            "schema_version": (
                "chemsmart.pyscf-dispersion-conformance.v1"
            ),
            "requested_method": "pbe",
            "requested_literal": "d4",
            "parsed_method": "pbe",
            "dispersion_version": "d4",
            "with_3body": True,
            "supported": True,
            "method_compatible": True,
            "status": "supported",
        },
        "solvent_ids": {"water", "methanol"},
    }

    assert preflight(settings, water, environment) == []
    assert solver.calls == 0


@pytest.mark.parametrize(
    ("charge", "multiplicity"),
    ((0, None), (None, 1)),
)
def test_preflight_rejects_partial_electronic_state_override(
    water, charge, multiplicity
):
    settings = {
        "functional": "pbe",
        "basis": "def2-svp",
        "charge": charge,
        "multiplicity": multiplicity,
        "jobtype": "sp",
        "engine": "cpu",
    }

    violations = preflight(settings, water, {})

    paired = [
        item
        for item in violations
        if item.rule_id == "pyscf.electrons.override_pair"
    ]
    assert len(paired) == 1
    assert paired[0].field == "charge_multiplicity"


def test_preflight_mapping_uses_canonical_density_fit_default(water):
    settings = {
        "functional": "pbe",
        "basis": "def2-svp",
        "charge": 0,
        "multiplicity": 1,
        "jobtype": "sp",
        "engine": "cpu",
        "aux_basis": "def2-universal-jkfit",
    }

    violations = preflight(settings, water, {})

    aux = [item for item in violations if item.field == "aux_basis"]
    assert len(aux) == 1
    assert aux[0].expected == "unset when density_fit is disabled"


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


@pytest.mark.parametrize("response_method", ["tda", "tddft"])
def test_td_preflight_accepts_only_bounded_cpu_preview(
    water, response_method
):
    settings = {
        "functional": "b3lyp",
        "basis": "def2-svp",
        "charge": 0,
        "multiplicity": 1,
        "jobtype": "td",
        "engine": "cpu",
        "response_method": response_method,
        "state_manifold": "singlet",
        "nstates": 5,
    }

    findings = preflight(settings, water, {})

    assert [finding.rule_id for finding in findings] == [
        "pyscf.td.preview_only_capability"
    ]


@pytest.mark.parametrize(
    ("change", "field"),
    [
        ({"response_method": "spin_flip"}, "response_method"),
        ({"state_manifold": "triplet"}, "state_manifold"),
        ({"nstates": 0}, "nstates"),
        ({"multiplicity": 2}, "multiplicity"),
        ({"solvent_model": "smd", "solvent_id": "water"}, "solvent_model"),
        ({"engine": "gpu"}, "engine"),
    ],
)
def test_td_preflight_rejects_out_of_scope_response_requests(
    water, change, field
):
    settings = {
        "functional": "b3lyp",
        "basis": "def2-svp",
        "charge": 0,
        "multiplicity": 1,
        "jobtype": "td",
        "engine": "cpu",
        "response_method": "tda",
        "state_manifold": "singlet",
        "nstates": 5,
    }
    settings.update(change)

    findings = preflight(settings, water, {})

    assert any(
        finding.rule_id == "pyscf.settings.invalid_value"
        and finding.field == field
        for finding in findings
    ) or any(
        finding.rule_id == "pyscf.td.gpu_preview_unsupported"
        and finding.field == field
        for finding in findings
    )


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
    assert receipt["stationary_point_classification"] == "unclassified"
    assert receipt["observed_imaginary_mode_count"] is None
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
        expected_imaginary_modes=0,
        imaginary_mode_cutoff_cm1=0.0,
    )

    rule_ids = {finding.rule_id for finding in receipt["findings"]}

    assert receipt["state"] == "failed"
    assert "pyscf.frequency.mode_count" in rule_ids
    assert "pyscf.frequency.nonfinite" in rule_ids
    assert "pyscf.frequency.imaginary_modes" in rule_ids


def _complete_result_spec(jobtype):
    stages = {
        "sp": ["scf"],
        "opt": ["scf", "opt"],
        "hess": ["scf", "hess"],
    }[jobtype]
    spec = {field: None for field in APPLIED_SPEC_FIELDS}
    spec.update(
        {
            "program": "pyscf",
            "result_contract_version": RESULT_CONTRACT_VERSION,
            "jobtype": jobtype,
            "engine": "cpu",
            "stages": stages,
            "symbols": ["O", "H", "H"],
            "positions": [
                [0.0, 0.0, 0.1174],
                [-0.757, 0.0, -0.4696],
                [0.757, 0.0, -0.4696],
            ],
            "unit": "Angstrom",
            "xc": "b3lypg",
            "reference_family": "rks",
            "method": "b3lyp",
            "basis": "def2-svp",
            "charge": 0,
            "spin": 0,
            "multiplicity": 1,
            "density_fit": False,
            "defgrid": "defgrid2",
            "scf_tol": 1.0e-9,
            "scf_maxiter": 100,
            "materializations": {
                "functional_definition": {
                    "schema_version": (
                        "chemsmart.pyscf-functional-definition.v3"
                    ),
                    "field": "xc",
                    "source": "pyscf.dft.libxc.parse_xc",
                    "source_key": "b3lypg",
                    "pyscf_version": "2.14.0",
                    "libxc_version": "7.0.0",
                    "environment_receipt_sha256": "e" * 64,
                    "parser_hybrid_decomposition": [0.0, 0.0, 0.0],
                    "exact_exchange_fraction": 0.2,
                    "range_separation_coefficients": [0.0, 0.2, 0.0],
                    "functional_ids": [402],
                    "functional_factors": [1.0],
                }
            },
            "requested_settings_sha256": "a" * 64,
            "settings_digest": "b" * 64,
        }
    )
    if jobtype == "opt":
        spec["opt_solver"] = "geometric"
        spec["opt_maxsteps"] = 100
    spec["applied_settings_sha256"] = PySCFScriptWriter.settings_digest(
        {field: spec.get(field) for field in APPLIED_SPEC_FIELDS}
    )
    spec["num_electrons"] = 10
    spec["nelec"] = [5, 5]
    return spec


def _write_complete_result(
    path,
    *,
    jobtype="sp",
    results=None,
    spec_updates=None,
    provenance_updates=None,
    status_override=None,
):
    spec = _complete_result_spec(jobtype)
    if spec_updates:
        spec.update(spec_updates)
        if "reference_family" not in spec_updates and (
            "multiplicity" in spec_updates
            or "xc" in spec_updates
            or "method" in spec_updates
        ):
            if spec.get("xc") is not None:
                spec["reference_family"] = (
                    "rks" if int(spec["multiplicity"]) == 1 else "uks"
                )
            elif str(spec.get("method") or "").lower() == "hf":
                electron_count = int(spec.get("num_electrons") or 10)
                if int(spec["multiplicity"]) == 1:
                    spec["reference_family"] = "rhf"
                elif electron_count == 1:
                    spec["reference_family"] = "rohf"
                else:
                    spec["reference_family"] = "uhf"
        spec["applied_settings_sha256"] = PySCFScriptWriter.settings_digest(
            {field: spec.get(field) for field in APPLIED_SPEC_FIELDS}
        )
    stages = {
        stage: {"converged": True} for stage in spec["stages"]
    }
    if "opt" in stages:
        stages["opt"].update(
            {
                "optimizer_converged": True,
                "final_scf_converged": True,
            }
        )
    values = {
        "energies": np.asarray([-76.3]),
        "positions": np.asarray(spec["positions"]),
        "atomic_numbers": np.asarray([8, 1, 1]),
        "mo_energy": np.asarray([-20.0, -2.0, -1.0, -0.8, -0.5, 0.5, 1.0]),
        "mo_occ": np.asarray([2.0, 2.0, 2.0, 2.0, 2.0, 0.0, 0.0]),
        "spin_square": np.asarray(0.0),
        "spin_square_effective_multiplicity": np.asarray(1.0),
    }
    if jobtype == "hess":
        matrix = np.eye(9, dtype=float) * 0.01
        frequencies, _ = _independent_mass_weighted_frequencies(
            matrix=matrix,
            symbols=tuple(spec["symbols"]),
            positions=np.asarray(spec["positions"], dtype=float),
        )
        values.update(
            {
                "hessian": matrix.reshape(3, 3, 3, 3).transpose(0, 2, 1, 3),
                "vibrational_frequencies": frequencies,
                "normal_modes": np.zeros((3, 3, 3)),
                "reduced_masses": np.ones(3),
                "force_constants": np.ones(3),
            }
        )
    if results:
        values.update(results)
    reference_class = {
        "rks": "pyscf.dft.rks.RKS",
        "uks": "pyscf.dft.uks.UKS",
        "rhf": "pyscf.scf.hf.RHF",
        "rohf": "pyscf.scf.rohf.HF1e",
        "uhf": "pyscf.scf.uhf.UHF",
    }[spec["reference_family"]]
    provenance = {
        "engine": "cpu",
        "pyscf_version": "2.14.0",
        "libxc_version": "7.0.0",
        "environment_receipt_sha256": "e" * 64,
        "settings_digest": spec["settings_digest"],
        "requested_settings_sha256": spec[
            "requested_settings_sha256"
        ],
        "applied_settings_sha256": spec["applied_settings_sha256"],
        "runtime": {
            "mean_field_class": reference_class,
            "mean_field_mro": [
                reference_class,
                "object",
            ],
        },
    }
    if provenance_updates:
        provenance.update(provenance_updates)
    status = (
        copy.deepcopy(status_override)
        if status_override is not None
        else {
            "engine_complete": True,
            "normal_termination": True,
            "stages": stages,
            "failure": None,
        }
    )
    status.setdefault(
        "properties", {"spin_square": {"status": "ok"}}
    )
    write_pyscf_h5(
        path,
        spec=spec,
        provenance=provenance,
        status=status,
        results=values,
    )
    return spec


def _result_expectation(jobtype):
    return {
        "functional": "b3lyp",
        "basis": "def2-svp",
        "charge": 0,
        "multiplicity": 1,
        "engine": "cpu",
        "density_fit": False,
        "defgrid": "defgrid2",
        "scf_tol": 1.0e-9,
        "scf_maxiter": 100,
        "jobtype": jobtype,
    }


def _validate_complete_result(
    path,
    *,
    jobtype="sp",
    settings=None,
    charge=0,
    multiplicity=1,
    symbols=("O", "H", "H"),
    expected_positions=None,
    stationary_point_policy=None,
):
    if expected_positions is None:
        expected_positions = _complete_result_spec(jobtype)["positions"]
    return validate_pyscf_result(
        path,
        settings=settings or _result_expectation(jobtype),
        expected_jobtype=jobtype,
        expected_charge=charge,
        expected_multiplicity=multiplicity,
        expected_symbols=symbols,
        expected_positions=expected_positions,
        expected_receipt={"require_applied_settings_sha256": True},
        stationary_point_policy=stationary_point_policy,
    )


def test_shared_result_validator_accepts_complete_sp(tmp_path):
    path = tmp_path / "water-sp.h5"
    _write_complete_result(path)

    receipt = validate_pyscf_result(
        path,
        settings=_result_expectation("sp"),
        expected_jobtype="sp",
        expected_charge=0,
        expected_multiplicity=1,
        expected_symbols=("O", "H", "H"),
        expected_positions=_complete_result_spec("sp")["positions"],
        expected_receipt={"require_applied_settings_sha256": True},
    )

    assert receipt["state"] == "validated"
    assert receipt["contract_validation"]["new_execution_admissible"] is True
    assert receipt["findings"] == []


def test_shared_result_validator_rejects_nonfinite_and_atom_reordering(
    tmp_path,
):
    path = tmp_path / "water-invalid.h5"
    _write_complete_result(path, results={"energies": [float("nan")]})

    receipt = validate_pyscf_result(
        path,
        settings=_result_expectation("sp"),
        expected_jobtype="sp",
        expected_charge=0,
        expected_multiplicity=1,
        expected_symbols=("H", "O", "H"),
        expected_positions=_complete_result_spec("sp")["positions"],
        expected_receipt={"require_applied_settings_sha256": True},
    )
    rules = {finding.rule_id for finding in receipt["findings"]}

    assert receipt["state"] == "failed"
    assert "pyscf.result.nonfinite" in rules
    assert "pyscf.result.atom_identity_mismatch" in rules


def test_shared_hessian_validator_uses_approved_imaginary_cutoff(tmp_path):
    path = tmp_path / "water-hess.h5"
    _write_complete_result(path, jobtype="hess")
    policy = {
        "policy_sha256": "c" * 64,
        "expected_imaginary_mode_count": 0,
        "imaginary_mode_cutoff_cm1": 20.0,
        "require_finite_modes": True,
        "require_symmetric_hessian": True,
    }

    receipt = validate_pyscf_result(
        path,
        settings=_result_expectation("hess"),
        expected_jobtype="hess",
        expected_charge=0,
        expected_multiplicity=1,
        expected_symbols=("O", "H", "H"),
        expected_positions=_complete_result_spec("hess")["positions"],
        expected_receipt={"require_applied_settings_sha256": True},
        stationary_point_policy=policy,
    )

    assert receipt["state"] == "validated"
    assert receipt["frequency_validation"][
        "stationary_point_classification"
    ] == "classified"
    assert receipt["frequency_validation"][
        "observed_imaginary_mode_count"
    ] == 0
    assert receipt["hessian_validation"]["symmetric"] is True


def test_hessian_without_policy_is_intact_but_not_scientifically_validated(
    tmp_path,
):
    path = tmp_path / "water-hess-unclassified.h5"
    _write_complete_result(path, jobtype="hess")

    receipt = _validate_complete_result(path, jobtype="hess")

    assert receipt["state"] == "unclassified"
    assert receipt["stationary_point_policy"]["classification_state"] == (
        "unclassified"
    )
    assert receipt["frequency_validation"][
        "stationary_point_classification"
    ] == "unclassified"
    assert receipt["hessian_validation"]["consistency"]["state"] == (
        "verified"
    )


@pytest.mark.parametrize(
    ("replacement", "expected_rule"),
    [
        (np.asarray([b"-76.3"]), "pyscf.result.nonnumeric_physical_array"),
        (np.asarray([True]), "pyscf.result.nonnumeric_physical_array"),
    ],
)
def test_result_validator_rejects_nonnumeric_native_physical_arrays(
    tmp_path, replacement, expected_rule
):
    path = tmp_path / "water-invalid-dtype.h5"
    _write_complete_result(path)
    with h5py.File(path, "a") as handle:
        del handle["results/energies"]
        dataset = handle["results"].create_dataset(
            "energies", data=replacement
        )
        dataset.attrs["unit"] = "Eh"

    receipt = _validate_complete_result(path)

    assert expected_rule in {finding.rule_id for finding in receipt["findings"]}


def test_result_validator_requires_units_on_physical_arrays(tmp_path):
    path = tmp_path / "water-missing-unit.h5"
    _write_complete_result(path)
    with h5py.File(path, "a") as handle:
        del handle["results/positions"].attrs["unit"]

    receipt = _validate_complete_result(path)

    assert "pyscf.result.unit_missing_or_mismatched" in {
        finding.rule_id for finding in receipt["findings"]
    }


def test_result_validator_accepts_unrestricted_alpha_beta_orbitals(tmp_path):
    path = tmp_path / "water-triplet.h5"
    mo_energy = np.asarray(
        [
            [-20.0, -2.0, -1.0, -0.8, -0.5, -0.2, 1.0],
            [-20.0, -2.0, -1.0, -0.8, 0.2, 0.5, 1.0],
        ]
    )
    mo_occ = np.asarray(
        [
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0],
            [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
        ]
    )
    _write_complete_result(
        path,
        spec_updates={"spin": 2, "multiplicity": 3, "nelec": [6, 4]},
        results={
            "mo_energy": mo_energy,
            "mo_occ": mo_occ,
            "spin_square": np.asarray(2.0),
            "spin_square_effective_multiplicity": np.asarray(3.0),
        },
    )
    settings = dict(_result_expectation("sp"), multiplicity=3)

    receipt = _validate_complete_result(
        path,
        settings=settings,
        multiplicity=3,
    )

    assert receipt["state"] == "validated"
    assert receipt["electronic_state_validation"][
        "orbital_representation"
    ] == "unrestricted"
    assert receipt["electronic_state_validation"][
        "observed_alpha_beta"
    ] == (6.0, 4.0)
    assert receipt["spin_diagnostic_validation"]["state"] == "available"
    assert receipt["spin_diagnostic_validation"]["observed_spin_square"] == 2.0


@pytest.mark.parametrize(
    ("spec_updates", "provenance_updates"),
    [
        ({"reference_family": "uks"}, None),
        (
            {},
            {
                "runtime": {
                    "mean_field_class": "pyscf.dft.uks.UKS",
                    "mean_field_mro": ["pyscf.dft.uks.UKS", "object"],
                }
            },
        ),
    ],
)
def test_result_validator_rejects_reference_substitution(
    tmp_path, spec_updates, provenance_updates
):
    path = tmp_path / "water-reference-substitution.h5"
    _write_complete_result(
        path,
        spec_updates=spec_updates,
        provenance_updates=provenance_updates,
    )

    receipt = _validate_complete_result(path)

    assert receipt["state"] == "failed"
    assert "pyscf.result.reference_family_mismatch" in {
        finding.rule_id for finding in receipt["findings"]
    }


def test_one_electron_hf_requires_hf1e_rohf_runtime_semantics(tmp_path):
    path = tmp_path / "hydrogen-doublet.h5"
    spec = _complete_result_spec("sp")
    spec.update(
        {
            "symbols": ["H"],
            "positions": [[0.0, 0.0, 0.0]],
            "xc": None,
            "ab_initio": "hf",
            "method": "hf",
            "reference_family": "rohf",
            "charge": 0,
            "multiplicity": 2,
            "spin": 1,
            "materializations": {},
            "num_electrons": 1,
            "nelec": [1, 0],
        }
    )
    spec["applied_settings_sha256"] = PySCFScriptWriter.settings_digest(
        {field: spec.get(field) for field in APPLIED_SPEC_FIELDS}
    )
    provenance = {
        "engine": "cpu",
        "pyscf_version": "2.14.0",
        "libxc_version": "7.0.0",
        "environment_receipt_sha256": "e" * 64,
        "settings_digest": spec["settings_digest"],
        "requested_settings_sha256": spec["requested_settings_sha256"],
        "applied_settings_sha256": spec["applied_settings_sha256"],
        "runtime": {
            "mean_field_class": "pyscf.scf.rohf.HF1e",
            "mean_field_mro": ["pyscf.scf.rohf.HF1e", "object"],
        },
    }
    write_pyscf_h5(
        path,
        spec=spec,
        provenance=provenance,
        status={
            "engine_complete": True,
            "normal_termination": True,
            "stages": {"scf": {"converged": True}},
            "failure": None,
            "properties": {"spin_square": {"status": "ok"}},
        },
        results={
            "energies": np.asarray([-0.5]),
            "positions": np.asarray([[0.0, 0.0, 0.0]]),
            "atomic_numbers": np.asarray([1]),
            "mo_energy": np.asarray([-0.5]),
            "mo_occ": np.asarray([1.0]),
            "spin_square": np.asarray(0.75),
            "spin_square_effective_multiplicity": np.asarray(2.0),
        },
    )

    receipt = validate_pyscf_result(
        path,
        settings={
            "jobtype": "sp",
            "engine": "cpu",
            "ab_initio": "hf",
            "basis": "def2-svp",
            "charge": 0,
            "multiplicity": 2,
        },
        expected_jobtype="sp",
        expected_charge=0,
        expected_multiplicity=2,
        expected_symbols=("H",),
        expected_positions=((0.0, 0.0, 0.0),),
        expected_receipt={"require_applied_settings_sha256": True},
    )

    assert receipt["state"] == "validated"
    assert receipt["reference_validation"]["runtime_family"] == "rohf"
    assert receipt["reference_validation"]["one_electron_semantics"] is True


def test_spin_contamination_is_advisory_without_universal_cutoff(tmp_path):
    path = tmp_path / "water-triplet-contaminated.h5"
    effective = float(np.sqrt(1.0 + 4.0 * 2.3))
    _write_complete_result(
        path,
        spec_updates={"spin": 2, "multiplicity": 3, "nelec": [6, 4]},
        results={
            "mo_energy": np.asarray(
                [
                    [-20.0, -2.0, -1.0, -0.8, -0.5, -0.2, 1.0],
                    [-20.0, -2.0, -1.0, -0.8, 0.2, 0.5, 1.0],
                ]
            ),
            "mo_occ": np.asarray(
                [
                    [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0],
                    [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
                ]
            ),
            "spin_square": np.asarray(2.3),
            "spin_square_effective_multiplicity": np.asarray(effective),
        },
    )

    receipt = _validate_complete_result(
        path,
        settings=dict(_result_expectation("sp"), multiplicity=3),
        multiplicity=3,
    )

    assert receipt["state"] == "validated"
    assert receipt["spin_diagnostic_validation"]["interpretation"] == (
        "deviation_observed_without_universal_threshold"
    )
    assert receipt["spin_diagnostic_validation"][
        "scientific_threshold_applied"
    ] is False
    assert {item["rule_id"] for item in receipt["advisories"]} >= {
        "pyscf.spin_diagnostic.deviation_observed"
    }


@pytest.mark.parametrize(
    ("spin_square", "effective_multiplicity"),
    [
        (float("nan"), 3.0),
        (2.0, 4.0),
        (-0.5, 1.0),
    ],
)
def test_spin_diagnostic_artifact_defects_fail(
    tmp_path, spin_square, effective_multiplicity
):
    path = tmp_path / "water-bad-spin-diagnostic.h5"
    _write_complete_result(
        path,
        results={
            "spin_square": np.asarray(spin_square),
            "spin_square_effective_multiplicity": np.asarray(
                effective_multiplicity
            ),
        },
    )

    receipt = _validate_complete_result(path)

    assert "pyscf.result.spin_diagnostic_invalid" in {
        finding.rule_id for finding in receipt["findings"]
    }


def test_current_open_shell_result_requires_spin_diagnostic_datasets(tmp_path):
    path = tmp_path / "water-triplet-missing-spin.h5"
    _write_complete_result(
        path,
        spec_updates={"spin": 2, "multiplicity": 3, "nelec": [6, 4]},
        results={
            "mo_energy": np.asarray(
                [
                    [-20.0, -2.0, -1.0, -0.8, -0.5, -0.2, 1.0],
                    [-20.0, -2.0, -1.0, -0.8, 0.2, 0.5, 1.0],
                ]
            ),
            "mo_occ": np.asarray(
                [
                    [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0],
                    [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
                ]
            ),
        },
        status_override={
            "engine_complete": True,
            "normal_termination": True,
            "stages": {"scf": {"converged": True}},
            "failure": None,
            "properties": {
                "spin_square": {
                    "status": "unavailable",
                    "failure": {"type": "NotImplementedError"},
                }
            },
        },
    )
    with h5py.File(path, "a") as handle:
        del handle["results/spin_square"]
        del handle["results/spin_square_effective_multiplicity"]

    receipt = _validate_complete_result(
        path,
        settings=dict(_result_expectation("sp"), multiplicity=3),
        multiplicity=3,
    )

    assert receipt["state"] == "failed"
    assert receipt["spin_diagnostic_validation"]["state"] == "unavailable"
    assert "pyscf.result.spin_diagnostic_invalid" in {
        finding.rule_id for finding in receipt["findings"]
    }
    assert {item["rule_id"] for item in receipt["advisories"]} >= {
        "pyscf.spin_diagnostic.unavailable"
    }


def test_result_validator_rejects_occupation_total_and_spin_mismatch(
    tmp_path,
):
    path = tmp_path / "water-bad-occupations.h5"
    _write_complete_result(
        path,
        results={"mo_occ": np.asarray([2.0, 2.0, 2.0, 2.0, 3.0, 0.0, 0.0])},
    )

    receipt = _validate_complete_result(path)
    fields = {finding.field for finding in receipt["findings"]}

    assert "results.mo_occ" in fields
    assert "results.mo_occ.electron_total" in fields


@pytest.mark.parametrize(
    "status",
    [
        {
            "engine_complete": True,
            "normal_termination": True,
            "stages": {"scf": {"converged": True}, "extra": {"converged": True}},
            "failure": None,
        },
        {
            "engine_complete": True,
            "normal_termination": True,
            "stages": {"scf": {"converged": True}},
            "failure": {"type": "RuntimeError", "message": "failed"},
        },
    ],
)
def test_result_validator_requires_exact_status_and_no_failure(
    tmp_path, status
):
    path = tmp_path / "water-status-invalid.h5"
    _write_complete_result(path, status_override=status)

    receipt = _validate_complete_result(path)

    assert receipt["state"] == "failed"
    assert "pyscf.result.stage_mismatch" in {
        finding.rule_id for finding in receipt["findings"]
    }


@pytest.mark.parametrize(
    "missing_path",
    [
        "status/failure",
        "status/properties",
        "status/properties/spin_square",
    ],
)
def test_current_result_contract_rejects_missing_status_fields(
    tmp_path, missing_path
):
    path = tmp_path / "water-missing-current-status.h5"
    _write_complete_result(path)
    with h5py.File(path, "a") as handle:
        del handle[missing_path]

    receipt = _validate_complete_result(path)

    assert receipt["state"] == "failed"
    assert receipt["contract_validation"]["new_execution_admissible"] is False
    assert "pyscf.result.contract_incomplete" in {
        finding.rule_id for finding in receipt["findings"]
    }


@pytest.mark.parametrize(
    "spec_updates",
    [{"charge": False}, {"multiplicity": True}],
)
def test_result_validator_rejects_boolean_electronic_state(
    tmp_path, spec_updates
):
    path = tmp_path / "water-boolean-state.h5"
    _write_complete_result(path, spec_updates=spec_updates)

    receipt = _validate_complete_result(path)

    assert receipt["state"] == "failed"
    assert "pyscf.result.electronic_state_mismatch" in {
        finding.rule_id for finding in receipt["findings"]
    }


def test_opt_requires_explicit_optimizer_and_final_scf_convergence(tmp_path):
    path = tmp_path / "water-opt-status-invalid.h5"
    _write_complete_result(
        path,
        jobtype="opt",
        status_override={
            "engine_complete": True,
            "normal_termination": True,
            "stages": {
                "scf": {"converged": True},
                "opt": {"converged": True},
            },
            "failure": None,
        },
    )

    receipt = _validate_complete_result(path, jobtype="opt")
    fields = {finding.field for finding in receipt["findings"]}

    assert "status.stages.opt.optimizer_converged" in fields
    assert "status.stages.opt.final_scf_converged" in fields
    assert "status.stages.opt.converged" in fields


def test_dft_requires_functional_materialization_and_hf_forbids_it(tmp_path):
    dft_path = tmp_path / "water-dft-no-materialization.h5"
    _write_complete_result(dft_path)
    with h5py.File(dft_path, "a") as handle:
        del handle["spec/materializations/functional_definition"]
    dft_receipt = _validate_complete_result(dft_path)

    assert "pyscf.provenance.materialization_invalid" in {
        finding.rule_id for finding in dft_receipt["findings"]
    }

    hf_path = tmp_path / "water-hf.h5"
    _write_complete_result(
        hf_path,
        spec_updates={
            "xc": None,
            "ab_initio": "hf",
            "method": "hf",
            "defgrid": None,
            "materializations": {},
        },
    )
    hf_settings = {
        "ab_initio": "hf",
        "basis": "def2-svp",
        "charge": 0,
        "multiplicity": 1,
        "engine": "cpu",
        "density_fit": False,
        "scf_tol": 1.0e-9,
        "scf_maxiter": 100,
        "jobtype": "sp",
    }
    hf_receipt = _validate_complete_result(hf_path, settings=hf_settings)

    assert hf_receipt["state"] == "validated"

    hf_bad_path = tmp_path / "water-hf-with-dft-materialization.h5"
    _write_complete_result(
        hf_bad_path,
        spec_updates={
            "xc": None,
            "ab_initio": "hf",
            "method": "hf",
            "defgrid": None,
        },
    )
    hf_bad_receipt = _validate_complete_result(
        hf_bad_path, settings=hf_settings
    )
    assert "pyscf.provenance.materialization_invalid" in {
        finding.rule_id for finding in hf_bad_receipt["findings"]
    }


def test_solvent_materialization_binds_target_environment_and_pyscf_version(
    tmp_path,
):
    environment_sha256 = "e" * 64
    base = _complete_result_spec("sp")
    materializations = dict(base["materializations"])
    materializations["solvent_dielectric"] = {
        "schema_version": "chemsmart.pyscf-solvent-materialization.v1",
        "field": "solvent_eps",
        "source": "pyscf.solvent.smd.solvent_db",
        "source_key": "water",
        "value": 78.3553,
        "unit": "dimensionless_relative_permittivity",
        "pyscf_version": "2.14.0",
        "environment_receipt_sha256": environment_sha256,
    }
    path = tmp_path / "water-smd.h5"
    _write_complete_result(
        path,
        spec_updates={
            "solvent_model": "smd",
            "solvent_call": "SMD",
            "solvent_method": "SMD",
            "solvent_id": "water",
            "solvent_eps": 78.3553,
            "materializations": materializations,
        },
        provenance_updates={
            "pyscf_version": "2.14.0",
            "environment_receipt_sha256": environment_sha256,
        },
    )
    settings = {
        **_result_expectation("sp"),
        "solvent_model": "smd",
        "solvent_id": "water",
    }

    receipt = _validate_complete_result(path, settings=settings)

    assert receipt["state"] == "validated"
    assert receipt["findings"] == []

    tampered = tmp_path / "water-smd-other-environment.h5"
    tampered_materializations = dict(materializations)
    tampered_materializations["solvent_dielectric"] = {
        **materializations["solvent_dielectric"],
        "environment_receipt_sha256": "f" * 64,
    }
    _write_complete_result(
        tampered,
        spec_updates={
            "solvent_model": "smd",
            "solvent_call": "SMD",
            "solvent_method": "SMD",
            "solvent_id": "water",
            "solvent_eps": 78.3553,
            "materializations": tampered_materializations,
        },
        provenance_updates={
            "pyscf_version": "2.14.0",
            "environment_receipt_sha256": environment_sha256,
        },
    )

    rejected = _validate_complete_result(tampered, settings=settings)

    assert "pyscf.provenance.materialization_invalid" in {
        finding.rule_id for finding in rejected["findings"]
    }


@pytest.mark.parametrize(
    ("provenance_field", "replacement"),
    [
        ("pyscf_version", "2.13.0"),
        ("libxc_version", "6.2.2"),
        ("environment_receipt_sha256", "f" * 64),
    ],
)
def test_functional_v3_rejects_target_environment_mismatch(
    tmp_path, provenance_field, replacement
):
    path = tmp_path / f"functional-{provenance_field}-mismatch.h5"
    _write_complete_result(
        path,
        provenance_updates={provenance_field: replacement},
    )

    receipt = _validate_complete_result(path)

    assert receipt["materialization_validation"]["functional_definition"][
        "state"
    ] == "invalid"
    assert "pyscf.provenance.materialization_invalid" in {
        finding.rule_id for finding in receipt["findings"]
    }


def test_legacy_functional_v2_is_readable_but_environment_unbound(tmp_path):
    path = tmp_path / "functional-v2-legacy.h5"
    legacy = copy.deepcopy(
        _complete_result_spec("sp")["materializations"]
    )
    definition = legacy["functional_definition"]
    definition["schema_version"] = "chemsmart.pyscf-functional-definition.v2"
    definition.pop("pyscf_version")
    definition.pop("environment_receipt_sha256")
    _write_complete_result(path, spec_updates={"materializations": legacy})

    receipt = _validate_complete_result(path)

    assert receipt["state"] == "qualified_legacy"
    assert receipt["contract_validation"]["new_execution_admissible"] is False
    assert receipt["materialization_validation"]["functional_definition"][
        "state"
    ] == "legacy_environment_unbound"
    assert {item["rule_id"] for item in receipt["advisories"]} >= {
        "pyscf.functional_definition.legacy_environment_unbound"
    }

def test_fixed_geometry_stages_match_input_but_opt_may_move(tmp_path):
    shifted = np.asarray(_complete_result_spec("sp")["positions"], dtype=float)
    shifted[0, 2] += 0.01
    sp_path = tmp_path / "water-sp-shifted.h5"
    _write_complete_result(sp_path, results={"positions": shifted})
    sp_receipt = _validate_complete_result(sp_path)
    assert any(
        finding.field == "results.positions"
        and finding.rule_id == "pyscf.result.geometry_invalid"
        for finding in sp_receipt["findings"]
    )

    opt_path = tmp_path / "water-opt-shifted.h5"
    _write_complete_result(
        opt_path,
        jobtype="opt",
        results={"positions": shifted},
    )
    opt_receipt = _validate_complete_result(opt_path, jobtype="opt")
    assert opt_receipt["state"] == "validated"


def test_result_validator_requires_pyscf_program_identity(tmp_path):
    path = tmp_path / "water-wrong-program.h5"
    _write_complete_result(path, spec_updates={"program": "orca"})

    receipt = _validate_complete_result(path)

    assert any(
        finding.field == "spec.program" for finding in receipt["findings"]
    )


def test_result_validator_binds_exact_expected_applied_digest(tmp_path):
    path = tmp_path / "water-applied-digest.h5"
    _write_complete_result(path)

    receipt = validate_pyscf_result(
        path,
        settings=_result_expectation("sp"),
        expected_jobtype="sp",
        expected_charge=0,
        expected_multiplicity=1,
        expected_symbols=("O", "H", "H"),
        expected_positions=_complete_result_spec("sp")["positions"],
        expected_receipt={
            "require_applied_settings_sha256": True,
            "applied_settings_sha256": "f" * 64,
        },
    )

    mismatches = [
        finding
        for finding in receipt["findings"]
        if finding.rule_id
        == "pyscf.provenance.applied_settings_digest_mismatch"
    ]
    assert {finding.field for finding in mismatches} >= {
        "spec.applied_settings_sha256",
        "provenance.applied_settings_sha256",
    }


def test_hessian_symmetry_and_reported_spectrum_are_unconditional(tmp_path):
    asymmetric_path = tmp_path / "water-asymmetric-hessian.h5"
    asymmetric = np.eye(9, dtype=float) * 0.01
    asymmetric[0, 1] = 0.1
    asymmetric_hessian = asymmetric.reshape(3, 3, 3, 3).transpose(0, 2, 1, 3)
    _write_complete_result(
        asymmetric_path,
        jobtype="hess",
        results={"hessian": asymmetric_hessian},
    )
    asymmetric_receipt = _validate_complete_result(
        asymmetric_path, jobtype="hess"
    )
    assert any(
        finding.field == "results.hessian.symmetry"
        for finding in asymmetric_receipt["findings"]
    )

    contradiction_path = tmp_path / "water-contradictory-hessian.h5"
    _write_complete_result(
        contradiction_path,
        jobtype="hess",
        results={"hessian": np.zeros((3, 3, 3, 3), dtype=float)},
    )
    contradiction_receipt = _validate_complete_result(
        contradiction_path, jobtype="hess"
    )
    assert "pyscf.result.hessian_frequency_inconsistent" in {
        finding.rule_id for finding in contradiction_receipt["findings"]
    }


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
    "xc": "b3lypg",
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
