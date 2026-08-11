"""Compute-free target-PySCF dispersion conformance tests."""

import ast
from types import SimpleNamespace

import pytest

from chemsmart.jobs.pyscf import environment as environment_module
from chemsmart.jobs.pyscf.validation import preflight


def _source_function(source, name):
    tree = ast.parse(source)
    node = next(
        item
        for item in tree.body
        if isinstance(item, ast.FunctionDef) and item.name == name
    )
    module = ast.Module(body=[node], type_ignores=[])
    ast.fix_missing_locations(module)
    scope = {}
    exec(compile(module, f"<{name}-only>", "exec"), scope)
    return scope[name]


@pytest.mark.parametrize(
    ("method", "literal", "status", "compatible"),
    [
        ("pbe", "d4", "supported", True),
        ("pbe", "not-a-dispersion", "invalid", False),
        ("pbe", "d3bj:b3lyp", "incompatible", False),
    ],
)
def test_target_probe_checks_exact_dispersion_literal(
    method, literal, status, compatible
):
    detail = _source_function(
        environment_module._PROBE_SCRIPT, "dispersion_detail"
    )

    observed = detail(method, literal)

    assert observed["requested_method"] == method
    assert observed["requested_literal"] == literal
    assert observed["status"] == status
    assert observed["method_compatible"] is compatible
    assert observed["supported"] is (status != "invalid")


def _settings():
    return {
        "functional": "pbe",
        "basis": "def2-svp",
        "charge": 0,
        "multiplicity": 1,
        "jobtype": "sp",
        "engine": "cpu",
        "dispersion": "d4",
    }


def _water():
    return SimpleNamespace(
        chemical_symbols=["O", "H", "H"], charge=0, multiplicity=1
    )


def _environment(conformance):
    return {
        "dependencies": {"pyscf-dispersion": True},
        "dispersion_conformance": conformance,
    }


def _conformance(**updates):
    value = {
        "schema_version": "chemsmart.pyscf-dispersion-conformance.v1",
        "requested_method": "pbe",
        "requested_literal": "d4",
        "parsed_method": "pbe",
        "dispersion_version": "d4",
        "with_3body": True,
        "supported": True,
        "method_compatible": True,
        "status": "supported",
    }
    value.update(updates)
    return value


@pytest.mark.parametrize(
    ("conformance", "expected_rule"),
    [
        (None, "pyscf.dispersion.literal_unverified"),
        (
            _conformance(
                supported=False,
                method_compatible=False,
                status="invalid",
                error_type="ValueError",
            ),
            "pyscf.dispersion.literal_invalid",
        ),
        (
            _conformance(
                parsed_method="b3lyp",
                method_compatible=False,
                status="incompatible",
            ),
            "pyscf.dispersion.method_incompatible",
        ),
        (
            _conformance(requested_literal="d3bj"),
            "pyscf.dispersion.literal_unverified",
        ),
    ],
)
def test_preflight_fails_closed_for_unbound_dispersion_evidence(
    conformance, expected_rule
):
    findings = preflight(_settings(), _water(), _environment(conformance))

    assert expected_rule in {finding.rule_id for finding in findings}


def test_preflight_accepts_exact_supported_dispersion_evidence():
    findings = preflight(
        _settings(), _water(), _environment(_conformance())
    )

    assert not [
        finding for finding in findings if finding.field == "dispersion"
    ]
