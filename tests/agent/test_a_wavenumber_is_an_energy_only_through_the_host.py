"""A wavenumber becomes an energy through one operation the host owns.

NOVEL-3 ino2 (2026-09-05): a session delivered both exchange couplings
under their declared ids in kcal/mol against declarations in cm^-1. The
completion gate returned the goal for a claim "of matching dimension"
without saying which dimension it had found, and the expectation row
printed ``not_comparable`` on a sixfold-consistent ferromagnetic sign
the declaration had predicted antiferromagnetic. The session asked four
times to convert hartree to cm^-1 and was refused each time, because
h*c*N_A lived only inside the zero-point kernel.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from chemsmart.analysis.quantity_expressions import (
    _OPERATIONS,
    CONVENTION_OPERATIONS,
    OPERATION_DESCRIPTIONS,
    OPERATION_INPUT_COUNTS,
    QuantityExpressionError,
    QuantityExpressionNodeV1,
    QuantityExpressionRequestV1,
    evaluate_quantity_expression,
)
from chemsmart.analysis.result_quantities import (
    ENERGY,
    FREQUENCY,
    make_quantity_value,
)

_HARTREE_IN_CM1 = 219474.6313705


def _value(quantity_id, value, unit, dimension):
    return make_quantity_value(
        quantity_id=quantity_id,
        source_value=value,
        source_unit=unit,
        value=value,
        unit=unit,
        dimension=dimension,
        evidence_ref=f"artifact:x#{quantity_id}",
    )


def _evaluate(inputs, node):
    return evaluate_quantity_expression(
        QuantityExpressionRequestV1(
            schema_version="chemsmart.quantity-expression-request.v1",
            expression_id="j",
            inputs=tuple(inputs),
            nodes=(node,),
            output_node_ids=(node.node_id,),
        )
    ).outputs[0]


@pytest.mark.capability("operation:wavenumber_to_energy")
@pytest.mark.capability("operation:energy_to_wavenumber")
def test_both_conversions_are_registered_as_conventions():
    for name in ("wavenumber_to_energy", "energy_to_wavenumber"):
        assert name in _OPERATIONS
        assert name in OPERATION_DESCRIPTIONS
        assert name in CONVENTION_OPERATIONS
        assert OPERATION_INPUT_COUNTS[name] == frozenset({1})


@pytest.mark.capability("operation:wavenumber_to_energy")
def test_a_wavenumber_restates_as_an_energy_in_the_unit_asked():
    """8065.544 cm^-1 is one electronvolt; the target unit is honoured
    and the canonical value is hartree."""

    out = _evaluate(
        [_value("j", 8065.544, "cm^-1", FREQUENCY)],
        QuantityExpressionNodeV1(
            node_id="j-ev",
            operation="wavenumber_to_energy",
            input_ids=("j",),
            target_unit="eV",
        ),
    )
    assert out.dimension == ENERGY
    assert out.source_unit == "eV"
    assert out.source_value == pytest.approx(1.0, abs=1e-4)
    assert out.unit == "hartree"
    assert out.value == pytest.approx(8065.544 / _HARTREE_IN_CM1, rel=1e-9)

    default = _evaluate(
        [_value("j", 1000.0, "cm^-1", FREQUENCY)],
        QuantityExpressionNodeV1(
            node_id="j-kj", operation="wavenumber_to_energy", input_ids=("j",)
        ),
    )
    assert default.source_unit == "kJ/mol"
    assert default.source_value == pytest.approx(11.9627, abs=1e-3)


@pytest.mark.capability("operation:energy_to_wavenumber")
def test_an_energy_restates_as_a_wavenumber_scalar_or_vector():
    out = _evaluate(
        [_value("e", 1.0, "hartree", ENERGY)],
        QuantityExpressionNodeV1(
            node_id="e-cm", operation="energy_to_wavenumber", input_ids=("e",)
        ),
    )
    assert out.dimension == FREQUENCY
    assert out.unit == "cm^-1"
    assert out.value == pytest.approx(_HARTREE_IN_CM1, rel=1e-12)

    vector = _evaluate(
        [_value("gaps", (0.001, 0.002), "hartree", ENERGY)],
        QuantityExpressionNodeV1(
            node_id="gaps-cm",
            operation="energy_to_wavenumber",
            input_ids=("gaps",),
        ),
    )
    assert vector.data_kind == "vector"
    assert list(vector.value) == pytest.approx(
        [0.001 * _HARTREE_IN_CM1, 0.002 * _HARTREE_IN_CM1], rel=1e-12
    )


def test_the_wrong_dimension_is_refused_never_silently_converted():
    with pytest.raises(QuantityExpressionError, match="takes a wavenumber"):
        _evaluate(
            [_value("e", 1.0, "hartree", ENERGY)],
            QuantityExpressionNodeV1(
                node_id="x",
                operation="wavenumber_to_energy",
                input_ids=("e",),
            ),
        )
    with pytest.raises(QuantityExpressionError, match="takes a molar energy"):
        _evaluate(
            [_value("j", 1.0, "cm^-1", FREQUENCY)],
            QuantityExpressionNodeV1(
                node_id="x",
                operation="energy_to_wavenumber",
                input_ids=("j",),
            ),
        )


_J = {
    "observable_id": "j_ohcl",
    "unit": "cm^-1",
    "dimension": FREQUENCY,
    "meaning": "exchange coupling J of the OH/Cl dimer",
    "expectation_basis": "superexchange through hydroxide",
    "expected_sign": "negative",
    "expected_low": -60.0,
    "expected_high": -1.0,
}


def _host(tmp_path):
    from tests.agent.test_a_guide_opens_when_something_asks import _host

    return _host(tmp_path, approved_requested_observable_declarations=[_J])


def _deliver(host, claim):
    host.analysis_claim_records["r1"] = SimpleNamespace(
        task_spec_sha256="a" * 64, claims=(claim,)
    )


@pytest.mark.capability("rule:plan.claim_carries_declared_id")
def test_an_id_joined_claim_of_the_wrong_dimension_names_the_route(
    tmp_path,
):
    """ino2's shape: J under its declared id, in kcal/mol against cm^-1.
    The miss names both dimensions and the operation that relates them;
    the expectation row compares the sign the unit cannot change and the
    band after the host's own conversion, and says diverged, because a
    ferromagnetic +4.8 kcal/mol is +1688 cm^-1 against a declared
    antiferromagnetic band."""

    host = _host(tmp_path)
    _deliver(
        host,
        SimpleNamespace(
            claim_id="j_ohcl",
            dimension=ENERGY,
            display_value=4.826,
            display_unit="kcal/mol",
        ),
    )
    misses, limitations = host._declared_observable_completion(
        task_spec_sha256="a" * 64
    )
    assert limitations == ("declared_observable:j_ohcl",)
    (miss,) = misses
    assert "'hartree'" in miss and "'cm^-1'" in miss
    assert "wavenumber_to_energy / energy_to_wavenumber" in miss

    (row,) = host._declared_observable_predictions(task_spec_sha256="a" * 64)
    assert row["delivered_claim_id"] == "j_ohcl"
    assert row["delivered_value_in_declared_unit"] == pytest.approx(
        4.826 * 349.755, rel=1e-3
    )
    assert row["agreement"] == "diverged"
    assert "band_untestable" not in row


def test_the_same_claim_restated_in_the_declared_unit_certifies(tmp_path):
    host = _host(tmp_path)
    _deliver(
        host,
        SimpleNamespace(
            claim_id="j_ohcl",
            dimension=FREQUENCY,
            display_value=-12.0,
            display_unit="cm^-1",
        ),
    )
    assert host._declared_observable_completion(task_spec_sha256="a" * 64) == (
        (),
        (),
    )
    (row,) = host._declared_observable_predictions(task_spec_sha256="a" * 64)
    assert row["agreement"] == "agreed"


def test_a_unit_the_host_cannot_relate_leaves_the_band_untestable(tmp_path):
    """A length against an energy has no conversion; the sign is still
    compared, the band is named untestable rather than silently agreed."""

    host = _host(tmp_path)
    _deliver(
        host,
        SimpleNamespace(
            claim_id="j_ohcl",
            dimension=(0, 1, 0, 0, 0, 0),
            display_value=-2.1,
            display_unit="angstrom",
        ),
    )
    (row,) = host._declared_observable_predictions(task_spec_sha256="a" * 64)
    assert "no conversion the host owns" in row["band_untestable"]
    assert row["agreement"] == "not_comparable"

    _deliver(
        host,
        SimpleNamespace(
            claim_id="j_ohcl",
            dimension=(0, 1, 0, 0, 0, 0),
            display_value=+2.1,
            display_unit="angstrom",
        ),
    )
    (row,) = host._declared_observable_predictions(task_spec_sha256="a" * 64)
    assert row["agreement"] == "diverged", "a wrong sign needs no band"
