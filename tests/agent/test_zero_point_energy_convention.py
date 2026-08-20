"""A zero-point energy is a convention ChemSmart owns; it must be selectable.

Found by executing a live plan rather than by reading one. Given only a
methodology section and an xyz, a session planned HF/6-31G(d) `opt` then `hess`
correctly and then built its scaled ZPVE out of primitives:

    sum -> scale 0.5 -> convert(cm^-1 -> kJ/mol) -> scale 0.9135

Two things are wrong with that, and both are the harness's fault.

The factor of one half is the definition of a zero-point level, and the model
typed it. `Thermochemistry.zero_point_energy` already owned the quantity, but
no expression operation exposed it, so the only route was arithmetic.

The `convert` node cannot run at all. A wavenumber is dimension FREQUENCY and a
molar energy is dimension ENERGY, so the unit engine refuses to relate them --
correctly, because `h*c*N_A` is a spectroscopic convention and not a unit
identity. The session predicted the refusal in its own report and proposed to
hand-multiply by `h*c*N_A`, which would have moved a physical constant into
model-authored arithmetic. That is the failure this project exists to prevent,
reached by a model that was reasoning carefully.

`harmonic_zero_point_energy` owns both numbers.
"""

import pytest

from chemsmart.analysis.aggregation import (
    AggregationError,
    harmonic_zero_point_energy,
)
from chemsmart.analysis.quantity_expressions import (
    _OPERATIONS,
    CONVENTION_OPERATIONS,
    OPERATION_DESCRIPTIONS,
    QuantityExpressionNodeV1,
    QuantityExpressionRequestV1,
    _unit_spec,
    evaluate_quantity_expression,
)
from chemsmart.analysis.result_quantities import (
    ENERGY,
    FREQUENCY,
    make_quantity_value,
)

#: HF/6-31G(d) water, from ChemSmart's own hess run of the sealed case.
_WATER_HF_631GD = (
    1826.0168416143915,
    4055.51188448563,
    4173.850868259163,
)


def test_the_operation_is_registered_described_and_a_convention():
    assert "harmonic_zero_point_energy" in _OPERATIONS
    assert "harmonic_zero_point_energy" in OPERATION_DESCRIPTIONS
    assert "harmonic_zero_point_energy" in CONVENTION_OPERATIONS


def test_it_reproduces_the_sealed_zero_point_energy():
    """The number a graded case actually needed, to its stated tolerance."""
    harmonic = harmonic_zero_point_energy(_WATER_HF_631GD)
    assert harmonic == pytest.approx(60.145, abs=0.5)
    assert 0.9135 * harmonic == pytest.approx(54.942, abs=0.5)
    assert harmonic_zero_point_energy(
        _WATER_HF_631GD, unit="hartree"
    ) == pytest.approx(0.022907841891569164, abs=1e-8)


def test_expression_keeps_zpe_source_display_but_uses_canonical_energy():
    """A ZPE produced in kJ/mol must add directly to a Hartree energy."""

    frequencies = make_quantity_value(
        quantity_id="frequencies",
        source_value=_WATER_HF_631GD,
        source_unit="cm^-1",
        value=_WATER_HF_631GD,
        unit="cm^-1",
        dimension=FREQUENCY,
        evidence_ref="artifact:water#frequencies",
    )
    energy = make_quantity_value(
        quantity_id="energy",
        source_value=-76.0,
        source_unit="hartree",
        value=-76.0,
        unit="hartree",
        dimension=ENERGY,
        evidence_ref="artifact:water#energy",
    )
    receipt = evaluate_quantity_expression(
        QuantityExpressionRequestV1(
            schema_version="chemsmart.quantity-expression-request.v1",
            expression_id="water-zpe",
            inputs=(frequencies, energy),
            nodes=(
                QuantityExpressionNodeV1(
                    node_id="zpe",
                    operation="harmonic_zero_point_energy",
                    input_ids=("frequencies",),
                ),
                QuantityExpressionNodeV1(
                    node_id="energy-plus-zpe",
                    operation="add",
                    input_ids=("energy", "zpe"),
                ),
            ),
            output_node_ids=("zpe", "energy-plus-zpe"),
        )
    )

    zpe, corrected = receipt.outputs
    assert zpe.source_unit == "kJ/mol"
    assert zpe.source_value == pytest.approx(60.145, abs=0.5)
    assert zpe.unit == "hartree"
    assert zpe.value == pytest.approx(0.022907841891569164, abs=1e-8)
    assert corrected.unit == "hartree"
    assert corrected.value == pytest.approx(-75.97709215810843, abs=1e-8)


def test_the_conversion_it_owns_is_one_the_unit_engine_refuses():
    """This is why the operation exists rather than a `convert` node."""
    frequency_dimension = _unit_spec("cm^-1")[0]
    energy_dimension = _unit_spec("kJ/mol")[0]
    assert frequency_dimension == FREQUENCY
    assert energy_dimension == ENERGY
    assert frequency_dimension != energy_dimension, (
        "if these ever became the same dimension the unit engine would relate "
        "a wavenumber to a molar energy silently, which is wrong"
    )


def test_an_imaginary_mode_is_skipped_rather_than_added():
    """An imaginary mode has no zero-point level; it must not contribute."""
    real_only = harmonic_zero_point_energy(_WATER_HF_631GD)
    with_imaginary = harmonic_zero_point_energy((-500.0, *_WATER_HF_631GD))
    assert with_imaginary == pytest.approx(real_only)


def test_a_structure_with_no_real_mode_is_refused_by_name():
    with pytest.raises(AggregationError, match="at least one real"):
        harmonic_zero_point_energy((-500.0, -120.0))


def test_an_unsupported_unit_lists_what_is_supported():
    with pytest.raises(AggregationError) as failure:
        harmonic_zero_point_energy(_WATER_HF_631GD, unit="rydberg")
    message = str(failure.value)
    assert "'rydberg'" in message
    assert "kj/mol" in message


def test_an_empty_or_non_finite_frequency_set_is_refused():
    with pytest.raises(AggregationError, match="at least one frequency"):
        harmonic_zero_point_energy(())
    with pytest.raises(AggregationError, match="finite"):
        harmonic_zero_point_energy((float("nan"), 1000.0))
