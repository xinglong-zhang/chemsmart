"""A barrier has a height and a position, and only one was reachable.

`min` and `max` reduce a series to a scalar and discard the index they reduced
over, so a session could report that a torsional barrier is 8.3 kcal/mol but
never that it falls at 180 degrees. That is not a gap composition can close:
no combination of the existing operations recovers an argument of an extremum,
and there is no interpolation either.

The surface used here is the real ORCA 6.1.1 relaxed scan of hydrogen peroxide
run through the ChemSmart CLI on this host, the same fixture the scan reader is
verified against, so the operation is exercised on values a program actually
produced rather than on a hand-made series.
"""

from __future__ import annotations

import pytest

from chemsmart.analysis.quantity_expressions import (
    QuantityExpressionError,
    QuantityExpressionNodeV1,
    QuantityExpressionRequestV1,
    evaluate_quantity_expression,
)
from chemsmart.analysis.result_quantities import (
    DIMENSIONLESS,
    ENERGY,
    make_quantity_value,
)
from chemsmart.io.orca.output import ORCAOutput

_SCAN = "tests/data/ORCATests/outputs/hooh_relaxed_scan_excerpt.out"


def _series(quantity_id, values, dimension, unit):
    return make_quantity_value(
        quantity_id=quantity_id,
        source_value=tuple(values),
        source_unit=unit,
        value=tuple(values),
        unit=unit,
        dimension=dimension,
        evidence_ref=f"receipt:{'a' * 64};quantity:{quantity_id}",
    )


def _locate(operation, values, coordinates):
    return evaluate_quantity_expression(
        QuantityExpressionRequestV1(
            schema_version="chemsmart.quantity-expression-request.v1",
            expression_id="torsion",
            inputs=(
                _series("e", values, ENERGY, "hartree"),
                _series("angle", coordinates, DIMENSIONLESS, "1"),
            ),
            nodes=(
                QuantityExpressionNodeV1(
                    node_id="where",
                    operation=operation,
                    input_ids=("e", "angle"),
                ),
            ),
            output_node_ids=("where",),
        )
    )


@pytest.fixture(scope="module")
def surface():
    profile = ORCAOutput(_SCAN).scan_profile
    return (
        [point["energy"] for point in profile],
        [point["coordinate"] for point in profile],
    )


def test_the_barrier_position_comes_out_of_a_real_surface(surface):
    energies, angles = surface

    receipt = _locate("coordinate_at_maximum", energies, angles)

    assert receipt.outputs[0].value == pytest.approx(180.0)


def test_the_well_position_comes_out_too(surface):
    energies, angles = surface

    receipt = _locate("coordinate_at_minimum", energies, angles)

    assert receipt.outputs[0].value == pytest.approx(60.0)


def test_the_answer_carries_the_coordinate_s_dimension(surface):
    """A position is not an energy, however it was found."""

    energies, angles = surface

    output = _locate("coordinate_at_maximum", energies, angles).outputs[0]

    assert output.dimension == DIMENSIONLESS


def test_height_and_position_are_answered_by_different_operations(surface):
    """The two halves of one question, and why max alone was not enough."""

    energies, angles = surface
    height = max(energies) - min(energies)
    position = _locate("coordinate_at_maximum", energies, angles)

    assert height * 627.5094740631 == pytest.approx(8.26, abs=0.05)
    assert position.outputs[0].value == pytest.approx(180.0)


def test_mismatched_series_lengths_are_refused(surface):
    energies, angles = surface

    with pytest.raises(QuantityExpressionError, match="one coordinate per"):
        _locate("coordinate_at_maximum", energies, angles[:-1])


def test_a_single_point_has_no_extremum(surface):
    energies, angles = surface

    with pytest.raises(QuantityExpressionError, match="has no extremum"):
        _locate("coordinate_at_maximum", energies[:1], angles[:1])


def test_a_wrong_input_count_names_the_order_expected():
    with pytest.raises(QuantityExpressionError, match="in order"):
        evaluate_quantity_expression(
            QuantityExpressionRequestV1(
                schema_version="chemsmart.quantity-expression-request.v1",
                expression_id="torsion",
                inputs=(_series("e", (1.0, 2.0), ENERGY, "hartree"),),
                nodes=(
                    QuantityExpressionNodeV1(
                        node_id="where",
                        operation="coordinate_at_maximum",
                        input_ids=("e",),
                    ),
                ),
                output_node_ids=("where",),
            )
        )
