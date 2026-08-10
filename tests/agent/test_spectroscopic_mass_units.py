"""Mass and spectroscopic frequency stay dimensional through expressions."""

import pytest

from chemsmart.analysis.quantity_expressions import (
    QuantityExpressionNodeV1,
    QuantityExpressionRequestV1,
    convert_normalized_value,
    evaluate_quantity_expression,
)


def test_atomic_mass_moment_and_rotational_constant_round_trip():
    receipt = evaluate_quantity_expression(
        QuantityExpressionRequestV1(
            schema_version="chemsmart.quantity-expression-request.v1",
            expression_id="rigid_rotor_units",
            inputs=(),
            nodes=(
                QuantityExpressionNodeV1(
                    node_id="mass",
                    operation="literal",
                    literal_value=1.00782503223,
                    literal_unit="u",
                ),
                QuantityExpressionNodeV1(
                    node_id="length",
                    operation="literal",
                    literal_value=1.0,
                    literal_unit="angstrom",
                ),
                QuantityExpressionNodeV1(
                    node_id="length_sq",
                    operation="multiply",
                    input_ids=("length", "length"),
                ),
                QuantityExpressionNodeV1(
                    node_id="moment",
                    operation="multiply",
                    input_ids=("mass", "length_sq"),
                ),
                QuantityExpressionNodeV1(
                    node_id="constant",
                    operation="literal",
                    literal_value=16.85763,
                    literal_unit="cm^-1*u*angstrom^2",
                ),
                QuantityExpressionNodeV1(
                    node_id="rotation",
                    operation="divide",
                    input_ids=("constant", "moment"),
                ),
            ),
            output_node_ids=("moment", "rotation"),
        )
    )

    moment, rotation = receipt.outputs
    assert moment.unit == "angstrom^2 u"
    assert rotation.unit == "cm^-1"
    assert rotation.value == pytest.approx(16.85763 / 1.00782503223)
    assert convert_normalized_value(
        rotation.value, rotation.dimension, "GHz"
    ) == pytest.approx(rotation.value * 29.9792458)


def test_geometry_to_rigid_rotor_conventions_own_the_full_physics():
    receipt = evaluate_quantity_expression(
        QuantityExpressionRequestV1(
            schema_version="chemsmart.quantity-expression-request.v1",
            expression_id="h2s_rigid_rotor",
            inputs=(),
            nodes=(
                QuantityExpressionNodeV1(
                    node_id="positions",
                    operation="literal",
                    literal_value=(
                        (0.0, 0.0, 0.000568),
                        (0.0, 0.966365, 0.929716),
                        (0.0, -0.966365, 0.929716),
                    ),
                    literal_unit="angstrom",
                ),
                QuantityExpressionNodeV1(
                    node_id="masses",
                    operation="literal",
                    literal_value=(
                        31.9720711744,
                        1.00782503223,
                        1.00782503223,
                    ),
                    literal_unit="u",
                ),
                QuantityExpressionNodeV1(
                    node_id="com",
                    operation="center_of_mass",
                    input_ids=("positions", "masses"),
                ),
                QuantityExpressionNodeV1(
                    node_id="moments",
                    operation="principal_moments_of_inertia",
                    input_ids=("positions", "masses"),
                ),
                QuantityExpressionNodeV1(
                    node_id="rotors",
                    operation="rigid_rotor_constants",
                    input_ids=("moments",),
                ),
            ),
            output_node_ids=("com", "moments", "rotors"),
        )
    )

    com, moments, rotors = receipt.outputs
    assert com.value == pytest.approx((0.0, 0.0, 0.05567134785115173))
    assert moments.unit == "angstrom^2 u"
    assert moments.value == pytest.approx(
        (1.6369433618117653, 1.8823376161986716, 3.519280978010437)
    )
    assert rotors.value == pytest.approx(
        (10.298237, 8.955689, 4.790078), rel=2.0e-7
    )
    assert receipt.output_dependencies[0].convention_operations == (
        "center_of_mass",
    )
    assert receipt.output_dependencies[1].convention_operations == (
        "principal_moments_of_inertia",
    )
    assert receipt.output_dependencies[2].convention_operations == (
        "principal_moments_of_inertia",
        "rigid_rotor_constants",
    )


def test_linear_rotor_constant_owns_the_spectroscopic_conversion():
    receipt = evaluate_quantity_expression(
        QuantityExpressionRequestV1(
            schema_version="chemsmart.quantity-expression-request.v1",
            expression_id="hcn_linear_rotor",
            inputs=(),
            nodes=(
                QuantityExpressionNodeV1(
                    node_id="positions",
                    operation="literal",
                    literal_value=(
                        (0.0, 0.0, -1.05886869),
                        (0.0, 0.0, 0.00639579967),
                        (0.0, 0.0, 1.15263877),
                    ),
                    literal_unit="angstrom",
                ),
                QuantityExpressionNodeV1(
                    node_id="masses",
                    operation="literal",
                    literal_value=(1.00782503223, 12.0, 14.00307400443),
                    literal_unit="u",
                ),
                QuantityExpressionNodeV1(
                    node_id="moments",
                    operation="principal_moments_of_inertia",
                    input_ids=("positions", "masses"),
                ),
                QuantityExpressionNodeV1(
                    node_id="rotation",
                    operation="linear_rotor_constant",
                    input_ids=("moments",),
                ),
            ),
            output_node_ids=("moments", "rotation"),
        )
    )

    moments, rotation = receipt.outputs
    assert moments.value == pytest.approx(
        (0.0, 11.237122838017903, 11.237122838017903)
    )
    assert rotation.value == pytest.approx(1.5001730803184132)
    assert {
        item.node_id
        for item in receipt.output_dependencies[1].model_authored_constants
    } == {"positions", "masses"}
    assert set(receipt.output_dependencies[1].convention_operations) == {
        "principal_moments_of_inertia",
        "linear_rotor_constant",
    }
