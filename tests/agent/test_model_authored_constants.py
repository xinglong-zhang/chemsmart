"""A result must say which of its digits the model supplied.

The receipt closure already proves which measurements a derived value depends
on.  It could not say whether the arithmetic on top of those measurements was
purely structural or whether the model also typed a number into it -- and for a
basis-set limit the typed number moves the answer further than the difference
between the calculations being extrapolated.  These tests pin that the two
cases are distinguishable from the receipt alone.
"""

from dataclasses import asdict

import pytest

from chemsmart.analysis.quantity_expressions import (
    ModelAuthoredConstantV1,
    QuantityContractError,
    QuantityExpressionNodeV1,
    QuantityExpressionOutputDependencyV1,
    QuantityExpressionRequestV1,
    evaluate_quantity_expression,
    quantity_expression_receipt_from_record,
)
from chemsmart.analysis.result_quantities import (
    ENERGY,
    canonical_quantity_sha256,
    make_quantity_value,
)

# Hartree-Fock total energies of water at the experimental equilibrium
# geometry in cc-pVDZ, cc-pVTZ and cc-pVQZ, computed by ChemSmart.
_SERIES = (-76.026798697, -76.057168515, -76.064835339)


def _measured(quantity_id, value):
    return make_quantity_value(
        quantity_id=quantity_id,
        source_value=value,
        source_unit="hartree",
        value=value,
        unit="hartree",
        dimension=ENERGY,
        evidence_ref=f"receipt:{'a' * 64};quantity:{quantity_id}",
    )


def _evaluate(expression_id, inputs, nodes, output):
    return evaluate_quantity_expression(
        QuantityExpressionRequestV1(
            schema_version="chemsmart.quantity-expression-request.v1",
            expression_id=expression_id,
            inputs=inputs,
            nodes=nodes,
            output_node_ids=(output,),
        )
    )


def test_a_parameter_free_limit_declares_no_model_constant():
    receipt = _evaluate(
        "parameter_free",
        (_measured("series", _SERIES),),
        (
            QuantityExpressionNodeV1(
                node_id="cbs",
                operation="exponential_cbs_limit",
                input_ids=("series",),
            ),
        ),
        "cbs",
    )
    dependency = receipt.output_dependencies[0]
    assert dependency.model_authored_constants == ()
    assert len(dependency.source_receipt_sha256s) == 1


def test_a_chosen_decay_exponent_is_named_in_the_output_closure():
    receipt = _evaluate(
        "model_parameterized",
        (_measured("e3", _SERIES[1]), _measured("e4", _SERIES[2])),
        (
            QuantityExpressionNodeV1(
                node_id="cbs",
                operation="scf_exponential_cbs_limit",
                input_ids=("e3", "e4"),
                cardinal_numbers=(3, 4),
                extrapolation_exponent=5.34,
            ),
        ),
        "cbs",
    )
    (constant,) = receipt.output_dependencies[0].model_authored_constants
    assert constant.node_id == "cbs"
    assert constant.role == "extrapolation_exponent"
    assert constant.value == "5.34"


def test_the_two_routes_disagree_by_more_than_the_last_basis_increment():
    """Without this attribution the two results look equally well founded."""

    parameter_free = _evaluate(
        "parameter_free",
        (_measured("series", _SERIES),),
        (
            QuantityExpressionNodeV1(
                node_id="cbs",
                operation="exponential_cbs_limit",
                input_ids=("series",),
            ),
        ),
        "cbs",
    ).outputs[0]
    parameterized = _evaluate(
        "model_parameterized",
        (_measured("e3", _SERIES[1]), _measured("e4", _SERIES[2])),
        (
            QuantityExpressionNodeV1(
                node_id="cbs",
                operation="scf_exponential_cbs_limit",
                input_ids=("e3", "e4"),
                cardinal_numbers=(3, 4),
                extrapolation_exponent=5.34,
            ),
        ),
        "cbs",
    ).outputs[0]
    assert abs(parameter_free.value - parameterized.value) > 1.0e-4


def test_a_constant_reaches_the_output_through_intermediate_nodes():
    """Attribution follows the DAG, so burying the literal does not hide it."""

    receipt = _evaluate(
        "buried",
        (_measured("e", _SERIES[2]),),
        (
            QuantityExpressionNodeV1(
                node_id="fudge",
                operation="literal",
                literal_value=0.0015,
                literal_unit="hartree",
            ),
            QuantityExpressionNodeV1(
                node_id="corrected",
                operation="add",
                input_ids=("e", "fudge"),
            ),
            QuantityExpressionNodeV1(
                node_id="reported",
                operation="ref",
                reference="corrected",
            ),
        ),
        "reported",
    )
    (constant,) = receipt.output_dependencies[0].model_authored_constants
    assert constant.node_id == "fudge"
    assert constant.role == "literal_value"


def test_an_output_that_never_touched_the_constant_stays_clean():
    """Attribution is per output, not per request."""

    receipt = evaluate_quantity_expression(
        QuantityExpressionRequestV1(
            schema_version="chemsmart.quantity-expression-request.v1",
            expression_id="mixed",
            inputs=(_measured("e3", _SERIES[1]), _measured("e4", _SERIES[2])),
            nodes=(
                QuantityExpressionNodeV1(
                    node_id="difference",
                    operation="subtract",
                    input_ids=("e4", "e3"),
                ),
                QuantityExpressionNodeV1(
                    node_id="scaled",
                    operation="scale",
                    input_ids=("difference",),
                    scale_factor=627.5094740631,
                ),
            ),
            output_node_ids=("difference", "scaled"),
        )
    )
    by_output = {
        item.output_id: item.model_authored_constants
        for item in receipt.output_dependencies
    }
    assert by_output["difference"] == ()
    assert [item.role for item in by_output["scaled"]] == ["scale_factor"]


def test_a_power_exponent_is_reported_as_such():
    receipt = _evaluate(
        "powered",
        (_measured("e", 2.0),),
        (
            QuantityExpressionNodeV1(
                node_id="cubed",
                operation="power",
                input_ids=("e",),
                literal_value=3.0,
            ),
        ),
        "cubed",
    )
    (constant,) = receipt.output_dependencies[0].model_authored_constants
    assert constant.role == "power_exponent"


def test_the_attribution_is_bound_into_the_receipt_digest():
    """An attribution that did not change the digest could be edited away."""

    clean = QuantityExpressionOutputDependencyV1(
        output_id="cbs", source_receipt_sha256s=("b" * 64,)
    )
    attributed = QuantityExpressionOutputDependencyV1(
        output_id="cbs",
        source_receipt_sha256s=("b" * 64,),
        model_authored_constants=(
            ModelAuthoredConstantV1(
                node_id="cbs", role="extrapolation_exponent", value="5.34"
            ),
        ),
    )
    assert canonical_quantity_sha256(clean) != canonical_quantity_sha256(
        attributed
    )


def test_a_persisted_receipt_rehydrates_its_attribution():
    receipt = _evaluate(
        "model_parameterized",
        (_measured("e3", _SERIES[1]), _measured("e4", _SERIES[2])),
        (
            QuantityExpressionNodeV1(
                node_id="cbs",
                operation="correlation_inverse_power_cbs_limit",
                input_ids=("e3", "e4"),
                cardinal_numbers=(3, 4),
                extrapolation_exponent=3.0,
            ),
        ),
        "cbs",
    )
    record = {
        key: value
        for key, value in asdict(receipt).items()
        if key != "receipt_sha256"
    }
    assert record["output_dependencies"][0]["model_authored_constants"], (
        "the persisted form must carry the attribution, not only the object"
    )
    rehydrated = quantity_expression_receipt_from_record(
        record, receipt_sha256=receipt.receipt_sha256
    )
    (constant,) = rehydrated.output_dependencies[0].model_authored_constants
    assert constant.role == "extrapolation_exponent"
    assert constant.value == "3.0"


def test_an_unknown_constant_role_is_refused_by_name():
    with pytest.raises(QuantityContractError, match="unsupported"):
        ModelAuthoredConstantV1(
            node_id="cbs", role="vibes", value="5.34"
        )


def test_constants_must_arrive_sorted_and_deduplicated():
    with pytest.raises(QuantityContractError, match="sorted and unique"):
        QuantityExpressionOutputDependencyV1(
            output_id="cbs",
            source_receipt_sha256s=(),
            model_authored_constants=(
                ModelAuthoredConstantV1(
                    node_id="z", role="scale_factor", value="2.0"
                ),
                ModelAuthoredConstantV1(
                    node_id="a", role="scale_factor", value="2.0"
                ),
            ),
        )


def test_the_profile_separates_toolkit_computation_from_model_assembly():
    """The discriminator must not depend on which paper is being reproduced.

    Both requests below return the same number from the same three measured
    energies.  One reached it through the operation that owns the convention;
    the other reassembled the closed form.  Nothing about water, Hartree-Fock,
    or basis sets appears in the check.
    """

    inputs = tuple(
        _measured(name, value)
        for name, value in zip(("e2", "e3", "e4"), _SERIES)
    )
    reassembled = evaluate_quantity_expression(
        QuantityExpressionRequestV1(
            schema_version="chemsmart.quantity-expression-request.v1",
            expression_id="reassembled",
            inputs=inputs,
            nodes=(
                QuantityExpressionNodeV1(
                    node_id="p24", operation="multiply", input_ids=("e2", "e4")
                ),
                QuantityExpressionNodeV1(
                    node_id="p33", operation="multiply", input_ids=("e3", "e3")
                ),
                QuantityExpressionNodeV1(
                    node_id="num", operation="subtract", input_ids=("p24", "p33")
                ),
                QuantityExpressionNodeV1(
                    node_id="s24", operation="add", input_ids=("e2", "e4")
                ),
                QuantityExpressionNodeV1(
                    node_id="t3",
                    operation="scale",
                    input_ids=("e3",),
                    scale_factor=2,
                ),
                QuantityExpressionNodeV1(
                    node_id="den", operation="subtract", input_ids=("s24", "t3")
                ),
                QuantityExpressionNodeV1(
                    node_id="limit", operation="divide", input_ids=("num", "den")
                ),
            ),
            output_node_ids=("limit",),
        )
    ).output_dependencies[0]
    computed = evaluate_quantity_expression(
        QuantityExpressionRequestV1(
            schema_version="chemsmart.quantity-expression-request.v1",
            expression_id="computed",
            inputs=inputs,
            nodes=(
                QuantityExpressionNodeV1(
                    node_id="limit",
                    operation="exponential_cbs_limit",
                    input_ids=("e2", "e3", "e4"),
                ),
            ),
            output_node_ids=("limit",),
        )
    ).output_dependencies[0]

    assert reassembled.convention_operations == ()
    assert reassembled.arithmetic_node_count == 7
    assert reassembled.model_authored_constants

    assert computed.convention_operations == ("exponential_cbs_limit",)
    assert computed.arithmetic_node_count == 0
    assert computed.model_authored_constants == ()


def test_a_convention_reached_through_arithmetic_reports_both():
    """Mixed derivations are common and must not read as either extreme."""

    dependency = evaluate_quantity_expression(
        QuantityExpressionRequestV1(
            schema_version="chemsmart.quantity-expression-request.v1",
            expression_id="mixed",
            inputs=(
                _measured("series", _SERIES),
                _measured("reference", _SERIES[2]),
            ),
            nodes=(
                QuantityExpressionNodeV1(
                    node_id="limit",
                    operation="exponential_cbs_limit",
                    input_ids=("series",),
                ),
                QuantityExpressionNodeV1(
                    node_id="relative",
                    operation="subtract",
                    input_ids=("limit", "reference"),
                ),
            ),
            output_node_ids=("relative",),
        )
    ).output_dependencies[0]
    assert dependency.convention_operations == ("exponential_cbs_limit",)
    assert dependency.arithmetic_node_count == 1


def test_an_unregistered_convention_name_is_refused():
    with pytest.raises(QuantityContractError, match="unregistered convention"):
        QuantityExpressionOutputDependencyV1(
            output_id="x",
            source_receipt_sha256s=(),
            convention_operations=("multiply",),
        )
