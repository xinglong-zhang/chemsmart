"""A ferrocene reference is a registered constant, never a typed literal.

NOVEL-3 ino3 (2026-09-05): asked to quote nickel thiolate potentials
against ferrocene in acetonitrile, a session typed the ferrocene
reference as a literal twice and delivered potentials that no reader
could trace to a source. The registry now holds the computed absolute
value with its stated accuracy and the experimental construction it was
benchmarked against, each with the entry it pairs with.
"""

import pytest

from chemsmart.analysis.literature_constants import literature_constant
from chemsmart.analysis.quantity_expressions import (
    QuantityExpressionNodeV1,
    QuantityExpressionRequestV1,
    evaluate_quantity_expression,
)


@pytest.mark.capability(
    "constant:ferrocene_absolute_reduction_potential_acetonitrile_namazian2010"
)
def test_the_computed_ferrocene_reference_carries_its_accuracy_and_pairing():
    computed = literature_constant(
        "ferrocene_absolute_reduction_potential_acetonitrile_namazian2010"
    )
    assert (computed.value, computed.unit) == (4.988, "V")
    assert "0.05-0.1 V" in computed.convention
    assert "same-level computed ferrocene pair" in computed.purpose
    assert computed.convention_family == (
        "namazian2010_g3mp2rad_cosmors_acetonitrile"
    )


def test_the_experimental_construction_is_one_family_that_adds_up():
    versus_sce = literature_constant(
        "ferrocene_vs_sce_acetonitrile_pavlishchuk2000"
    )
    sce = literature_constant(
        "sce_absolute_potential_acetonitrile_as_used_namazian2010"
    )
    assert versus_sce.convention_family == sce.convention_family
    assert versus_sce.value + sce.value == pytest.approx(4.980)
    assert "not registered" in versus_sce.note


def test_a_constant_node_resolves_each_ferrocene_entry():
    for name in (
        "ferrocene_absolute_reduction_potential_acetonitrile_namazian2010",
        "ferrocene_vs_sce_acetonitrile_pavlishchuk2000",
        "sce_absolute_potential_acetonitrile_as_used_namazian2010",
    ):
        receipt = evaluate_quantity_expression(
            QuantityExpressionRequestV1(
                schema_version="chemsmart.quantity-expression-request.v1",
                expression_id="fc",
                inputs=(),
                nodes=(
                    QuantityExpressionNodeV1(
                        node_id="ref",
                        operation="constant",
                        input_ids=(),
                        constant_name=name,
                    ),
                ),
                output_node_ids=("ref",),
            )
        )
        (value,) = receipt.outputs
        assert value.source_unit == "V"
        assert value.source_value == literature_constant(name).value
