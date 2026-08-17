"""An omitted semantic role must be derived, not demanded.

The commonest thing asked of this tool chain is a comparison, and a comparison
repeats a source quantity name by construction: four thermochemistry receipts
each yield a quantity called ``gibbs_free_energy``.  For five cycles the host
asked the model to disambiguate that repetition by hand-authoring a role per
input, a requirement stated in terms of receipt internals the model never sees.
Each attempt to explain it better halved the failure without clearing it.

The host now falls back to the input's own id, which the expression contract
already requires to be unique.  These tests execute that path rather than
reading the wording that describes it: the first proves an expression with a
repeated source quantity and no declared roles is *accepted*, and the second
proves the enforcement that remains -- two inputs the model gave the *same*
explicit role are still refused, because that would make the evidence
reference ambiguous.
"""

from __future__ import annotations

import pytest

from chemsmart.analysis.quantity_expressions import (
    QuantityExpressionError,
    QuantityExpressionNodeV1,
    QuantityExpressionRequestV1,
    evaluate_quantity_expression,
    quantity_expression_semantic_signature,
)
from chemsmart.analysis.result_quantities import (
    ENERGY,
    QuantityContractError,
    make_quantity_value,
)

#: Gibbs free energies of two rotamers, in hartree.  The values are immaterial
#: here; what matters is that both arrive under one source quantity name.
_GAUCHE = -229.123456
_ANTI = -229.120987


def _thermochemistry_input(input_id, value, *, role=""):
    """One input drawn from a thermochemistry receipt.

    ``quantity:gibbs_free_energy`` is the *source* name inside the producing
    receipt and is identical for every species; ``input_id`` is the model's own
    handle and is unique by contract.  That difference is the whole subject of
    these tests.
    """

    evidence = f"receipt:{'a' * 64};quantity:gibbs_free_energy"
    if role:
        evidence = f"{evidence};semantic-role:{role}"
    return make_quantity_value(
        quantity_id=input_id,
        source_value=value,
        source_unit="hartree",
        value=value,
        unit="hartree",
        dimension=ENERGY,
        evidence_ref=evidence,
    )


def _difference(inputs):
    return QuantityExpressionRequestV1(
        schema_version="chemsmart.quantity-expression-request.v1",
        expression_id="rotamer-preference",
        inputs=inputs,
        nodes=(
            QuantityExpressionNodeV1(
                node_id="dg",
                operation="subtract",
                input_ids=tuple(item.quantity_id for item in inputs),
            ),
        ),
        output_node_ids=("dg",),
    )


def test_a_repeated_source_quantity_needs_no_declared_role():
    """The case that failed for five cycles must now simply work."""

    request = _difference(
        (
            _thermochemistry_input("g_gauche", _GAUCHE),
            _thermochemistry_input("g_anti", _ANTI),
        )
    )

    receipt = evaluate_quantity_expression(request)

    assert receipt.outputs[0].value == pytest.approx(_GAUCHE - _ANTI)
    # The signature is what the role feeds; it must resolve without one.
    assert quantity_expression_semantic_signature(request)


def test_the_derived_roles_keep_the_two_occurrences_apart():
    """Deriving a role must distinguish, not merely avoid raising."""

    swapped = _difference(
        (
            _thermochemistry_input("g_anti", _GAUCHE),
            _thermochemistry_input("g_gauche", _ANTI),
        )
    )
    original = _difference(
        (
            _thermochemistry_input("g_gauche", _GAUCHE),
            _thermochemistry_input("g_anti", _ANTI),
        )
    )

    assert quantity_expression_semantic_signature(
        original
    ) != quantity_expression_semantic_signature(swapped)


def test_two_inputs_given_the_same_explicit_role_are_still_refused():
    """The retained enforcement: a shared role is genuinely ambiguous.

    Cycle 029's reason for keeping this refusal stands -- two inputs sharing a
    role collapse onto one slot, so the evidence reference no longer names a
    single occurrence.  Deriving an absent role never required relaxing it.
    """

    request = _difference(
        (
            _thermochemistry_input("g_gauche", _GAUCHE, role="rotamer"),
            _thermochemistry_input("g_anti", _ANTI, role="rotamer"),
        )
    )

    with pytest.raises(
        (QuantityContractError, QuantityExpressionError),
        match="distinct semantic_role",
    ):
        quantity_expression_semantic_signature(request)
