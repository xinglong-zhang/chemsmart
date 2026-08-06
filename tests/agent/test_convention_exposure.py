"""A convention ChemSmart owns must be one the model can select.

The water Hartree-Fock episode is one instance of a class: when a convention
lives in ``chemsmart.analysis`` but is not reachable as an expression
operation, the model has no way to use it and rebuilds it out of arithmetic,
supplying the physical constants itself.  Auditing that per paper would find
one gap at a time and only after a run.

This module closes it structurally.  It enumerates what the analysis layer
owns and requires each entry to be either exposed to the model or explicitly
excluded with a reason, so a convention added later cannot silently become
unreachable.
"""

import inspect

import pytest

from chemsmart.analysis import aggregation
from chemsmart.analysis.quantity_expressions import (
    ARITHMETIC_OPERATIONS,
    CONVENTION_OPERATIONS,
    OPERATION_DESCRIPTIONS,
    QuantityExpressionNodeV1,
    QuantityExpressionRequestV1,
    evaluate_quantity_expression,
)
from chemsmart.analysis.result_quantities import (
    ENERGY,
    TEMPERATURE,
    make_quantity_value,
)

#: Analysis helpers that are deliberately not model-selectable, and why.
#: An entry here is a decision, not an oversight; adding one requires a reason
#: that survives review.
_NOT_MODEL_SELECTABLE = {
    "count_imaginary_modes": "exposed as imaginary_mode_count",
    "convert_energy": (
        "unit handling; the evaluator normalizes units itself and the claim "
        "contract enforces display units"
    ),
    "dataclass": "re-exported from the standard library, not a convention",
    "cbs_extrapolation": (
        "composes the SCF and correlation limits; the model composes them as "
        "two expression nodes so each half is separately attributable"
    ),
    "extrapolate_correlation_helgaker": (
        "the fixed-exponent special case of "
        "correlation_inverse_power_cbs_limit, which is exposed"
    ),
    "extrapolate_correlation_inverse_power": (
        "exposed as correlation_inverse_power_cbs_limit"
    ),
    "extrapolate_scf_exponential": "exposed as scf_exponential_cbs_limit",
    "extrapolate_exponential_three_point": "exposed as exponential_cbs_limit",
    "state_energy_difference": (
        "a typed record around a subtraction the evaluator already performs "
        "dimensionally, with no convention of its own"
    ),
}


def _owned_conventions():
    return {
        name
        for name, value in vars(aggregation).items()
        if inspect.isfunction(value) and not name.startswith("_")
    }


def test_every_owned_convention_is_exposed_or_explicitly_excluded():
    owned = _owned_conventions()
    unaccounted = sorted(
        name
        for name in owned
        if name not in OPERATION_DESCRIPTIONS
        and name not in _NOT_MODEL_SELECTABLE
    )
    assert not unaccounted, (
        "chemsmart.analysis.aggregation owns these conventions but the model "
        f"cannot select them: {unaccounted}. Expose them as expression "
        "operations or record why not."
    )


def test_the_exclusion_list_does_not_outlive_what_it_excludes():
    """A stale exclusion hides the fact that a convention went away."""

    owned = _owned_conventions()
    stale = sorted(set(_NOT_MODEL_SELECTABLE) - owned)
    assert not stale, f"excluded names no longer exist: {stale}"


def test_an_exclusion_that_claims_another_name_must_point_at_a_real_one():
    """"Exposed as X" is the commonest reason and the easiest to get wrong."""

    for name, reason in _NOT_MODEL_SELECTABLE.items():
        assert reason.strip(), name
        for word in reason.replace(",", " ").split():
            if word.endswith("_limit") or word in OPERATION_DESCRIPTIONS:
                assert word in OPERATION_DESCRIPTIONS, (
                    f"{name} is excluded as though {word} were exposed, "
                    "but no such operation exists"
                )


def test_the_operation_classes_partition_the_vocabulary():
    plumbing = set(OPERATION_DESCRIPTIONS) - CONVENTION_OPERATIONS - (
        ARITHMETIC_OPERATIONS
    )
    assert CONVENTION_OPERATIONS & ARITHMETIC_OPERATIONS == set()
    assert (
        CONVENTION_OPERATIONS | ARITHMETIC_OPERATIONS | plumbing
        == set(OPERATION_DESCRIPTIONS)
    )


def _quantity(quantity_id, value, unit, dimension):
    return make_quantity_value(
        quantity_id=quantity_id,
        source_value=value,
        source_unit=unit,
        value=value,
        unit=unit,
        dimension=dimension,
        evidence_ref=f"receipt:{'a' * 64};quantity:{quantity_id}",
    )


_KCAL = 0.0015936010974


def _populations(nodes, inputs, output):
    return evaluate_quantity_expression(
        QuantityExpressionRequestV1(
            schema_version="chemsmart.quantity-expression-request.v1",
            expression_id="boltzmann",
            inputs=inputs,
            nodes=nodes,
            output_node_ids=(output,),
        )
    )


def test_conformer_populations_are_reachable_in_one_node():
    """The commonest post-processing step in the literature, not rebuilt."""

    receipt = _populations(
        (
            QuantityExpressionNodeV1(
                node_id="p",
                operation="boltzmann_populations",
                input_ids=("g", "t"),
            ),
        ),
        (
            _quantity("g", (0.0, 0.5 * _KCAL, 1.2 * _KCAL), "hartree", ENERGY),
            _quantity("t", 298.15, "K", TEMPERATURE),
        ),
        "p",
    )
    populations = receipt.outputs[0].value
    assert sum(populations) == pytest.approx(1.0)
    assert populations[0] == pytest.approx(0.6402, abs=5e-4)
    dependency = receipt.output_dependencies[0]
    assert dependency.convention_operations == ("boltzmann_populations",)
    assert dependency.arithmetic_node_count == 0
    assert dependency.model_authored_constants == ()


def test_a_boltzmann_average_weights_a_property_by_those_populations():
    receipt = _populations(
        (
            QuantityExpressionNodeV1(
                node_id="avg",
                operation="boltzmann_average",
                input_ids=("v", "g", "t"),
            ),
        ),
        (
            _quantity("v", (1.0, 2.0, 3.0), "angstrom", (0, 1, 0, 0, 0, 0)),
            _quantity("g", (0.0, 0.5 * _KCAL, 1.2 * _KCAL), "hartree", ENERGY),
            _quantity("t", 298.15, "K", TEMPERATURE),
        ),
        "avg",
    )
    assert receipt.outputs[0].value == pytest.approx(1.4443, abs=5e-4)
    assert receipt.outputs[0].unit == "angstrom"


def test_the_temperature_must_be_a_temperature_not_a_bare_number():
    from chemsmart.analysis.quantity_expressions import QuantityExpressionError

    with pytest.raises(QuantityExpressionError, match="must be a temperature"):
        _populations(
            (
                QuantityExpressionNodeV1(
                    node_id="p",
                    operation="boltzmann_populations",
                    input_ids=("g", "t"),
                ),
            ),
            (
                _quantity("g", (0.0, _KCAL), "hartree", ENERGY),
                _quantity("t", 298.15, "1", (0, 0, 0, 0, 0, 0)),
            ),
            "p",
        )


def test_mismatched_values_and_states_are_refused_with_both_counts():
    from chemsmart.analysis.quantity_expressions import QuantityExpressionError

    with pytest.raises(QuantityExpressionError) as failure:
        _populations(
            (
                QuantityExpressionNodeV1(
                    node_id="avg",
                    operation="boltzmann_average",
                    input_ids=("v", "g", "t"),
                ),
            ),
            (
                _quantity("v", (1.0, 2.0), "angstrom", (0, 1, 0, 0, 0, 0)),
                _quantity(
                    "g", (0.0, _KCAL, 2 * _KCAL), "hartree", ENERGY
                ),
                _quantity("t", 298.15, "K", TEMPERATURE),
            ),
            "avg",
        )
    assert "2 values and 3 energies" in str(failure.value)


def test_degeneracies_move_the_answer_and_belong_to_the_owned_weighting():
    """n-Butane: the two gauche forms are enantiomers, so gauche counts twice.

    Gibbs free energies computed by ChemSmart at B3LYP/6-31G*.  Folding the
    multiplicity in by hand would put it in the model's arithmetic; declaring
    it keeps the weighting owned and the multiplicity visible as an input.
    """

    anti, gauche = -158.35045218500127, -158.34900918428787
    inputs = (
        _quantity("g", (anti, gauche), "hartree", ENERGY),
        _quantity("t", 298.15, "K", TEMPERATURE),
        _quantity("d", (1.0, 2.0), "1", (0, 0, 0, 0, 0, 0)),
    )
    single = _populations(
        (
            QuantityExpressionNodeV1(
                node_id="p",
                operation="boltzmann_populations",
                input_ids=("g", "t"),
            ),
        ),
        inputs[:2],
        "p",
    ).outputs[0]
    weighted = _populations(
        (
            QuantityExpressionNodeV1(
                node_id="p",
                operation="boltzmann_populations",
                input_ids=("g", "t", "d"),
            ),
        ),
        inputs,
        "p",
    )
    assert single.value[0] == pytest.approx(0.8218, abs=5e-4)
    assert weighted.outputs[0].value[0] == pytest.approx(0.6974, abs=5e-4)
    # The measured ~68% anti is only reachable with the multiplicity.
    assert abs(single.value[0] - weighted.outputs[0].value[0]) > 0.1
    dependency = weighted.output_dependencies[0]
    assert dependency.convention_operations == ("boltzmann_populations",)
    assert dependency.arithmetic_node_count == 0


def test_a_degeneracy_vector_must_be_dimensionless():
    from chemsmart.analysis.quantity_expressions import QuantityExpressionError

    with pytest.raises(QuantityExpressionError, match="dimensionless"):
        _populations(
            (
                QuantityExpressionNodeV1(
                    node_id="p",
                    operation="boltzmann_populations",
                    input_ids=("g", "t", "d"),
                ),
            ),
            (
                _quantity("g", (0.0, _KCAL), "hartree", ENERGY),
                _quantity("t", 298.15, "K", TEMPERATURE),
                _quantity("d", (1.0, 2.0), "hartree", ENERGY),
            ),
            "p",
        )


def test_a_wrong_degeneracy_count_names_both_counts():
    from chemsmart.analysis.quantity_expressions import QuantityExpressionError

    with pytest.raises(QuantityExpressionError) as failure:
        _populations(
            (
                QuantityExpressionNodeV1(
                    node_id="p",
                    operation="boltzmann_populations",
                    input_ids=("g", "t", "d"),
                ),
            ),
            (
                _quantity("g", (0.0, _KCAL, 2 * _KCAL), "hartree", ENERGY),
                _quantity("t", 298.15, "K", TEMPERATURE),
                _quantity("d", (1.0, 2.0), "1", (0, 0, 0, 0, 0, 0)),
            ),
            "p",
        )
    assert "2 degeneracies for 3 states" in str(failure.value)


_FREQ = (0, 0, 0, 0, 1, 0)


def _count(freqs, cutoff=None):
    inputs = (_quantity("f", freqs, "cm^-1", _FREQ),)
    ids = ("f",)
    if cutoff is not None:
        inputs += (_quantity("c", cutoff, "cm^-1", _FREQ),)
        ids += ("c",)
    return _populations(
        (
            QuantityExpressionNodeV1(
                node_id="n", operation="imaginary_mode_count", input_ids=ids
            ),
        ),
        inputs,
        "n",
    )


def test_counting_imaginary_modes_is_one_node_not_twenty_two():
    """Observed live: a session built 22 arithmetic nodes to count negatives.

    Every optimization must show zero and every transition state exactly one,
    so this is among the most universal post-processing steps in the field.
    Left unexposed, the definition of "imaginary" lived in the model's
    arithmetic rather than in the toolkit.
    """

    receipt = _count((121.6, 300.0, 900.0))
    assert receipt.outputs[0].value == 0.0
    dependency = receipt.output_dependencies[0]
    assert dependency.convention_operations == ("imaginary_mode_count",)
    assert dependency.arithmetic_node_count == 0
    assert dependency.model_authored_constants == ()


def test_a_transition_state_shows_exactly_one():
    assert _count((-450.0, 300.0, 900.0)).outputs[0].value == 1.0


def test_a_cutoff_ignores_a_numerically_noisy_near_zero_mode():
    """The treatment the PySCF stationary-point policy already applies."""

    noisy = (-3.2, 300.0, 900.0)
    assert _count(noisy).outputs[0].value == 1.0
    assert _count(noisy, cutoff=10.0).outputs[0].value == 0.0
    assert _count((-450.0, 300.0), cutoff=10.0).outputs[0].value == 1.0


def test_the_count_is_dimensionless_whatever_the_frequencies_were():
    assert _count((121.6, 300.0)).outputs[0].unit == "1"


def test_a_non_frequency_input_is_refused_by_name():
    from chemsmart.analysis.quantity_expressions import QuantityExpressionError

    with pytest.raises(QuantityExpressionError, match="frequency vector"):
        _populations(
            (
                QuantityExpressionNodeV1(
                    node_id="n",
                    operation="imaginary_mode_count",
                    input_ids=("e",),
                ),
            ),
            (_quantity("e", (-1.0, 2.0), "hartree", ENERGY),),
            "n",
        )
