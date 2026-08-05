"""Post-processing aggregation: analytic fixtures, not regression snapshots."""

import math

import pytest

from chemsmart.analysis.aggregation import (
    AggregationError,
    boltzmann_average,
    boltzmann_populations,
    cbs_extrapolation,
    convert_energy,
    extrapolate_correlation_helgaker,
    extrapolate_scf_exponential,
    state_energy_difference,
)


def test_correlation_extrapolation_recovers_an_exact_x_cubed_series():
    limit, coefficient = -1.0, 1.0
    recovered = extrapolate_correlation_helgaker(
        smaller_cardinal=4,
        larger_cardinal=5,
        smaller_correlation_energy=limit + coefficient * 4**-3,
        larger_correlation_energy=limit + coefficient * 5**-3,
    )
    assert recovered == pytest.approx(limit, abs=1e-12)


def test_scf_extrapolation_recovers_an_exact_exponential_series():
    limit, amplitude, alpha = -76.0, 0.5, 3.9
    recovered = extrapolate_scf_exponential(
        smaller_cardinal=4,
        larger_cardinal=5,
        smaller_scf_energy=limit + amplitude * math.exp(-alpha * math.sqrt(4)),
        larger_scf_energy=limit + amplitude * math.exp(-alpha * math.sqrt(5)),
        alpha=alpha,
    )
    assert recovered == pytest.approx(limit, abs=1e-9)


def test_components_are_extrapolated_separately_and_reported_separately():
    """SCF and correlation converge by different laws and must not be merged."""

    alpha = 3.9
    scf_limit, scf_amplitude = -76.0, 0.5
    corr_limit, corr_coefficient = -1.0, 1.0
    result = cbs_extrapolation(
        smaller_cardinal=4,
        larger_cardinal=5,
        smaller_scf_energy=scf_limit
        + scf_amplitude * math.exp(-alpha * math.sqrt(4)),
        larger_scf_energy=scf_limit
        + scf_amplitude * math.exp(-alpha * math.sqrt(5)),
        smaller_correlation_energy=corr_limit + corr_coefficient * 4**-3,
        larger_correlation_energy=corr_limit + corr_coefficient * 5**-3,
        scf_alpha=alpha,
    )
    assert result.scf_energy == pytest.approx(scf_limit, abs=1e-9)
    assert result.correlation_energy == pytest.approx(corr_limit, abs=1e-12)
    assert result.total_energy == pytest.approx(
        scf_limit + corr_limit, abs=1e-9
    )
    assert (result.smaller_cardinal, result.larger_cardinal) == (4, 5)


def test_cardinal_ordering_and_range_are_enforced():
    with pytest.raises(AggregationError, match="must exceed"):
        extrapolate_correlation_helgaker(
            smaller_cardinal=5,
            larger_cardinal=4,
            smaller_correlation_energy=-1.0,
            larger_correlation_energy=-1.0,
        )
    with pytest.raises(AggregationError, match="at least 2"):
        extrapolate_correlation_helgaker(
            smaller_cardinal=1,
            larger_cardinal=3,
            smaller_correlation_energy=-1.0,
            larger_correlation_energy=-1.0,
        )


def test_state_difference_is_signed_and_states_its_direction():
    result = state_energy_difference(
        upper_state_energy=-39.14999528,
        lower_state_energy=-39.16874454,
        upper_state_label="singlet",
        lower_state_label="triplet",
    )
    assert result.value == pytest.approx(11.765, abs=1e-3)
    assert result.convention == "E(upper) - E(lower)"
    assert (result.upper_state_label, result.lower_state_label) == (
        "singlet",
        "triplet",
    )


def test_energy_unit_conversion_round_trips_known_factors():
    assert convert_energy(1.0, "kcal/mol") == pytest.approx(
        627.509474, abs=1e-6
    )
    assert convert_energy(1.0, "eV") == pytest.approx(27.211386, abs=1e-6)
    with pytest.raises(AggregationError, match="unsupported energy unit"):
        convert_energy(1.0, "furlongs")


def test_boltzmann_populations_normalise_and_order_correctly():
    populations = boltzmann_populations((0.0, 1.0, 2.0))
    assert sum(populations) == pytest.approx(1.0, abs=1e-12)
    assert populations[0] > populations[1] > populations[2]


def test_boltzmann_populations_are_invariant_to_an_energy_offset():
    shifted = boltzmann_populations((10.0, 11.0, 12.0))
    base = boltzmann_populations((0.0, 1.0, 2.0))
    for left, right in zip(base, shifted):
        assert left == pytest.approx(right, abs=1e-12)


def test_boltzmann_average_weights_values_by_population():
    assert boltzmann_average((1.0, 2.0), (0.0, 0.0)) == pytest.approx(1.5)
    assert boltzmann_average((1.0, 2.0), (0.0, 100.0)) == pytest.approx(
        1.0, abs=1e-6
    )
    with pytest.raises(AggregationError, match="same length"):
        boltzmann_average((1.0,), (0.0, 1.0))


def test_three_point_exponential_recovers_an_exact_series():
    """A live run reported this gap: HF CBS is a three-parameter fit.

    The engine exposed only linear fits, so a paper stating "extrapolated with
    a three-parameter exponential formula" had no expressible producer for its
    complete-basis limit.
    """

    from chemsmart.analysis.aggregation import (
        extrapolate_exponential_three_point,
    )

    limit, amplitude, ratio = -76.06, 0.85, 0.28
    series = tuple(limit + amplitude * ratio**n for n in (2, 3, 4))
    assert extrapolate_exponential_three_point(series) == pytest.approx(
        limit, abs=1e-10
    )


def test_a_series_without_exponential_curvature_is_refused():
    from chemsmart.analysis.aggregation import (
        extrapolate_exponential_three_point,
    )

    with pytest.raises(AggregationError, match="no exponential curvature"):
        extrapolate_exponential_three_point((1.0, 1.0, 1.0))
    with pytest.raises(AggregationError, match="shrink monotonically"):
        extrapolate_exponential_three_point((0.0, 1.0, 3.0))
    with pytest.raises(AggregationError, match="exactly three"):
        extrapolate_exponential_three_point((0.0, 1.0))


def test_the_expression_engine_exposes_the_exponential_limit():
    from chemsmart.analysis.quantity_expressions import (
        QuantityExpressionNodeV1,
        QuantityExpressionRequestV1,
        evaluate_quantity_expression,
    )
    from chemsmart.analysis.result_quantities import (
        ENERGY,
        make_quantity_value,
    )

    limit, amplitude, ratio = -76.06, 0.85, 0.28
    series = [limit + amplitude * ratio**n for n in (2, 3, 4)]
    quantity = make_quantity_value(
        quantity_id="hf_series",
        source_value=series,
        source_unit="hartree",
        value=series,
        unit="hartree",
        dimension=ENERGY,
        evidence_ref="artifact:x#" + "0" * 64,
    )
    receipt = evaluate_quantity_expression(
        QuantityExpressionRequestV1(
            schema_version="chemsmart.quantity-expression-request.v1",
            expression_id="hf.cbs",
            inputs=(quantity,),
            nodes=(
                QuantityExpressionNodeV1(
                    node_id="e_cbs",
                    operation="exponential_cbs_limit",
                    input_ids=("hf_series",),
                ),
            ),
            output_node_ids=("e_cbs",),
        )
    )
    assert receipt.outputs[0].value == pytest.approx(limit, abs=1e-10)
    assert receipt.outputs[0].unit == "hartree"
