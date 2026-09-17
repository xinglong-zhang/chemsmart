"""An electrode potential is energy per charge, not a number labelled volts.

Before this, ``normalize_numeric_value(1.0, "V")`` raised ``unsupported unit``.
There were eight dimension bases -- energy, length, temperature, angle,
frequency, pressure, dipole moment, mass -- and none was electric, so a redox
potential could not be a typed quantity, a registered literature constant, or
a claim. Every workflow family that ends in a potential stopped there.

Charge is added as the ninth base rather than potential, and that is the whole
point: with charge present, potential is *derived* as energy per charge, so
dG = -nFE typechecks instead of being a convention nobody can verify, and the
Faraday constant dissolves into the unit system -- it is a definition, not a
value taken from the literature, and it does not belong in a registry of
measured constants.

One volt acting on one elementary charge is one electronvolt, by definition of
the electronvolt, so the volt conversion introduces no new physical constant
either.
"""

from __future__ import annotations

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.analysis_claims import AnalysisReportedQuantityV1
from chemsmart.analysis.quantity_expressions import (
    QuantityExpressionNodeV1,
    QuantityExpressionRequestV1,
    canonical_unit_for_dimension,
    evaluate_quantity_expression,
    normalize_numeric_value,
)
from chemsmart.analysis.result_quantities import (
    CHARGE,
    DIPOLE_MOMENT,
    ELECTRIC_POTENTIAL,
    ENERGY,
    QuantityContractError,
    make_quantity_value,
)

_SCHEMA = "chemsmart.quantity-expression-request.v1"


def _quantity(quantity_id, value, unit):
    canonical, canonical_unit, dimension = normalize_numeric_value(value, unit)
    return make_quantity_value(
        quantity_id=quantity_id,
        source_value=value,
        source_unit=unit,
        value=canonical,
        unit=canonical_unit,
        dimension=dimension,
        evidence_ref=f"quantity:{quantity_id}",
    )


@pytest.mark.parametrize("spelling", ["V", "volt", " v "])
def test_a_volt_is_an_electric_potential(spelling):
    _, unit, dimension = normalize_numeric_value(1.0, spelling)
    assert dimension == ELECTRIC_POTENTIAL
    assert unit == "hartree e^-1"


def test_the_atomic_unit_of_potential_is_the_expected_number_of_volts():
    """A wrong factor here would be invisible in every downstream number."""

    one_volt_in_atomic_units, _, _ = normalize_numeric_value(1.0, "V")
    assert 1.0 / one_volt_in_atomic_units == pytest.approx(27.211386, abs=1e-4)


def test_a_millivolt_is_a_thousandth_of_a_volt():
    volt, _, _ = normalize_numeric_value(1.0, "V")
    millivolt, _, _ = normalize_numeric_value(1.0, "mV")
    assert millivolt == pytest.approx(volt * 1e-3)


def test_the_elementary_charge_is_the_canonical_charge_unit():
    value, unit, dimension = normalize_numeric_value(1.0, "e")
    assert (value, unit, dimension) == (1.0, "e", CHARGE)
    assert canonical_unit_for_dimension(CHARGE) == "e"
    assert canonical_unit_for_dimension(ELECTRIC_POTENTIAL) == "hartree e^-1"


def test_charge_times_potential_is_energy():
    """The relation the ninth base exists to make checkable.

    A one-electron half reaction referenced to an absolute electrode
    potential of 4.28 V costs 98.70 kcal/mol, and nothing in the expression
    asserts that -- it falls out of the dimensions.
    """

    nodes = (
        QuantityExpressionNodeV1(
            node_id="potential", operation="ref", reference="E_abs"
        ),
        QuantityExpressionNodeV1(
            node_id="charge", operation="ref", reference="n_e"
        ),
        QuantityExpressionNodeV1(
            node_id="work",
            operation="multiply",
            input_ids=("potential", "charge"),
        ),
        QuantityExpressionNodeV1(
            node_id="work_kcal",
            operation="convert",
            input_ids=("work",),
            target_unit="kcal/mol",
        ),
    )
    receipt = evaluate_quantity_expression(
        QuantityExpressionRequestV1(
            schema_version=_SCHEMA,
            expression_id="nFE",
            inputs=(
                _quantity("E_abs", 4.28, "V"),
                _quantity("n_e", 1.0, "e"),
            ),
            nodes=nodes,
            output_node_ids=("work", "work_kcal"),
        )
    )
    outputs = {item.quantity_id: item for item in receipt.outputs}
    assert outputs["work"].dimension == ENERGY
    assert outputs["work"].unit == "hartree"
    assert outputs["work_kcal"].source_value == pytest.approx(98.70, abs=0.01)


def test_energy_over_charge_is_a_potential():
    """The inverse direction, which is how a computed dG becomes a potential."""

    nodes = (
        QuantityExpressionNodeV1(
            node_id="dg", operation="ref", reference="delta_g"
        ),
        QuantityExpressionNodeV1(
            node_id="charge", operation="ref", reference="n_e"
        ),
        QuantityExpressionNodeV1(
            node_id="potential",
            operation="divide",
            input_ids=("dg", "charge"),
        ),
        QuantityExpressionNodeV1(
            node_id="in_volts",
            operation="convert",
            input_ids=("potential",),
            target_unit="V",
        ),
    )
    receipt = evaluate_quantity_expression(
        QuantityExpressionRequestV1(
            schema_version=_SCHEMA,
            expression_id="dg_to_e",
            inputs=(
                _quantity("delta_g", -98.70, "kcal/mol"),
                _quantity("n_e", 1.0, "e"),
            ),
            nodes=nodes,
            output_node_ids=("potential", "in_volts"),
        )
    )
    outputs = {item.quantity_id: item for item in receipt.outputs}
    assert outputs["potential"].dimension == ELECTRIC_POTENTIAL
    assert outputs["in_volts"].source_value == pytest.approx(-4.28, abs=0.01)


def test_a_potential_can_become_a_claim():
    """Nine-integer dimensions are admitted by the claim record."""

    claim = AnalysisReportedQuantityV1(
        claim_id="e_half",
        source_kind="quantity_expression",
        source_receipt_sha256="a" * 64,
        quantity_id="potential",
        quantity_value_sha256="b" * 64,
        display_value=0.79,
        display_unit="V",
        canonical_value=0.029,
        canonical_unit="hartree e^-1",
        dimension=ELECTRIC_POTENTIAL,
        data_kind="scalar",
    )
    assert len(claim.dimension) == 9


@pytest.mark.parametrize("length", [6, 7, 8, 9, 10])
def test_supported_dimension_widths_stay_valid(length):
    """Widening must not invalidate a receipt already written."""

    claim = AnalysisReportedQuantityV1(
        claim_id="legacy",
        source_kind="quantity_extraction",
        source_receipt_sha256="a" * 64,
        quantity_id="energy",
        quantity_value_sha256="b" * 64,
        display_value=1.0,
        display_unit="hartree",
        canonical_value=1.0,
        canonical_unit="hartree",
        dimension=(1,) + (0,) * (length - 1),
        data_kind="scalar",
    )
    assert len(claim.dimension) == length


def test_an_eleven_integer_dimension_is_still_refused():
    with pytest.raises(ContractError):
        AnalysisReportedQuantityV1(
            claim_id="bogus",
            source_kind="quantity_extraction",
            source_receipt_sha256="a" * 64,
            quantity_id="energy",
            quantity_value_sha256="b" * 64,
            display_value=1.0,
            display_unit="hartree",
            canonical_value=1.0,
            canonical_unit="hartree",
            dimension=(1,) + (0,) * 10,
            data_kind="scalar",
        )
    with pytest.raises(QuantityContractError):
        make_quantity_value(
            quantity_id="bogus",
            source_value=1.0,
            source_unit="1",
            value=1.0,
            unit="1",
            dimension=(0,) * 11,
            evidence_ref="quantity:bogus",
        )


def test_the_dipole_moment_base_is_untouched():
    """Charge times length is not silently a dipole, and must not become one.

    Dipole moment stays its own base because rewriting it would change the
    dimension stored in every dipole receipt already written, and a receipt's
    dimension is part of its identity.
    """

    assert normalize_numeric_value(1.0, "debye")[2] == DIPOLE_MOMENT
    assert normalize_numeric_value(1.0, "e bohr")[2] == DIPOLE_MOMENT
    assert normalize_numeric_value(1.0, "e a0")[2] == DIPOLE_MOMENT
    assert DIPOLE_MOMENT != CHARGE


def _redox(delta_g_kcal, electrons, reference_volts=None):
    inputs = [
        _quantity("in_dg", delta_g_kcal, "kcal/mol"),
        _quantity("in_n", float(electrons), "1"),
    ]
    nodes = [
        QuantityExpressionNodeV1(
            node_id="dg", operation="ref", reference="in_dg"
        ),
        QuantityExpressionNodeV1(
            node_id="n", operation="ref", reference="in_n"
        ),
        QuantityExpressionNodeV1(
            node_id="e_abs",
            operation="gibbs_to_redox_potential",
            input_ids=("dg", "n"),
        ),
        QuantityExpressionNodeV1(
            node_id="e_abs_volts",
            operation="convert",
            input_ids=("e_abs",),
            target_unit="V",
        ),
    ]
    outputs = ["e_abs_volts"]
    if reference_volts is not None:
        inputs.append(_quantity("in_ref", reference_volts, "V"))
        nodes += [
            QuantityExpressionNodeV1(
                node_id="ref", operation="ref", reference="in_ref"
            ),
            QuantityExpressionNodeV1(
                node_id="e_vs_ref",
                operation="subtract",
                input_ids=("e_abs", "ref"),
            ),
            QuantityExpressionNodeV1(
                node_id="e_vs_ref_volts",
                operation="convert",
                input_ids=("e_vs_ref",),
                target_unit="V",
            ),
        ]
        outputs.append("e_vs_ref_volts")
    receipt = evaluate_quantity_expression(
        QuantityExpressionRequestV1(
            schema_version=_SCHEMA,
            expression_id="redox",
            inputs=tuple(inputs),
            nodes=tuple(nodes),
            output_node_ids=tuple(outputs),
        )
    )
    return {item.quantity_id: item for item in receipt.outputs}


def test_the_redox_operation_owns_the_iupac_sign():
    """A favourable reduction has a negative dG and a positive potential.

    Getting this backwards is the classic redox error and it is silent: the
    magnitude is right and only the sign says the cell runs the other way.
    The operation owns it so no expression has to restate it.
    """

    outputs = _redox(-98.70, 1)
    assert outputs["e_abs_volts"].source_value == pytest.approx(4.28, abs=0.01)
    assert outputs["e_abs_volts"].dimension == ELECTRIC_POTENTIAL

    unfavourable = _redox(+98.70, 1)
    assert unfavourable["e_abs_volts"].source_value == pytest.approx(
        -4.28, abs=0.01
    )


def test_referencing_an_electrode_is_ordinary_subtraction():
    """Which electrode was referenced stays visible in the expression.

    A half reaction whose absolute potential equals the absolute potential of
    the standard hydrogen electrode is zero on the SHE scale, by definition.
    """

    outputs = _redox(-98.70, 1, reference_volts=4.28)
    assert outputs["e_vs_ref_volts"].source_value == pytest.approx(
        0.0, abs=0.01
    )


def test_a_two_electron_process_halves_the_potential():
    one = _redox(-98.70, 1)["e_abs_volts"].source_value
    two = _redox(-98.70, 2)["e_abs_volts"].source_value
    assert two == pytest.approx(one / 2.0, abs=1e-6)


@pytest.mark.parametrize("electrons", [0.0, -1.0])
def test_a_non_positive_electron_count_is_refused(electrons):
    with pytest.raises(Exception, match="positive electron count"):
        _redox(-98.70, electrons)


def test_the_redox_operation_refuses_wrong_input_dimensions():
    """A temperature where an electron count belongs must not silently work."""

    nodes = (
        QuantityExpressionNodeV1(
            node_id="dg", operation="ref", reference="in_dg"
        ),
        QuantityExpressionNodeV1(
            node_id="t", operation="ref", reference="in_t"
        ),
        QuantityExpressionNodeV1(
            node_id="e",
            operation="gibbs_to_redox_potential",
            input_ids=("dg", "t"),
        ),
    )
    with pytest.raises(Exception, match="electron count"):
        evaluate_quantity_expression(
            QuantityExpressionRequestV1(
                schema_version=_SCHEMA,
                expression_id="bad",
                inputs=(
                    _quantity("in_dg", -1.0, "kcal/mol"),
                    _quantity("in_t", 298.15, "K"),
                ),
                nodes=nodes,
                output_node_ids=("e",),
            )
        )
