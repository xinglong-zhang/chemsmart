"""A refused selector must name the tool that does provide the quantity.

Preparing a case built on a published frequency scaling factor exposed this:
`zero_point_energy` is not a selector, and the refusal said only

    unsupported quantity selector: 'zero_point_energy'

The quantity is not missing.  ChemSmart's RRHO engine computes it, and the
`derive_thermochemistry` tool returns it -- one tool away.  A caller told only
that the selector is unsupported goes looking for a selector that will never
exist, or rebuilds the quantity from frequencies by hand.

This is the same courtesy the host registries already extend when an ID is
unknown: say what would provide it.
"""

import pytest

from chemsmart.analysis.result_quantities import (
    QUANTITIES_FROM_ANOTHER_TOOL,
    SUPPORTED_SELECTORS,
    QuantityContractError,
    QuantitySelectorV1,
)


@pytest.mark.parametrize("selector", sorted(QUANTITIES_FROM_ANOTHER_TOOL))
def test_a_quantity_another_tool_owns_names_that_tool(selector):
    with pytest.raises(QuantityContractError) as failure:
        QuantitySelectorV1(quantity_id="q", selector=selector)
    message = str(failure.value)
    assert QUANTITIES_FROM_ANOTHER_TOOL[selector] in message
    assert "not by result extraction" in message


def test_an_unknown_selector_still_lists_what_is_supported():
    """Redirection must not replace the plain listing for a true typo."""

    with pytest.raises(QuantityContractError) as failure:
        QuantitySelectorV1(quantity_id="q", selector="bogus_quantity")
    message = str(failure.value)
    assert "supported selectors" in message
    assert "energy" in message


def test_the_redirect_table_never_shadows_a_real_selector():
    """A quantity that is selectable must not be redirected elsewhere."""

    overlap = sorted(set(QUANTITIES_FROM_ANOTHER_TOOL) & SUPPORTED_SELECTORS)
    assert (
        not overlap
    ), f"{overlap} are selectable here and would be wrongly redirected"


def test_every_redirect_target_is_a_tool_the_model_can_call():
    from chemsmart.agent.tool_specs import build_command_compiled_tool_surface

    exposed = {
        item["function"]["name"]
        for item in build_command_compiled_tool_surface().tool_definitions
    }
    missing = sorted(set(QUANTITIES_FROM_ANOTHER_TOOL.values()) - exposed)
    assert (
        not missing
    ), f"the refusal points at {missing}, which the model cannot call"


def test_the_redirected_quantities_are_the_ones_that_tool_returns():
    """A redirect to a tool that does not produce the quantity is worse than
    none, so the table is checked against the engine's own output names."""

    import inspect

    from chemsmart.analysis import result_quantities

    source = inspect.getsource(result_quantities)
    body = source[source.index("energy_values = {") :]
    for selector in QUANTITIES_FROM_ANOTHER_TOOL:
        if selector == "entropy":
            assert 'quantity_id="entropy"' in source, selector
        else:
            assert f'"{selector}"' in body, selector
