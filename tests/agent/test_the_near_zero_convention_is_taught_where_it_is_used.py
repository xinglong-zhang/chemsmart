"""The 20 cm^-1 convention existed and reached 1 of 20 live sessions.

The near-zero-frequency convention was documented in the thermochemistry
module, exposed as the typed quantity near_zero_mode_count, and explained
in docs -- and in the four-cycle real-ability protocol exactly one of
twenty sessions applied it in a validation rule. Everyone else wrote a
strict >= 0 rule, which would have failed a physically acceptable
benchmark minimum at -17.5 cm^-1. Knowledge that is not visible at the
decision point does not exist; the tool descriptions the model actually
reads now teach the default and the user-override rule.
"""

from __future__ import annotations

from chemsmart.agent.skills import resolve_skill
from chemsmart.agent.tool_specs import build_command_compiled_tool_surface


def _description(surface, tool_name, *path):
    spec = next(
        item["function"]
        for item in surface.tool_definitions
        if item["function"]["name"] == tool_name
    )
    node = spec["parameters"]["properties"]
    for key in path:
        node = node[key]
    return node


def test_the_validation_rule_guidance_teaches_both_questions():
    surface = build_command_compiled_tool_surface()
    spec = next(
        item["function"]
        for item in surface.tool_definitions
        if item["function"]["name"] == "plan_scientific_workflow"
    )
    text = str(spec["parameters"])

    assert "20 cm-1" in text
    assert "near_zero_mode_count" in text
    assert "overrides" in text
    assert "say which" in text.lower() or "which you asked" in text.lower()


def test_the_thermochemistry_tool_names_the_convention_and_the_override():
    surface = build_command_compiled_tool_surface()
    spec = next(
        item["function"]
        for item in surface.tool_definitions
        if item["function"]["name"] == "derive_thermochemistry"
    )
    text = spec["description"]

    assert "near_zero_mode_count" in text
    assert "20 cm-1" in text
    assert "overrides" in text


def test_the_conventions_skill_states_the_default_and_the_override():
    resolved = resolve_skill("scientific-conventions")
    body = resolved.body

    assert "20 cm" in body
    assert "near_zero_mode_count" in body
    assert "default" in body.lower()
    assert "no imaginary modes, never" in body
