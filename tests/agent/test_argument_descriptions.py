"""A required argument must say what it wants.

Measured before the shared registry existed: 7 of 93 required arguments across
the model-visible surface carried a description, so a model had to infer 86
mandatory fields from their spelling alone.  Two live sessions failed on
``record_scientific_decision`` for exactly that reason -- one omitted
``alternatives`` entirely, another wrote the sentence
``extract_result_quantities (waiting on td_initial.td_result)`` into
``stage_order``, which is an identifier field.

Neither model was being careless.  Nine mandatory scientific-provenance fields
were presented with nothing but their names.
"""

import pytest

from chemsmart.agent.tool_specs import (
    ARGUMENT_DESCRIPTIONS,
    build_approved_execution_tool_surface,
    build_command_compiled_tool_surface,
)

_SURFACES = (
    build_command_compiled_tool_surface,
    build_approved_execution_tool_surface,
)


def _required_fields(surface):
    for item in surface.tool_definitions:
        function = item["function"]
        parameters = function.get("parameters") or {}
        properties = parameters.get("properties") or {}
        for name in parameters.get("required", ()):
            yield function["name"], name, properties.get(name) or {}


@pytest.mark.parametrize("build", _SURFACES)
def test_every_required_argument_on_every_surface_is_described(build):
    undescribed = sorted(
        f"{tool}.{field}"
        for tool, field, schema in _required_fields(build())
        if not str(schema.get("description", "")).strip()
    )
    assert not undescribed, (
        f"{len(undescribed)} required arguments have no description: "
        f"{undescribed[:12]}. A model must infer them from their spelling."
    )


def test_the_field_that_broke_two_live_sessions_now_explains_itself():
    surface = build_command_compiled_tool_surface()
    definition = next(
        item["function"]
        for item in surface.tool_definitions
        if item["function"]["name"] == "record_scientific_decision"
    )
    properties = definition["parameters"]["properties"]
    assert "not a decision" in properties["alternatives"]["description"]
    stage_order = properties["stage_order"]["description"]
    assert "lower-case identifier" in stage_order
    assert "dependency prose" in stage_order, (
        "the model pasted a frontier dependency sentence here; the field must "
        "say not to"
    )


def test_a_shared_name_carries_one_meaning_wherever_it_appears():
    """One meaning per argument name is the point of keying by name.

    The full text may still differ: a field constrained to a public identifier
    also states that format rule, and a field without the pattern does not.
    What must not differ is what the argument is *for*.
    """

    divergent = []
    for build in _SURFACES:
        for tool, field, schema in _required_fields(build()):
            meaning = ARGUMENT_DESCRIPTIONS.get(field)
            if meaning and meaning not in str(schema.get("description", "")):
                divergent.append(f"{tool}.{field}")
    assert not divergent, sorted(set(divergent))


def test_a_tool_specific_description_is_not_overwritten_by_the_shared_one():
    """The shared text is a floor, not a ceiling."""

    surface = build_command_compiled_tool_surface()
    definition = next(
        item["function"]
        for item in surface.tool_definitions
        if item["function"]["name"] == "render_project_yaml"
    )
    sections = definition["parameters"]["properties"]["sections"]
    assert sections["description"], "sections must be described"


#: Arguments described for a tool that is built but deliberately withheld.
#: Keeping the text means re-exposing the tool is a one-line surface change.
#:
#: ``repair_command`` reads a counterexample registry with no producer.  The
#: rest belong to the two tools withheld for the same reason once the
#: reachability rule was checked against every late-bound registry rather than
#: the three someone had remembered to list: ``inspect_calculation_artifact``
#: reads run receipts and settings objects, ``assess_program_candidate`` reads
#: claim evidence, and Runtime V2 binds none of the three.
_WITHHELD = {
    "counterexample_id",
    # inspect_calculation_artifact
    "run_receipt_id",
    "settings_id",
    # assess_program_candidate
    "basis_mode",
    "constraint_kinds",
    "job_families",
    "method_family",
    "method_name",
    "request_id",
    "requested_engine",
    "requested_program",
    "requires_double_hybrid",
    "requires_post_hf",
    "selected_engine",
    "selected_program",
    "source_claim_sha256s",
}


def _property_names(schema, into):
    """Collect property names at every depth the describer reaches.

    `_describe_tool_definitions` descends through nested object properties
    and array items, so a name that appears only inside an item schema --
    `blocked_reason` inside `support_repairs`, for example -- is genuinely
    described and is not stale. A guard that looked one level deep flagged
    exactly such an entry.
    """

    if not isinstance(schema, dict):
        return
    properties = schema.get("properties")
    if isinstance(properties, dict):
        for name, value in properties.items():
            into.add(name)
            _property_names(value, into)
    items = schema.get("items")
    if isinstance(items, dict):
        _property_names(items, into)


def test_the_registry_carries_no_entry_no_surface_uses():
    """A stale entry describes a field that no longer exists."""

    used = set(_WITHHELD)
    for build in _SURFACES:
        for item in build().tool_definitions:
            _property_names(
                item["function"].get("parameters") or {}, used
            )
    stale = sorted(set(ARGUMENT_DESCRIPTIONS) - used)
    assert not stale, f"described but never exposed: {stale}"
