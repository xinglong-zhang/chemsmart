"""A rejected argument must carry the path, the value, and the rule.

Observed live: a DeepSeek session reproducing a published Hartree-Fock
basis-set series was told

    tool argument analysis_nodes[].outputs[].output_id does not match its
    pattern

which names neither which node, nor which output, nor what the value was, nor
what the pattern requires.  The model had eleven analysis nodes.  An earlier
repair added those three facts to one dataclass validator; this pins them at
the generic argument validator, where they cover every tool at once.
"""

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.tool_runtime import _validate_json_value

_NODE_SCHEMA = {
    "type": "object",
    "additionalProperties": False,
    "required": ["node_id", "outputs"],
    "properties": {
        "node_id": {"type": "string", "pattern": "^[a-z][a-z0-9_]{0,63}$"},
        "outputs": {
            "type": "array",
            "minItems": 1,
            "maxItems": 4,
            "items": {
                "type": "object",
                "additionalProperties": False,
                "required": ["output_id"],
                "properties": {
                    "output_id": {
                        "type": "string",
                        "pattern": "^[a-z][a-z0-9_]{0,63}$",
                    },
                    "kind": {"type": "string", "enum": ["energy", "geometry"]},
                    "weight": {"type": "number", "exclusiveMinimum": 0},
                },
            },
        },
    },
}


def _reject(value, schema=_NODE_SCHEMA, name="analysis_nodes[3]"):
    with pytest.raises(ContractError) as failure:
        _validate_json_value(name, value, schema)
    return str(failure.value)


def test_a_failing_pattern_names_the_index_the_value_and_the_pattern():
    message = _reject(
        {
            "node_id": "cbs_exponential_234",
            "outputs": [
                {"output_id": "e_hf_cbs"},
                {"output_id": "cbs-exponential-234"},
            ],
        }
    )
    assert "analysis_nodes[3].outputs[1].output_id" in message, message
    assert "'cbs-exponential-234'" in message, message
    assert "^[a-z][a-z0-9_]{0,63}$" in message, message


def test_a_failing_enum_lists_what_was_allowed():
    message = _reject(
        {
            "node_id": "cbs",
            "outputs": [{"output_id": "e", "kind": "wavefunction"}],
        }
    )
    assert "'wavefunction'" in message
    assert "['energy', 'geometry']" in message


def test_a_wrong_type_names_both_the_expected_and_the_received_type():
    message = _reject({"node_id": 4, "outputs": [{"output_id": "e"}]})
    assert "must be string" in message
    assert "int" in message


def test_a_missing_key_reports_what_was_supplied_alongside_it():
    message = _reject({"outputs": [{"output_id": "e"}]})
    assert "['node_id']" in message
    assert "['outputs']" in message


def test_an_unexpected_key_lists_the_accepted_ones():
    message = _reject(
        {
            "node_id": "cbs",
            "outputs": [{"output_id": "e", "units": "hartree"}],
        }
    )
    assert "['units']" in message
    assert "'output_id'" in message


def test_a_bound_violation_reports_the_number_and_the_bound():
    message = _reject(
        {"node_id": "cbs", "outputs": [{"output_id": "e", "weight": 0}]}
    )
    assert "must be greater than 0" in message


def test_an_item_count_violation_reports_both_counts():
    empty = _reject({"node_id": "cbs", "outputs": []})
    assert "0 items" in empty and "required 1" in empty
    crowded = _reject(
        {
            "node_id": "cbs",
            "outputs": [{"output_id": f"e{index}"} for index in range(5)],
        }
    )
    assert "5 items" in crowded and "allowed 4" in crowded


def test_a_rejection_never_echoes_a_whole_payload_back():
    """The message is for the caller, not a channel for arbitrary bytes."""

    message = _reject(
        {
            "node_id": "cbs",
            "outputs": [{"output_id": "x" * 5000}],
        }
    )
    assert len(message) < 400, "a rejection must stay bounded"
    assert "..." in message


def test_a_nested_object_of_the_wrong_shape_is_summarized_not_dumped():
    message = _reject(
        {"node_id": "cbs", "outputs": [{"output_id": {"a": 1, "b": 2}}]}
    )
    assert "dict of 2 entries" in message
    assert "must be string" in message


def test_an_empty_registry_names_what_would_fill_it():
    """Observed live: a session asked twice for a counterexample that no
    failure had produced.  Being told the registry was empty did not say what
    fills it, so the retry was identical."""

    from chemsmart.agent.tool_runtime import _REGISTRY_PRODUCERS

    class _Host:
        _get = staticmethod(
            __import__(
                "chemsmart.agent.tool_runtime", fromlist=["x"]
            ).CommandCompiledToolHostV1._get
        )

    with pytest.raises(ContractError) as failure:
        _Host._get({}, "2f8a651d", "counterexample")
    message = str(failure.value)
    assert "'2f8a651d'" in message
    assert "no counterexample is bound yet" in message
    assert "fails inspection, safe preview, or program validation" in message

    # A non-empty registry lists what exists instead; the hint would be noise.
    with pytest.raises(ContractError) as populated:
        _Host._get({"abc": 1}, "xyz", "counterexample")
    assert "bound counterexample IDs: ['abc']" in str(populated.value)
    assert "one is bound" not in str(populated.value)

    assert "trusted artifact" in _REGISTRY_PRODUCERS


def test_every_registry_label_the_host_looks_up_has_a_producer_hint():
    """A label added later without a hint reopens the unactionable case."""

    import re as _re
    from pathlib import Path

    from chemsmart.agent.tool_runtime import _REGISTRY_PRODUCERS

    source = Path("chemsmart/agent/tool_runtime.py").read_text()
    labels = set(
        _re.findall(r'self\._get\(\s*[^,]+,\s*[^,]+,\s*"([a-z ]+)"', source)
    )
    labels |= set(
        _re.findall(r'_get\(\s*self\.[a-z_]+,\s*[^,]+,\s*\n\s*"([a-z ]+)"', source)
    )
    missing = sorted(labels - set(_REGISTRY_PRODUCERS))
    assert not missing, (
        f"these host registries have no producer hint: {missing}; a caller "
        "told one is empty would not know what fills it"
    )


def test_a_field_that_may_not_apply_accepts_an_explicit_null():
    """Observed live: a result-extraction node has no temperature.

    The node contract types temperature_k as ``float | None``, but the schema
    said ``number``, so omitting the key succeeded and saying null failed --
    the explicit statement was the rejected one.
    """

    from chemsmart.agent.tool_specs import build_command_compiled_tool_surface

    definition = next(
        item["function"]
        for item in build_command_compiled_tool_surface().tool_definitions
        if item["function"]["name"] == "plan_scientific_workflow"
    )
    node = definition["parameters"]["properties"]["analysis_nodes"]["items"]
    schema = node["properties"]["temperature_k"]
    assert schema["type"] == ["number", "null"]

    _validate_json_value("analysis_nodes[0].temperature_k", None, schema)
    _validate_json_value("analysis_nodes[0].temperature_k", 298.15, schema)

    with pytest.raises(ContractError, match="must be greater than 0"):
        _validate_json_value("analysis_nodes[0].temperature_k", 0, schema)
    with pytest.raises(ContractError) as wrong:
        _validate_json_value("analysis_nodes[0].temperature_k", "warm", schema)
    assert "['number', 'null']" in str(wrong.value)
    assert "'warm'" in str(wrong.value)


def test_a_union_type_is_not_treated_as_an_unknown_type():
    """Before the union was handled, a list type accepted anything at all."""

    schema = {"type": ["number", "null"]}
    with pytest.raises(ContractError, match="must be one of"):
        _validate_json_value("x", {"not": "a number"}, schema)
