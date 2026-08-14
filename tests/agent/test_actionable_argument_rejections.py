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

    from chemsmart.agent.tool_specs import REGISTRY_PRODUCERS

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

    assert "trusted artifact" in REGISTRY_PRODUCERS


def test_every_registry_label_the_host_looks_up_has_a_producer_hint():
    """A label added later without a hint reopens the unactionable case."""

    import re as _re
    from pathlib import Path

    from chemsmart.agent.tool_specs import REGISTRY_PRODUCERS

    source = Path("chemsmart/agent/tool_runtime.py").read_text()
    labels = set(
        _re.findall(r'self\._get\(\s*[^,]+,\s*[^,]+,\s*"([a-z ]+)"', source)
    )
    labels |= set(
        _re.findall(
            r'_get\(\s*self\.[a-z_]+,\s*[^,]+,\s*\n\s*"([a-z ]+)"', source
        )
    )
    missing = sorted(labels - set(REGISTRY_PRODUCERS))
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


def test_a_withheld_tool_is_not_merely_described_differently():
    """This test replaces one that asserted repair_command's wording.

    That tool is no longer exposed: nothing in the runtime fills the registry
    it reads, so it could never succeed and three wording changes could not
    help.  What survives here is the argument description, which stays correct
    for when a producer is wired -- see test_reachable_tool_surface.
    """

    from chemsmart.agent.tool_specs import ARGUMENT_DESCRIPTIONS

    assert "only after such a failure" in (
        ARGUMENT_DESCRIPTIONS["counterexample_id"]
    )


def test_every_late_bound_tool_states_its_precondition():
    """The failure is a pattern, not one tool.

    repair_command was called with no counterexample bound in four sessions;
    assess_program_candidate was called with no claim evidence bound in
    another.  Both take an argument that indexes a host registry which only
    something else can fill, and in both cases the rule appeared only in the
    rejection.
    """

    from chemsmart.agent.tool_specs import (
        LATE_BOUND_ARGUMENTS,
        REGISTRY_PRODUCERS,
        build_command_compiled_tool_surface,
    )

    silent = []
    for item in build_command_compiled_tool_surface().tool_definitions:
        function = item["function"]
        properties = (function.get("parameters") or {}).get("properties") or {}
        if any(name in LATE_BOUND_ARGUMENTS for name in properties):
            if "PRECONDITION" not in function["description"]:
                silent.append(function["name"])
    assert (
        not silent
    ), f"these tools take a late-bound argument and do not say so: {silent}"

    # Every late-bound argument must name a registry the producers table knows,
    # so the sentence before the call and the rejection after it agree.
    unknown = sorted(
        set(LATE_BOUND_ARGUMENTS.values()) - set(REGISTRY_PRODUCERS)
    )
    assert not unknown, unknown


def test_the_precondition_and_the_rejection_come_from_one_table():
    """Two texts that could drift would eventually contradict each other."""

    from chemsmart.agent import tool_runtime
    from chemsmart.agent.tool_specs import (
        REGISTRY_PRODUCERS,
        build_command_compiled_tool_surface,
    )

    assert tool_runtime.REGISTRY_PRODUCERS is REGISTRY_PRODUCERS

    definition = next(
        item["function"]
        for item in build_command_compiled_tool_surface().tool_definitions
        if item["function"]["name"] == "preview_command"
    )
    producer = REGISTRY_PRODUCERS["canonical invocation"]
    assert producer in definition["description"]

    with pytest.raises(ContractError) as failure:
        tool_runtime.CommandCompiledToolHostV1._get(
            {}, "abc", "canonical invocation"
        )
    assert producer in str(failure.value)


def _count_expression_intent(unit):
    from chemsmart.agent.scientific_toolchain import (
        AnalysisNodeIntentV1,
        AnalysisOutputIntentV1,
    )

    return AnalysisNodeIntentV1(
        node_id="derive-count",
        analysis_kind="quantity_expression",
        dependencies=(),
        inputs=(),
        selectors=(),
        outputs=(
            AnalysisOutputIntentV1(
                output_id="n_imaginary",
                quantity_kind="count",
                unit=unit,
            ),
        ),
        expression_nodes=(
            {
                "node_id": "n_imaginary",
                "operation": "literal",
                "literal_value": 0,
                "literal_unit": "1",
            },
        ),
        expression_output_node_ids=("n_imaginary",),
        temperature_k=None,
        pressure_atm=None,
        support_state="planned",
        blocked_reason="",
    )


def test_a_blank_unit_names_the_output_and_the_value_to_use():
    """Seen in three separate live cases, always on a count.

    "analysis output unit must not be empty" named neither which output was
    blank nor what a dimensionless quantity should carry, so a model that had
    written a count with no unit had to guess between "", "none", "count" and
    "1".
    """

    from chemsmart.agent.analysis_nodes import AnalysisContractError
    from chemsmart.agent.scientific_toolchain import (
        ScientificToolchainContractError,
    )

    with pytest.raises(ScientificToolchainContractError) as failure:
        _count_expression_intent("")
    message = str(failure.value)
    assert "'n_imaginary'" in message
    assert "count" in message
    assert "'1'" in message

    assert _count_expression_intent("1").outputs[0].unit == "1"

    # The other plane must answer the same mistake the same way. Both now read
    # one contract-level constant, so this asserts a shared object rather than
    # two copies that happen to be equal today.
    from chemsmart.agent import analysis_nodes, scientific_toolchain
    from chemsmart.agent._contracts import DIMENSIONLESS_UNIT_HINT

    assert analysis_nodes.DIMENSIONLESS_UNIT_HINT is DIMENSIONLESS_UNIT_HINT
    assert (
        scientific_toolchain.DIMENSIONLESS_UNIT_HINT is DIMENSIONLESS_UNIT_HINT
    )

    from chemsmart.agent.analysis_nodes import AnalysisOutputSpecV1

    with pytest.raises(AnalysisContractError) as other_plane:
        AnalysisOutputSpecV1(
            quantity_id="n_imaginary", unit="", dimension=(0,) * 6
        )
    assert DIMENSIONLESS_UNIT_HINT in str(other_plane.value)
    assert AnalysisContractError is not None


def test_an_unsupported_count_unit_is_rejected_before_plan_sealing():
    from chemsmart.agent.scientific_toolchain import (
        ScientificToolchainContractError,
    )

    with pytest.raises(ScientificToolchainContractError) as failure:
        _count_expression_intent("count")

    message = str(failure.value)
    assert "analysis output 'n_imaginary' (count)" in message
    assert "unsupported unit 'count'" in message
    assert "counts use physical unit '1'" in message
