"""Independent slips in one payload must cost one round trip, not several.

A workflow DAG is authored in a single ~13 KB payload.  The validator used to
raise on the first violation it found, so a session that mistyped an identifier
in node seven, a rule id in node eight and a unit in node nine paid three full
resubmissions of the entire graph to learn three unrelated facts.  Across the
recorded sessions, most multi-failure streaks on ``plan_scientific_workflow``
carried violations on entirely different fields -- independent, not a cascade
that had to be unwound in order.

These tests execute the validator against real tool schemas rather than reading
the code: the first proves several independent slips arrive together, the
second proves a lone violation keeps its exact previous wording, and the third
proves descent stops at a wrongly-typed value instead of piling on rules that
value could never satisfy.
"""

from __future__ import annotations

import pytest

from chemsmart.agent._contracts import ContractError
from chemsmart.agent.tool_runtime import (
    _collect_json_violations,
    _validate_json_value,
    _validate_tool_arguments,
)
from chemsmart.agent.tool_specs import build_command_compiled_tool_surface

_SURFACE = build_command_compiled_tool_surface()


def _plan_arguments(**overrides):
    """A minimal analysis-only DAG, shaped like the ones sessions author."""

    node = {
        "node_id": "extract-a",
        "analysis_kind": "result_extraction",
        "dependencies": [],
        "inputs": [],
        "selectors": [],
        "outputs": [
            {
                "output_id": "e-total",
                "quantity_kind": "electronic_energy",
                "unit": "hartree",
            }
        ],
        "expression_nodes": [],
        "expression_output_node_ids": [],
        "support_state": "planned",
        "blocked_reason": "",
    }
    node.update(overrides)
    return {
        "plan_id": "acidity",
        "workflow_id": "acidity",
        "calculation_nodes": [],
        "analysis_nodes": [node],
        "required_output_ids": ["dg"],
    }


def test_three_independent_slips_arrive_in_one_rejection():
    """The streak this change exists for, reproduced in one payload."""

    arguments = _plan_arguments(
        node_id="Extract-A",  # upper case: pattern
        outputs=[
            {
                "output_id": "dS",  # mixed case: pattern
                "quantity_kind": "entropy",
                "unit": "cal/(mol*K)",
            }
        ],
        support_state="resolvable",  # not an analysis support state: enum
    )

    with pytest.raises(ContractError) as excinfo:
        _validate_tool_arguments(
            _SURFACE, "plan_scientific_workflow", arguments
        )

    message = str(excinfo.value)
    assert message.startswith("3 arguments are invalid:")
    # Every slip is named with its own path, so one revision can fix all three.
    assert "analysis_nodes[0].node_id" in message
    assert "analysis_nodes[0].outputs[0].output_id" in message
    assert "analysis_nodes[0].support_state" in message


def test_one_slip_keeps_its_exact_previous_wording():
    """Nothing that reads these messages learns a second shape."""

    arguments = _plan_arguments(node_id="Extract-A")

    with pytest.raises(ContractError) as excinfo:
        _validate_tool_arguments(
            _SURFACE, "plan_scientific_workflow", arguments
        )

    message = str(excinfo.value)
    assert not message.startswith("1 arguments are invalid")
    assert message.startswith("tool argument analysis_nodes[0].node_id is ")
    assert "does not match the required pattern" in message


def test_a_wrongly_typed_value_is_reported_once_not_piled_on():
    """A value being replaced wholesale needs one reason, not four."""

    findings: list[str] = []
    schema = {
        "type": "string",
        "pattern": "^[a-z]+$",
        "enum": ["alpha", "beta"],
    }
    _collect_json_violations("node.name", 17, schema, findings)

    assert len(findings) == 1
    assert "must be string" in findings[0]


def test_the_single_raising_entry_point_still_raises_first():
    """Direct callers of the raising helper are unaffected."""

    schema = {"type": "string", "pattern": "^[a-z][a-z0-9_.-]*$"}
    _validate_json_value("x", "start-xyz", schema)
    with pytest.raises(ContractError, match="required pattern"):
        _validate_json_value("x", "Start.XYZ", schema)


def test_an_unknown_argument_no_longer_hides_the_rest():
    """The unknown name and the real slips are reported together."""

    arguments = _plan_arguments(node_id="Extract-A")
    arguments["not_a_field"] = 1

    with pytest.raises(ContractError) as excinfo:
        _validate_tool_arguments(
            _SURFACE, "plan_scientific_workflow", arguments
        )

    message = str(excinfo.value)
    assert "does not accept ['not_a_field']" in message
    assert "analysis_nodes[0].node_id" in message
