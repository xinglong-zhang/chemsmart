"""The host must state the workflow structure it already knows.

A model that has lost track of its own DAG re-plans, and re-planning is where a
stage carrying the answer gets dropped.  These tests pin that the harness says
what each node waits for -- and that the statement is derived, never asserted
by the model.
"""

from dataclasses import dataclass

import pytest

from chemsmart.agent._contracts import ContractError, canonical_data
from chemsmart.agent.feedback import CAUSAL_FEEDBACK_V1, project_tool_feedback
from chemsmart.agent.workflow_context import (
    WorkflowNodeContextV1,
    project_workflow_context,
)


@dataclass(frozen=True)
class _Input:
    binding_id: str
    artifact_class: str = "geometry_xyz"
    artifact_id: str = ""
    producer_node_id: str = ""
    producer_output_id: str = ""


@dataclass(frozen=True)
class _Output:
    output_id: str
    artifact_class: str = "geometry_xyz"


@dataclass(frozen=True)
class _Node:
    node_id: str
    dependencies: tuple = ()
    inputs: tuple = ()
    expected_outputs: tuple = ()


def _fan_out_fan_in():
    """optimize -> two single points -> one extrapolation."""

    return [
        _Node(
            "opt",
            (),
            (_Input("geom", artifact_id="start.xyz"),),
            (_Output("opt_geom"),),
        ),
        _Node(
            "sp_qz",
            ("opt",),
            (
                _Input(
                    "geom",
                    producer_node_id="opt",
                    producer_output_id="opt_geom",
                ),
            ),
            (_Output("e_qz", "energy"),),
        ),
        _Node(
            "sp_5z",
            ("opt",),
            (
                _Input(
                    "geom",
                    producer_node_id="opt",
                    producer_output_id="opt_geom",
                ),
            ),
            (_Output("e_5z", "energy"),),
        ),
        _Node(
            "cbs",
            ("sp_qz", "sp_5z"),
            (
                _Input(
                    "qz",
                    "energy",
                    producer_node_id="sp_qz",
                    producer_output_id="e_qz",
                ),
                _Input(
                    "fz",
                    "energy",
                    producer_node_id="sp_5z",
                    producer_output_id="e_5z",
                ),
            ),
            (_Output("e_cbs", "energy"),),
        ),
    ]


def _project(materialized=("start.xyz",)):
    return project_workflow_context(
        workflow_id="wf.probe",
        nodes=_fan_out_fan_in(),
        materialized_artifact_ids=set(materialized),
    )


def test_only_the_root_is_ready_before_anything_has_run():
    projection = _project()
    assert projection.ready_node_ids == ("opt",)
    assert projection.waiting_node_ids == ("cbs", "sp_5z", "sp_qz")


def test_a_waiting_node_names_the_exact_upstream_output_it_needs():
    """'Unresolved' is not information; naming the producer output is."""

    node = _project().node("sp_qz")
    assert node.state == "waiting"
    assert node.waiting_on == ("opt",)
    assert node.unsatisfied_inputs[0].producer_node_id == "opt"
    assert node.unsatisfied_inputs[0].producer_output_id == "opt_geom"
    assert "opt.opt_geom" in node.reason


def test_fan_out_is_visible_from_the_producer():
    """Branching is a property of the DAG the model should not have to recall."""

    assert _project().node("opt").dependents == ("sp_5z", "sp_qz")


def test_the_frontier_advances_as_artifacts_materialize():
    """The projection is a function of state, so it moves without re-planning."""

    nodes = _fan_out_fan_in()
    # Rebind the two single points to an artifact the optimization produced.
    nodes[1] = _Node(
        "sp_qz",
        ("opt",),
        (
            _Input(
                "geom",
                artifact_id="opt.xyz",
                producer_node_id="opt",
                producer_output_id="opt_geom",
            ),
        ),
        (_Output("e_qz", "energy"),),
    )
    projection = project_workflow_context(
        workflow_id="wf.probe",
        nodes=nodes,
        materialized_artifact_ids={"start.xyz", "opt.xyz"},
    )
    assert "sp_qz" in projection.ready_node_ids
    assert "sp_5z" in projection.waiting_node_ids


def test_a_validated_handoff_advances_the_unchanged_producer_edge():
    """Re-planning must not erase an OPT geometry already handed to HESS."""

    projection = project_workflow_context(
        workflow_id="wf.probe",
        nodes=_fan_out_fan_in(),
        materialized_artifact_ids={"start.xyz", "opt.xyz"},
        materialized_producer_inputs={
            ("sp_qz", "geom", "opt", "opt_geom"),
        },
        completed_node_ids={"opt"},
    )

    assert projection.completed_node_ids == ("opt",)
    assert projection.node("opt").state == "completed"
    assert projection.node("sp_qz").state == "ready"
    assert projection.node("sp_5z").state == "waiting"


def test_control_dependency_waits_for_completion_without_a_data_edge():
    nodes = (
        _Node("first", (), (), (_Output("done", "status"),)),
        _Node("second", ("first",), (), (_Output("next", "status"),)),
    )

    waiting = project_workflow_context(
        workflow_id="wf.control",
        nodes=nodes,
    )
    assert waiting.ready_node_ids == ("first",)
    assert waiting.node("second").waiting_on == ("first",)
    assert waiting.node("second").reason == "waiting for first completion"

    advanced = project_workflow_context(
        workflow_id="wf.control",
        nodes=nodes,
        completed_node_ids={"first"},
    )
    assert advanced.ready_node_ids == ("second",)


def test_an_input_with_no_producer_and_no_artifact_is_blocked_not_waiting():
    """Waiting resolves itself; a dangling input never will, so say so."""

    projection = project_workflow_context(
        workflow_id="wf.probe",
        nodes=[_Node("orphan", (), (_Input("geom"),), (_Output("out"),))],
    )
    node = projection.node("orphan")
    assert node.state == "blocked"
    assert "no artifact and no producer" in node.reason


def test_an_explicitly_blocked_node_keeps_its_reason():
    projection = project_workflow_context(
        workflow_id="wf.probe",
        nodes=_fan_out_fan_in(),
        materialized_artifact_ids={"start.xyz"},
        blocked_reasons={"cbs": "no extrapolation executor exists"},
    )
    assert projection.blocked_node_ids == ("cbs",)
    assert projection.node("cbs").reason == "no extrapolation executor exists"


def test_a_ready_node_cannot_carry_unsatisfied_inputs():
    with pytest.raises(ContractError, match="ready node"):
        WorkflowNodeContextV1(
            node_id="x",
            state="ready",
            waiting_on=(),
            unsatisfied_inputs=(
                type(
                    "I",
                    (),
                    {
                        "binding_id": "a",
                        "artifact_class": "",
                        "producer_node_id": "",
                        "producer_output_id": "",
                    },
                )(),
            ),
            produces=(),
            dependents=(),
            reason="",
        )


def test_a_waiting_node_must_say_what_it_waits_for():
    with pytest.raises(ContractError, match="must name what it waits for"):
        WorkflowNodeContextV1(
            node_id="x",
            state="waiting",
            waiting_on=(),
            unsatisfied_inputs=(),
            produces=(),
            dependents=(),
            reason="",
        )


def test_the_projection_survives_causal_feedback_with_its_reasons_intact():
    """A causal arm that silently dropped this would measure nothing."""

    envelope = {
        "tool": "plan_command_workflow",
        "status": "ok",
        "result": canonical_data(
            {
                "workflow_context": _project(),
                # The model's own submitted DAG is echoed back and dropped.
                "nodes": [{"node_id": "echo"}],
            }
        ),
    }
    projected = project_tool_feedback(
        tool="plan_command_workflow",
        result=envelope,
        mode=CAUSAL_FEEDBACK_V1,
    ).content["result"]

    assert "nodes" not in projected, "the echo must still be dropped"
    context = projected["workflow_context"]
    assert len(context["nodes"]) == 4
    waiting = next(
        item for item in context["nodes"] if item["node_id"] == "sp_qz"
    )
    assert "opt.opt_geom" in waiting["reason"]


def test_a_rejected_node_names_the_node_and_the_offending_value():
    """A live run hit the old message, retried unchanged, and stalled.

    The runtime sorts these lists before building the node, so duplication was
    the only reachable cause -- yet the message said only "must be sorted and
    unique", naming neither the node nor the value.
    """

    from chemsmart.agent.workflows import (
        ArtifactInputIntentV1,
        ArtifactOutputIntentV1,
        CommandNodeIntentV1,
    )

    def _input(binding_id):
        return ArtifactInputIntentV1(
            binding_id=binding_id,
            artifact_class="geometry_xyz",
            artifact_id="start.xyz",
            producer_node_id="",
            producer_output_id="",
        )

    common = dict(
        node_id="sp_ccpvqz",
        program="pyscf",
        jobtype="sp",
        project_role="sp",
        dependencies=(),
        expected_outputs=(
            ArtifactOutputIntentV1(output_id="e", artifact_class="energy"),
        ),
        unresolved_fields=(),
    )

    with pytest.raises(ContractError) as duplicate:
        CommandNodeIntentV1(
            **common, inputs=(_input("geometry"), _input("geometry"))
        )
    message = str(duplicate.value)
    assert "sp_ccpvqz" in message, "the failing node must be named"
    assert "geometry" in message, "the offending value must be named"
    assert "repeats" in message

    with pytest.raises(ContractError, match="out of order"):
        CommandNodeIntentV1(**common, inputs=(_input("zz"), _input("aa")))

    with pytest.raises(ContractError, match="declares no expected output"):
        CommandNodeIntentV1(
            **{**common, "expected_outputs": ()}, inputs=(_input("g"),)
        )
